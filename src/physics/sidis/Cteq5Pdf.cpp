// src/physics/sidis/Cteq5Pdf.cpp — cteq5/Ctq5Pdf.f (the CTEQ5 parton distributions, table + polynomial interpolation
// in x and Q). See include/simc/physics/sidis/Cteq5Pdf.h.
// Operation order follows the Fortran source (D-024). All arithmetic is REAL(8) (the Implicit Double Precision lines of
// Ctq5Pdf.f are commented out; -fdefault-real-8).
#include "simc/physics/sidis/Cteq5Pdf.h"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "simc/core/Errors.h"

namespace simc {
namespace sidis {

const std::array<const char*, 9> Cteq5Pdf::kFiles = {"cteq5m.tbl",  "cteq5d.tbl",  "cteq5l.tbl",
                                                     "cteq5hj.tbl", "cteq5hq.tbl", "cteq5f3.tbl",
                                                     "cteq5f4.tbl", "cteq5m1.tbl", "ctq5hq1.tbl"};

namespace {

// Fortran list-directed numbers (Read(Nu,*)): blank/comma separated items; the remainder of the last record touched is
// discarded. The CTEQ tables use no repeat counts (r*c), slashes or null values; such items are reported as errors.
class ListReader {
public:
    explicit ListReader(const std::vector<std::string>& lines) : lines_(lines) {}
    /// Fortran Read(Nu,'(A)') Line: skip one record.
    bool skipRecord() {
        if (line_ >= lines_.size()) return false;
        ++line_;
        pos_ = 0;
        return true;
    }
    /// Next list-directed item (crossing records). False at end of file.
    bool next(std::string& tok) {
        for (;;) {
            if (!started_) {
                if (line_ >= lines_.size()) return false;
                started_ = true;
                pos_ = 0;
            }
            const std::string& l = lines_[line_];
            while (pos_ < l.size() && (l[pos_] == ' ' || l[pos_] == '\t' || l[pos_] == ',' || l[pos_] == '\r')) ++pos_;
            if (pos_ < l.size()) {
                const std::size_t b = pos_;
                while (pos_ < l.size() && l[pos_] != ' ' && l[pos_] != '\t' && l[pos_] != ',' && l[pos_] != '\r') ++pos_;
                tok = l.substr(b, pos_ - b);
                return true;
            }
            ++line_;   // continue on the next record
            pos_ = 0;
            if (line_ >= lines_.size()) return false;
        }
    }
    /// End of one READ statement: the rest of the current record is discarded.
    void endRead() {
        if (started_) {
            ++line_;
            pos_ = 0;
            started_ = false;
        }
    }
    /// Start of a READ statement.
    void beginRead() { started_ = false; }

private:
    const std::vector<std::string>& lines_;
    std::size_t line_{0}, pos_{0};
    bool started_{false};
};

bool parseReal(const std::string& t, double& v) {
    if (t.find_first_of("*/") != std::string::npos) return false;
    std::string s = t;
    for (char& c : s)
        if (c == 'd' || c == 'D') c = 'E';   // Fortran D exponent
    errno = 0;
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end && *end == '\0' && end != s.c_str();
}

bool parseInt(const std::string& t, int& v) {
    char* end = nullptr;
    const long l = std::strtol(t.c_str(), &end, 10);
    if (!end || *end != '\0' || end == t.c_str()) return false;
    v = static_cast<int>(l);
    return true;
}

bool fail(std::string* err, const std::string& m) {
    if (err) *err = m;
    return false;
}

}  // namespace

// Fortran: Ctq5Pdf.f:POLINT (l.297-334)
double polint(const double* xa, const double* ya, int n, double x, double& dy) {
    constexpr int NMAX = 10;
    double c[NMAX + 1], d[NMAX + 1];   // C(1:NMAX), D(1:NMAX), 1-based
    auto XA = [&](int i) { return xa[i - 1]; };
    auto YA = [&](int i) { return ya[i - 1]; };
    int ns = 1;
    double dif = std::abs(x - XA(1));
    for (int i = 1; i <= n; ++i) {
        const double dift = std::abs(x - XA(i));
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = YA(i);
        d[i] = YA(i);
    }
    double y = YA(ns);
    ns = ns - 1;
    for (int m = 1; m <= n - 1; ++m) {
        for (int i = 1; i <= n - m; ++i) {
            const double ho = XA(i) - x;
            const double hp = XA(i + m) - x;
            const double w = c[i + 1] - d[i];
            double den = ho - hp;
            if (den == 0.) throw FortranStop("Ctq5Pdf.f:POLINT l.322: PAUSE (two equal abscissae in the CTEQ5 table)");
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        if (2 * ns < n - m) {
            dy = c[ns + 1];
        } else {
            dy = d[ns];
            ns = ns - 1;
        }
        y = y + dy;
    }
    return y;
}

// Fortran: Ctq5Pdf.f:SetCtq5 (l.199-237) + ReadTbl (l.239-276)
bool Cteq5Pdf::load(const std::string& dir, int iset, std::string* err) {
    if (iset < 1 || iset > static_cast<int>(kFiles.size()))
        return fail(err, "Invalid Iset number in SetCtq5 : " + std::to_string(iset) + " (Ctq5Pdf.f:221 STOP)");
    const std::string path = dir + "/" + kFiles[static_cast<std::size_t>(iset - 1)];
    std::ifstream in(path);
    if (!in) return fail(err, " Data file " + path + " cannot be opened in SetCtq5!! (Ctq5Pdf.f:234 STOP)");
    std::vector<std::string> lines;
    for (std::string l; std::getline(in, l);) lines.push_back(l);

    ListReader r(lines);
    std::string tok;
    auto readReals = [&](std::vector<double*> dst) -> bool {
        r.beginRead();
        for (double* p : dst) {
            if (!r.next(tok) || !parseReal(tok, *p)) return false;
        }
        r.endRead();
        return true;
    };
    // l.250-256
    if (!r.skipRecord() || !r.skipRecord()) return fail(err, path + ": short header");
    double Dr = 0, Fl = 0;
    {
        std::vector<double*> v{&Dr, &Fl, &Al()};
        for (auto& a : Amass) v.push_back(&a);
        if (!readReals(v)) return fail(err, path + ": cannot read Dr, Fl, Al, Amass(1:6)");
    }
    Iorder = static_cast<int>(std::lround(Dr));   // Nint
    Nfl = static_cast<int>(std::lround(Fl));
    Alambda = Al();
    // l.258-259
    if (!r.skipRecord()) return fail(err, path + ": short header");
    r.beginRead();
    if (!r.next(tok) || !parseInt(tok, Nx) || !r.next(tok) || !parseInt(tok, Nt) || !r.next(tok) || !parseInt(tok, NfMx))
        return fail(err, path + ": cannot read NX, NT, NfMx");
    r.endRead();
    if (Nx < 0 || Nx > MXX || Nt < 0 || Nt > MXQ || NfMx < 0 || NfMx > MXF)
        return fail(err, path + ": NX, NT, NfMx outside the Fortran dimensions MXX, MXQ, MXF");
    // l.261-262
    if (!r.skipRecord()) return fail(err, path + ": short header");
    {
        std::vector<double*> v{&Qini, &Qmax};
        for (int i = 0; i <= Nt; ++i) v.push_back(&QL(i));
        if (!readReals(v)) return fail(err, path + ": cannot read QINI, QMAX, QL(0:NT)");
    }
    // l.264-265
    if (!r.skipRecord()) return fail(err, path + ": short header");
    {
        std::vector<double*> v{&Xmin};
        for (int i = 0; i <= Nx; ++i) v.push_back(&XV(i));
        if (!readReals(v)) return fail(err, path + ": cannot read XMIN, XV(0:NX)");
    }
    // l.267-269
    for (int iq = 0; iq <= Nt; ++iq) QL(iq) = std::log(QL(iq) / Al());
    // l.274-279: Read (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts): an end of file or a bad item ends the read silently;
    // the values read so far are kept, the others keep their COMMON value (0).
    const int nblk = (Nx + 1) * (Nt + 1);
    npts_ = nblk * (NfMx + 3);
    if (npts_ > MXPQX) return fail(err, path + ": Npts > MXPQX");
    for (int i = 1; i <= MXPQX; ++i) UPD(i) = 0.0;
    if (!r.skipRecord()) return fail(err, path + ": no parton table");
    r.beginRead();
    npts_read_ = 0;
    for (int i = 1; i <= npts_; ++i) {
        double v = 0;
        if (!r.next(tok) || !parseReal(tok, v)) break;
        UPD(i) = v;
        ++npts_read_;
    }
    iset_ = iset;
    warn = true;
    first = true;
    return true;
}

// Fortran: Ctq5Pdf.f:Ctq5Pdf (l.61-95)
double Cteq5Pdf::pdf(int iparton, double x, double& q, std::ostream* log) {
    if (!loaded()) throw std::logic_error("Cteq5Pdf::pdf before load (SetCtq5)");
    if (x < 0. || x > 1.) {   // X .lt. 0D0 .or. X .gt. 1D0 (binary128 compare of a REAL(8): exact)
        if (log) *log << " X out of range in Ctq5Pdf: " << x << "\n";
        throw FortranStop("Ctq5Pdf.f:74 STOP: X out of range in Ctq5Pdf: " + std::to_string(x));
    }
    if (q < Alambda) {
        if (log) *log << " Q out of range in Ctq5Pdf: " << q << "\n Setting to Alambda\n";
        q = Alambda;   // the caller's argument is changed (l.79)
    }
    if (iparton < -NfMx || iparton > NfMx) {
        if (warn) {
            warn = false;
            if (log) *log << " Warning: Iparton out of range in Ctq5Pdf: " << iparton << "\n";
        }
        return 0.;
    }
    double r = partonX(iparton, x, q, log);
    if (r < 0.) r = 0.;
    return r;
}

// Fortran: Ctq5Pdf.f:PartonX (l.97-196)
double Cteq5Pdf::partonX(int iprtn, double x, double q, std::ostream* log) {
    const double qg = std::log(q / Al());   // Work with Log (Q)
    // Find lower end of interval containing X
    int jl = -1, ju = Nx + 1;
    while (ju - jl > 1) {
        const int jm = (ju + jl) / 2;
        if (x > XV(jm))
            jl = jm;
        else
            ju = jm;
    }
    int jx = jl - (M - 1) / 2;
    if (x < Xmin && first) {
        first = false;
        if (log) {
            char b[160];
            std::snprintf(b, sizeof b, " WARNING: X << Xmin, extrapolation used; X, Xmin =%12.4E%12.4E", x, Xmin);
            *log << b << "\n";
        }
        if (jx < 0) jx = 0;
    } else if (jx > Nx - M) {
        jx = Nx - M;
    }
    // Find the interval where Q lies
    jl = -1;
    ju = Nt + 1;
    while (ju - jl > 1) {
        const int jm = (ju + jl) / 2;
        if (qg > QL(jm))
            jl = jm;
        else
            ju = jm;
    }
    int jq = jl - (M - 1) / 2;
    if (jq < 0) {
        jq = 0;
    } else if (jq > Nt - M) {
        jq = Nt - M;
        if (q > Qmax && log) {
            char b[160];
            std::snprintf(b, sizeof b, " WARNING: Q > Qmax, extrapolation used; Q, Qmax =%12.4E%12.4E", q, Qmax);
            *log << b << "\n";
        }
    }
    // X <= XV(0) after the first call leaves Jx = -1 (x = 0, or a NaN x): the Fortran then reads XV(-1) = Al and, for
    // JFL = Jq = 0, UPD(0) = QL(MXQ) - the neighbouring words of COMMON /CtqPar1/, reproduced by the storage layout.
    const int ip = iprtn >= 3 ? -iprtn : iprtn;
    // Find the off-set in the linear array Upd
    const int jfl = ip + NfMx;
    const int j0 = (jfl * (Nt + 1) + jq) * (Nx + 1) + jx;
    // Now interpolate in x for M1 Q's
    double fq[M1], df = 0;
    for (int iq = 1; iq <= M1; ++iq) {
        const int j1 = j0 + (Nx + 1) * (iq - 1) + 1;
        fq[iq - 1] = polint(&XV(jx), &UPD(j1), M1, x, df);
    }
    // Finish off by interpolating in Q
    double ddf = 0;
    return polint(&QL(jq), fq, M1, qg, ddf);
}

}  // namespace sidis
}  // namespace simc
