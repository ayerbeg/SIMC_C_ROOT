// include/simc/core/Errors.h
//
// The SIMC++ exception hierarchy (task P7.4: one hierarchy instead of a FortranStop/NotPorted per module).
//
//   simc::Error         base of every error that SIMC++ reports as a run error (std::runtime_error)
//   simc::FortranStop   a Fortran STOP statement was reached in a ported routine; the message names the routine
//                       (file:line) and carries the Fortran text. The run cannot continue, exactly as the Fortran.
//   simc::NotPorted     a feature that the Fortran supports but SIMC++ does not (yet): reaction families, arms,
//                       options. The message names the feature and the planned task.
//
// Programming errors (a stage used before init(), a missing collaborator) are std::logic_error; bad user input in the
// JSON pipeline is std::invalid_argument. The simc++ CLI turns every std::exception into exit status 3.
#ifndef SIMC_CORE_ERRORS_H
#define SIMC_CORE_ERRORS_H

#include <stdexcept>
#include <string>

namespace simc {

class Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class FortranStop : public Error {
public:
    using Error::Error;
};

class NotPorted : public Error {
public:
    using Error::Error;
};

}  // namespace simc

#endif  // SIMC_CORE_ERRORS_H
