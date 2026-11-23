% docs/physics/latex/simcnote.cls - document class of the SIMC++ physics notes (one short article per reaction plus the
% common note 00_common). Same look as the SIMC++ books (docs/common/latex/simcpp.cls): the palette and TikZ styles of
% simcfig.sty, the typefaces of simcfonts.sty, the \code macro, the tables and the callout boxes; an article layout
% with a title band instead of a title page.
%   Compile: cmake --build <build> --target physics_pdf, or
%     python3 -I docs/physics/build_physics_pdf.py --source-root . --build-dir DIR   (see docs/physics/README.md).
%   TEXINPUTS must contain docs/physics/latex, docs/common/latex// and docs/common/figures (the script sets it).
% Engine: XeLaTeX. Packages (TeX Live): fontspec, unicode-math (TeX Gyre Pagella Math), amsmath, mathtools, xcolor,
%   tikz, geometry, microtype, etoolbox, environ, array, longtable, booktabs, tabularx, enumitem, caption, titlesec,
%   fancyhdr, graphicx, xr-hyper, hyperref, bookmark.
%
% Conventions of the notes (checked by build_physics_pdf.py):
%   * every displayed equation is written with the environment  feq  (one line) or  feqa  (an align block), whose
%     mandatory argument is the code citation printed under the equation, e.g.
%         \begin{feq}{event.f:518}  E' = ... \label{eq:x} \end{feq}
%     plain equation/align/\[ ... \] environments are rejected by the build script;
%   * \fsrc{file:line} cites the Fortran (simc_gfortran upstream 59ccdf6), \csrc{file:function} the C++;
%   * cross references: \eqr{eq:x} -> Eq. (n), \secr{sec:x} -> section, \tabr{tab:x} -> Table n; the reaction notes
%     refer to the common note 00 with \ceqr, \csecr, \ctabr (xr-hyper, \externaldocument[C-]{00_common}).
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{simcnote}[2026/11/23 SIMC++ physics notes]
\LoadClass[11pt,a4paper]{article}

\RequirePackage{amsmath,mathtools}
\RequirePackage[table]{xcolor}
\RequirePackage{simcfonts}
\RequirePackage{unicode-math}
\setmathfont{texgyrepagella-math.otf}
\RequirePackage{simcfig}
\RequirePackage[a4paper,inner=24mm,outer=24mm,top=26mm,bottom=26mm,headheight=15pt,headsep=8mm,footskip=11mm]{geometry}
\RequirePackage{microtype}
\RequirePackage{etoolbox,environ,array,longtable,booktabs,tabularx,enumitem,caption,titlesec,fancyhdr,graphicx}
\RequirePackage[hyphens]{url}
\RequirePackage{xr-hyper}
\RequirePackage[colorlinks=true,linkcolor=simcblue,urlcolor=simcblue,citecolor=simcblue,
  bookmarksnumbered=true,bookmarksopen=true,pdfdisplaydoctitle=true]{hyperref}
\RequirePackage{bookmark}

\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5\baselineskip plus 2pt minus 1pt}
\linespread{1.04}
\setcounter{tocdepth}{2}
\setcounter{secnumdepth}{3}
\color{simcnavy!92!black}
\urlstyle{tt}
\emergencystretch=3em
\raggedbottom
\numberwithin{equation}{section}
\allowdisplaybreaks[2]

% ---- metadata ----------------------------------------------------------------------------------------------------
% \simcnotemeta{number}{short title (header)}{title}{subtitle}
\newcommand\simcnotenumber{}\newcommand\simcnoteshort{}\newcommand\simcnotetitle{}\newcommand\simcnotesubtitle{}
\newcommand\simcnotemeta[4]{\renewcommand\simcnotenumber{#1}\renewcommand\simcnoteshort{#2}%
  \renewcommand\simcnotetitle{#3}\renewcommand\simcnotesubtitle{#4}%
  \hypersetup{pdftitle={SIMC++ physics note #1: #2},pdfsubject={#4},pdfauthor={SIMC++ project},
    pdfkeywords={SIMC++, SIMC, Monte Carlo, electron scattering, physics}}}
% version and date: simcnote_version.tex is written at build time by build_physics_pdf.py (git describe, commit date)
\newcommand\simcnoteversion{unversioned}\newcommand\simcnotedate{\today}
\InputIfFileExists{simcnote_version.tex}{}{}
\newcommand\simclogo{simcpp_logo}   % docs/common/figures/simcpp_logo.png (TEXINPUTS)

% ---- inline elements ---------------------------------------------------------------------------------------------
\newcommand\simccodeinner[1]{{\ttfamily\hyphenchar\font=-1\color{simcblue!70!simcnavy}#1}}
\DeclareRobustCommand\code[1]{\ifmmode\text{\simccodeinner{#1}}\else\leavevmode\simccodeinner{#1}\fi}
% Fortran citation (file:line of simc_gfortran upstream 59ccdf6) and C++ citation (file:function)
\DeclareRobustCommand\fsrc[1]{\leavevmode{\ttfamily\small\color{simcgray!80!simcnavy}#1}}
\DeclareRobustCommand\csrc[1]{\leavevmode{\ttfamily\small\color{simcblue!80!simcnavy}#1}}
% a finding of docs/translation/markdown/FORTRAN_DIFFERENCES.md
\DeclareRobustCommand\finding[1]{\leavevmode{\sffamily\small\bfseries\color{simcred!85!black}#1}}
\newcommand\simcth[1]{{\sffamily\bfseries\color{white}\renewcommand\code[1]{{\ttfamily\bfseries##1}}#1}}
\renewcommand\footnoterule{\noindent{\color{simcrule}\rule{0.3\linewidth}{0.5pt}}\par\vspace{2pt}}

% ---- equations with their code citation -------------------------------------------------------------------------
% The citation is set flush right in small grey type directly under the display.
\newcommand\simceqsrc[1]{\par\nopagebreak\vspace{-0.55\baselineskip}\nopagebreak
  {\noindent\hfill{\sffamily\scriptsize\color{simcgray}\textsc{code}\enspace\ttfamily#1}\par}\vspace{0.15\baselineskip}}
\NewEnviron{feq}[1]{\postdisplaypenalty=10000 \begin{equation}\BODY\end{equation}\simceqsrc{#1}}
\NewEnviron{feqa}[1]{\postdisplaypenalty=10000 \begin{align}\BODY\end{align}\simceqsrc{#1}}

% ---- headings ----------------------------------------------------------------------------------------------------
\titleformat{\section}{\normalfont\sffamily\Large\bfseries\color{simcnavy}\RaggedRight}
  {\color{simcblue}\thesection}{0.7em}{}
\titleformat{\subsection}{\normalfont\sffamily\large\bfseries\color{simcnavy}\RaggedRight}
  {\color{simcblue}\thesubsection}{0.6em}{}
\titleformat{\subsubsection}{\normalfont\sffamily\normalsize\bfseries\color{simcnavy}}
  {\color{simcblue}\thesubsubsection}{0.5em}{}
\titlespacing*{\section}{0pt}{16pt plus 4pt minus 2pt}{4pt}
\titlespacing*{\subsection}{0pt}{11pt plus 3pt minus 2pt}{2pt}
\titlespacing*{\subsubsection}{0pt}{9pt plus 2pt}{1pt}
\RequirePackage{ragged2e}

% ---- title block ------------------------------------------------------------------------------------------------
\newcommand\simcnotemaketitle{%
  \thispagestyle{simcfirst}%
  \noindent
  \begin{tikzpicture}[remember picture, overlay]
    \fill[simcblue] (current page.north west) rectangle ([yshift=-12mm]current page.north east);
    \node[anchor=west, text=white, font=\sffamily\bfseries\small] at ([xshift=24mm,yshift=-6mm]current page.north west)
      {SIMC++ PHYSICS NOTE \simcnotenumber};
    \node[anchor=east] at ([xshift=-24mm,yshift=-27mm]current page.north east)
      {\includegraphics[width=34mm]{\simclogo}};
  \end{tikzpicture}%
  \vspace*{2mm}
  \parbox[t]{\linewidth-40mm}{\sffamily\bfseries\fontsize{24}{29}\selectfont\color{simcnavy}\RaggedRight\simcnotetitle}\par
  \vspace{2mm}{\color{simcblue}\rule{38mm}{1.3pt}}\par\vspace{1mm}
  {\sffamily\large\color{simcgray}\RaggedRight\simcnotesubtitle\par}
  \vspace{2mm}
  {\sffamily\footnotesize\color{simcgray}SIMC++ source \texttt{\simcnoteversion} of \simcnotedate\enspace\textbullet\enspace
   specification: Fortran SIMC (\texttt{simc\_gfortran}) upstream commit \texttt{59ccdf6}\par}
  \vspace{3mm}}

% chapter-like abstract with the blue bar of the books
\newcommand\simcabstract[1]{%
  \par\noindent
  \begin{tikzpicture}
    \node[inner sep=0pt, text width=\linewidth-7mm, anchor=north west, align=left,
          font=\normalsize\itshape, text=simcnavy!85] (a) at (7mm,0) {#1\par};
    \fill[simcblue] ([xshift=-7mm]a.north west) rectangle ([xshift=-5.6mm]a.south west);
  \end{tikzpicture}\par\vspace{4pt}}

% compact table of contents
\newcommand\simcnotetoc{%
  \begingroup\hypersetup{linkcolor=simcnavy}\small
  \setlength{\parskip}{0pt}%
  \renewcommand\contentsname{Contents}%
  \tableofcontents\endgroup\vspace{2mm}{\color{simcrule}\hrule height 0.5pt}\vspace{2mm}}
\RequirePackage{titletoc}
\titlecontents{section}[1.8em]{\addvspace{2pt}\sffamily\bfseries\color{simcnavy}}
  {\contentslabel{1.8em}}{\hspace*{-1.8em}}{\hfill\contentspage}
\titlecontents{subsection}[4.2em]{\normalfont\small}
  {\contentslabel{2.4em}}{\hspace*{-2.4em}}{\titlerule*[0.6pc]{.}\contentspage}
\titleformat{name=\section,numberless}{\normalfont\sffamily\Large\bfseries\color{simcnavy}}{}{0pt}{}

% ---- page styles -------------------------------------------------------------------------------------------------
\newcommand\simcheadfont{\sffamily\footnotesize\color{simcgray}}
\fancypagestyle{simcmain}{%
  \fancyhf{}%
  \fancyhead[L]{\simcheadfont SIMC++ physics note \simcnotenumber\enspace\textbullet\enspace\simcnoteshort}%
  \fancyhead[R]{\raisebox{-1.2pt}{\includegraphics[height=9pt]{\simclogo}}}%
  \fancyfoot[L]{\sffamily\scriptsize\color{simcgray}\simcnoteversion}%
  \fancyfoot[R]{\sffamily\small\bfseries\color{simcblue}\thepage}%
  \renewcommand\headrule{{\color{simcrule}\hrule height 0.5pt}}%
  \renewcommand\footrulewidth{0pt}}
\fancypagestyle{simcfirst}{%
  \fancyhf{}%
  \fancyfoot[L]{\sffamily\scriptsize\color{simcgray}\simcnoteversion}%
  \fancyfoot[R]{\sffamily\small\bfseries\color{simcblue}\thepage}%
  \renewcommand\headrulewidth{0pt}\renewcommand\footrulewidth{0pt}}
\pagestyle{simcmain}

% ---- lists -------------------------------------------------------------------------------------------------------
\setlist{itemsep=1pt, parsep=1pt, topsep=2pt, partopsep=0pt}
\setlist[itemize,1]{label={\color{simcblue}\small\textbullet}, leftmargin=1.4em}
\setlist[itemize,2]{label={\color{simcgray}\scriptsize\textendash}, leftmargin=1.3em}
\setlist[enumerate,1]{label={\sffamily\color{simcblue}\arabic*.}, leftmargin=1.7em}
\setlist[description]{font=\sffamily\bfseries\color{simcnavy}, leftmargin=1.2em}

% ---- captions and tables -----------------------------------------------------------------------------------------
\DeclareCaptionFont{simccap}{\sffamily\small\color{simcnavy}}
\DeclareCaptionLabelFormat{simclabel}{{\color{simcblue}\bfseries#1~#2}}
\captionsetup{font=simccap, labelformat=simclabel, labelsep=quad, justification=raggedright,
  singlelinecheck=false, skip=5pt}
\captionsetup[longtable]{position=above, skip=4pt}
\colorlet{simchead}{simcnavy}
\colorlet{simczebra}{simcblue!5}
\setlength{\LTpre}{4pt}\setlength{\LTpost}{4pt}
\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt plus 1fill}
% \begin{simctable}{size}{colspec} ... \end{simctable}: a zebra longtable with a navy header row
\newenvironment{simctable}[2]{%
  \par\begingroup\csname#1\endcsname\setlength{\parskip}{0pt}\renewcommand\arraystretch{1.2}\setlength{\tabcolsep}{4pt}%
  \arrayrulecolor{simcnavy}\rowcolors{2}{white}{simczebra}%
  \begin{longtable}{@{}#2@{}}}%
  {\end{longtable}\endgroup\par}
\newcolumntype{L}[1]{>{\RaggedRight\arraybackslash}p{#1}}
\aboverulesep=0pt \belowrulesep=0pt
\heavyrulewidth=0.9pt \lightrulewidth=0.5pt

% ---- callouts ----------------------------------------------------------------------------------------------------
\newcommand\simccallout[4]{% colour, symbol, title, body
  \par\addvspace{4pt}\noindent
  \begin{tikzpicture}
    \node[inner xsep=10pt, inner ysep=6pt, text width=\linewidth-23pt, anchor=north west, align=left,
          fill=#1!6, font=\normalsize] (b) at (3pt,0)
      {{\sffamily\bfseries\color{#1}#3}\par\smallskip #4\par};
    \fill[#1] ([xshift=-3pt]b.north west) rectangle (b.south west);
    \node[overlay, circle, fill=#1, text=white, font=\sffamily\bfseries\footnotesize, inner sep=0pt, minimum size=5.2mm]
      at ([xshift=-1.5pt,yshift=-3.6mm]b.north west) {#2};
  \end{tikzpicture}\par\addvspace{4pt}}
\NewEnviron{simcnote}{\simccallout{simcblue}{i}{Note}{\BODY}}
\NewEnviron{simcwarning}{\simccallout{simcred}{!}{Warning}{\BODY}}
\NewEnviron{simclegacy}{\simccallout{simcnavy}{F}{Fortran behaviour}{\BODY}}
\NewEnviron{simcfinding}{\simccallout{simcred}{?}{Candidate finding (not in FORTRAN\_DIFFERENCES.md)}{\BODY}}

% ---- cross references (labels sec:, eq:, tab:; prefix C- for the common note 00 via xr-hyper) ----------------
\newcommand\eqr[1]{Eq.~(\ref{#1})}
\newcommand\secr[1]{\S\ref{#1}}
\newcommand\tabr[1]{Table~\ref{#1}}
\newcommand\ceqr[1]{Eq.~(\ref{C-#1}) of note~00}
\newcommand\csecr[1]{\S\ref{C-#1} of note~00}
\newcommand\ctabr[1]{Table~\ref{C-#1} of note~00}
% \simccommon: load the labels of the common note (built first, in the same directory)
\newcommand\simccommon{\externaldocument[C-]{00_common}[SIMCpp_Physics_00_common.pdf]}

% ---- references --------------------------------------------------------------------------------------------------
\renewcommand\refname{References}
\endinput
