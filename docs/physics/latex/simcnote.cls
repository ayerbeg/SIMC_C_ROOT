% docs/physics/latex/simcnote.cls - document class of the SIMC++ physics notes and of the bound book
% "SIMC++ Physics Reference" built from the same sources.
%
%   \documentclass{simcnote}         a standalone note (article layout with a title band; one PDF per reaction)
%   \documentclass[book]{simcnote}   the book (docs/physics/latex/physics_reference.tex): the class of the SIMC++
%                                    manuals (docs/common/latex/simcpp.cls) with the note macros; every note becomes a
%                                    chapter through \simcinputnote{<file>} (note 00 forms Part I, the reactions Part II)
%
% Build: cmake --build <build> --target physics_pdf (publishes), or
%   python3 -I docs/physics/build_physics_pdf.py --source-root . --build-dir DIR   (see docs/physics/README.md).
% TEXINPUTS must contain docs/physics/latex, docs/common/latex// and docs/common/figures (the script sets it).
% Engine: XeLaTeX. Packages (TeX Live): fontspec, unicode-math (TeX Gyre Pagella Math), amsmath, mathtools, xcolor,
%   tikz, geometry, microtype, etoolbox, environ, array, longtable, booktabs, enumitem, caption, titlesec, titletoc,
%   fancyhdr, graphicx, ragged2e, xr-hyper (standalone notes), hyperref, bookmark; the book mode needs those of
%   simcpp.cls in addition.
%
% Writing rules (checked by build_physics_pdf.py; see docs/physics/README.md and template_reaction.tex):
%   * the note file is  \documentclass{simcnote} [\simccommon] \begin{document} \simcnotemeta{NN}{short}{title}{subtitle}
%     \simcnotemaketitle \simcabstract{...} \simcnotetoc  <sections>  \end{document};  everything a note needs is
%     after \begin{document} (the book inputs the file without its preamble);
%   * displayed equations: \begin{feq}{<C++ file: function>}{<Fortran file:lines>} ... \end{feq} (one line) or feqa
%     (an align block). Both citations are mandatory; the C++ one is printed first. Write "none" as the Fortran
%     citation of a formula that exists only in SIMC++ (and say why in the text);
%   * labels with \lbl{eq:x} (never \label): the note number is prepended, so the labels of all notes can live in
%     one book; references with \eqr{eq:x}, \secr{sec:x}, \tabr{tab:x}; to the common note 00 with \ceqr, \csecr,
%     \ctabr (in a standalone note \simccommon in the preamble loads the labels of 00_common.aux);
%   * literature: \begin{simcrefs} \simcref{key} text ... \end{simcrefs}, cited with \citer{key};
%   * the reaction prime: \pr (H(e,e\pr p)) in titles and text, never $'$ in a title;
%   * internal finding ids (\finding{...}) only inside the simcfindings table of the last section.
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{simcnote}[2026/11/24 SIMC++ physics notes and Physics Reference]
\newif\ifsimcbook
\DeclareOption{book}{\simcbooktrue}
\ProcessOptions\relax

\ifsimcbook
  % ---------------------------------------------------------------------------------------------------------------
  % book mode: the class of the manuals (title page, part pages, chapter openings, running headers, tables, boxes)
  \LoadClass{simcpp}
\else
  % ---------------------------------------------------------------------------------------------------------------
  % standalone note: article layout with the look of the books
  \LoadClass[11pt,a4paper]{article}
  \RequirePackage[table]{xcolor}
  \RequirePackage{simcfonts}
  \RequirePackage{simcfig}
  \RequirePackage[a4paper,inner=24mm,outer=24mm,top=26mm,bottom=26mm,headheight=15pt,headsep=8mm,footskip=11mm]{geometry}
  \RequirePackage{microtype}
  \RequirePackage{etoolbox,environ,array,longtable,booktabs,enumitem,caption,titlesec,titletoc,fancyhdr,graphicx,
    ragged2e}
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
\fi

% ---- mathematics (both modes) ---------------------------------------------------------------------------------------
\RequirePackage{amsmath,mathtools}
\RequirePackage{unicode-math}
\setmathfont{texgyrepagella-math.otf}
\numberwithin{equation}{section}
\allowdisplaybreaks[2]
\RequirePackage{environ,array}

% ---- version and date: simcnote_version.tex is written at build time by build_physics_pdf.py --------------------
\newcommand\simcnoteversion{unversioned}\newcommand\simcnotedate{\today}
\InputIfFileExists{simcnote_version.tex}{}{}
\providecommand\simclogo{simcpp_logo}   % docs/common/figures/simcpp_logo.png (TEXINPUTS)

% ---- the note being typeset: number, labels prefix ---------------------------------------------------------------
\newcommand\simcnotenumber{}\newcommand\simcnoteshort{}\newcommand\simcnotetitle{}\newcommand\simcnotesubtitle{}
\newcommand\simcnoteprefix{N--}
% \lbl{key}: a label of the current note (stored as N<number>-key)
\newcommand\lbl[1]{\label{\simcnoteprefix#1}}
\newcommand\eqr[1]{Eq.~(\ref{\simcnoteprefix#1})}
\newcommand\eqsr[2]{Eqs.~(\ref{\simcnoteprefix#1})--(\ref{\simcnoteprefix#2})}
\newcommand\secr[1]{\S\ref{\simcnoteprefix#1}}
\newcommand\tabr[1]{Table~\ref{\simcnoteprefix#1}}
% references to the common note 00 (labels N00-...)
\ifsimcbook
  \newcommand\ceqr[1]{Eq.~(\ref{N00-#1})}
  \newcommand\ceqsr[2]{Eqs.~(\ref{N00-#1})--(\ref{N00-#2})}
  \newcommand\csecr[1]{\S\ref{N00-#1}}
  \newcommand\ctabr[1]{Table~\ref{N00-#1}}
  \newcommand\simccommon{}
\else
  \newcommand\ceqr[1]{Eq.~(\ref{N00-#1}) of note~00}
  \newcommand\ceqsr[2]{Eqs.~(\ref{N00-#1})--(\ref{N00-#2}) of note~00}
  \newcommand\csecr[1]{\S\ref{N00-#1} of note~00}
  \newcommand\ctabr[1]{Table~\ref{N00-#1} of note~00}
  % in the preamble of a reaction note: the labels of the common note (built first, in the same directory)
  \newcommand\simccommon{\externaldocument{00_common}[SIMCpp_Physics_00_common.pdf]}
\fi

% ---- inline elements ------------------------------------------------------------------------------------------------
% code text may break after . / _ , : (long file and function names in narrow table cells)
\ExplSyntaxOn
\tl_new:N \l__simcnote_code_tl
\cs_new_protected:Npn \simcbreakable #1
  {
    \tl_set:Nn \l__simcnote_code_tl {#1}
    \regex_replace_all:nnN { (\.|/|\c{_}|,|:) } { \0 \c{allowbreak} } \l__simcnote_code_tl
    \tl_use:N \l__simcnote_code_tl
  }
\ExplSyntaxOff
% a file or a Fortran location in plain type, breakable after / . _ : ,
\DeclareRobustCommand\file[1]{\leavevmode{\simcbreakable{#1}}}
\newcommand\simccodeinner[1]{{\ttfamily\hyphenchar\font=-1\color{simcblue!70!simcnavy}\simcbreakable{#1}}}
\DeclareRobustCommand\code[1]{\ifmmode\text{\simccodeinner{#1}}\else\leavevmode\simccodeinner{#1}\fi}
% Fortran citation (file:lines of simc_gfortran upstream 59ccdf6) and C++ citation (file: function)
\DeclareRobustCommand\fsrc[1]{\leavevmode{\ttfamily\small\color{simcgray!80!simcnavy}\simcbreakable{#1}}}
\DeclareRobustCommand\csrc[1]{\leavevmode{\ttfamily\small\color{simcblue!80!simcnavy}\simcbreakable{#1}}}
% an id of docs/translation/markdown/FORTRAN_DIFFERENCES.md (only in the findings table, simcfindings)
\DeclareRobustCommand\finding[1]{\leavevmode{\sffamily\small\bfseries\color{simcred!85!black}#1}}
% the prime of a reaction, e.g. H(e,e\pr p): a text apostrophe of the current font (no detached math prime in titles)
\DeclareRobustCommand\pr{\ifmmode{}^{\prime}\else\textquoteright\fi}
\providecommand\simcth[1]{{\sffamily\bfseries\color{white}\renewcommand\code[1]{{\ttfamily\bfseries##1}}#1}}

% ---- equations with their two code citations --------------------------------------------------------------------------
% \begin{feq}{C++ file: function}{Fortran file:lines} ... \end{feq}; the citations are printed flush right under the
% display, the C++ one first; no page break between the display and its citations.
\newcommand\simceqsrc[2]{\par\nopagebreak\vspace{-0.55\baselineskip}\nopagebreak
  {\noindent\hfill{\sffamily\scriptsize\color{simcgray}%
   \textsc{simc++}\enspace{\ttfamily\color{simcblue!70!simcnavy}#1}\enspace\textbullet\enspace
   \textsc{fortran}\enspace{\ttfamily#2}}\par}\vspace{0.15\baselineskip}}
\NewEnviron{feq}[2]{\postdisplaypenalty=10000 \begin{equation}\BODY\end{equation}\simceqsrc{#1}{#2}}
\NewEnviron{feqa}[2]{\postdisplaypenalty=10000 \begin{align}\BODY\end{align}\simceqsrc{#1}{#2}}

% ---- tables: \begin{simctable}{size}{colspec} (one item per row) ------------------------------------------------------
\colorlet{simchead}{simcnavy}
\colorlet{simczebra}{simcblue!5}
\setlength{\LTpre}{4pt}\setlength{\LTpost}{4pt}
\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt plus 1fill}
\ifdefined\simctable \let\simctable\relax \let\endsimctable\relax \fi
\newenvironment{simctable}[2]{%
  \par\begingroup\csname#1\endcsname\setlength{\parskip}{0pt}\renewcommand\arraystretch{1.2}\setlength{\tabcolsep}{4pt}%
  \arrayrulecolor{simcnavy}\rowcolors{2}{white}{simczebra}%
  \begin{longtable}{@{}#2@{}}}%
  {\end{longtable}\endgroup\par}
% the findings table of section 9 (the only place for finding ids); same layout as simctable
\newenvironment{simcfindings}[2]{\begin{simctable}{#1}{#2}}{\end{simctable}}
% L{w}: a ragged-right column of total width w INCLUDING the column separation, so that widths summing to
% \linewidth fill the line exactly
\newcolumntype{L}[1]{>{\RaggedRight\arraybackslash}p{\dimexpr#1-2\tabcolsep\relax}}

% ---- references: \begin{simcrefs} \simcref{key} text ... \end{simcrefs}; \citer{key} -------------------------------------
\newenvironment{simcrefs}{\section*{References}%
  \ifsimcbook\else\addcontentsline{toc}{section}{References}\fi
  \begin{enumerate}[label={[\arabic*]}, leftmargin=2.2em, itemsep=2pt]\small}{\end{enumerate}}
\newcommand\simcref[1]{\item\lbl{ref:#1}}
\newcommand\citer[1]{[\ref{\simcnoteprefix ref:#1}]}

\ifsimcbook
  % =============================================================================================================
  % book mode
  % \simcnotemeta{number}{short}{title}{subtitle}: starts the chapter of the note
  \newcommand\simcnotemeta[4]{\renewcommand\simcnotenumber{#1}\renewcommand\simcnoteshort{#2}%
    \renewcommand\simcnotetitle{#3}\renewcommand\simcnotesubtitle{#4}\renewcommand\simcnoteprefix{N#1-}%
    \chapter[#3]{#3}\label{N#1-note}%
    {\sffamily\color{simcgray}\RaggedRight #4\par}\medskip
    {\sffamily\footnotesize\color{simcgray}Physics note #1. Standalone PDF: \texttt{docs/physics/pdf/}; source:
     \texttt{docs/physics/latex/}.\par}\medskip}
  \newcommand\simcnotemaketitle{}
  \newcommand\simcnotetoc{}
  \NewEnviron{simcfindingbox}{\simccallout{simcred}{?}{Candidate finding (not yet in the list of Fortran differences)}{\BODY}}
  % \simcinputnote{file}: input a complete note file without its preamble (\documentclass ... \begin{document})
  \newcommand\simcinputnote[1]{\begingroup
    \let\documentclass\@gobble
    \renewenvironment{document}{}{}%
    \input{#1}\endgroup\clearpage}
  % front matter of the book (no list of figures)
  \newcommand\simcphysfrontmatter{%
    \frontmatter\pagestyle{simcmain}%
    \simctitlepage
    \thispagestyle{empty}\null\vfill
    {\small\sffamily\color{simcgray}
     \textbf{\color{simcnavy}\simctitle}\\[2pt]\simcsubtitle\par\medskip
     Version \texttt{\simcversion} of \simcdate. This book binds the SIMC++ physics notes; each chapter is also
     available as a standalone PDF (\texttt{docs/physics/pdf/}), with copies next to the plug-in code
     (\texttt{src/physics/<reaction>/}). Single source: \texttt{docs/physics/latex/}. Build:
     \texttt{cmake --build build --target physics\_pdf}.\par\medskip
     SIMC++ is the C++ translation of SIMC (Fortran, \texttt{simc\_gfortran}), the Monte Carlo simulation of
     electron-scattering experiments with magnetic spectrometers developed at Jefferson Lab.\par}
    \clearpage
    \pdfbookmark[0]{Contents}{simc-toc}%
    \begingroup\hypersetup{linkcolor=simcnavy}\tableofcontents\clearpage\listoftables\clearpage\endgroup}
\else
  % =============================================================================================================
  % standalone mode
  \newcommand\simcnotemeta[4]{\renewcommand\simcnotenumber{#1}\renewcommand\simcnoteshort{#2}%
    \renewcommand\simcnotetitle{#3}\renewcommand\simcnotesubtitle{#4}\renewcommand\simcnoteprefix{N#1-}%
    \hypersetup{pdftitle={SIMC++ physics note #1: #2},pdfsubject={#4},pdfauthor={SIMC++ project},
      pdfkeywords={SIMC++, SIMC, Monte Carlo, electron scattering, physics}}}
  \NewEnviron{simcfindingbox}{\simccallout{simcred}{?}{Candidate finding (not yet in the list of Fortran differences)}{\BODY}}

  % ---- headings
  \titleformat{\section}{\normalfont\sffamily\Large\bfseries\color{simcnavy}\RaggedRight}
    {\color{simcblue}\thesection}{0.7em}{}
  \titleformat{\subsection}{\normalfont\sffamily\large\bfseries\color{simcnavy}\RaggedRight}
    {\color{simcblue}\thesubsection}{0.6em}{}
  \titleformat{\subsubsection}{\normalfont\sffamily\normalsize\bfseries\color{simcnavy}}
    {\color{simcblue}\thesubsubsection}{0.5em}{}
  \titleformat{name=\section,numberless}{\normalfont\sffamily\Large\bfseries\color{simcnavy}}{}{0pt}{}
  \titlespacing*{\section}{0pt}{16pt plus 4pt minus 2pt}{4pt}
  \titlespacing*{\subsection}{0pt}{11pt plus 3pt minus 2pt}{2pt}
  \titlespacing*{\subsubsection}{0pt}{9pt plus 2pt}{1pt}

  % ---- title band
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
     \\part of the \emph{SIMC++ Physics Reference} (\code{docs/physics/pdf/})\par}
    \vspace{3mm}}

  % ---- abstract with the blue bar of the books
  \newcommand\simcabstract[1]{%
    \par\noindent
    \begin{tikzpicture}
      \node[inner sep=0pt, text width=\linewidth-7mm, anchor=north west, align=left,
            font=\normalsize\itshape, text=simcnavy!85] (a) at (7mm,0) {#1\par};
      \fill[simcblue] ([xshift=-7mm]a.north west) rectangle ([xshift=-5.6mm]a.south west);
    \end{tikzpicture}\par\vspace{4pt}}

  % ---- compact table of contents
  \newcommand\simcnotetoc{%
    \begingroup\hypersetup{linkcolor=simcnavy}\small
    \setlength{\parskip}{0pt}%
    \renewcommand\contentsname{Contents}%
    \tableofcontents\endgroup\vspace{2mm}{\color{simcrule}\hrule height 0.5pt}\vspace{2mm}}
  \titlecontents{section}[1.8em]{\addvspace{2pt}\sffamily\bfseries\color{simcnavy}}
    {\contentslabel{1.8em}}{\hspace*{-1.8em}}{\hfill\contentspage}
  \titlecontents{subsection}[4.2em]{\normalfont\small}
    {\contentslabel{2.4em}}{\hspace*{-2.4em}}{\titlerule*[0.6pc]{.}\contentspage}

  % ---- page styles
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

  % ---- lists
  \setlist{itemsep=1pt, parsep=1pt, topsep=2pt, partopsep=0pt}
  \setlist[itemize,1]{label={\color{simcblue}\small\textbullet}, leftmargin=1.4em}
  \setlist[itemize,2]{label={\color{simcgray}\scriptsize\textendash}, leftmargin=1.3em}
  \setlist[enumerate,1]{label={\sffamily\color{simcblue}\arabic*.}, leftmargin=1.7em}
  \setlist[description]{font=\sffamily\bfseries\color{simcnavy}, leftmargin=1.2em}

  % ---- captions
  \DeclareCaptionFont{simccap}{\sffamily\small\color{simcnavy}}
  \DeclareCaptionLabelFormat{simclabel}{{\color{simcblue}\bfseries#1~#2}}
  \captionsetup{font=simccap, labelformat=simclabel, labelsep=quad, justification=raggedright,
    singlelinecheck=false, skip=5pt}
  \captionsetup[longtable]{position=above, skip=4pt}
  \aboverulesep=0pt \belowrulesep=0pt
  \heavyrulewidth=0.9pt \lightrulewidth=0.5pt
  \renewcommand\footnoterule{\noindent{\color{simcrule}\rule{0.3\linewidth}{0.5pt}}\par\vspace{2pt}}

  % ---- callouts (the boxes of the books)
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
  \NewEnviron{simclegacy}{\simccallout{simcnavy}{F}{Legacy SIMC}{\BODY}}
  \RequirePackage{bookmark}
\fi
\endinput
