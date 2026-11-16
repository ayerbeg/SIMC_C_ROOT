% docs/common/latex/simcpp.cls - document class of the three SIMC++ books (user manual, developer manual, translation
% report): title page with the SIMC++ logo, part pages, chapter openings with an abstract, running headers, callout
% boxes, command/code listings, reference tables, the parameter reference, the index, the dashboard page.
% Shared by the LaTeX editions in docs/manual/latex/ and docs/translation/latex/ (see docs/README.md).
% Compile a book: cmake --build <build> --target pdf, or python3 -I tools/build_docs_pdf.py --source-root . --build-dir DIR
%   (figures first, then xelatex, makeindex -s docs/common/latex/simcpp.ist, xelatex until the references settle).
%   By hand: TEXINPUTS must contain docs/common/latex//, docs/common/figures and the directory of the figure PDFs.
% Engine: XeLaTeX. Packages: fontspec, xcolor, tikz, geometry, microtype, etoolbox, calc, xparse, environ, array,
%   longtable, booktabs, pdflscape, graphicx, enumitem, listings, caption, titlesec, titletoc, fancyhdr, makeidx,
%   placeins, ragged2e, url, hyperref, bookmark (all in TeX Live). Fonts: TeX Gyre Pagella, TeX Gyre Heros,
%   DejaVu Sans Mono (simcfonts.sty). Palette and TikZ styles: simcfig.sty.
% Figures: the books include the PDFs built from docs/{manual,translation}/figures/*.tex (preview.sty), so these must
%   be compiled first (tools/build_docs_pdf.py does it).
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{simcpp}[2026/11/16 SIMC++ documentation class]
\LoadClass[11pt,a4paper,oneside,openany]{book}

\RequirePackage[table]{xcolor}
\RequirePackage{simcfonts}
\RequirePackage{simcfig}
\RequirePackage[a4paper,inner=25mm,outer=25mm,top=27mm,bottom=27mm,headheight=15pt,headsep=9mm,footskip=12mm]{geometry}
\RequirePackage{microtype}
\RequirePackage{etoolbox,calc,xparse,environ,array,longtable,booktabs,pdflscape,graphicx,enumitem,listings,
  caption,titlesec,titletoc,fancyhdr,makeidx,placeins,ragged2e}
\RequirePackage[hyphens]{url}
\RequirePackage[colorlinks=true,linkcolor=simcblue,urlcolor=simcblue,citecolor=simcblue,
  bookmarksnumbered=true,bookmarksopen=true,bookmarksopenlevel=0,pdfdisplaydoctitle=true]{hyperref}
\RequirePackage{bookmark}

\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55\baselineskip plus 2pt minus 1pt}
\linespread{1.04}
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{3}
\color{simcnavy!92!black}
\hypersetup{linkcolor=simcblue}
\urlstyle{tt}
\emergencystretch=3em
\raggedbottom
\makeindex

% ---- metadata: \simcsetmeta in the book's main .tex; version and date from simcpp_version.tex, written at build time
% by tools/md2latex.py --generated (git describe and commit date); without it "unversioned" and today
\newcommand\simctitle{}\newcommand\simcsubtitle{}\newcommand\simcversion{}\newcommand\simcdate{}
\newcommand\simckicker{}
\newcommand\simcsetmeta[5]{\renewcommand\simctitle{#1}\renewcommand\simcsubtitle{#2}%
  \renewcommand\simcversion{#3}\renewcommand\simcdate{#4}\renewcommand\simckicker{#5}%
  \hypersetup{pdftitle={#1},pdfsubject={#2},pdfauthor={SIMC++ project},pdfkeywords={SIMC++, SIMC, Monte Carlo, electron scattering, spectrometers}}}
\newcommand\simclogo{simcpp_logo}   % docs/common/figures/simcpp_logo.png (TEXINPUTS)
\newcommand\simcbuildversion{unversioned}\newcommand\simcbuilddate{\today}
\InputIfFileExists{simcpp_version.tex}{}{}
% the dashboard page of the translation report is data, written at build time; placeholder when it is missing
\newcommand\simcnodashboard{\chapter*{The translation at a glance}\addcontentsline{toc}{chapter}{The translation at a glance}%
  \label{dashboard}This page is generated at build time from the catalogue (tools/md2latex.py --generated); it was not
  found when this PDF was compiled.\clearpage}

% ---- inline elements ------------------------------------------------------------------------------------------------
\DeclareRobustCommand\code[1]{\leavevmode{\ttfamily\hyphenchar\font=-1\color{simcblue!70!simcnavy}#1}}
\newcommand\simcth[1]{{\sffamily\bfseries\color{white}\renewcommand\code[1]{{\ttfamily\bfseries##1}}#1}}
\newcommand\simccontinued{\sffamily\scriptsize\color{simcgray}(continued)}
\newcommand\simccontinues{\sffamily\scriptsize\color{simcgray}continued on the next page}
\newcommand\simcrule{\par\noindent{\color{simcrule}\rule{\linewidth}{0.5pt}}\par}
\renewcommand\footnoterule{\noindent{\color{simcrule}\rule{0.3\linewidth}{0.5pt}}\par\vspace{2pt}}
\renewcommand\@makefnmark{\hbox{\textsuperscript{\sffamily\color{simcblue}\@thefnmark}}}

% ---- headings ------------------------------------------------------------------------------------------------------
\titleformat{\chapter}[display]
  {\normalfont\sffamily}
  {\begin{tikzpicture}[baseline]
     \node[anchor=base west, inner sep=0, font=\sffamily\small\bfseries, text=simcgray]
       {\MakeUppercase{\chaptertitlename}};
     \node[anchor=base east, inner sep=0, font=\sffamily\bfseries\fontsize{54}{54}\selectfont, text=simcblue]
       at (\linewidth,0) {\thechapter};
   \end{tikzpicture}}
  {2pt}
  {\color{simcnavy}\fontsize{25}{30}\selectfont\bfseries\RaggedRight}
  [{\vspace{5pt}\color{simcblue}\titlerule[1.4pt]}]
\titleformat{name=\chapter,numberless}[display]
  {\normalfont\sffamily}{}{0pt}
  {\color{simcnavy}\fontsize{25}{30}\selectfont\bfseries\RaggedRight}
  [{\vspace{5pt}\color{simcblue}\titlerule[1.4pt]}]
\titlespacing*{\chapter}{0pt}{-6mm}{8mm}
\titleformat{\section}{\normalfont\sffamily\Large\bfseries\color{simcnavy}\RaggedRight}
  {\color{simcblue}\thesection}{0.75em}{}
\titleformat{\subsection}{\normalfont\sffamily\large\bfseries\color{simcnavy}\RaggedRight}
  {\color{simcblue}\thesubsection}{0.7em}{}
\titleformat{\subsubsection}{\normalfont\sffamily\normalsize\bfseries\color{simcnavy}}
  {\color{simcblue}\thesubsubsection}{0.6em}{}
\titlespacing*{\section}{0pt}{16pt plus 4pt minus 2pt}{5pt}
\titlespacing*{\subsection}{0pt}{12pt plus 3pt minus 2pt}{3pt}
\titlespacing*{\subsubsection}{0pt}{10pt plus 2pt}{2pt}

% chapter abstract (the first paragraph of a chapter, or <!-- abstract: ... -->)
\newcommand\simcabstract[1]{%
  \par\vspace{-2pt}\noindent
  \begin{tikzpicture}
    \node[inner sep=0pt, text width=\linewidth-7mm, anchor=north west, align=left,
          font=\normalsize\itshape, text=simcnavy!85] (a) at (7mm,0) {#1\par};
    \fill[simcblue] ([xshift=-7mm]a.north west) rectangle ([xshift=-5.6mm]a.south west);
  \end{tikzpicture}\par\vspace{6pt}}

% ---- table of contents ---------------------------------------------------------------------------------------------
\titlecontents{part}[0pt]{\addvspace{14pt}\sffamily\bfseries\color{simcblue}}
  {\MakeUppercase{\partname}~\thecontentslabel\quad}{}{\hfill}[\addvspace{2pt}]
\titlecontents{chapter}[2.4em]{\addvspace{6pt}\sffamily\bfseries\color{simcnavy}}
  {\contentslabel{2.4em}}{\hspace*{-2.4em}}{\hfill\contentspage}
\titlecontents{section}[5.3em]{\normalfont\small}
  {\contentslabel{2.9em}}{\hspace*{-2.9em}}{\titlerule*[0.6pc]{.}\contentspage}
\titlecontents{figure}[2.8em]{\normalfont\small}{\contentslabel{2.8em}}{}{\titlerule*[0.6pc]{.}\contentspage}
\titlecontents{table}[2.8em]{\normalfont\small}{\contentslabel{2.8em}}{}{\titlerule*[0.6pc]{.}\contentspage}

% ---- page styles ---------------------------------------------------------------------------------------------------
\newcommand\simcheadfont{\sffamily\footnotesize\color{simcgray}}
\fancypagestyle{simcmain}{%
  \fancyhf{}%
  \fancyhead[L]{\simcheadfont\nouppercase{\leftmark}}%
  \fancyhead[R]{\raisebox{-1.2pt}{\includegraphics[height=9pt]{\simclogo}}\hspace{0.5em}\simcheadfont\simctitle}%
  \fancyfoot[L]{\sffamily\scriptsize\color{simcgray}\simcversion}%
  \fancyfoot[R]{\sffamily\small\bfseries\color{simcblue}\thepage}%
  \renewcommand\headrule{{\color{simcrule}\hrule height 0.5pt}}%
  \renewcommand\footrulewidth{0pt}}
\fancypagestyle{plain}{%
  \fancyhf{}\fancyfoot[R]{\sffamily\small\bfseries\color{simcblue}\thepage}%
  \fancyfoot[L]{\sffamily\scriptsize\color{simcgray}\simcversion}%
  \renewcommand\headrulewidth{0pt}\renewcommand\footrulewidth{0pt}}
\renewcommand\chaptermark[1]{\markboth{\ifnum\c@secnumdepth>\m@ne\if@mainmatter\thechapter\quad\fi\fi#1}{}}
\renewcommand\sectionmark[1]{}
\pagestyle{simcmain}

% ---- title page ----------------------------------------------------------------------------------------------------
\newcommand\simctitlepage{%
  \begin{titlepage}
  \thispagestyle{empty}
  \begin{tikzpicture}[remember picture, overlay]
    \fill[simcblue] (current page.north west) rectangle ([yshift=-14mm]current page.north east);
    \node[anchor=west, text=white, font=\sffamily\bfseries\small] at ([xshift=25mm,yshift=-7mm]current page.north west)
      {\MakeUppercase{\simckicker}};
    \node[anchor=north] at ([yshift=-36mm]current page.north)
      {\includegraphics[width=128mm]{\simclogo}};
    \node[anchor=north west, text width=160mm, align=left, text=simcnavy,
          font=\sffamily\bfseries\fontsize{34}{40}\selectfont]
      (t) at ([xshift=25mm,yshift=-162mm]current page.north west) {\simctitle};
    \fill[simcblue] ([yshift=-5mm]t.south west) rectangle ([yshift=-6.4mm,xshift=40mm]t.south west);
    \node[anchor=north west, text width=150mm, align=left, text=simcgray,
          font=\sffamily\fontsize{15}{19}\selectfont] at ([yshift=-11mm]t.south west) {\simcsubtitle};
    \fill[simcnavy] (current page.south west) rectangle ([yshift=34mm]current page.south east);
    \fill[simcgreen] ([yshift=34mm]current page.south west) rectangle ([yshift=35.2mm,xshift=70mm]current page.south west);
    \fill[simcred] ([yshift=34mm,xshift=70mm]current page.south west) rectangle ([yshift=35.2mm,xshift=90mm]current page.south west);
    \node[anchor=south west, text=white, align=left, font=\sffamily\small] at ([xshift=25mm,yshift=11mm]current page.south west)
      {\textbf{Version} \texttt{\simcversion}\\[2pt]\textbf{Date} \simcdate};
    \node[anchor=south east, text=white!80!simcblue, align=right, font=\sffamily\footnotesize]
      at ([xshift=-25mm,yshift=11mm]current page.south east)
      {LaTeX edition of the SIMC++ documentation\\[2pt] (repository directory \texttt{docs/})};
  \end{tikzpicture}
  \end{titlepage}}

% colophon on the back of the title page
\newcommand\simccolophon{%
  \thispagestyle{empty}\null\vfill
  {\small\sffamily\color{simcgray}
  \textbf{\color{simcnavy}\simctitle}\\[2pt]\simcsubtitle\par\medskip
  Version \texttt{\simcversion} of \simcdate. This is the LaTeX edition of the SIMC++ documentation; a Markdown
  edition of the same text, readable in the repository browser, is kept next to it (\texttt{docs/README.md}).
  Build: \texttt{cmake --build build --target pdf}.\par\medskip
  SIMC++ is the C++ translation of SIMC (Fortran, \texttt{simc\_gfortran}), the Monte Carlo simulation of
  electron-scattering experiments with magnetic spectrometers developed at Jefferson Lab.\par}
  \clearpage}

% ---- book structure ------------------------------------------------------------------------------------------------
\newcommand\simcfrontmatter{%
  \frontmatter\pagestyle{simcmain}%
  \simctitlepage\simccolophon
  \pdfbookmark[0]{Contents}{simc-toc}%
  \begingroup\hypersetup{linkcolor=simcnavy}%
  \tableofcontents\clearpage
  \listoffigures\clearpage
  \listoftables\clearpage\endgroup}
\newcommand\simcmainmatter{\mainmatter\pagestyle{simcmain}}
\newcommand\simcpart[2]{%
  \clearpage\thispagestyle{empty}%
  \refstepcounter{part}\phantomsection
  \addcontentsline{toc}{part}{\protect\numberline{\thepart}#1}%
  \markboth{}{}%
  \begin{tikzpicture}[remember picture, overlay]
    \fill[simcnavy] (current page.north west) rectangle ([yshift=-118mm]current page.north east);
    \fill[simcblue] ([yshift=-118mm]current page.north west) rectangle ([yshift=-121mm]current page.north east);
    \node[anchor=south west, text=simcblue!45!white, font=\sffamily\bfseries\large]
      at ([xshift=25mm,yshift=-62mm]current page.north west) {\MakeUppercase{\partname}~\Roman{part}};
    \node[anchor=north west, text width=160mm, align=left, text=white, font=\sffamily\bfseries\fontsize{32}{38}\selectfont]
      at ([xshift=25mm,yshift=-66mm]current page.north west) {#1};
    \node[anchor=north west, text width=130mm, align=left, text=simcnavy, font=\large\itshape]
      at ([xshift=25mm,yshift=-136mm]current page.north west) {#2};
    \node[anchor=south east, opacity=0.9] at ([xshift=-25mm,yshift=20mm]current page.south east)
      {\includegraphics[width=38mm]{\simclogo}};
  \end{tikzpicture}%
  \clearpage}
\newcommand\simcappendix{%
  \appendix
  \clearpage\thispagestyle{empty}\phantomsection
  \addcontentsline{toc}{part}{Appendices}%
  \begin{tikzpicture}[remember picture, overlay]
    \fill[simcnavy] (current page.north west) rectangle ([yshift=-118mm]current page.north east);
    \fill[simcblue] ([yshift=-118mm]current page.north west) rectangle ([yshift=-121mm]current page.north east);
    \node[anchor=north west, text=white, font=\sffamily\bfseries\fontsize{32}{38}\selectfont]
      at ([xshift=25mm,yshift=-66mm]current page.north west) {Appendices};
  \end{tikzpicture}\clearpage}
\newcommand\simcbackmatter{%
  \backmatter\clearpage\phantomsection
  \addcontentsline{toc}{chapter}{Index}%
  \printindex}

% ---- index ---------------------------------------------------------------------------------------------------------
\renewenvironment{theindex}{%
  \chapter*{Index}\markboth{Index}{}%
  {\small\color{simcgray}Input parameters, command-line options and classes are grouped under their own headings.
   The page numbers point to the definition or to a mention; in the user manual the definition of a parameter is
   its entry in the parameter reference.\par}\medskip
  \small\setlength{\parskip}{0pt}\setlength{\columnsep}{8mm}%
  \begin{multicols}{2}\raggedright
  \def\item{\par\addvspace{7pt}\hangindent 40pt\sffamily\bfseries\color{simcblue}}%
  \def\subitem{\par\hangindent 40pt\hspace*{12pt}\normalfont\color{simcnavy}}%
  \def\subsubitem{\par\hangindent 40pt\hspace*{24pt}\normalfont\color{simcnavy}}}{\end{multicols}}
\RequirePackage{multicol}
\newcommand\simcindexgroup[1]{\par\addvspace{8pt}{\sffamily\bfseries\color{simcblue}#1}\par\nopagebreak\smallskip}

% ---- lists ---------------------------------------------------------------------------------------------------------
\setlist{itemsep=2pt, parsep=1pt, topsep=3pt, partopsep=0pt}
\setlist[itemize,1]{label={\color{simcblue}\small\textbullet}, leftmargin=1.4em}
\setlist[itemize,2]{label={\color{simcgray}\scriptsize\textendash}, leftmargin=1.3em}
\setlist[enumerate,1]{label={\sffamily\color{simcblue}\arabic*.}, leftmargin=1.7em}
\setlist[enumerate,2]{label={\sffamily\color{simcblue}\alph*)}, leftmargin=1.6em}

% ---- captions and figures ------------------------------------------------------------------------------------------
\DeclareCaptionFont{simccap}{\sffamily\small\color{simcnavy}}
\DeclareCaptionLabelFormat{simclabel}{{\color{simcblue}\bfseries#1~#2}}
\captionsetup{font=simccap, labelformat=simclabel, labelsep=quad, justification=raggedright,
  singlelinecheck=false, skip=6pt}
\captionsetup[longtable]{position=above, skip=4pt}
\newsavebox\simcfigbox
\newcommand\simcfigure[1]{%
  \sbox\simcfigbox{\includegraphics{#1}}%
  \centering
  \ifdim\wd\simcfigbox>\linewidth
    \includegraphics[width=\linewidth,height=0.72\textheight,keepaspectratio]{#1}%
  \else\ifdim\ht\simcfigbox>0.72\textheight
    \includegraphics[height=0.72\textheight,keepaspectratio]{#1}%
  \else\usebox\simcfigbox\fi\fi\par}
\setcounter{topnumber}{3}\renewcommand\topfraction{0.9}\renewcommand\textfraction{0.08}
\renewcommand\floatpagefraction{0.75}

% ---- tables (md2latex: \begin{simctable}{size}{ncol}{colspec}) -----------------------------------------------------
\colorlet{simchead}{simcnavy}
\colorlet{simczebra}{simcblue!5}
\setlength{\LTpre}{6pt}\setlength{\LTpost}{6pt}
\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt plus 1fill}
\newenvironment{simctable}[3]{%
  \par\begingroup\csname#1\endcsname\setlength{\parskip}{0pt}\renewcommand\arraystretch{1.22}\setlength{\tabcolsep}{4pt}%
  \sffamily\simctablefont
  \arrayrulecolor{simcnavy}\rowcolors{2}{white}{simczebra}%
  \begin{longtable}{@{}#3@{}}}%
  {\end{longtable}\endgroup\par}
\newcommand\simctablefont{\rmfamily\simcliningserif}
\aboverulesep=0pt \belowrulesep=0pt
\heavyrulewidth=0.9pt \lightrulewidth=0.5pt

% ---- callouts (> **Note:** / Tip / Warning / Legacy SIMC) ---------------------------------------------------------
\newcommand\simccallout[4]{% colour, symbol, title, body
  \par\addvspace{6pt}\noindent
  \begin{tikzpicture}
    \node[inner xsep=10pt, inner ysep=7pt, text width=\linewidth-23pt, anchor=north west, align=left,
          fill=#1!6, font=\normalsize] (b) at (3pt,0)
      {{\sffamily\bfseries\color{#1}#3}\par\smallskip #4\par};
    \fill[#1] ([xshift=-3pt]b.north west) rectangle (b.south west);
    \node[overlay, circle, fill=#1, text=white, font=\sffamily\bfseries\footnotesize, inner sep=0pt, minimum size=5.2mm]
      at ([xshift=-1.5pt,yshift=-3.6mm]b.north west) {#2};
  \end{tikzpicture}\par\addvspace{6pt}}
\NewEnviron{simcnote}{\simccallout{simcblue}{i}{Note}{\BODY}}
\NewEnviron{simctip}{\simccallout{simcgreen}{\checkmark}{Tip}{\BODY}}
\NewEnviron{simcwarning}{\simccallout{simcred}{!}{Warning}{\BODY}}
\NewEnviron{simclegacy}{\simccallout{simcnavy}{F}{Legacy SIMC}{\BODY}}
\providecommand\checkmark{\ensuremath{\surd}}
\renewenvironment{quote}{\list{}{\leftmargin=1.5em\rightmargin=1.5em}\item\relax\itshape}{\endlist}

% ---- listings ------------------------------------------------------------------------------------------------------
\definecolor{simcshellbg}{HTML}{F1F4F8}
\definecolor{simccodebg}{HTML}{FAFBFC}
\lstdefinelanguage{simcjson}{morestring=[b]", morecomment=[l]{//}, literate=
  {:}{{{\color{simcnavy}:}}}1 {,}{{{\color{simcnavy},}}}1 {\{}{{{\color{simcblue}\{}}}1 {\}}{{{\color{simcblue}\}}}}1
  {[}{{{\color{simcblue}[}}}1 {]}{{{\color{simcblue}]}}}1}
\lstdefinelanguage{simccmake}{sensitive=false, morecomment=[l]\#, morestring=[b]",
  morekeywords={cmake_minimum_required,project,find_package,add_executable,add_library,target_link_libraries,
    target_compile_options,target_include_directories,set,add_subdirectory,add_test,install,if,else,endif,
    option,set_tests_properties,REQUIRED,PRIVATE,PUBLIC,INTERFACE,OBJECT}}
\lstdefinestyle{simcbase}{basicstyle=\ttfamily\footnotesize, columns=fullflexible, keepspaces=true,
  breaklines=true, breakatwhitespace=false, showstringspaces=false, upquote=true, tabsize=4,
  postbreak=\mbox{\textcolor{simcblue}{$\hookrightarrow$}\space}, aboveskip=6pt, belowskip=8pt,
  xleftmargin=6pt, xrightmargin=4pt, framexleftmargin=6pt, framexrightmargin=4pt, framextopmargin=3pt,
  framexbottommargin=3pt, extendedchars=true, texcl=false, mathescape=false, escapechar={}}
\lstdefinestyle{simccpp}{language=C++, keywordstyle=\color{simcblue}\bfseries,
  commentstyle=\color{simcgray}\itshape, stringstyle=\color{simcgreen!80!black},
  morekeywords={override,final,nullptr,constexpr,noexcept,std,unique_ptr,make_unique,size_t,uint64_t}}
\lstdefinestyle{simcjson}{language=simcjson, stringstyle=\color{simcblue!80!black}}
\lstdefinestyle{simccmake}{language=simccmake, keywordstyle=\color{simcblue}\bfseries,
  commentstyle=\color{simcgray}\itshape, stringstyle=\color{simcgreen!80!black}}
\lstdefinestyle{simcbash}{language=bash, commentstyle=\color{simcgray}\itshape, keywordstyle={},
  morekeywords={}, deletekeywords={cd,export,for,do,done,if,then,fi,echo,test,in,set,source,exit,read,cat,time,
  type,ls,wait,sed,diff,head,tail,local}}
% command lines (```bash): a shaded box with a blue bar
\lstnewenvironment{simcshell}[1][]{\lstset{style=simcbase, style=simcbash, backgroundcolor=\color{simcshellbg},
  frame=leftline, framerule=2.4pt, rulecolor=\color{simcblue}, #1}}{}
% source code (```cpp, ```json, ```cmake): white box with thin rules, coloured syntax
\lstnewenvironment{simccode}[1][]{\lstset{style=simcbase, backgroundcolor=\color{simccodebg},
  frame=tb, framerule=0.6pt, rulecolor=\color{simcnavy!45}, #1}}{}
% text, program output, diagrams (``` without language)
\lstnewenvironment{simcplain}[1][]{\lstset{style=simcbase, backgroundcolor=\color{simcgray!6},
  frame=single, framerule=0.4pt, rulecolor=\color{simcgray!35}, #1}}{}
\newcommand\simcsource[1]{\par\addvspace{4pt}\noindent{\sffamily\scriptsize\color{simcgray}%
  \textbf{source file}\enspace\code{#1}}\par\nopagebreak\vspace{-4pt}}
\newcommand\simccodesize[1]{\ifstrequal{#1}{normal}{\footnotesize}{\ifstrequal{#1}{small}{\scriptsize}{\tiny}}}

% ---- parameter reference --------------------------------------------------------------------------------------------
% one entry per parameter: name (and the metadata on the right), then the description; md2latex writes the entries
\newenvironment{simcparams}{\par\addvspace{2pt}\begingroup\setlength{\parskip}{0pt}}{\par\endgroup\addvspace{4pt}}
% \simcparam{name}{label}{units}{default}{used by}{meaning}{index key}
\newcommand\simcparam[7]{%
  \par\addvspace{5pt}\needspacebox{3\baselineskip}\noindent
  \phantomsection\label{#2}\simcidx{Parameters!#7}%
  {\ttfamily\bfseries\small\color{simcnavy}#1}\hfill
  {\sffamily\scriptsize\color{simcgray}#3\ifstrempty{#4}{}{\enspace\textbullet\enspace default \textbf{#4}}%
   \ifstrempty{#5}{}{\enspace\textbullet\enspace used by \textcolor{simcblue}{#5}}}\par\nopagebreak
  \vspace{1pt}{\leftskip=1.2em\small #6\par}\vspace{3pt}%
  {\color{simcrule}\rule{\linewidth}{0.3pt}}\par}
\newcommand\needspacebox[1]{\penalty-100\vskip 0pt plus #1\penalty-200\vskip 0pt plus -#1}
% every index entry goes through one macro, so that the same entry is written identically from any context
\newcommand\simcidx[1]{\index{#1}}
\newcommand\simcindextop[1]{\par\addvspace{6pt}{\sffamily\bfseries\color{simcblue}#1}}

% ---- translation report: the dashboard page -----------------------------------------------------------------------
\newcommand\simcdashtile[3]{% value, label, colour
  \begin{tikzpicture}
    \node[fill=#3!8, draw=#3!70, line width=0.6pt, rounded corners=3pt, minimum width=34mm, minimum height=23mm,
          inner sep=4pt] (t) {};
    \node[font=\sffamily\bfseries\fontsize{26}{28}\selectfont, text=#3] at ([yshift=3mm]t.center) {#1};
    \node[font=\sffamily\scriptsize, text=simcnavy, text width=31mm, align=center] at ([yshift=-6mm]t.center) {#2};
  \end{tikzpicture}}
\newcommand\simcdashtiles[5]{% findings, routines, cited, uncited, levels
  \noindent\simcdashtile{#2}{Fortran routines linked into \texttt{simc}}{simcnavy}\hfill
  \simcdashtile{#3}{routines cited by the C++ code (#4 not ported, with a reason)}{simcblue}\hfill
  \simcdashtile{#1}{findings catalogued (bugs, quirks, deviations)}{simcred}\hfill
  \simcdashtile{#5}{validation levels against the Fortran}{simcgreen}\par}
\newenvironment{simcdashbars}{\par\medskip\noindent\begin{tikzpicture}[x=1mm,y=1mm]}{\end{tikzpicture}\par}
% \simcdashbar{index}{Category}{meaning}{count}{fraction of the largest}{colour}
\newcommand\simcdashbar[6]{%
  \node[anchor=east, font=\sffamily\bfseries\small, text=simcnavy] at (12,-#1*7) {#2};
  \node[anchor=west, font=\sffamily\scriptsize, text=simcgray, text width=58mm] at (14,-#1*7) {#3};
  \fill[simcblue!15] (75,-#1*7-2.3) rectangle (155,-#1*7+2.3);
  \fill[#6] (75,-#1*7-2.3) rectangle ({75+80*#5},-#1*7+2.3);
  \node[anchor=west, font=\sffamily\bfseries\small, text=simcnavy] at ({75+80*#5+1},-#1*7) {#4};}
\newenvironment{simcdashlevels}{\par\medskip\noindent\begingroup\small\setlength{\parskip}{0pt}\renewcommand\arraystretch{1.3}%
  \rowcolors{2}{white}{simczebra}%
  \begin{longtable}{@{}>{\sffamily\bfseries\color{simcblue}}p{0.24\linewidth}>{\RaggedRight\arraybackslash}p{0.72\linewidth}@{}}
  \toprule\rowcolor{simchead}\simcth{level} & \simcth{the question it answers}\\\midrule\endhead
  \bottomrule\endlastfoot}{\end{longtable}\endgroup}
\newcommand\simcdashlevel[2]{#1 & #2\\}
\endinput
