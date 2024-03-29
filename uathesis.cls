\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{uathesis}[2004/11/04 v1.4 UA Thesis]
%
% Updated on website on 4/7/14 by blinhart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% University of Akron Thesis Class (D. Dane Quinn)
%    thanks also to: Hasnain Rangwalla
%                    Curtis Clemons
%
% Comments and suggestions can be directed to quinn@uakron.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Options for this class:
%    ms       M.S. format
%    phd      Ph.D. format
%
%    cpyr     produces a Copyright page
%
%    lot      produces a ``List of Tables''
%    lof      produces a ``List of Figures''
%
%%%%%%%%%%%%%%%
%
% To create the UA signature page, add the following:
%    \advisor{Name}
%    \chair{Name}
%    \collegedean{Name}
%    \gradschdean{Name}
%
%  Optional for the M.S. thesis:
%    \coadvisor{Name}
%    \facreader{Name}
%      add faculty readers as necessary
%
%  Optional for the Ph.D. dissertation:
%    \coadvisor{Name}
%    \committee{Name}
%      add committee members as necessary
%
%%%%%%%%%%%%%%%
%
% Appendicies:
%   Appendicies can be added with the command:
%     \appendix{n}
%
%   where ``n'' is the total number of appendicies in the thesis.
%   Each appendix is included as a chapter (i.e., \chapter{Appendix Name}).
%
%   Following, UA thesis guidelines, the text is formatted differently for
%   n=1 than for n>1.
%
%%%%%%%%%%%%%%%
%
% All other formating *should* be generated automatically
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newif\if@cpyr
\newif\if@lot
\newif\if@lof
\newif\if@ms
\newif\if@phd

\newcommand{\@fntsz}{}
%\DeclareOption{8pt}{\renewcommand{\@fntsz}{8pt}}
%\DeclareOption{9pt}{\renewcommand{\@fntsz}{9pt}}
%\DeclareOption{10pt}{\renewcommand{\@fntsz}{10pt}}
%\DeclareOption{11pt}{\renewcommand{\@fntsz}{11pt}}
\DeclareOption{12pt}{\renewcommand{\@fntsz}{12pt}}

\DeclareOption{cpyr}{\@cpyrtrue}
\DeclareOption{lot}{\@lottrue}
\DeclareOption{lof}{\@loftrue}

\@mstrue
\@phdfalse%
\newcommand{\@uatype}{Thesis}
\newcommand{\@uadegree}{Master of Science}

\DeclareOption{ms}{%
  \@mstrue%
  \@phdfalse%
  \renewcommand{\@uatype}{Thesis}
  \renewcommand{\@uadegree}{Master of Science}
}
\DeclareOption{phd}{%
  \@phdtrue%
  \@msfalse%
  \renewcommand{\@uatype}{Dissertation}
  \renewcommand{\@uadegree}{Doctor of Philosophy}
}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ExecuteOptions{12pt}
\ProcessOptions

\LoadClass[12pt]{report}

\renewcommand{\baselinestretch}{2.00}

\input{thesis\@fntsz.clo}

\normalsize

%% Common page layout variables

\setlength{\textwidth}{6.00in}

\setlength{\headheight}{0.00em}
\setlength{\headsep}{0.00em}
\setlength{\parindent}{8.00ex}
\setlength{\baselineskip}{12pt plus2pt minus1pt}
\setlength{\tabcolsep}{1.0pt}

\setlength{\oddsidemargin}{0.50in}
\setlength{\evensidemargin}{\oddsidemargin}
\setlength{\topmargin}{0.00in}
\addtolength{\topmargin}{0.00em}

\setlength{\smallskipamount}{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength{\medskipamount}{6\p@ \@plus 2\p@ \@minus 2\p@}
\setlength{\bigskipamount}{12\p@ \@plus 4\p@ \@minus 4\p@}
\setlength{\topskip}{1.00em}
\setlength{\footskip}{2.00em}

%% Pagestyle definitions

\newcommand{\ps@thesis}{%
  \setlength{\textheight}{9.00in}
  \addtolength{\textheight}{-2.00em}

  \renewcommand{\@oddhead}{\@empty}%
  \renewcommand{\@evenhead}{\@oddhead}%
  \renewcommand{\@oddfoot}{\hfill {\thepage} \hfill}%
  \renewcommand{\@evenfoot}{\@oddfoot}}

\newcommand{\ps@emthesis}{%
  \setlength{\textheight}{9.00in}
  \addtolength{\textheight}{0.00em}

  \renewcommand{\@oddhead}{\@empty}%
  \renewcommand{\@evenhead}{\@oddhead}%
  \renewcommand{\@oddfoot}{\@empty}%
  \renewcommand{\@evenfoot}{\@oddfoot}%
}


\renewcommand\@makefntext[1]{%
  \setlength{\footnotesep}{\baselineskip}
  \renewcommand{\baselinestretch}{1.00}\normalsize\footnotesize%
  \parindent 1em%
  \noindent
  \hb@xt@1.8em{\hss\@makefnmark}#1}

\@settopoint{\textwidth}
\@settopoint{\textheight}
\@settopoint{\oddsidemargin}
\@settopoint{\marginparwidth}
\@settopoint{\evensidemargin}

\newcommand{\@secfont}{\rmfamily}
\newcommand{\@append}{chapter}

\clubpenalty=8000
\widowpenalty=8000

%% Thesis Preamble

\renewcommand{\maketitle}{\par
  \begingroup
  \global\@topnum\z@   % Prevents figures from going at top of page.
  \@maketitle
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
  \pagestyle{thesis}
}

  % Preamble Defaults

\renewcommand{\@title}{Create with %
  {\texttt{$\backslash$\MakeLowercase{title}\{\}}} command}
\renewcommand{\@author}{Create with %
  {\texttt{$\backslash$\MakeLowercase{author}\{\}}} command}
\newcommand{\@conferralmonth}{Create with %
  {\texttt{$\backslash$\MakeLowercase{conferraldate\{month\}\{year\}}}}
  command}
\newcommand{\@conferralyear}{}
\newcommand{\@advisor}{\@uasigline{Advisor}{Create with %
    {\texttt{$\backslash$\MakeLowercase{advisor}\{\}}} command}}
\newcommand{\@chair}{\@uasigline{Department Chair}{Create with %
    {\texttt{$\backslash$\MakeLowercase{chair}\{\}}} command}}
\newcommand{\@collegedean}{\@uasigline{Dean of the College}{Create with %
    {\texttt{$\backslash$\MakeLowercase{collegedean}\{\}}} command}}
\newcommand{\@gradschdean}{\@uasigline{Interim Dean of the Graduate School}{%
    Create with {\texttt{$\backslash$\MakeLowercase{gradschdean}\{\}}} command}}


\renewcommand{\title}[1]{\renewcommand{\@title}{#1}}
\renewcommand{\author}[1]{\renewcommand{\@author}{#1}}
\newcommand{\conferraldate}[2]%
           {\renewcommand{\@conferralmonth}{#1}%
             \renewcommand{\@conferralyear}{#2}%
           }

%% Signature Page
 
\newcommand{\@uasigline}[2]{%
  \parbox[t]{\linewidth}{%
    \renewcommand{\baselinestretch}{1.00}\large\normalsize\vspace{3.00em}%
    \underline{\makebox[0.975\linewidth]{}}\vspace{0.00em}\\%
    #1\\%
    #2
  }
}

\newcommand{\advisor}[1]{%
  \renewcommand{\@advisor}{\@uasigline{Advisor}{#1}}}
\newcommand{\chair}[1]{%
  \renewcommand{\@chair}{\@uasigline{Department Chair}{#1}}}
\newcommand{\collegedean}[1]{%
  \renewcommand{\@collegedean}{\@uasigline{Dean of the College}{#1}}}
\newcommand{\gradschdean}[1]{%
  \renewcommand{\@gradschdean}{\@uasigline{Interim Dean of the Graduate School}{#1}}}

\newif\if@coadvflag
\@coadvflagfalse
\newcommand{\coadvisor}[1]{%
  \@coadvflagtrue
  \newcommand{\@coadv}{\@uasigline{Co-Advisor}{#1}}
}


\newcount\@facreadcnt
\newcount\@tmpcnt\@tmpcnt\z@
\newif\if@facreadflag
\@facreadflagfalse
\newcommand{\facreader}[1]{%
  \@facreadflagtrue
  \global\advance\@facreadcnt1
  \expandafter\gdef\csname @facread\romannumeral\the\@facreadcnt\endcsname{%
    \@uasigline{Faculty Reader}{#1}}
}
\newcommand{\@facreader}{%
  \ifnum\@tmpcnt<\@facreadcnt
  \global\advance\@tmpcnt1
  \csname @facread\romannumeral\the\@tmpcnt\endcsname
  \let\next\@facreader
  \else
  \let\next\relax
  \fi
  \next
}


\newcount\@commcnt
\newcount\@tmpcnt\@tmpcnt\z@
\newif\if@commflag
\@commflagfalse
\newcommand{\committee}[1]{%
  \@commflagtrue
  \global\advance\@commcnt1
  \expandafter\gdef\csname @comm\romannumeral\the\@commcnt\endcsname{%
    \@uasigline{Committee Member}{#1}}
}
\newcommand{\@committee}{%
  \ifnum\@tmpcnt<\@commcnt
  \global\advance\@tmpcnt1
  \csname @comm\romannumeral\the\@tmpcnt\endcsname
  \let\next\@committee
  \else
  \let\next\relax
  \fi
  \next
}


\newcommand{\@uaapproved}{%
  \if@ms
      {\@advisor}\\
      \if@coadvflag
          {\@coadv}\\
      \fi
      \if@facreadflag
          {\@facreader}\\
      \fi
      {\@chair}
  \fi
  \if@phd
      {\@advisor}\\
      \if@coadvflag
          {\@coadv}\\
      \fi
      \if@commflag
          {\@committee}
      \fi
  \fi
}


\newcommand{\@uaaccepted}{%
  \if@ms
      {\@collegedean}\\
      {\@gradschdean}\\
      \@uasigline{Date}{\null}
  \fi
  \if@phd
      {\@chair}\\
      {\@collegedean}\\
      {\@gradschdean}\\
      \@uasigline{Date}{\null}
  \fi
}


\newcommand{\@maketitle}{%
  \begingroup

  % Copyright page
  \if@cpyr
  \thispagestyle{emthesis}%
  \centering \par\null \vfill%
  \copyright {\@conferralyear} \par%
  \MakeUppercase{\@author} \par%
  \MakeUppercase{All Rights Reserved} \par%
  \clearpage
  \fi

  % Title page
  \pagenumbering{roman}%
  \setcounter{page}{1}
  \thispagestyle{emthesis}
  \centering \@secfont%
  \vspace*{72\p@}%
  \vspace*{-\baselineskip}%
  \MakeUppercase{\@title} \null%
  \normalsize \vfill \vfill \vfill %
  A {\@uatype}\\Presented to\\The Graduate Faculty of %
  The University of Akron%
  \par \vfill \vfill%
  In Partial Fulfillment\\of the Requirements for the Degree\\%
  {\@uadegree}\\%
  \vfill \vfill \@author\\%
  {\@conferralmonth,}~{\@conferralyear}%

  \endgroup

  \pagestyle{thesis}
  \clearpage

  % Signature Page

  \begingroup
  \centering \@secfont%
  \vspace*{72\p@}%
  \vspace*{-\baselineskip}%
  \MakeUppercase{\@title} \null%
  \normalsize \vfill \vfill \vfill \@author\\
  \vfill \vfill \vfill {\@uatype} \vfill%
  \begin{minipage}[t]{0.46875\linewidth} \raggedright
    Approved:\\
    {\@uaapproved}
  \end{minipage} \hfill
  \begin{minipage}[t]{0.46875\linewidth} \raggedright
    Accepted:\\
    {\@uaaccepted}
  \end{minipage}\par \vfill%
  \endgroup

  \include{abstract}
  \include{acknowledgements}

  \tableofcontents
  \if@lot \listoftables \fi
  \if@lof \listoffigures \fi

  \addtocontents{toc}{\MakeUppercase{\@chapapp} \newline}

  \clearpage
  \pagenumbering{arabic}
}

\renewcommand{\abstract}{\chapter*{Abstract} \suppressfloats}
\renewcommand{\endabstract}{}

\def\acknowledgements{\chapter*{Acknowledgements}}
\def\endacknowledgements{}



%\renewcommand{\acknowledgements}{\chapter*{\Acknowledgements} \suppressfloats}
%\renewcommand{\endacknowledgements}{}

%% Table of Contents

\renewcommand{\contentsname}{Table of Contents}
\renewcommand{\tableofcontents}{%
  \if@twocolumn
  \@restonecoltrue\onecolumn
  \else
  \@restonecolfalse
  \fi
  \chapter*{\normalsize \contentsname
    \@mkboth{%
      \MakeUppercase{\contentsname}}{\MakeUppercase{\contentsname}}}%
\vskip -\baselineskip
  \hfill Page \relax \newline%
  \@starttoc{toc}%
  \if@restonecol\twocolumn\fi
}

%% List of Tables

\renewcommand{\listoftables}{%
  \if@twocolumn
  \@restonecoltrue\onecolumn
  \else
  \@restonecolfalse
  \fi
  \chapter*{\normalsize \listtablename
    \@mkboth{\MakeUppercase{\listtablename}}%
            {\MakeUppercase{\listtablename}}}%
  \addcontentsline{toc}{part}{\MakeUppercase{\listtablename}}%
\vskip -\baselineskip
  Table \hfill Page \relax \newline%
  \@starttoc{lot}%
  \if@restonecol\twocolumn\fi
}

%% List of Figures

\renewcommand{\listoffigures}{%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{\normalsize \listfigurename
    \@mkboth{\MakeUppercase{\listfigurename}}%
            {\MakeUppercase{\listfigurename}}}%
  \addcontentsline{toc}{part}{\MakeUppercase{\listfigurename}}%
\vskip -\baselineskip
  Figure \hfill Page \relax \newline%
  \@starttoc{lof}%
  \if@restonecol\twocolumn\fi
}

%% dotted toc lines

\newlength{\@appindent}

\renewcommand*{\l@part}{\@dottedtocline{0}{0.00em}{0.00em}}

\newlength{\@chtocnum}
\settowidth{\@chtocnum}{III.}
%\addtolength{\@chtocnum}{1.00em}
\addtolength{\@chtocnum}{0.50em}
%\renewcommand*{\l@chapter}{\@dottedtocline{1}{1.50em}{\@chtocnum}}
\renewcommand*{\l@chapter}{\@dottedtocline{0}{0.0em}{\@chtocnum}}
\newlength{\@apptocnum}
\settowidth{\@apptocnum}{\MakeUppercase{\@chapapp}\hspace{1.00ex}{A.}}
\@settopoint{\@apptocnum}
\addtolength{\@apptocnum}{1.00em}
\@settopoint{\@apptocnum}
\newcommand*{\l@appendix}{\@dottedtocline{0}{0.00em}{\@apptocnum}}
%\newcommand*{\l@appendices}{\@dottedtocline{2}{2.50em}{\@apptocnum}}
\newcommand*{\l@appendices}{\@dottedtocline{2}{2.0em}{\@apptocnum}}

%\newcommand*{\l@appendices}{\@dottedtocline{2}{1.50em}{\@apptocnum}}
\newlength{\@sectocnum}
\settowidth{\@sectocnum}{77}
\addtolength{\@sectocnum}{1.00em}
%\renewcommand*{\l@section}{\@dottedtocline{2}{3.00em}{\@sectocnum}}
%\renewcommand*{\l@section}{\@dottedtocline{2}{2.50em}{\@sectocnum}}
\renewcommand*{\l@section}{\@dottedtocline{2}{2.00em}{\@sectocnum}}
\renewcommand*{\l@subsection}{\@dottedtocline{3}{4.50em}{3.50em}}
\renewcommand*{\l@subsubsection}{\@dottedtocline{4}{6.00em}{4.50em}}
\renewcommand*{\l@paragraph}{\@dottedtocline{5}{7.50em}{5.00em}}
\renewcommand*{\l@subparagraph}{\@dottedtocline{6}{9.00em}{6.00em}}

\newlength{\@figtocnum}
\settowidth{\@figtocnum}{9.99}
\addtolength{\@figtocnum}{1.00em}
\renewcommand*{\l@figure}{\@dottedtocline{0}{0.00em}{\@figtocnum}}

\newlength{\@tabtocnum}
\settowidth{\@tabtocnum}{9.99}
\addtolength{\@tabtocnum}{1.00em}
\renewcommand*{\l@table}{\@dottedtocline{0}{0.00em}{\@figtocnum}}

\renewcommand{\numberline}[1]{\makebox[\@tempdima][l]{{#1}\hfill}}
\renewcommand{\@dottedtocline}[5]{
  \ifnum #1>\c@tocdepth
  \else
      {\noindent \begin{minipage}[t]{\linewidth}%
          \renewcommand{\baselinestretch}{1.00}\large\normalsize%
          \setlength{\@tempdima}{#3}%
          \setlength{\rightskip}{\@tocrmarg} \addtolength{\rightskip}{24\p@}%
          \setlength{\parfillskip}{-\rightskip}
          \setlength{\leftskip}{#2}
          \setlength{\parindent}{\leftskip}%
          \addtolength{\parindent}{-0.50ex}%
          \addtolength{\leftskip}{0.50ex}%
          \addtolength{\leftskip}{\@tempdima} \@afterindenttrue
          \interlinepenalty\@M
          \leavevmode
          \null \nobreak \hspace{-\leftskip}
              {{#4}\nobreak%
                \leaders\hbox{$\m@th \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
                  mu$}\hfil \nobreak%
                \makebox[\@pnumwidth][r]{\hfil{#5}}\par\null}%
        \end{minipage} \newline}%
      \fi}

%% Appendix

\newcommand{\@appendnum}{0}
\renewcommand{\appendix}[1]{
	\renewcommand{\@appendnum}{#1}%
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \renewcommand{\@chapapp}{\appendixname}%
  \renewcommand{\@chptoc}{}%
  \renewcommand{\@chplist}{}%
  \renewcommand{\@chpnum}{\@Alph}%
  \ifnum \@appendnum>1
    \renewcommand{\@append}{appendices}
    \renewcommand{\@chptoc}{\MakeUppercase{\@chapapp}\space}%
    \if@openright
      \cleardoublepage
    \else
      \clearpage
    \fi
    \thispagestyle{thesis}%
    \null \vfill \vfill \vfill \vfill
    \begin{center}
      \MakeUppercase{Appendices}
    \end{center}
    \vfill \vfill \vfill \vfill \vfill \vfill \null
    \addcontentsline{toc}{part}{\MakeUppercase{Appendices}}%
  \else
  \renewcommand{\@chptoc}{\MakeUppercase{\@chapapp}\space}%
    \renewcommand{\@append}{appendix}
%    \addtocontents{toc}{\noindent \MakeUppercase{\@chapapp} \hfill%
%      \relax \newline} \null
  \fi
}

%% Bibliography

\def\thebibliography#1{\chapter*{\bibname\@mkboth
  {\bibname}{\bibname}}\list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
  \advance\leftmargin\labelsep
  \setlength{\parsep}{0pt}
  \setlength{\itemsep}{\baselineskip}
  \addtolength{\itemsep}{-14pt}
  \usecounter{enumi}}
  \def\newblock{\hskip .11em plus .33em minus .07em}
  \fussy\clubpenalty8000\widowpenalty8000
  \sfcode`\.=1000\relax
  \addcontentsline{toc}{part}{\MakeUppercase{\bibname}}%
}

\def\bibitem{\@ifnextchar[\@lbibitem\@bibitem}
\def\@lbibitem[#1]#2{\item[\@biblabel{#1}\hfill]\if@filesw
		{\let\protect\noexpand
    \immediate
    \write\@auxout{\string\bibcite{#2}{#1}}}\fi\ignorespaces
 	\renewcommand{\baselinestretch}{1.00} \large \normalsize
}
\def\@bibitem#1{\item\if@filesw \immediate\write\@auxout
  {\string\bibcite{#1}{\the\value{\@listctr}}}\fi\ignorespaces
 	\renewcommand{\baselinestretch}{1.00} \large \normalsize
  \fussy\interlinepenalty=10000
}

%% Figures

\renewcommand{\fnum@figure}{{\@secfont Figure \thefigure}}
\renewcommand{\thefigure}
     {\ifnum\c@chapter>\z@\thechapter.\fi\@arabic\c@figure}

\newtoks\@stfigure
\def\subfigures{\refstepcounter{figure}%
  \edef\@savedfigure{\the\c@figure}%
  \@stfigure=\expandafter{\thefigure}%   %only want \thefigure
  \edef\@savedthefigure{\the\@stfigure}% % expanded once
  \edef\oldthefigure{\thefigure}%
  \setcounter{figure}{0}%
  \def\thefigure{\oldthefigure.\arabic{figure}}}

\def\endsubfigures{\setcounter{figure}{\@savedfigure}%
  \@stfigure=\expandafter{\@savedthefigure}%
  \edef\thefigure{\the\@stfigure}\global\@ignoretrue}

\setlength\belowcaptionskip{2.00em}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    \renewcommand{\baselinestretch}{1.00} \small \normalsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

%% Parts

\renewcommand{\chapter}{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \pagestyle{thesis}% 
  \global\@topnum\z@%
  \@afterindentfalse%
  \secdef\@chapter\@schapter}

\newcommand{\@chptoc}{}
\newcommand{\@chplist}{\@Roman}
\newcommand{\@chpnum}{\@arabic}
\renewcommand{\thechapter}{\@chpnum\c@chapter}

\def\@chapter[#1]#2{%
  \ifnum \c@secnumdepth > \m@ne
    \refstepcounter{chapter}%
    \ifnum \@appendnum = 1%
      \typeout{\@chapapp.}%
      \addcontentsline{toc}{part}{\MakeUppercase{\appendixname}}%
    \else
      \typeout{\@chapapp\space\thechapter.}%
      \addcontentsline{toc}{\@append}{%
        \protect\numberline{\@chptoc\@chplist\thechapter.}%
          {\MakeUppercase{#1}}}%
    \fi
  \else
    \addcontentsline{toc}{\@append}{#1}%
  \fi
  \chaptermark{#1}%
  \if@twocolumn
    \@topnewpage[\@makechapterhead{#2}]%
  \else
    \@makechapterhead{#2}%
    \@afterheading
  \fi}

\def\@makechapterhead#1{%
  \vspace*{72\p@}%
  \vspace*{-\baselineskip}%
  {\parindent \z@ \centering \normalfont%
    \ifnum \c@secnumdepth > \m@ne%
      \ifnum \@appendnum = 1%
        \MakeUppercase{\@chapapp}%
      \else
        \MakeUppercase{\@chapapp}\space{\@chplist\thechapter}%
      \fi
      \par\nobreak%
    \fi%
    \interlinepenalty\@M%
    \MakeUppercase {#1}\par\nobreak%
    \vspace{1.00em}%
  }
}

\def\@schapter#1{%
  \if@twocolumn%
    \@topnewpage[\@makeschapterhead{#1}]%
  \else%
    \@makeschapterhead{#1}%
    \@afterheading%
  \fi%
}

\def\@makeschapterhead#1{%
  \vspace*{72\p@}%
  \vspace*{-\baselineskip}%
  {\parindent \z@ \centering%
    \normalfont%
    \interlinepenalty\@M%
    \MakeUppercase{#1} \par\nobreak%
    \vskip 40\p@
  }
}

\renewcommand{\section}{\@startsection{section}{1}{\z@}%
  {-2.0em plus -1.0em minus -0.2em}{1.0em plus 0.20em}%
  {\normalsize \@secfont}}
\renewcommand{\subsection}{\@startsection{subsection}{2}{\z@}%
  {-2.0ex plus -1ex minus -.2ex}{0.5ex plus 0.2ex}%
  {\noindent \@secfont}}
\renewcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}%
  {-2.0ex plus -1ex minus -.2ex}{0.5ex plus 0.2ex}%
  {\noindent \@secfont}}
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
  {3.25ex plus 1ex minus .2ex}{-1em}%
  {\@secfont}}
\renewcommand{\subparagraph}{\@startsection{subparagraph}{4}{\parindent}%
  {3.25ex plus 1ex minus .2ex}{-1em}%
  {\@secfont}}

\renewcommand{\@listi}{\leftmargin\leftmargini
    \parsep 4\p@ \@plus2\p@ \@minus\p@
    \topsep 8\p@ \@plus2\p@ \@minus4\p@
    \itemsep4\p@ \@plus2\p@ \@minus\p@}
\let\@listI\@listi
\@listi
\renewcommand{\@listii}{\leftmargin\leftmarginii
    \labelwidth\leftmarginii
    \advance\labelwidth-\labelsep
    \topsep 4\p@ \@plus2\p@ \@minus\p@
    \parsep 2\p@ \@plus\p@  \@minus\p@
    \itemsep \parsep}
\renewcommand{\@listiii}{\leftmargin\leftmarginiii
    \labelwidth\leftmarginiii
    \advance\labelwidth-\labelsep
    \topsep 2\p@ \@plus\p@\@minus\p@
    \parsep \z@
    \partopsep \p@ \@plus\z@ \@minus\p@
    \itemsep \topsep}
\renewcommand{\@listiv}{\leftmargin\leftmarginiv
    \labelwidth\leftmarginiv
    \advance\labelwidth-\labelsep}
\renewcommand{\@listv}{\leftmargin\leftmarginv
    \labelwidth\leftmarginv
    \advance\labelwidth-\labelsep}
\renewcommand{\@listvi}{\leftmargin\leftmarginvi
    \labelwidth\leftmarginvi
    \advance\labelwidth-\labelsep}


\normalfont
\pagestyle{emthesis}

\endinput

