\documentclass[ignorenonframetext,]{beamer}
\setbeamertemplate{caption}[numbered]
\setbeamertemplate{caption label separator}{:}
\setbeamercolor{caption name}{fg=normal text.fg}
\usepackage{amssymb,amsmath}
\usepackage{ifxetex,ifluatex}
\usepackage{fixltx2e} % provides \textsubscript
\usepackage{lmodern}
\ifxetex
  \usepackage{fontspec,xltxtra,xunicode}
  \defaultfontfeatures{Mapping=tex-text,Scale=MatchLowercase}
  \newcommand{\euro}{€}
\else
  \ifluatex
    \usepackage{fontspec}
    \defaultfontfeatures{Mapping=tex-text,Scale=MatchLowercase}
    \newcommand{\euro}{€}
  \else
    \usepackage[T1]{fontenc}
    \usepackage[utf8]{inputenc}
      \fi
\fi
% use upquote if available, for straight quotes in verbatim environments
\IfFileExists{upquote.sty}{\usepackage{upquote}}{}
% use microtype if available
\IfFileExists{microtype.sty}{\usepackage{microtype}}{}
\usepackage{longtable,booktabs}
\usepackage{caption}
% These lines are needed to make table captions work with longtable:
\makeatletter
\def\fnum@table{\tablename~\thetable}
\makeatother

% Comment these out if you don't want a slide with just the
% part/section/subsection/subsubsection title:
\AtBeginPart{
  \let\insertpartnumber\relax
  \let\partname\relax
  \frame{\partpage}
}
\AtBeginSection{
  \let\insertsectionnumber\relax
  \let\sectionname\relax
  \frame{\sectionpage}
}
\AtBeginSubsection{
  \let\insertsubsectionnumber\relax
  \let\subsectionname\relax
  \frame{\subsectionpage}
}

\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt plus 2pt minus 1pt}
\setlength{\emergencystretch}{3em}  % prevent overfull lines
\setcounter{secnumdepth}{0}

\date{}

\begin{document}

\begin{frame}

\thispagestyle{empty}

\vspace{3cm}

\vspace{4cm}

\Huge\textbf{Gaussian-Process-Model-based}

\vspace{4mm}

\Huge\textbf{System-Identification}

\vspace{4mm}

\Huge\textbf{Toolbox for Matlab}

\vspace{4mm}

\Huge\textbf{Version 1.2.2}

\vspace{3mm}

\vspace{3cm}

{\LargeMartin Stepančič and Juš Kocijan}

\vfill

\Large\today

\pagebreak\setcounter{page}{1} \newpage

\newpage\thispagestyle{empty} .\newpage\pagenumbering{arabic}
\setcounter{page}{1}

\end{frame}

\begin{frame}{Introduction}

The idea of this toolbox is to facilitate dynamic systems identification
with Gaussian-process (GP) models. The presented toolbox is continuously
developing and is put together with hope to be useful as a springboard
for the modelling of dynamic systems with GP models.

The GP model belongs to the class of black-box models. GP modelling
differs from most other black-box identification approaches in that it
does not try to approximate the modelled system by fitting the
parameters of the selected basis functions, but rather it searches for
the relationship among the measured data. The model is composed of
input-output data that describes the behaviour of the modelled system
and the covariance function that describes the relation with respect to
the input-output data. The prediction of the GP model output is given as
a normal distribution, expressed in terms of the mean and the variance.
The mean value represents the most likely output, and the variance can
be interpreted as a measure of its confidence.

System identification is composed of methods to build mathematical
models of dynamic systems from measured data. It is one of the
scientific pillars used for dynamic-systems analysis and control design.
The identification of a dynamic system means that we are looking for a
relationship between past observations and future outputs.
Identification can be interpreted as the concatenation of a mapping from
measured data to a regression vector, followed by a nonlinear mapping
from the regression vector to the output space. Various machine-learning
methods and statistical methods are employed to determine the nonlinear
mapping from the regression vector to the output space. One of the
possible methods for a description of the nonlinear mapping used in
identification is GP models. It is straightforward to employ GP models
for the discrete-time modelling of dynamic systems within the
prediction-error framework.

Many dynamic systems are often considered as complex; however,
simplified input-output behaviour representations are sufficient for
certain purposes, e.g., feedback control design, prediction models for
supervisory control, etc.

More on the topic of system identification with GP models and the use of
this models for control design can be found in the book:\\Juš Kocijan
(2016) Modelling and Control of Dynamic Systems Using Gaussian Process
Models, Springer.

\end{frame}

\begin{frame}{GP-Model-based System-Identification Toolbox for Matlab}

\begin{block}{Prerequisites}

As this toolbox is intended to use within Matlab environment the user
should have Matlab installed. It works on Matlab 7 and later, but there
should be no problems using the toolbox on previous versions of Matlab,
e.g., 6 or 5.

It is also assumed that the GPML toolbox\footnote<.->{It can be obtained
  from \emph{http://www.gaussianprocess.org/gpml}.}, general purpose GP
modelling toolbox for Matlab, is installed. The GP-model-based
system-identification toolbox serves as upgrade to GPML toolbox.

The user should posses some familiarity with the Matlab structure and
programming.

\end{block}

\begin{block}{Installing GPdyn toolbox}

Unzip the file GPdyn into chosen directory and add path, with
subdirectories, to Matlab path.

\end{block}

\begin{block}{Overview of the GPdyn toolbox}

GPdyn files are contained in several directories, depending on their
purpose:

\begin{description}[<+->]
\item[training functions,]
used for training GP models of dynamic systems;
\item[GP-model evaluation functions,]
used for simulating the dynamic GP model;
\item[LMGP-model evaluation functions,]
which are used when modelling and simulating the system with a GP model
with incorporated local models (LMGP model);
\item[utilities functions,]
that are various support functions;
\item[demo functions,]
which demonstrate the use of the toolbox for identification of dynamic
systems.
\end{description}

\clearpage

The list of included functions, demos and one model is given in
following tables.

\renewcommand{\arraystretch}{1.1}

%%%%%%%%%%%%%%%%%%%%%%%%%%  VSE TABELE SKUPAJ

\begin{longtable}[c]{@{}ll@{}}
\toprule
\multicolumn{2}{|l|}{\tabelarule \bf GP-model training functions}
&\tabularnewline
\midrule
\endhead
\fun{trainGParx} & GP-model training of ARX model\tabularnewline
\fun{trainGPoe} & GP-model training of OE model\tabularnewline
\fun{gp\_initial} & - finding initial values of hyperparameters with
random search\tabularnewline
\fun{minimizeDE} & minimize a multivariate function using differential
evolution\tabularnewline
\multicolumn{2}{l}{\tabelarule \bf } &\tabularnewline
\multicolumn{2}{|l|}{\tabelarule \bf Covariance functions}
&\tabularnewline
\multicolumn{2}{|l|}{\tabelarule included and explained in enclosed GPML
toolbox} &\tabularnewline
\multicolumn{2}{l}{\tabelarule \bf } &\tabularnewline
\multicolumn{2}{|l|}{\tabelarule \bf GP-model evaluation}
&\tabularnewline
\fun{simulGPnaive} & GP model simulation without the propagation of
uncertainty\tabularnewline
\fun{simulGPmcmc} & GP model simulation with Monte Carlo
approximation\tabularnewline
\fun{simulGPtaylorSE} & GP model simulation with analytical
approximation of statistical\tabularnewline
& moments with a Taylor expansion for the squared
exponential\tabularnewline
& covariance function\tabularnewline
\fun{simulGPexactSE} & GP model simulation with exact matching of
statistical moments\tabularnewline
& for the squared exponential covariance function\tabularnewline
\fun{simulGPexactLIN} & GP model simulation with exact matching of
statistical moments\tabularnewline
& for the linear covariance function\tabularnewline
\fun{predGPnaive} & multi-step-ahead prediction of GP model
without\tabularnewline
& the propagation of uncertainty\tabularnewline
\fun{gpx} & modified version of GP rutine from the GPML
toolbox\tabularnewline
\fun{gmx\_sample} & creates samples of mixture components\tabularnewline
\fun{gpTaylorSEard} & GP model prediction with stochastic inputs
for\tabularnewline
& the squared exponential covariance function with Taylor
expansion\tabularnewline
\fun{gpExactLINard} & GP model prediction with stochastic inputs
for\tabularnewline
& the linear covariance function\tabularnewline
\fun{gpExactSEard} & GP model prediction with stochastic inputs
for\tabularnewline
& the squared exponential covariance function\tabularnewline
\bottomrule
\end{longtable}

\vspace{3mm}

\begin{tabular}{|l|l|c|}
%


\multicolumn{2}{l}{\tabelarule \bf } \\ \hline
%
 \multicolumn{2}{|l|}{\tabelarule \bf LMGP-model evaluation} \\
  \hline \fun{simulLMGPnaive} & LMGP model simulation without the propagation of uncertainty \\
 \hline \fun{simulLMGPmcmc} & LMGP model simulation with Monte Carlo approximation\\
  \hline \fun{trainLMGP} & LMGP model training \\
  \hline \fun{gpSD00} & - LMGP model prediction \\
  & - data likelihood and its derivatives \\
\hline
%

\multicolumn{2}{l}{\tabelarule \bf } \\


 \hline \multicolumn{2}{|l|}{\tabelarule \bf Supporting functions}\\
 \hline \fun{add\_noise\_to\_vector} & adding white noise to noise-free simulation results\\
 \hline \fun{construct} & construction of the input regressors\\
  & from system's input signals\\
 \hline \fun{eval\_func} & method to evaluate covariance, mean and likelihood functions\\
 \hline \fun{likelihood} & calculates negative log marginal likelihood\\
 \hline \fun{lipschitz} & the method for the lag-space selection, based on Lipschitz quotients\\
 \hline \fun{validate} & checking of the parameters match \\
 \hline \fun{loss} & performance measures \\
 \hline \fun{mcmc\_test\_pdfs} & testing sampled probability distributions\\
 \hline \fun{plotgp} & plot results (output and error) of the GP model prediction \\
 \hline \fun{plotgpe} & plot error of the GP model prediction \\
 \hline \fun{plotgpy} & plot output of the GP model prediction \\
 \hline \fun{preNorm} & preprocessing of data \\
 \hline \fun{postNorm} & postprocessing of data \\
 \hline \fun{postNormVar} & postprocessing of predicted variance\\
 \hline \fun{sig\_prbs} & generating pseudo-random binary signal \\
 \hline \fun{sig\_prs\_minmax} & generating pseudo-random signal \\ \hline


\end{tabular}

\begin{tabular}{|l|l|}


 \multicolumn{2}{l}{\tabelarule \bf } \\

%
 \hline \multicolumn{2}{|l|}{\tabelarule \bf Demos}\\

 \hline \fun{demo\_example\_present} & present the system used in demos  \\
 \hline \fun{demo\_example\_gp\_data} & generate data for the identification and validation  \\
 & of the GP model \\
 \hline \fun{demo\_example\_gp\_norm} & normalization of input and output data \\
 \hline \fun{demo\_example\_gp\_training} & training of the GP model \\
 \hline \fun{demo\_example\_gp\_simulation} & validation with simulation of the GP model \\
 \hline \fun{demo\_example\_lmgp\_data} & generate data for the identification and validation \\
 &  of the LMGP model \\
 \hline \fun{demo\_example\_lmgp\_training} & training of the LMGP model \\
 \hline \fun{demo\_example\_lmgp\_simulation} & simulation of the LMGP model\\
 \hline \fun{demo\_example} & system simulation\\
 \hline \fun{demo\_example\_derivative} & obtaining system's derivatives\\
 \hline \fun{demo\_example\_LM\_ident} & identification of system's local models \\ \hline

\end{tabular}

\clearpage

\end{block}

\begin{block}{How to use this toolbox}

\begin{block}{Demos}

A simple nonlinear dynamic system is used to demonstrate the
identification and simulation of the GP models:
\[y(k+1) = \frac{y(k)}{1+y^2(k)} + u^3(k) \label{eq:narendra}\] The
system was used as an example of dynamic system identification with
artificial neural networks in:\\K.S. Narendra and K. Parthasarathy.
Identification and Control of Dynamical Systems Using Neural Networks,
IEEE Transactions on Neural Networks, Vol.1 No. 1, 4--27, 1990.

\begin{description}[<+->]
\item[demo\_example\_present,]
presents this system.
\end{description}

Following three demos present the identification of dynamic systems with
the GP model:

\begin{description}[<+->]
\item[demo\_example\_gp\_data,]
which presents how to obtain and assemble data for identification;
\item[demo\_example\_gp\_norm,]
which shows how to normalise input and output data for training;
\item[demo\_example\_gp\_training,]
which demonstrates the identification with a GP model;
\item[demo\_example\_gp\_simulation,]
which shows how to simulate the GP model.
\end{description}

The use of the GP model with incorporated local models is presented with
demos:

\begin{description}[<+->]
\item[demo\_example\_lmgp\_data,]
which presents how to obtain and assemble data for identification;
\item[demo\_example\_lmgp\_training,]
which demonstrates the training (=identifying) the LMGP model;
\item[demo\_example\_lmgp\_simulation,]
which shows how to simulate the LMGP model.
\end{description}

\end{block}

\begin{block}{Acknowledgements}

We would like to thank all past, present and future contributors to this
toolbox.

\end{block}

\end{block}

\end{frame}

\end{document}
