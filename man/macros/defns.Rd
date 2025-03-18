%% macro definitions for spatstat man pages
%% Authors
\newcommand{\adrian}{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}}
\newcommand{\rolf}{Rolf Turner \email{rolfturner@posteo.net}}
\newcommand{\ege}{Ege Rubak \email{rubak@math.aau.dk}}
\newcommand{\spatstatAuthors}{\adrian, \rolf and \ege}
\newcommand{\spatstatAuthorsComma}{\adrian, \rolf, \ege}
%% Contributors with emails
\newcommand{\pavel}{Pavel Grabarnik \email{pavel.grabar@issp.serpukhov.su}}
\newcommand{\dominic}{Dominic Schuhmacher \email{dominic.schuhmacher@mathematik.uni-goettingen.de}, URL \code{http://dominic.schuhmacher.name/}}
\newcommand{\wei}{Ang Qi Wei \email{aqw07398@hotmail.com}}
\newcommand{\colette}{Marie-Colette van Lieshout \email{Marie-Colette.van.Lieshout@cwi.nl}}
\newcommand{\rasmus}{Rasmus Plenge Waagepetersen \email{rw@math.auc.dk}}
\newcommand{\abdollah}{Abdollah Jalilian \email{jalilian@razi.ac.ir}}
\newcommand{\ottmar}{Ottmar Cronie \email{ottmar@chalmers.se}}
\newcommand{\stephenEglen}{Stephen Eglen \email{S.J.Eglen@damtp.cam.ac.uk}}
\newcommand{\mehdi}{Mehdi Moradi \email{m2.moradi@yahoo.com}}
\newcommand{\yamei}{Ya-Mei Chang \email{yamei628@gmail.com}}
\newcommand{\martinH}{Martin Hazelton \email{Martin.Hazelton@otago.ac.nz}}
\newcommand{\tilman}{Tilman Davies \email{Tilman.Davies@otago.ac.nz}}
% Names with accents
\newcommand{\Bogsted}{\ifelse{latex}{\out{B\o gsted}}{Bogsted}}
\newcommand{\Cramer}{\ifelse{latex}{\out{Cram\'er}}{Cramer}}
\newcommand{\Francois}{\ifelse{latex}{\out{Fran\c{c}ois}}{Francois}}
\newcommand{\Frederic}{\ifelse{latex}{\out{Fr{\'e}d{\'e}ric}}{Frederic}}
\newcommand{\Hogmander}{\ifelse{latex}{\out{H{\"o}gmander}}{Hogmander}}
\newcommand{\Jyvaskyla}{\ifelse{latex}{\out{Jyv\"askyl\"a}}{Jyvaskyla}}
\newcommand{\Lucia}{\ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia}}
\newcommand{\Matern}{\ifelse{latex}{\out{Mat\'ern}}{Matern}}
\newcommand{\Moller}{\ifelse{latex}{\out{M\o ller}}{Moller}}
\newcommand{\Oehlschlaegel}{\ifelse{latex}{\out{Oehlschl\"{a}gel}}{Oehlschlaegel}}
\newcommand{\Prokesova}{\ifelse{latex}{\out{Proke\u{s}ov{\'{a}}}}{Prokesova}}
\newcommand{\Sarkka}{\ifelse{latex}{\out{S\"{a}rkk\"{a}}}{Sarkka}}
\newcommand{\Sanchez}{\ifelse{latex}{\out{S\'{a}nchez}}{Sanchez}}
\newcommand{\Martin}{\ifelse{latex}{\out{Mart\'{\i}n}}{Martin}}
\newcommand{\Dominguez}{\ifelse{latex}{\out{Dom\'{\i}nguez}}{Dominguez}}
\newcommand{\Rodriguez}{\ifelse{latex}{\out{Rodr\'{\i}guez}}{Rodriguez}}
\newcommand{\Gonzalez}{\ifelse{latex}{\out{Gonz\'{a}lez}}{Gonzalez}}
%% List of all Gibbs interactions 
\newcommand{\GibbsInteractionsList}{\code{\link[MPKG]{AreaInter}}, \code{\link[MPKG]{BadGey}}, \code{\link[MPKG]{Concom}}, \code{\link[MPKG]{DiggleGatesStibbard}}, \code{\link[MPKG]{DiggleGratton}}, \code{\link[MPKG]{Fiksel}}, \code{\link[MPKG]{Geyer}}, \code{\link[MPKG]{Hardcore}}, \code{\link[MPKG]{HierHard}}, \code{\link[MPKG]{HierStrauss}}, \code{\link[MPKG]{HierStraussHard}}, \code{\link[MPKG]{Hybrid}}, \code{\link[MPKG]{LennardJones}}, \code{\link[MPKG]{MultiHard}}, \code{\link[MPKG]{MultiStrauss}}, \code{\link[MPKG]{MultiStraussHard}}, \code{\link[MPKG]{OrdThresh}}, \code{\link[MPKG]{Ord}}, \code{\link[MPKG]{Pairwise}}, \code{\link[MPKG]{PairPiece}}, \code{\link[MPKG]{Penttinen}}, \code{\link[MPKG]{Poisson}}, \code{\link[MPKG]{Saturated}}, \code{\link[MPKG]{SatPiece}}, \code{\link[MPKG]{Softcore}}, \code{\link[MPKG]{Strauss}}, \code{\link[MPKG]{StraussHard}} and \code{\link[MPKG]{Triplets}}}
%% List of interactions recognised by RMH code
\newcommand{\rmhInteractionsList}{\code{\link[MPKG]{AreaInter}}, \code{\link[MPKG]{BadGey}}, \code{\link[MPKG]{DiggleGatesStibbard}}, \code{\link[MPKG]{DiggleGratton}}, \code{\link[MPKG]{Fiksel}}, \code{\link[MPKG]{Geyer}}, \code{\link[MPKG]{Hardcore}}, \code{\link[MPKG]{Hybrid}}, \code{\link[MPKG]{LennardJones}}, \code{\link[MPKG]{MultiStrauss}}, \code{\link[MPKG]{MultiStraussHard}}, \code{\link[MPKG]{PairPiece}}, \code{\link[MPKG]{Penttinen}}, \code{\link[MPKG]{Poisson}}, \code{\link[MPKG]{Softcore}}, \code{\link[MPKG]{Strauss}}, \code{\link[MPKG]{StraussHard}} and \code{\link[MPKG]{Triplets}}}
%% Frequent references
\newcommand{\baddrubaturnbook}{Baddeley, A., Rubak, E. and Turner, R. (2015) \emph{Spatial Point Patterns: Methodology and Applications with R}. Chapman and Hall/CRC Press. }
%% Citations of recent articles that will change rapidly
\newcommand{\baddchangclustersim}{Baddeley, A. and Chang, Y.-M. (2023) Robust algorithms for simulating cluster point processes. \emph{Journal of Statistical Computation and Simulation} \bold{93}, 1950--1975.}
\newcommand{\smoothpcfpaper}{Baddeley, A., Davies, T.M. and Hazelton, M.L. (2025) An improved estimator of the pair correlation function of a spatial point process. \emph{Biometrika}, to appear.}
\newcommand{\smoothpcfpapercite}{Baddeley, Davies and Hazelton (2025)}


