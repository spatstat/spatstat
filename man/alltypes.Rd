\name{alltypes}
\alias{alltypes}
\title{Calculate Summary Statistic for All Types in a Multitype Point Pattern}
\description{
  Given a marked point pattern, this computes the estimates of
  a selected summary function (\eqn{F},\eqn{G}, \eqn{J}, \eqn{K} etc)
  of the pattern, for all possible combinations of marks,
  and returns these functions in an array.
}

\usage{
  alltypes(X, fun="K", \dots,
           dataname=NULL,verb=FALSE,envelope=FALSE,reuse=TRUE)
}

\arguments{
  \item{X}{The observed point pattern, for which summary function
    estimates are required.  An object of class \code{"ppp"} or \code{"lpp"}.
  }

  \item{fun}{The summary function. Either an \R function,
    or a character string indicating the summary function
    required.  Options for strings are
    \code{"F"}, \code{"G"}, \code{"J"}, \code{"K"}, \code{"L"}, \code{"pcf"},
    \code{"Gcross"}, \code{"Jcross"}, \code{"Kcross"}, \code{"Lcross"},
    \code{"Gdot"}, \code{"Jdot"}, \code{"Kdot"}, \code{"Ldot"}.
  }
  \item{\dots}{
    Arguments passed to the summary function
    (and to the function \code{\link{envelope}} if appropriate)
  }

  \item{dataname}{Character string giving an optional (alternative)
    name to the point pattern, different from what is given
    in the call.  This name, if supplied, may be used by
    \code{\link{plot.fasp}()} in forming the title of the plot.
    If not supplied it defaults to the parsing of the argument
    supplied as \code{X} in the call.
  }

  \item{verb}{
    Logical value.  If \code{verb} is
    true then terse ``progress reports'' (just the values of the
    mark indices) are printed out when the calculations for that
    combination of marks are completed. 
  }
  \item{envelope}{
    Logical value. If \code{envelope} is true, then simulation envelopes
    of the summary function will also be computed. See Details.
  }
  \item{reuse}{
    Logical value indicating whether the envelopes in each panel
    should be based on the same set of simulated patterns
    (\code{reuse=TRUE}) or on different, independent sets of simulated
    patterns (\code{reuse=FALSE}).
  }
}

\details{
  This routine is a convenient way to analyse the dependence between
  types in a multitype point pattern.
  It computes the estimates of a selected summary function of the
  pattern, for all possible combinations of marks.
  It returns these functions in an array
  (an object of class \code{"fasp"}) amenable to plotting
  by \code{\link{plot.fasp}()}.

  The argument \code{fun} specifies the summary function that will
  be evaluated for each type of point, or for each pair of types.
  It may be either an \R function or a character string.
  
  Suppose that the points have possible types \eqn{1,2,\ldots,m}
  and let \eqn{X_i}{X[i]} denote the pattern of points of type \eqn{i} only.

  If \code{fun="F"} then this routine
  calculates, for each possible type \eqn{i},
  an estimate of the Empty Space Function \eqn{F_i(r)}{F[i](r)} of
  \eqn{X_i}{X[i]}. See \code{\link{Fest}}
  for explanation of the empty space function.
  The estimate is computed by applying \code{\link{Fest}}
  to \eqn{X_i}{X[i]} with the optional arguments \code{\dots}.

  If \code{fun} is
  \code{"Gcross"}, \code{"Jcross"}, \code{"Kcross"} or \code{"Lcross"},
  the routine calculates, for each pair of types \eqn{(i,j)},
  an estimate of the ``\code{i}-to\code{j}'' cross-type function
  \eqn{G_{ij}(r)}{G[i,j](r)},
  \eqn{J_{ij}(r)}{J[i,j](r)},
  \eqn{K_{ij}(r)}{K[i,j](r)} or
  \eqn{L_{ij}(r)}{L[i,j](r)} respectively describing the
  dependence between 
  \eqn{X_i}{X[i]} and \eqn{X_j}{X[j]}.
  See \code{\link{Gcross}}, \code{\link{Jcross}}, \code{\link{Kcross}}
  or \code{\link{Lcross}} respectively for explanation of these
  functions.
  The estimate is computed by applying the relevant function
  (\code{\link{Gcross}} etc)
  to \code{X} using each possible value of the arguments \code{i,j},
  together with the optional arguments \code{\dots}.
  
  If \code{fun} is \code{"pcf"} the routine calculates
  the cross-type pair correlation function \code{\link{pcfcross}}
  between each pair of types.

  If \code{fun} is 
  \code{"Gdot"}, \code{"Jdot"}, \code{"Kdot"} or \code{"Ldot"},
  the routine calculates, for each type \eqn{i},
  an estimate of the ``\code{i}-to-any'' dot-type function
  \eqn{G_{i\bullet}(r)}{G[i.](r)},
  \eqn{J_{i\bullet}(r)}{J[i.](r)} or
  \eqn{K_{i\bullet}(r)}{K[i.](r)} or
  \eqn{L_{i\bullet}(r)}{L[i.](r)} respectively describing the
  dependence between \eqn{X_i}{X[i]} and \eqn{X}{X}.
  See \code{\link{Gdot}}, \code{\link{Jdot}}, \code{\link{Kdot}}
  or \code{\link{Ldot}} respectively for explanation of these functions.
  The estimate is computed by applying the relevant function
  (\code{\link{Gdot}} etc)
  to \code{X} using each possible value of the argument \code{i},
  together with the optional arguments \code{\dots}.

  The letters \code{"G"}, \code{"J"}, \code{"K"} and \code{"L"}
  are interpreted as abbreviations for \code{\link{Gcross}},
  \code{\link{Jcross}}, \code{\link{Kcross}} and \code{\link{Lcross}}
  respectively, assuming the point pattern is
  marked. If the point pattern is unmarked, the appropriate
  function \code{\link{Fest}}, \code{\link{Jest}},
  \code{\link{Kest}} or \code{\link{Lest}} is invoked instead.

  If \code{envelope=TRUE}, then as well as computing the value of the
  summary function for each combination of types, the algorithm also
  computes simulation envelopes of the summary function for each
  combination of types. The arguments \code{\dots} are passed to the function
  \code{\link{envelope}} to control the number of
  simulations, the random process generating the simulations,
  the construction of envelopes, and so on. 

  When \code{envelope=TRUE} it is possible that errors could occur
  because the simulated point patterns do not satisfy the requirements
  of the summary function (for example, because the simulated pattern
  is empty and \code{fun} requires at least one point). If the number
  of such errors exceeds the maximum permitted number \code{maxnerr},
  then the envelope algorithm will give up, and will return
  the empirical summary function for the data point pattern,
  \code{fun(X)}, in place of the envelope.
}
\value{
  A function array (an object of class \code{"fasp"},
  see \code{\link{fasp.object}}). This can be plotted
  using \code{\link{plot.fasp}}.

  If the pattern is not marked, the resulting ``array'' has dimensions
  \eqn{1 \times 1}{1 x 1}. Otherwise the following is true:

  If \code{fun="F"},
  the function array has dimensions \eqn{m \times 1}{m * 1}
  where \eqn{m} is the number of different marks in the point pattern.
  The entry at position \code{[i,1]} in this array
  is the result of applying \code{\link{Fest}} to the
  points of type \code{i} only.

  If \code{fun} is \code{"Gdot"}, \code{"Jdot"}, \code{"Kdot"}
  or \code{"Ldot"}, the function array
  again has dimensions \eqn{m \times 1}{m * 1}.
  The entry at position \code{[i,1]} in this array
  is the result of \code{Gdot(X, i)}, \code{Jdot(X, i)}
  \code{Kdot(X, i)} or \code{Ldot(X, i)} respectively.

  If \code{fun} is \code{"Gcross"}, \code{"Jcross"}, \code{"Kcross"}
  or \code{"Lcross"} 
  (or their abbreviations \code{"G"}, \code{"J"}, \code{"K"} or \code{"L"}),
  the function array has dimensions \eqn{m \times m}{m * m}.
  The \code{[i,j]} entry of the function array
  (for \eqn{i \neq j}{i != j}) is the
  result of applying the function \code{\link{Gcross}},
  \code{\link{Jcross}}, \code{\link{Kcross}} or\code{\link{Lcross}} to
  the pair of types \code{(i,j)}. The diagonal
  \code{[i,i]} entry of the function array is the result of
  applying the univariate function \code{\link{Gest}},
  \code{\link{Jest}}, \code{\link{Kest}} or \code{\link{Lest}} to the
  points of type \code{i} only.

  If \code{envelope=FALSE}, then
  each function entry \code{fns[[i]]} retains the format
  of the output of the relevant estimating routine
  \code{\link{Fest}}, \code{\link{Gest}}, \code{\link{Jest}},
  \code{\link{Kest}},  \code{\link{Lest}}, \code{\link{Gcross}},
  \code{\link{Jcross}} ,\code{\link{Kcross}}, \code{\link{Lcross}},
  \code{\link{Gdot}}, \code{\link{Jdot}}, \code{\link{Kdot}} or
  \code{\link{Ldot}}
  The default formulae for plotting these functions are 
  \code{cbind(km,theo) ~ r} for F, G, and J functions, and
  \code{cbind(trans,theo) ~ r} for K and L functions.

  If \code{envelope=TRUE}, then each function entry \code{fns[[i]]}
  has the same format as the output of the \code{\link{envelope}} command.
}
\note{
  Sizeable amounts of memory may be needed during the calculation.
}

\seealso{
  \code{\link{plot.fasp}},
  \code{\link{fasp.object}},
  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Jest}},
  \code{\link{Kest}},
  \code{\link{Lest}},
  \code{\link{Gcross}},
  \code{\link{Jcross}},
  \code{\link{Kcross}},
  \code{\link{Lcross}},
  \code{\link{Gdot}},
  \code{\link{Jdot}},
  \code{\link{Kdot}},
  \code{\link{envelope}}.
}
\examples{
   # bramblecanes (3 marks).
   bram <- bramblecanes
   \testonly{
      bram <- bram[c(seq(1, 744, by=20), seq(745, 823, by=4))]
   }
   bF <- alltypes(bram,"F",verb=TRUE)
   plot(bF)
   if(interactive()) {
     plot(alltypes(bram,"G"))
     plot(alltypes(bram,"Gdot"))
   }
   
   # Swedishpines (unmarked).
  swed <- swedishpines
   \testonly{
     swed <- swed[1:25]
   }
   plot(alltypes(swed,"K"))

   plot(alltypes(amacrine, "pcf"), ylim=c(0,1.3))

   # A setting where you might REALLY want to use dataname:
   \dontrun{
   xxx <- alltypes(ppp(Melvin$x,Melvin$y,
                window=as.owin(c(5,20,15,50)),marks=clyde),
                fun="F",verb=TRUE,dataname="Melvin")
   }

   # envelopes
   bKE <- alltypes(bram,"K",envelope=TRUE,nsim=19)
   \dontrun{
   bFE <- alltypes(bram,"F",envelope=TRUE,nsim=19,global=TRUE)
   }

   # extract one entry
   as.fv(bKE[1,1])
   
}
\author{\adrian
  and \rolf.
}
\keyword{spatial}
\keyword{nonparametric}
