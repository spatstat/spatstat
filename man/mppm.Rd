\name{mppm}
\alias{mppm}
\title{Fit Point Process Model to Several Point Patterns}
\description{
  Fits a Gibbs point process model to several point patterns
  simultaneously. 
}
\usage{
   mppm(formula, data, interaction=Poisson(), ...,
            iformula=NULL,
%#%^!ifdef RANDOMEFFECTS
            random=NULL,
%#%^!endif
        use.gam = FALSE)
}
\arguments{
  \item{formula}{
    A formula describing the systematic part of the model.
    Variables in the formula are names of columns in \code{data}.
  }
  \item{data}{
    A hyperframe (object of class \code{"hyperframe"},
    see \code{\link{hyperframe}}) containing the
    point pattern responses and the explanatory variables. 
  }
  \item{interaction}{
    Interpoint interaction(s) appearing in the model.
    Either an object of class \code{"interact"}
    describing the point process interaction
    structure, or a hyperframe (with the same number of
    rows as \code{data}) whose entries are objects of class
    \code{"interact"}.     
  }
  \item{\dots}{Arguments passed to \code{\link{ppm}} controlling
    the fitting procedure.
  }
  \item{iformula}{
    Optional. A formula (with no left hand side)
    describing the interaction to be applied to each case.
    Each variable name in the formula should either be the name of a column
    in the hyperframe \code{interaction}, or the name of a column
    in the hyperframe \code{data} that is a vector or factor.
  }
%#%^!ifdef RANDOMEFFECTS  
  \item{random}{
    Optional. A formula (with no left hand side)
    describing a random effect. Variable names in the formula
    may be any of the column names of \code{data} and \code{interaction}.
    The formula must be recognisable to \code{\link{lme}}. 
  }
%#%^!endif  
  \item{use.gam}{Logical flag indicating whether to fit the model
    using \code{\link{gam}} or \code{\link{glm}}.
  }
}
\details{
  This function fits a common point process model to a dataset
  containing several different point patterns.
  
  It extends the capabilities of the function \code{\link{ppm}}
  to deal with data such as
  \itemize{
    \item replicated observations of spatial point patterns
    \item two groups of spatial point patterns
    \item a designed experiment in which the response from each unit
    is a point pattern.
  }
  
  The syntax of this function is similar to that of
  standard \R model-fitting functions like \code{\link{lm}} and
  \code{\link{glm}}. The first argument \code{formula}  is an \R formula
  describing the systematic part of the model. The second argument
  \code{data} contains the responses and the explanatory variables.
  Other arguments determine the stochastic structure of the model.

  Schematically,
  the data are regarded as the results of a designed experiment
  involving \eqn{n} experimental units.  Each unit has a
  \sQuote{response}, and optionally some \sQuote{explanatory variables}
  (covariates) describing the experimental conditions for that unit.
  In this context,
  \emph{the response from each unit is a point pattern}.
  The value of a particular covariate for each unit can be
  either a single value (numerical, logical or factor),
  or a spatial covariate.
  A \sQuote{spatial} covariate is a quantity that depends on spatial location,
  for example, the soil acidity or altitude at each location.
  For the purposes of \code{mppm}, a spatial covariate must be stored
  as a pixel image (object of class \code{"im"}) which gives the values
  of the covariate at a fine grid of locations.

  The argument \code{data} is a hyperframe (a generalisation of
  a data frame, see \code{\link{hyperframe}}). This is like a data frame
  except that the entries can be objects of any class.
  The hyperframe has one row for each experimental unit,
  and one column for each variable (response or explanatory variable).

  The \code{formula} should be an \R formula.
  The left hand side of \code{formula} determines the \sQuote{response}
  variable. This should be a single name, which
  should correspond to a column in \code{data}.

  The right hand side of \code{formula} determines the 
  spatial trend of the model. It specifies the linear predictor,
  and effectively represents the \bold{logarithm}
  of the spatial trend.
  Variables in the formula must be the names of columns of
  \code{data}, or one of the reserved names
  \describe{
    \item{x,y}{Cartesian coordinates of location}
    \item{marks}{Mark attached to point}
    \item{id}{which is a factor representing the
      serial number (\eqn{1} to \eqn{n}) of the point pattern,
      i.e. the row number in the data hyperframe.
    }
  }

  The column of responses in \code{data}
  must consist of point patterns (objects of class \code{"ppp"}).
  The individual point pattern responses 
  can be defined in different spatial windows.
  If some of the point patterns are marked, then they must all be
  marked, and must have the same type of marks.

  The scope of models that can be fitted to each pattern is the same as the
  scope of \code{\link{ppm}}, that is, Gibbs point processes with
  interaction terms that belong to a specified list, including
  for example the Poisson process, Strauss process, Geyer's saturation
  model, and piecewise constant pairwise interaction models.
%#%^!ifdef RANDOMEFFECTS  
  Additionally, it is possible to include random effects
  as explained in the section on Random Effects below.
%#%^!endif
  
  The stochastic part of the model is determined by
  the arguments \code{interaction} and (optionally) \code{iformula}.
  \itemize{
    \item 
    In the simplest case, \code{interaction} is 
    an object of class \code{"interact"},
    determining the interpoint interaction structure of the point
    process model, for all experimental units.
    \item
    Alternatively, \code{interaction} may be a hyperframe,
    whose entries are objects of class \code{"interact"}.
    It should have the same number of rows as \code{data}.
    \itemize{
      \item
      If \code{interaction} consists of only one column,
      then the entry in row \code{i} is taken to be the
      interpoint interaction for the \code{i}th experimental unit
      (corresponding to the \code{i}th row of \code{data}).
      \item
      If \code{interaction} has more than one column,
      then the argument \code{iformula} is also required.
      Each row of \code{interaction} determines
      several interpoint interaction structures that might be applied
      to the corresponding row of \code{data}.
      The choice of interaction is determined by \code{iformula};
      this should be an \R formula,
      without a left hand side.
      For example if \code{interaction} has two columns called
      \code{A} and \code{B} then \code{iformula = ~B} indicates that the
      interpoint interactions are taken from the second column.
    }
  }
  Variables in \code{iformula}
  typically refer to column names of \code{interaction}.
  They can also be names of columns in
  \code{data}, but only for columns of numeric, logical or factor
  values. For example \code{iformula = ~B * group} (where \code{group}
  is a column of \code{data} that contains a factor) causes the
  model with interpoint interaction \code{B} to be fitted
  with different interaction parameters for each level of \code{group}.
}
%#%^!ifdef RANDOMEFFECTS
\section{Random Effects}{
  It is also possible to include random effects in the
  trend term. The argument \code{random} is a formula,
  with no left-hand side, that specifies the structure of the
  random effects. The formula should be recognisable to
  \code{\link{lme}} (see the description of the argument \code{random}
  for \code{\link{lme}}).

  The names in the formula \code{random} may be any of the covariates
  supplied by \code{data}.
  Additionally the formula may involve the name
  \code{id}, which is a factor representing the
  serial number (\eqn{1} to \eqn{n}) of the point pattern in the
  list \code{X}.
}
%#%^!endif
\value{
  An object of class \code{"mppm"} representing the
  fitted model.

  There are methods for
  \code{print}, \code{summary} and \code{coef} for this class.
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Baddeley, A., Bischof, L., Sintorn, I.-M., Haggarty, S.,
  Bell, M. and Turner, R. 
  Analysis of a designed experiment where the response is a spatial
  point pattern. In preparation.

  Baddeley, A., Rubak, E. and Turner, R. (2015)
  \emph{Spatial Point Patterns: Methodology and Applications with R}.
  London: Chapman and Hall/CRC Press. 
}
\author{
  Adrian Baddeley, Ida-Maria Sintorn and Leanne Bischoff.
  Implemented by 
  Adrian Baddeley
  \email{Adrian.Baddeley@uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  Rolf Turner
  \email{r.turner@auckland.ac.nz}
  and Ege Rubak
  \email{rubak@math.aau.dk}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{print.mppm}},
  \code{\link{summary.mppm}},
  \code{\link{coef.mppm}},
}
\examples{
# Waterstriders data
 H <- hyperframe(Y = waterstriders)
 mppm(Y ~ 1,  data=H)
 mppm(Y ~ 1,  data=H, Strauss(7))
 mppm(Y ~ id, data=H)
 mppm(Y ~ x,  data=H)

# Synthetic data from known model
n <- 10
H <- hyperframe(V=1:n,
                U=runif(n, min=-1, max=1),
                M=factor(letters[1 + (1:n) \%\% 3]))
H$Z <- setcov(square(1))
H$U <- with(H, as.im(U, as.rectangle(Z)))
H$Y <- with(H, rpoispp(eval.im(exp(2+3*Z))))

fit <- mppm(Y ~Z + U + V, data=H)
}
\keyword{spatial}
\keyword{models}

