\name{model.depends}
\alias{model.depends}
\alias{model.is.additive}
\alias{model.covariates}
\alias{has.offset.term}
\alias{has.offset}
\title{
  Identify Covariates Involved in each Model Term
}
\description{
  Given a fitted model (of any kind), identify which of the covariates
  is involved in each term of the model.
}
\usage{
model.depends(object)
model.is.additive(object)
model.covariates(object, fitted=TRUE, offset=TRUE)
has.offset.term(object)
has.offset(object)
}
\arguments{
  \item{object}{
    A fitted model of any kind.
  }
  \item{fitted,offset}{
    Logical values determining which type of covariates to include.
  }
}
\details{
  The \code{object} can be a fitted model of any kind,
  including models of the classes \code{\link{lm}}, \code{\link{glm}}
  and \code{\link{ppm}}.

  To be precise,
  \code{object} must belong to a class for which there are methods
  for \code{\link{formula}}, \code{\link{terms}}
  and \code{\link{model.matrix}}.
  
  The command \code{model.depends} determines the relationship between
  the original covariates (the data supplied when \code{object} was
  fitted) and the canonical covariates (the columns of the design matrix).
  It returns a logical matrix, with one row for each canonical
  covariate, and one column for each of the original covariates,
  with the \code{i,j} entry equal to \code{TRUE} if the
  \code{i}th canonical covariate depends on the \code{j}th
  original covariate.

  If the model formula of \code{object} includes offset terms
  (see \code{\link{offset}}), then the return value of \code{model.depends}
  also has an attribute \code{"offset"}. This is a logical value or
  matrix with one row for each offset term and one column for each of
  the original covariates, with the \code{i,j} entry equal to \code{TRUE} if the
  \code{i}th offset term depends on the \code{j}th
  original covariate.

  The command \code{model.covariates} returns a character vector
  containing the names of all (original) covariates that were actually
  used to fit the model. By default, this includes all covariates that
  appear in the model formula, including offset terms as well as 
  canonical covariate terms. To omit the offset terms, set
  \code{offset=FALSE}. To omit the canonical covariate terms,
  set \code{fitted=FALSE}.

  The command \code{model.is.additive} determines whether the model
  is additive, in the sense that there is no canonical covariate that
  depends on two or more original covariates. It returns a logical value.

  The command \code{has.offset.term} is a faster way to determine whether the
  model \emph{formula} includes an \code{offset} term.

  The functions \code{model.depends} and \code{has.offset.term}
  only detect \code{offset} terms which are present
  in the model formula. They do not detect numerical offsets in the
  model object, that were inserted using the \code{offset} argument
  in \code{lm}, \code{glm} etc. To detect the presence of offsets
  of both kinds, use \code{has.offset}.
}
\value{
  A logical value or matrix.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{model.matrix}}
}
\examples{
   x <- 1:10
   y <- 3*x + 2
   z <- rep(c(-1,1), 5)
   fit <- lm(y ~ poly(x,2) + sin(z))
   model.depends(fit)
   model.covariates(fit)
   model.is.additive(fit)

   fitoff1 <- lm(y ~ x + offset(z))
   fitoff2 <- lm(y ~ x, offset=z)
   has.offset.term(fitoff1)
   has.offset(fitoff1)
   has.offset.term(fitoff2)
   has.offset(fitoff2)
}
\author{Adrian Baddeley \email{Adrian.Baddeley@curtin.edu.au}
  
  
  and Rolf Turner \email{r.turner@auckland.ac.nz}
  
}
\keyword{spatial}
\keyword{models}

