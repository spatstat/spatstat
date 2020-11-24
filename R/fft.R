#'
#'           fft.R
#' 
#'    choose code for computing Discrete Fourier Transform
#' 
#'    $Revision: 1.1 $  $Date: 2020/11/24 01:10:13 $

fft2D <- function(z, inverse=FALSE, west=fftwAvailable()) {
  if(west) return(fftwtools::fftw2d(data=z, inverse=inverse))
  return(stats::fft(z=z, inverse=inverse))
}

fftwAvailable <- function() {
  ok <- requireNamespace("fftwtools", quietly=TRUE)
  return(ok)
}

