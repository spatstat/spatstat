#
#  primefactors.R
#
#  $Revision: 1.4 $   $Date: 2012/04/18 09:16:32 $
#

primesbelow <- local({

  # all primes below 1000
  p1000 <- c(
   2,   3,   5,   7,  11,  13,  17,  19,  23,
  29,  31,  37,  41,  43,  47,  53,  59,  61,
  67,  71,  73,  79,  83,  89,  97, 101, 103,
 107, 109, 113, 127, 131, 137, 139, 149, 151,
 157, 163, 167, 173, 179, 181, 191, 193, 197,
 199, 211, 223, 227, 229, 233, 239, 241, 251,
 257, 263, 269, 271, 277, 281, 283, 293, 307,
 311, 313, 317, 331, 337, 347, 349, 353, 359,
 367, 373, 379, 383, 389, 397, 401, 409, 419,
 421, 431, 433, 439, 443, 449, 457, 461, 463,
 467, 479, 487, 491, 499, 503, 509, 521, 523,
 541, 547, 557, 563, 569, 571, 577, 587, 593,
 599, 601, 607, 613, 617, 619, 631, 641, 643,
 647, 653, 659, 661, 673, 677, 683, 691, 701,
 709, 719, 727, 733, 739, 743, 751, 757, 761,
 769, 773, 787, 797, 809, 811, 821, 823, 827,
 829, 839, 853, 857, 859, 863, 877, 881, 883,
 887, 907, 911, 919, 929, 937, 941, 947, 953,
 967, 971, 977, 983, 991, 997)

  primesbelow <- function(nmax) {
    if(nmax <= 1000) return(p1000[p1000 <=  nmax])
    eratosthenes(nmax, c(p1000, 1001:nmax))
  }
  primesbelow
})

eratosthenes <- function(nmax, startset=2:nmax) {
  # The Sieve of Eratosthenes
  if(nmax < 2) return(numeric(0))
  numbers <- startset
  prime <- 2
  repeat{
    retain <-  (numbers <= prime) | (numbers %% prime != 0)
    numbers <- numbers[retain]
    remaining <- (numbers > prime)
    if(!any(remaining))
      break
    prime <- min(numbers[remaining])
  }
  return(numbers)
}
  
primefactors <- function(n, prmax) {
  if(missing(prmax)) prmax <- floor(sqrt(n))
  primes <- primesbelow(prmax)
  divides.n <- (n %% primes == 0)
  if(!any(divides.n)) 
    return(n)
  else {
    divisors <- primes[divides.n]
    prmax <- max(divisors)
    m <- n/prod(divisors)
    if(m == 1) return(divisors)
    else {
      mfactors <- primefactors(m, prmax=prmax)
      return(sort(c(divisors, mfactors)))
    }
  }
}

is.prime <- function(n) { length(primefactors(n)) == 1 }

least.common.multiple <- function(n, m) {
  nf <- primefactors(n)
  mf <- primefactors(m)
  p <- sort(unique(c(nf,mf)))
  nfac <- table(factor(nf, levels=p))
  mfac <- table(factor(mf, levels=p))
  prod(p^pmax(nfac,mfac))
}

greatest.common.divisor <- function(n, m) {
  nf <- primefactors(n)
  mf <- primefactors(m)
  p <- sort(unique(c(nf,mf)))
  nfac <- table(factor(nf, levels=p))
  mfac <- table(factor(mf, levels=p))
  prod(p^pmin(nfac,mfac))
}
  
divisors <- function(n) {
  p <- primefactors(n)

  up <- sort(unique(p))
  k <- table(factor(p, levels=up))

  rest <- function(kk, uu) {
    powers <- uu[1]^(0:(kk[1]))
    if(length(uu) == 1)
      return(powers)
    rr <- rest(kk[-1], uu[-1])
    products <- as.vector(outer(powers, rr, "*"))
    return(sort(unique(products)))
    }

  return(rest(k, up))
}

    
