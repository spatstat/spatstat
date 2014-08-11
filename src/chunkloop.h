/*
  chunkloop.h

  Divide a loop into chunks 

  Convenient for divide-and-recombine,
  and reducing calls to R_CheckUserInterrupt, etc.

  $Revision: 1.1 $  $Date: 2012/03/27 01:34:42 $
  
*/

#define OUTERCHUNKLOOP(IVAR, LOOPLENGTH, ICHUNK, CHUNKSIZE) \
  IVAR = 0; \
  ICHUNK = 0; \
  while(IVAR < LOOPLENGTH) 

#define INNERCHUNKLOOP(IVAR, LOOPLENGTH, ICHUNK, CHUNKSIZE) \
    ICHUNK += CHUNKSIZE; \
    if(ICHUNK > LOOPLENGTH) ICHUNK = LOOPLENGTH; \
    for(; IVAR < ICHUNK; IVAR++) 

#define XOUTERCHUNKLOOP(IVAR, ISTART, IEND, ICHUNK, CHUNKSIZE) \
  IVAR = ISTART; \
  ICHUNK = 0; \
  while(IVAR <= IEND) 

#define XINNERCHUNKLOOP(IVAR, ISTART, IEND, ICHUNK, CHUNKSIZE)	\
    ICHUNK += CHUNKSIZE; \
    if(ICHUNK > IEND) ICHUNK = IEND; \
    for(; IVAR <= IEND; IVAR++) 

#define CHUNKLOOP_H




