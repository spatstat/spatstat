#
# interactive analysis of point patterns
#
#   $Revision: 1.17 $   $Date: 2014/08/27 09:48:23 $
#
#

istat <- function(x, xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  verifyclass(x, "ppp")
  # generate simulations of CSR for use in envelopes
  simx <- envelope(x, fun=NULL, nsim=39,
                   internal=list(csr=TRUE, eject="patterns"))
  # initial value of smoothing parameter
  sigma0 <- with(x$window, min(diff(xrange),diff(yrange)))/8
  # create panel
  require(rpanel)
  p <- rpanel::rp.control(paste("istat(", xname, ")", sep=""),
                          x=x,           # point pattern
                          xname=xname,   # name of point pattern
                          simx=simx,   # simulated realisations of CSR
                          stat="data",
                          envel="none",
                          sigma=sigma0,
                          size=c(600, 400))
# Split panel into two halves  
# Left half of panel: display
# Right half of panel: controls
  rpanel::rp.grid(p, "gdisplay",
                  pos=list(row=0,column=0), width=400, height=400)
  rpanel::rp.grid(p, "gcontrols",
                  pos=list(row=0,column=1), width=200, height=400)

#----- Display side ------------

  # This line is to placate the package checker
  mytkr2 <- NULL
  
  rpanel::rp.tkrplot(p, mytkr2, do.istat,
                     pos=list(row=0,column=0,grid="gdisplay"))

  redraw <- function(panel) {
    rpanel::rp.tkrreplot(panel, mytkr2)
    panel
  }
  
#----- Control side ------------
  nextrow <- 0
  pozzie <- function(n=nextrow,s='w')
    list(row=n,column=0,grid="gcontrols",sticky=s)
  
# choice of summary statistic
  ftable <- c(data="data",
              density="kernel smoothed",
              Kest="K-function",
              Lest="L-function",
              pcf="pair correlation",
              Kinhom="inhomogeneous K",
              Linhom="inhomogeneous L",
              Fest="empty space function F",
              Gest="nearest neighbour function G",
              Jest="J-function")
  fvals <- names(ftable)
  flabs <- as.character(ftable)
  stat <- NULL
  rpanel::rp.radiogroup(p, stat, vals=fvals, labels=flabs,
                        title="statistic", action=redraw,
                        pos=pozzie(0))
  nextrow <- 1
# envelopes?
  envel <- NULL
  evals <- c("none", "pointwise", "simultaneous")
  elabs <- c("No simulation envelopes",
             "Pointwise envelopes under CSR",
             "Simultaneous envelopes under CSR")
  rpanel::rp.radiogroup(p, envel, vals=evals, labels=elabs,
                        title="Simulation envelopes", action=redraw,
                        pos=pozzie(nextrow))
  nextrow <- nextrow + 1
# smoothing parameters
  sigma <- NULL
  rect <- as.rectangle(x$window)
  winwid  <- min(abs(diff(rect$xrange)), abs(diff(rect$yrange)))
  rpanel::rp.slider(p, sigma, winwid/80, winwid/2, action=redraw, 
                    title="sigma",
                    initval=winwid/8, showvalue=TRUE, pos=pozzie(nextrow, ''))
  nextrow <- nextrow + 1
  pcfbw <- pcfbwinit <- 0.15/sqrt(5 * x$n/area.owin(x$window))
  rpanel::rp.slider(p, pcfbw, pcfbwinit/10, 4 * pcfbwinit, action=redraw, 
                    title="bw", initval=pcfbwinit,
                    showvalue=TRUE, pos=pozzie(nextrow, ''))
  nextrow <- nextrow + 1
# button to print a summary at console
  rpanel::rp.button(p, title="Print summary information",
                    action=function(panel) { print(summary(panel$x)); panel},
                    pos=pozzie(nextrow))
  nextrow <- nextrow + 1
# quit button 
  rpanel::rp.button(p, title="Quit", quitbutton=TRUE,
                    action= function(panel) { panel }, pos=pozzie(nextrow))

  invisible(NULL)
}

# function that updates the plot when the control panel is operated

do.istat <- function(panel) { 
  x     <- panel$x
  xname <- panel$xname
  envel <- panel$envel
  stat  <- panel$stat
  sigma <- panel$sigma
  simx  <- panel$simx
  if(stat=="data") {
    plot(x, main=xname)
    return(panel)
  }
  out <- 
    switch(envel,
           none=switch(stat,
             density=density(x, sigma=sigma),
             Kest=Kest(x),
             Lest=Lest(x), 
             pcf=pcf(x, bw=panel$pcfbw),
             Kinhom=Kinhom(x, sigma=sigma),
             Linhom=Linhom(x, sigma=sigma),
             Fest=Fest(x),
             Gest=Gest(x),
             Jest=Jest(x)),
           pointwise=switch(stat,
             density=density(x, sigma=sigma),
             Kest=envelope(x, Kest, nsim=39, simulate=simx),
             Lest=envelope(x, Lest, nsim=39, simulate=simx),
             pcf=envelope(x, pcf, bw=panel$pcfbw, nsim=39, simulate=simx),
             Kinhom=envelope(x, Kinhom, nsim=39, sigma=sigma, simulate=simx),
             Linhom=envelope(x, Linhom, nsim=39, sigma=sigma, simulate=simx),
             Fest=envelope(x, Fest, nsim=39, simulate=simx),
             Gest=envelope(x, Gest, nsim=39, simulate=simx),
             Jest=envelope(x, Jest, nsim=39, simulate=simx)),
           simultaneous=switch(stat,
             density=density(x, sigma=sigma),
             Kest=envelope(x, Kest, nsim=19, global=TRUE, simulate=simx),
             Lest=envelope(x, Lest, nsim=19, global=TRUE, simulate=simx),
             pcf=envelope(x, pcf, bw=panel$pcfbw, nsim=19, global=TRUE, simulate=simx),
             Kinhom=envelope(x, Kinhom, nsim=19, sigma=sigma, global=TRUE, simulate=simx),
             Linhom=envelope(x, Linhom, nsim=19, sigma=sigma, global=TRUE, simulate=simx),
             Fest=envelope(x, Fest, nsim=19, global=TRUE, simulate=simx),
             Gest=envelope(x, Gest, nsim=19, global=TRUE, simulate=simx),
             Jest=envelope(x, Jest, nsim=19, global=TRUE, simulate=simx))
           )
  # plot it
  if(stat %in% c("density", "Kinhom", "Linhom")) {
    plot(out, main=paste(stat, "(", xname, ", sigma)", sep=""))
    if(stat == "density")
      points(x)
  } else if(stat == "pcf")
    plot(out, main=paste("pcf(", xname, ", bw)", sep=""))
  else 
    plot(out, main=paste(stat, "(", xname, ")", sep=""))

  return(panel)
}

