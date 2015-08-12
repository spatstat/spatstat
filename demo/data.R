if(dev.cur() <= 1) {
  dd <- getOption("device")
  if(is.character(dd)) dd <- get(dd)
  dd()
}

oldpar <- par(ask = interactive() && dev.interactive(orNone=TRUE))
oldoptions <- options(warn=-1)

plot(amacrine)

plot(anemones, markscale=1)

ants.extra$plotit()

plot(austates)

plot(bei.extra$elev, main="Beilschmiedia")
plot(bei, add=TRUE, pch=16, cex=0.3)

plot(betacells)

plot(bramblecanes, cols=1:3)
plot(split(bramblecanes))

plot(bronzefilter,markscale=2)

plot(cells)

plot(chicago, main="Chicago Street Crimes", col="grey",
     cols=c("red", "blue", "black", "blue", "red", "blue", "blue"),
     chars=c(16,2,22,17,24,15,6), leg.side="left", show.window=FALSE)

chorley.extra$plotit()

plot(clmfires, which.marks="cause", cols=2:5, cex=0.25,
     main="Castilla-La Mancha forest fires")
plot(clmfires.extra$clmcov200, main="Covariates for forest fires")

plot(copper$Points, main="Copper")
plot(copper$Lines, add=TRUE)

plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }),
      parargs=list(mar=rep(1,4)))

plot(dendrite, leg.side="bottom", main="", cex=0.75, cols=2:4)

plot(demopat)

plot(finpines, main="Finnish pines")

wildM1 <- with(flu, virustype == "wt" & stain == "M2-M1")
plot(flu[wildM1, 1, drop=TRUE],
     main=c("flu data", "wild type virus, M2-M1 stain"),
     chars=c(16,3), cex=0.4, cols=2:3)

plot(gordon, main="People in Gordon Square", pch=16)

plot(gorillas, which.marks=1, chars=c(1,3), cols=2:3, main="Gorilla nest sites")

plot(hamster, cols=c(2,4))

plot(heather)

plot(humberside)

plot(hyytiala, cols=2:5)

plot(japanesepines)

plot(lansing)
plot(split(lansing))

plot(longleaf)

plot(mucosa, chars=c(1,3), cols=c("red", "green"))
plot(mucosa.subwin, add=TRUE, lty=3)

plot(murchison, main="Murchison data")

plot(murchison$greenstone, main="Murchison data", col="lightgreen")
plot(murchison$gold, add=TRUE, pch=3, col="blue")
plot(murchison$faults, add=TRUE, col="red")

plot(nbfires, use.marks=FALSE, pch=".")
plot(split(nbfires), use.marks=FALSE, chars=".")
plot(split(nbfires)$"2000", which.marks="fire.type",
     main=c("New Brunswick fires 2000", "by fire type"),
     cols=c("blue", "green", "red", "cyan"),
     leg.side="left")

plot(nztrees)
plot(trim.rectangle(as.owin(nztrees), c(0,5), 0), add=TRUE, lty=3)

plot(osteo[1:10,], tick.marks=FALSE, xlab="", ylab="", zlab="")

plot(paracou, cols=2:3, chars=c(16,3))

ponderosa.extra$plotit()

pyr <- pyramidal
pyr$grp <- abbreviate(pyramidal$group, minlength=7)
plot(pyr, quote(plot(Neurons, pch=16, main=grp)), main="Pyramidal Neurons")
rm(pyr)

plot(redwood)

redwoodfull.extra$plotit()

plot(residualspaper$Fig1)
plot(residualspaper$Fig4a)
plot(residualspaper$Fig4b)
plot(residualspaper$Fig4c)

shapley.extra$plotit(main="Shapley")

plot(simdat)

plot(spiders, pch=16, show.window=FALSE)

plot(sporophores, chars=c(16,1,2), cex=0.6)
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)

plot(spruces, maxsize=min(nndist(spruces)))

plot(swedishpines)

plot(urkiola, cex=0.5, cols=2:3)

plot(waka, markscale=0.04, main=c("Waka national park", "tree diameters"))

plot(waterstriders)

par(oldpar)
options(oldoptions)
