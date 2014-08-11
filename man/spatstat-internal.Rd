\name{spatstat-internal}
\title{Internal spatstat functions}
\alias{[.hyperframe}
\alias{[<-.hyperframe}
\alias{$.hyperframe}
\alias{$<-.hyperframe}
\alias{[.localpcfmatrix}
\alias{[.splitppx}
\alias{[<-.splitppx}
\alias{acedist.show}
\alias{acedist.noshow}
\alias{adjustthinrange}
\alias{affinexy}
\alias{affinexypolygon}
\alias{anycrossing.psp}
\alias{apply23sum}
\alias{applytolayers}
\alias{area.xypolygon}       
\alias{areaGain.diri}       
\alias{areaGain.grid}       
\alias{areaLoss.diri}       
\alias{areaLoss.grid}       
\alias{assemble.plot.objects}       
\alias{AsymmDistance.psp}
\alias{as.breakpts}
\alias{as.data.frame.bw.optim}
\alias{as.data.frame.fv}
\alias{as.fv}
\alias{as.im.scan.test}
\alias{as.im.linim}
\alias{as.list.hyperframe}
\alias{as.listof}
\alias{as.units}
\alias{badprobability}
\alias{beachcolours}
\alias{bermantestEngine}
\alias{bdrylength.xypolygon} 
\alias{bdry.mask}
\alias{bind.ratfv}
\alias{blankcoefnames}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{bt.frame}
\alias{bw.optim}
\alias{cartesian}
\alias{cat.factor}
\alias{cellmiddles}
\alias{censtimeCDFest}
\alias{change.default.expand}          
\alias{checkfields}          
\alias{check.finite}
\alias{check.hist.lengths}
\alias{check.named.list}
\alias{check.named.thing}
\alias{check.named.vector}
\alias{check.nvector}
\alias{check.nmatrix}
\alias{check.1.real}
\alias{check.range}
\alias{check.testfun}
\alias{clarkevansCalc}
\alias{clear.simplepanel}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closepairs}
\alias{coef.summary.ppm}
\alias{compatible.rat}
\alias{complaining}
\alias{compileCDF}
\alias{compileK}
\alias{compilepcf}
\alias{conform.ratfv}
\alias{crosspairs}
\alias{crosspaircounts}
\alias{cobble.xy}
\alias{codetime}
\alias{commasep}
\alias{conform.imagelist}
\alias{countingweights}
\alias{damaged.ppm}
\alias{datagen.runifpointOnLines}
\alias{datagen.runifpoisppOnLines}
\alias{datagen.rpoisppOnLines}
\alias{default.clipwindow}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{dflt.redraw}
\alias{densitypointsEngine}
\alias{diagnose.ppm.engine}
\alias{digital.volume}
\alias{dilate.owin}
\alias{dim.fasp}               
\alias{dim.hyperframe}               
\alias{dim.im}               
\alias{dim.msr}
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{dimnames.msr}
\alias{distpl}               
\alias{distppl}
\alias{distppll}
\alias{distppllmin}
\alias{distributecbind}
\alias{dist2dpath}
\alias{divisors}
\alias{do.as.im}
\alias{do.call.matched}
\alias{do.istat}
\alias{dotexpr.to.call}
\alias{dropifsingle}
\alias{edge.Ripley}
\alias{edge.Trans}
\alias{edge.Trans.modif}
\alias{emptywindow}
\alias{ensure2vector}
\alias{envelopeEngine}
\alias{envelopeTest}
\alias{envelope.matrix}
\alias{equalpairs}          
\alias{equalpairs.quad}          
\alias{equals.quad}          
\alias{equalsfun.quad}          
\alias{eratosthenes}
\alias{erodemask}
\alias{erode.owin}
\alias{evalCovar}
\alias{evalCovar.ppm}
\alias{evalCovar.lppm}
\alias{evalCovariate}
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{even.breaks.owin}
\alias{evenly.spaced}
\alias{exactdt}              
\alias{exactPdt}
\alias{explain.ifnot}
\alias{extractAIC.slrm}
\alias{extractAtomicQtests}
\alias{fave.order}
\alias{f3engine}
\alias{f3Cengine}
\alias{fasp}
\alias{findbestlegendpos}
\alias{findCovariate}
\alias{findcbind}
\alias{firstfactor}
\alias{fii}
\alias{fillNA}
\alias{flipxypolygon}
\alias{forbidNA}
\alias{forbid.logi}
\alias{FormatFaspFormulae}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{fvnames}
\alias{fvnames<-}
\alias{g3engine}
\alias{g3Cengine}
\alias{gpcmethod}
\alias{greatest.common.divisor}
\alias{getdataname}
\alias{getfields}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getlambda.lpp}
\alias{getSumFun}
\alias{geyercounts}
\alias{GLMpredict}
\alias{good.names}
\alias{gridindex}            
\alias{grid1index}
\alias{grow.rectangle}
\alias{grow.simplepanel}
\alias{handle.r.b.args}
\alias{handle.rshift.args}
\alias{ho.engine}
\alias{hsvNA}
\alias{identical.formulae}
\alias{idorempty}
\alias{implemented.for.K}
\alias{inpoint}
\alias{inside.range}     
\alias{inside.triangle}     
\alias{inside.xypolygon}     
\alias{instantiate.interact}
\alias{intersect.ranges}     
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}
\alias{is.atomicQtest}
\alias{is.cadlag}
\alias{is.data}
\alias{is.expandable}
\alias{is.expandable.ppm}
\alias{is.expandable.rmhmodel}
\alias{is.fv}
\alias{is.hole.xypolygon}
\alias{is.hyperframe}
\alias{is.infline}
\alias{is.interact}
\alias{is.marked.default}    
\alias{is.marked.psp}    
\alias{is.marked.quad}    
\alias{is.multitype.quad}    
\alias{is.multitype.default}    
\alias{is.parseable}
\alias{is.pp3}
\alias{is.ppx}
\alias{is.prime}
\alias{is.psp}
\alias{is.slrm}
\alias{is.tess}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Knone.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{Kpcf.kppm}               
\alias{killinteraction}
\alias{km.rs.opt}
\alias{kppmMinCon}
\alias{kppmComLik}
\alias{labels.ppm}
\alias{layout.boxes}
\alias{least.common.multiple}
\alias{levels.im}
\alias{levels<-.im}
\alias{lhs.of.formula}
\alias{linequad}
\alias{linearKengine}
\alias{linearpcfengine}
\alias{listof}
\alias{localKengine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{local2lpp}
\alias{logi.engine}
\alias{passthrough}
\alias{paste.expr}
\alias{prettyinside}
\alias{lookup.im}
\alias{majorminorversion}
\alias{make.even.breaks}
\alias{make.parseable}
\alias{makeunits}
\alias{markappend}
\alias{markcbind}
\alias{markformat}
\alias{markformat.ppp}
\alias{markformat.ppx}
\alias{markformat.psp}
\alias{markformat.default}
\alias{mark.scale.default}
\alias{markspace.integral}
\alias{marks.default}           
\alias{marks.quad}           
\alias{\%mapp\%}
\alias{markappendop}
\alias{marksubset}
\alias{markreplicateop}
\alias{\%mrep\%}
\alias{marksubsetop}
\alias{\%msub\%}
\alias{mask2df}
\alias{matcolall}
\alias{matcolany}
\alias{matcolsum}
\alias{matrixsample}
\alias{matrowall}
\alias{matrowany}
\alias{matrowsum}
\alias{maxflow}
\alias{meanlistfv}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{model.se.image}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{MultiPair.checkmatrix}
\alias{multiply.only.finite.entries}
\alias{na.handle.im}
\alias{names.hyperframe}
\alias{names<-.hyperframe}
\alias{nearest.pixel}
\alias{nearest.valid.pixel}
\alias{newstyle.coeff.handling}
\alias{niceround}
\alias{nncleanEngine}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{nobjects}
\alias{nobjects.ppp}
\alias{nobjects.ppx}
\alias{nobjects.psp}
\alias{n.quad}
\alias{numalign}
\alias{numeric.columns}
\alias{nzpaste}
\alias{offsetsinformula}
\alias{onecolumn}
\alias{optimStatus}
\alias{ordinal}
\alias{outdated.interact}
\alias{overlap.owin}
\alias{overlap.trapezium}
\alias{overlap.xypolygon}    
\alias{oversize.quad}    
\alias{owinpolycheck}
\alias{owinpoly2mask}
\alias{pairs.listof}
\alias{param.quad}
\alias{paren}
\alias{partialModelMatrix}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{pickoption}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.addvar}
\alias{plot.barplotdata}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.fii}
\alias{plot.localpcfmatrix}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.plotpairsim}
\alias{plot.pppmatching}
\alias{plot.profilepl}
\alias{plot.qqppm}
\alias{plot.quadratcount}
\alias{plot.quadrattest}
\alias{plot.scan.test}
\alias{polynom}
\alias{ppllengine}
\alias{ppmCovariates}
\alias{ppm.influence}
\alias{pppdist.mat}
\alias{pppdist.prohorov}
\alias{ppsubset}
\alias{prange}
\alias{prefixfv}
\alias{primefactors}       
\alias{primesbelow}
\alias{printStatus}
\alias{print.addvar}    
\alias{print.bt.frame}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.colourmap}
\alias{print.diagppm}
\alias{print.distfun}
\alias{print.envelope}
\alias{print.fasp}       
\alias{print.fii}       
\alias{print.fv}       
\alias{print.funxy}       
\alias{print.hyperframe}
\alias{print.interact}       
\alias{print.isf}
\alias{print.layered}
\alias{print.linim}
\alias{print.localpcfmatrix}
\alias{print.logippm}
\alias{print.lut}
\alias{print.minconfit}
\alias{print.msr}
\alias{print.nnfun}
\alias{print.parres}
\alias{print.plotpairsim}
\alias{print.plotppm}
\alias{print.pppmatching}
\alias{print.profilepl}
\alias{print.quadrattest}
\alias{print.qqppm}
\alias{print.rat}
\alias{print.rmhcontrol}
\alias{print.rmhexpand}
\alias{print.rmhmodel}
\alias{print.rmhstart}
\alias{print.rmhInfoList}
\alias{print.splitppp}
\alias{print.splitppx}
\alias{print.summary.fii}
\alias{print.summary.hyperframe}
\alias{print.summary.listof}
\alias{print.summary.lut}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{print.summary.rmhexpand}
\alias{print.summary.splitppp}
\alias{print.summary.splitppx}
\alias{print.summary.units}
\alias{print.tess}
\alias{print.timed}
\alias{prolongseq}
\alias{quad}
\alias{RandomFieldsSafe}
\alias{ratfv}
\alias{rectquadrat.breaks}
\alias{rectquadrat.countEngine}
\alias{redraw.simplepanel}
\alias{repair.image.xycoords}
\alias{rgbNA}
\alias{rhs.of.formula}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{RmhExpandRule}
\alias{rmhsnoop}
\alias{run.simplepanel}
\alias{quadrat.testEngine}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{pointgrid}
\alias{rastersample}
\alias{rasterx.im}
\alias{rastery.im}
\alias{rasterxy.im}
\alias{reach.fii}
\alias{rebadge.fv}
\alias{rebound}
\alias{rebound.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{reconcile.fv}
\alias{repair.old.factor.image}
\alias{reincarnate.interact}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.defaults}
\alias{resolve.1.default}
\alias{resolve.2D.kernel}
\alias{restrict.mask}
\alias{reverse.xypolygon}    
\alias{revcumsum}
\alias{rmax.rule}
\alias{rotxy}
\alias{rotxypolygon}
\alias{row.names.hyperframe}
\alias{row.names<-.hyperframe}
\alias{runifpoispp}          
\alias{runifpoisppOnLines}          
\alias{runifrect}
\alias{rmhResolveControl}
\alias{rmhResolveExpansion}
\alias{rmhResolveTypes}
\alias{rmhSnoopEnv}
\alias{rmhcontrol.rmhcontrol}
\alias{rmhcontrol.list}
\alias{rmhEngine}
\alias{rmhmodel.rmhmodel}
\alias{rmhstart.rmhstart}
\alias{rmhstart.list}
\alias{rmpoint.I.allim}
\alias{rpoint.multi}
\alias{safelookup}
\alias{samefunction}
\alias{scanmeasure}
\alias{scanmeasure.ppp}
\alias{scanmeasure.im}
\alias{scanBinomLRTS}
\alias{scanPoisLRTS}
\alias{scanLRTS}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{sensiblevarname}
\alias{sewpcf}
\alias{sewsmod}
\alias{shiftxy}              
\alias{shiftxypolygon}              
\alias{short.deparse}
\alias{signalStatus}
\alias{simplepanel}
\alias{simplify.xypolygon}              
\alias{simulrecipe}              
\alias{singlestring}              
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{smoothpointsEngine}              
\alias{sort.im}
\alias{spatstat.rawdata.location}
\alias{spatstatClusterModelInfo}
\alias{spatstatRmhInfo}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{sphere.volume}
\alias{store.versionstring.spatstat}
\alias{str.hyperframe}
\alias{stratrand}
\alias{strausscounts}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summarise.trend}
\alias{summary.envelope}
\alias{summary.fii}
\alias{summary.hyperframe}
\alias{summary.logippm}
\alias{summary.lut}
\alias{summary.profilepl}
\alias{summary.pppmatching}
\alias{summary.ppx}
\alias{summary.splitppx}
\alias{summary.rmhexpand}
\alias{sumsymouter}
\alias{superimposeMarks}
\alias{sympoly}
\alias{termsinformula}
\alias{test.crossing.psp}
\alias{test.selfcrossing.psp}
\alias{tilecentroids}        
\alias{tilenames}        
\alias{trap.extra.arguments}
\alias{trianglediameters}
\alias{trim.mask}        
\alias{tweak.fv.entry}
\alias{unitname.default}
\alias{unitname<-.default}
\alias{update.interact}
\alias{update.rmhstart}
\alias{update.slrm}
\alias{validradius}        
\alias{validate.mask}        
\alias{validate.quad}        
\alias{validposint}        
\alias{vanilla.fv}        
\alias{variablesinformula}
\alias{verbalogic}
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{versionstring.spatstat}
\alias{verifyclass}
\alias{verify.xypolygon}
\alias{warn.ignored.args}
\alias{with.msr}
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xy.grid}
\alias{X2testEngine}
\alias{xtfrm.im}
\alias{xypolygon2psp}
\alias{xypolyselfint}
\description{
  Internal spatstat functions.
}
\usage{
\method{[}{hyperframe}(x, i, j, drop=FALSE, ...)
\method{[}{hyperframe}(x, i, j) <- value
\method{$}{hyperframe}(x, name)
\method{$}{hyperframe}(x, i) <- value
\method{[}{splitppx}(x, ...)
\method{[}{splitppx}(x, ...) <- value
acedist.show(X, Y, n, d, timelag = 0)
acedist.noshow(X, Y, n, d)
adjustthinrange(ur,vstep,vr)
affinexy(X, mat = diag(c(1, 1)), vec = c(0, 0), invert=FALSE)
affinexypolygon(p, mat, vec, detmat)
anycrossing.psp(A,B)
apply23sum(x)
applytolayers(L, FUN, ...)
area.xypolygon(polly)
areaGain.diri(u, X, r, ..., W=as.owin(X))
areaGain.grid(u, X, r, ..., W=NULL, ngrid=spatstat.options("ngrid.disc"))
areaLoss.diri(X, r, ..., W=as.owin(X), subset=NULL)
areaLoss.grid(X, r, ..., W=as.owin(X), subset=NULL,
                         method = c("count", "distmap"),
                         ngrid = spatstat.options("ngrid.disc"),
                         exact = FALSE)
assemble.plot.objects(xlim, ylim, ..., lines, polygon)
AsymmDistance.psp(X, Y, metric, method)
as.breakpts(\dots)
\method{as.data.frame}{fv}(x, ...)
\method{as.data.frame}{bw.optim}(x, ...)
as.fv(x)
\method{as.list}{hyperframe}(x, ...)
\method{as.im}{scan.test}(X, \dots)
\method{as.im}{linim}(X, \dots)
as.listof(x)
as.units(s)
badprobability(x, NAvalue)
beachcolours(heightrange, sealevel, monochrome, ncolours)
bermantestEngine(model, covariate, which, alternative, ..., modelname, covname, dataname)
bdrylength.xypolygon(polly)
bdry.mask(W)
bind.ratfv(x, numerator, denominator, labl, desc, preferred, ratio)
blankcoefnames(x)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
bt.frame(Q, trend=~1, interaction=NULL, ..., covariates=NULL,
         correction="border", rbord=0, use.gam=FALSE, allcovar=FALSE)
bw.optim(cv, h, iopt, \dots, cvname, hname, criterion)
cartesian(pp, markset, fac = TRUE)
cat.factor(\dots, recursive=FALSE)
cellmiddles(W, nx, ny, npix, gi)
censtimeCDFest(o, cc, d, breaks, \dots, KM, RS, HAN, RAW, han.denom,
pmax)
change.default.expand(x, newdefault)
checkfields(X,L)          
check.finite(x, context, xname, fatal, usergiven)
check.hist.lengths(hist,breaks)
check.named.list(x, nam, context, namopt)
check.named.vector(x, nam, context, namopt)
check.named.thing(x, nam, namopt, xtitle, valid, type, context, fatal)
check.nvector(v, npoints, fatal=TRUE, things="data points", naok=FALSE)
check.nmatrix(m, npoints, fatal=TRUE, things="data points", naok=FALSE, squarematrix=TRUE)
check.1.real(x, context, fatal)
check.range(x, fatal)
check.testfun(f, f1, X)
clarkevansCalc(X, correction, clipregion, working)
clear.simplepanel(P)
clip.psp(x, window, check=TRUE)
cliprect.psp(x, window)
clippoly.psp(s, window)
closepairs(X,rmax,ordered,what)
\method{coef}{summary.ppm}(object, \dots)
\method{compatible}{rat}(A, B, \dots)
complaining(whinge, fatal, value)
compileCDF(D, B, r, \dots, han.denom, check)
compileK(D, r, weights, denom, check, ratio)
compilepcf(D, r, weights, denom, check, endcorrect, \dots)
conform.ratfv(x)
crosspairs(X,Y,rmax,what)
crosspaircounts(X,Y,r)
cobble.xy(x, y, f, fatal, \dots)
codetime(x, hms)
commasep(x, join)
conform.imagelist(X, Zlist)
countingweights(id, areas, check = TRUE)
damaged.ppm(object)
datagen.runifpointOnLines(n, L)
datagen.runifpoisppOnLines(lambda, L)
datagen.rpoisppOnLines(lambda, L, lmax, \dots, check=TRUE)
default.clipwindow(object, epsilon)
default.n.tiling(X, nd, ntile, npix, eps, verbose)
default.ntile(X)
dflt.redraw(button, name, env)
densitypointsEngine(x, sigma, \dots,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted)
diagnose.ppm.engine(object, \dots, type="eem", typename, opt,
                         sigma=NULL, rbord = reach(object), compute.sd=TRUE,
                         compute.cts=TRUE, rv=NULL, oldstyle=FALSE)
digital.volume(range, nval, vside)
dilate.owin(\dots)
\method{dim}{fasp}(x)
\method{dim}{hyperframe}(x)
\method{dim}{im}(x)
\method{dim}{msr}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
\method{dimnames}{msr}(x)
distpl(p, l)               
distppl(p, l)
distppll(p, l, mintype, method, listit)
distppllmin(p, l, big)
distributecbind(x)
dist2dpath(dist, method="C")
divisors(n)
do.as.im(x, action, \dots, W, eps, dimyx, xy, na.replace)
do.call.matched(fun, arglist, funargs, extrargs, sieve)
do.istat(panel)
dotexpr.to.call(expr, dot, evaluator)
dropifsingle(x)
edge.Ripley(X, r, W, method, maxweight)
edge.Trans(X, Y, W, exact, paired, trim)
edge.Trans.modif(X, Y, WX, WY, exact, paired, trim)
emptywindow(w)
ensure2vector(x)
envelopeEngine(X, fun, simul,
           nsim=99, nrank=1, \dots, 
           verbose=TRUE, clipdata=TRUE, 
           transform=NULL, global=FALSE, ginterval=NULL,
           savefuns=FALSE, savepatterns=FALSE, nsim2=nsim,
           VARIANCE=FALSE, nSD=2,
           Yname=NULL, maxnerr=nsim, internal=NULL, cl=NULL,
           envir.user=envir.user, expected.arg="r", do.pwrong=FALSE)
envelopeTest(X, \dots,
            power, rinterval, use.theo, tie.rule, save.envelope,
            Xname, verbose, internal)
\method{envelope}{matrix}(Y, ..., rvals, observed, theory, funX,
      nsim, nsim2, jsim, jsim.mean,
      type, csr, use.theory, nrank, ginterval, nSD,
      savefuns, check, Yname, do.pwrong)
equalpairs(U, X, marked=FALSE)
equalpairs.quad(Q)
equals.quad(Q)          
equalsfun.quad(Q)          
eratosthenes(nmax, startset)
erodemask(w,r,strict)
erode.owin(\dots)
evalCovar(model, covariate, \dots)
\method{evalCovar}{ppm}(model, covariate, \dots,
          dimyx, eps, jitter, modelname, covname, dataname)
\method{evalCovar}{lppm}(model, covariate, \dots,
          eps, nd, jitter, modelname, covname, dataname)
evalCovariate(covariate, locations)
evalInteraction(X,P,E,interaction,correction,\dots,precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,\dots,precomputed,savecomputed)
even.breaks.owin(w)
evenly.spaced(x, tol)
exactdt(X, \dots)              
exactPdt(w)
explain.ifnot(expr, context)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, ...)
extractAtomicQtests(x)
fave.order(x)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
findbestlegendpos(...)
findCovariate(covname, scope, scopename=NULL)
findcbind(root, depth, maxdepth)
firstfactor(x)
fii(interaction, coefs, Vnames, IsOffset)
fillNA(x, value)
flipxypolygon(p)
forbidNA(x, context, xname, fatal, usergiven)
forbid.logi(object)
FormatFaspFormulae(f, argname)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
fvnames(X, a)
fvnames(X, a) <- value
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
gpcmethod(fname, signature)
greatest.common.divisor(n,m)
getdataname(defaultvalue, \dots, dataname)
getfields(X, L, fatal = TRUE)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getlambda.lpp(lambda, X, \dots)
getSumFun(abbreviation, classname, ismarked, fatal)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
GLMpredict(fit, data, coefs, changecoef)
good.names(nama, defaults, suffices)
gridindex(x, y, xrange, yrange, nx, ny)            
grid1index(x, xrange, nx)
grow.rectangle(W, xmargin=0, ymargin=xmargin)
grow.simplepanel(P, side, len, new.clicks, new.redraws, \dots, aspect)
handle.r.b.args(r = NULL, breaks = NULL, window, eps = NULL, rmaxdefault)
handle.rshift.args(W, \dots, radius, width, height, edge, clip, edgedefault)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
hsvNA(h, s, v, \dots)
identical.formulae(x,y)
idorempty(w, r, caller)
implemented.for.K(correction, windowtype, explicit)
inpoint(W)
inside.range(x, r)
inside.triangle(x, y, xx, yy)
inside.xypolygon(pts, polly, test01, method)
instantiate.interact(x, par)
intersect.ranges(a,b,fatal)
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
is.atomicQtest(x)
is.cadlag(s)
is.data(Q)
is.expandable(x)
\method{is.expandable}{ppm}(x)
\method{is.expandable}{rmhmodel}(x)
is.fv(x)
is.hole.xypolygon(polly)
is.hyperframe(x)
is.infline(x)
is.interact(x)
\method{is.marked}{default}(\dots)  
\method{is.marked}{psp}(X, \dots)
\method{is.marked}{quad}(X, na.action="warn", \dots)
\method{is.multitype}{default}(X, \dots)  
\method{is.multitype}{quad}(X, na.action="warn", \dots)
is.parseable(x)
is.pp3(x)
is.ppx(x)
is.prime(n)
is.psp(x)
is.slrm(x)
is.tess(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks)
Kpcf.kppm(model, what)
killinteraction(model)
km.rs.opt(o, cc, d, breaks, KM, RS)
kppmMinCon(X, Xname, po, clusters, statistic, statargs, \dots)
kppmComLik(X, Xname, po, clusters, control, weightfun, rmax, \dots)
\method{labels}{ppm}(object, \dots)
layout.boxes(B, n, horizontal, aspect, usefrac)
least.common.multiple(n,m)
\method{levels}{im}(x)
\method{levels}{im}(x) <- value
lhs.of.formula(x)
linequad(X, Y, ..., eps, nd)
linearKengine(X, \dots, r, reweight, denom, correction, showworking)
linearpcfengine(X, \dots, r, reweight, denom, correction)
listof(\dots)
localKengine(X, \dots, wantL, lambda, correction, verbose, rvalue)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
local2lpp(L, seg, tp, X)
logi.engine(\dots)
passthrough(.Fun, \dots, .Fname)
paste.expr(x)
prettyinside(x, \dots)
\method{print}{localpcfmatrix}(x, \dots)
\method{print}{logippm}(x, \dots)
\method{plot}{localpcfmatrix}(x, \dots)
\method{[}{localpcfmatrix}(x, i, ...)
lookup.im(Z, x, y, naok, strict)
majorminorversion(v)
make.even.breaks(bmax, npos, bstep)
make.parseable(x)
makeunits(sing, plur, mul)
markappend(\dots)
markcbind(\dots)
markformat(x)
\method{markformat}{ppp}(x)
\method{markformat}{ppx}(x)
\method{markformat}{psp}(x)
\method{markformat}{default}(x)
mark.scale.default(marx, w, markscale=NULL, maxsize=NULL)
markspace.integral(X)
\method{marks}{default}(x, \dots)
\method{marks}{quad}(x, dfok=FALSE, \dots)
markappendop(x, y)
marksubset(x, index, format)
marksubsetop(x, i)
markreplicateop(x, n)
mask2df(w)
matcolall(x)
matcolany(x)
matcolsum(x)            
matrixsample(mat, newdim, phase, scale, na.value)
matrowall(x)
matrowany(x)
matrowsum(x)
maxflow(costm)
meanlistfv(z)
meanX.owin(w)            
meanY.owin(w)
model.se.image(fit, W, \dots, what)
mpl.engine(Q, trend, interaction, \dots, covariates, covfunargs, correction,
	 rbord, use.gam, gcontrol, famille,
         forcefit, nd, eps, allcovar, callstring,
         precomputed, savecomputed, preponly,
         rename.intercept, justQ)
mpl.get.covariates(covariates, locations, type, covfunargs)
mpl.prepare(Q, X, P, trend, interaction, covariates, 
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots,
            covfunargs, allcovar, precomputed, savecomputed,
            vnamebase, vnameprefix, warn.illegal, warn.unidentifiable)
MultiPair.checkmatrix(mat, n, matname, naok, zerook)
multiply.only.finite.entries(x, a)
na.handle.im(X, na.replace)
\method{names}{hyperframe}(x)
\method{names}{hyperframe}(x) <- value
nearest.pixel(x, y, im)
nearest.valid.pixel(x, y, im)
newstyle.coeff.handling(object)
niceround(x, m)
nncleanEngine(kthNND, k, d, \dots, tol, plothist, verbose, maxit)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
nobjects(x)
\method{nobjects}{ppp}(x)
\method{nobjects}{ppx}(x)
\method{nobjects}{psp}(x)
n.quad(Q)
numalign(i, nmax, zero)
numeric.columns(M, logical, others)
nzpaste(\dots, sep, collapse)
offsetsinformula(x)
onecolumn(m)
optimStatus(x, call)
printStatus(x, errors.only)
redraw.simplepanel(P, verbose)
run.simplepanel(P, verbose)
signalStatus(x, errors.only)
simplepanel(title, B, boxes, clicks, redraws, exit, env)
ordinal(k)
outdated.interact(object)
overlap.owin(A,B)
overlap.trapezium(xa, ya, xb, yb, verb = FALSE)
overlap.xypolygon(P, Q)
oversize.quad(Q, \dots, nU, nX)
owinpolycheck(W, verbose=TRUE)
owinpoly2mask(w, rasta, check=TRUE)
\method{pairs}{listof}(\dots, plot=TRUE)
param.quad(Q)
paren(x, type)
partialModelMatrix(X,D,model,callstring,\dots)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, ...,
               r = NULL, breaks = NULL,
               kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
               correction = c("translate", "Ripley"),
               sigma = NULL, varcov = NULL,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J")
pickoption(what="option", key, keymap, \dots,
           exact=FALSE, list.on.err=TRUE, die=TRUE, multi=FALSE)
ploterodewin(W1, W2, col.edge, col.inside, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{barplotdata}(x, \dots)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{fii}(x, \dots)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{pppmatching}(x, addmatch = NULL, main = NULL, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{profilepl}(x, \dots, add=FALSE, main=NULL, tag=TRUE, coeff=NULL, xvariable=NULL)
\method{plot}{qqppm}(x, \dots, limits=TRUE, monochrome=FALSE,
           limcol=if(monochrome) "black" else "red")
\method{plot}{quadratcount}(x, \dots, add, entries, dx, dy, show.tiles)
\method{plot}{quadrattest}(x, \dots)
\method{plot}{scan.test}(x, \dots, do.window)
polynom(x, \dots)
ppllengine(X, Y, action="project", check=FALSE)
ppmCovariates(model)
ppm.influence(fit,
   what = c("leverage","influence","dfbetas","derivatives"),
   ..., iScore = NULL, iHessian = NULL, iArgs=list(), drop = FALSE,
   method=c("C", "interpreted"))
pppdist.mat(X, Y, cutoff = 1, q = 1, matching = TRUE, precision = 9, approximation = 10)
pppdist.prohorov(X, Y, n, dfix, type, cutoff = 1, matching = TRUE, ccode = TRUE, precision = 9, approximation = 10) 
ppsubset(X, I)
prange(x)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
primefactors(n, prmax)
primesbelow(nmax)
\method{print}{addvar}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{colourmap}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{distfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{fii}(x, \dots)
\method{print}{funxy}(x, \dots)
\method{print}{fv}(x, \dots)
\method{print}{hyperframe}(x, ...)
\method{print}{interact}(x, \dots, family=TRUE, brief=FALSE)       
\method{print}{isf}(x, \dots)
\method{print}{layered}(x, \dots)
\method{print}{linim}(x, \dots)
\method{print}{lut}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{nnfun}(x, \dots)
\method{print}{parres}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{pppmatching}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rat}(x, \dots)
\method{print}{rmhcontrol}(x, \dots)
\method{print}{rmhexpand}(x, \dots, prefix=TRUE)
\method{print}{rmhmodel}(x, \dots)
\method{print}{rmhstart}(x, \dots)
\method{print}{rmhInfoList}(x, \dots)
\method{print}{splitppp}(x, \dots)
\method{print}{splitppx}(x, \dots)
\method{print}{summary.fii}(x, \dots)
\method{print}{summary.hyperframe}(x, ...)
\method{print}{summary.listof}(x, \dots)
\method{print}{summary.lut}(x, \dots)
\method{print}{summary.owin}(x, \dots)
\method{print}{summary.ppp}(x, \dots, dp=3)
\method{print}{summary.psp}(x, \dots)
\method{print}{summary.rmhexpand}(x, \dots)
\method{print}{summary.splitppp}(x, \dots)
\method{print}{summary.splitppx}(x, \dots)
\method{print}{summary.units}(x, \dots)
\method{print}{tess}(x, \dots, brief=FALSE)
\method{print}{timed}(x, \dots)
prolongseq(x, newrange)
quad(data, dummy, w, param)
RandomFieldsSafe()
ratfv(df, numer, denom, \dots, ratio)
rectquadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
rectquadrat.countEngine(x, y, xbreaks, ybreaks, weights)
repair.image.xycoords(x)
rgbNA(red, green, blue, \dots)
rhs.of.formula(x)
rhohatEngine(model, covariate, reference, volume, \dots,
               method, smoother, resolution, 
               n, bw, adjust, from, to, 
               bwref, covname, covunits, confidence, modelcall, callstring)
rhohatCalc(ZX, Zvalues, lambda, baseline, \dots,
           method, smoother,
           n, bw, adjust, from, to, 
           bwref, covname, confidence,
           covunits, modelcall, callstring, savestuff)
RmhExpandRule(nama)
rmhsnoop(\dots, Wsim, Wclip, R, xcoords, ycoords, mlevels, mcodes, irep, itype, 
     proptype, proplocn, propmark, propindx, numerator, denominator)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, Xname, fitname)
quadscheme.replicated(data, dummy, method = "grid", \dots)
quadscheme.spatial(data, dummy, method = "grid", \dots)
pointgrid(W, ngrid)
rastersample(X, Y)
rasterx.im(x)
rastery.im(x)
rasterxy.im(x, drop)
\method{reach}{fii}(x, \dots, epsilon)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
new.dotnames, new.preferred, new.formula, new.tags)
rebound(x, rect)
\method{rebound}{im}(x, rect)
\method{rebound}{ppp}(x, rect)
\method{rebound}{psp}(x, rect)
\method{rebound}{owin}(x, rect)
reconcile.fv(\dots)
repair.old.factor.image(x)
reincarnate.interact(object)
resid4plot(RES, plot.neg="image", plot.smooth="imagecontour",
           spacing=0.1, srange=NULL,monochrome=FALSE, main=NULL, \dots)
resid1plot(RES, opt, plot.neg="image", plot.smooth="imagecontour",
              srange=NULL, monochrome=FALSE, main=NULL, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD, xlab,ylab, \dots)
resolve.defaults(\dots, .MatchNull=TRUE, .StripNull=FALSE)
resolve.1.default(.A, \dots)
resolve.2D.kernel(\dots, sigma, varcov, x, mindist, adjust, bwfun)
restrict.mask(M, W)
reverse.xypolygon(p, adjust=FALSE)
revcumsum(x)
rmax.rule(fun, W, lambda)
rotxy(X, angle = pi/2)
rotxypolygon(p, angle = pi/2)
rmhResolveControl(control, model)
rmhResolveExpansion(win, control, imagelist, itype)
rmhResolveTypes(model, start, control)
rmhSnoopEnv(Xinit, Wclip, R)
\method{rmhcontrol}{rmhcontrol}(\dots)
\method{rmhcontrol}{list}(\dots)
rmhEngine(InfoList, \dots, verbose, kitchensink, preponly, snoop) 
\method{rmhmodel}{rmhmodel}(model, \dots)
\method{rmhstart}{rmhstart}(start, \dots)
\method{rmhstart}{list}(start, \dots)
rmpoint.I.allim(n, f, types)
\method{row.names}{hyperframe}(x)
\method{row.names}{hyperframe}(x) <- value
rpoint.multi(n, f, fmax, marks, win, giveup, verbose, warn)
runifpoispp(lambda, win = owin(c(0, 1), c(0, 1)))
runifpoisppOnLines(lambda, L)
runifrect(n, win = owin(c(0, 1), c(0, 1)))
safelookup(Z, X, factor, warn)
samefunction(f, g)
scanmeasure(X, \dots)
\method{scanmeasure}{ppp}(X, r, \dots, method)
\method{scanmeasure}{im}(X, r, \dots)
scanPoisLRTS(nZ, nG, muZ, muG, alternative)
scanBinomLRTS(nZ, nG, muZ, muG, alternative)
scanLRTS(X, r, \dots, method, baseline, case, alternative)
second.moment.calc(x, sigma=NULL, edge=TRUE, what="Kmeasure", debug=FALSE,
\dots, varcov=NULL, expand=FALSE)
second.moment.engine(x, sigma=NULL, edge=TRUE, what="Kmeasure", debug=FALSE,
\dots, obswin=as.owin(x), varcov=NULL, npts=NULL)
sensiblevarname(guess, fallback, maxlen)
sewpcf(d, w, denargs, lambda2area)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
shiftxy(X, vec = c(0, 0))
shiftxypolygon(p, vec = c(0, 0))
short.deparse(x, maxlen)
simplify.xypolygon(p, dmin)
simulrecipe(type, expr, envir, csr, pois)
singlestring(s, coll)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
smoothpointsEngine(x, values, sigma, \dots,
                    weights, varcov, leaveoneout, sorted)
\method{sort}{im}(x, \dots)
spatstat.rawdata.location(\dots)
spatstatClusterModelInfo(name)
spatstatRmhInfo(cifname)
spatialCDFframe(model, covariate, \dots)
spatialCDFtest(model, covariate, test, \dots, dimyx, eps, jitter, modelname, covname, dataname)
sphere.volume(range, nval = 10)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
store.versionstring.spatstat()
\method{str}{hyperframe}(object, \dots)
stratrand(window, nx, ny, k = 1)
strausscounts(U,X,r,EqualPairs)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
summarise.trend(trend, w, a)
\method{summary}{envelope}(object,\dots)
\method{summary}{fii}(object,\dots)
\method{summary}{hyperframe}(object, ..., brief=FALSE)
\method{summary}{logippm}(object, \dots)
\method{summary}{lut}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{pppmatching}(object, \dots)
\method{summary}{ppx}(object, \dots)
\method{summary}{rmhexpand}(object, \dots)
\method{summary}{splitppx}(object, \dots)
sumsymouter(x, w)
superimposeMarks(arglist, nobj)
sympoly(x, y, n)
termsinformula(x)
test.crossing.psp(A,B)
test.selfcrossing.psp(A)
tilecentroids(W, nx, ny)
tilenames(x)
trap.extra.arguments(\dots, .Context, .Fatal)
trianglediameters(iedge, jedge, edgelength, \dots, nvert, check)
trim.mask(M, R, tolerant)
tweak.fv.entry(x, current.tag, new.labl=NULL, new.desc=NULL, new.tag=NULL)
\method{unitname}{default}(x)
\method{unitname}{default}(x) <- value
\method{update}{interact}(object, \dots)
\method{update}{rmhstart}(object, \dots)
\method{update}{slrm}(object, ..., evaluate = TRUE, env = parent.frame())
validradius(r, caller)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
validposint(n, caller, fatal)
vanilla.fv(x)
variablesinformula(x)
verbalogic(x, op)
versionstring.interact(object)
versionstring.ppm(object)
versionstring.spatstat()
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
verify.xypolygon(p, fatal=TRUE)
warn.ignored.args(\dots, context)
\method{with}{msr}(data, expr, \dots)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
xy.grid(xr, yr, nx, ny, dx, dy)
X2testEngine(OBS, EXP, \dots, method, df, nsim, conditional, alternative, testname, dataname)
\method{xtfrm}{im}(x)
xypolyselfint(p, eps, proper, yesorno, checkinternal)
xypolygon2psp(p, w, check)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
