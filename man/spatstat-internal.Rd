\name{spatstat-internal} 
\title{Internal spatstat functions}
\alias{[.pp3}
\alias{[.localpcfmatrix}
\alias{[.rat}
\alias{[.splitppx}
\alias{[.diagramobj}
\alias{[<-.splitppx}
\alias{acedist.show}
\alias{acedist.noshow}
\alias{active.interactions}
\alias{adaptcoef}
\alias{adjust.ratfv}
\alias{affinexy}
\alias{affinexypolygon}
\alias{ang2rad}
\alias{anycrossing.psp}
\alias{ApplyConnected}
\alias{applytolayers}
\alias{applyPolyclipArgs}
\alias{areadelta2}
\alias{areaGain.diri}       
\alias{areaGain.grid}       
\alias{areaLoss.diri}       
\alias{areaLoss.grid}
\alias{areaLoss.poly}
\alias{assemble.plot.objects}       
\alias{AsymmDistance.psp}
\alias{as.breakpts}
\alias{as.character.unitname}
\alias{as.data.frame.bw.optim}
\alias{as.data.frame.fv}
\alias{as.double.im}
\alias{as.linfun.linfun}
\alias{as.list.hyperframe}
\alias{as.listof}
\alias{as.owin.lintess}
\alias{as.unitname}
\alias{augment.msr}
\alias{avenndist}
\alias{BartCalc}
\alias{bbEngine}
\alias{bermantestCalc}
\alias{bermantestEngine}
\alias{bdry.mask}
\alias{bind.ratfv}
\alias{blankcoefnames}
\alias{bounding.box3}
\alias{break.holes}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{bt.frame}
\alias{bw.optim}
\alias{calc.DR}
\alias{calc.NNIR}
\alias{calc.SAVE}
\alias{calc.SIR}
\alias{calc.TSE}
\alias{cannot.update}
\alias{cartesian}
\alias{cellmiddles}
\alias{censtimeCDFest}
\alias{change.default.expand}          
\alias{checkbigmatrix}          
\alias{checkfields}          
\alias{checksolve}          
\alias{check.finespacing}
\alias{check.hist.lengths}
\alias{check.mat.mul}
\alias{check.testfun}
\alias{circticks}
\alias{clarkevansCalc}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closethresh}
\alias{coef.summary.kppm}
\alias{coef.summary.ppm}
\alias{coef.vblogit} 
\alias{coerce.marks.numeric}
\alias{compatible.rat}
\alias{compileCDF}
\alias{conform.ratfv}
\alias{crosspairquad}
\alias{cobble.xy}
\alias{codetime}
\alias{col.args.to.grey}
\alias{commonPolyclipArgs}
\alias{conform.imagelist}
\alias{countingweights}
\alias{CressieReadStatistic}
\alias{CressieReadSymbol}
\alias{CressieReadName}
\alias{CVforPCF}
\alias{damaged.ppm}
\alias{data.mppm}
\alias{datagen.runifpointOnLines}
\alias{datagen.runifpoisppOnLines}
\alias{datagen.rpoisppOnLines}
\alias{default.clipwindow}
\alias{default.linnet.tolerance}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{deltasuffstat}
\alias{Deviation}
\alias{dfbetas.ppmInfluence}
\alias{dflt.redraw}
\alias{densitycrossEngine}
\alias{densitypointsEngine}
\alias{diagnose.ppm.engine}
\alias{diagramobj}
\alias{digestCovariates}
\alias{digital.volume}
\alias{dilate.owin}
\alias{dim.fasp}               
\alias{dim.hyperframe}               
\alias{dim.im}               
\alias{dim.msr}
\alias{dim.owin}
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{dimnames.msr}
\alias{distributecbind}
\alias{dist2dpath}
\alias{do.as.im}
\alias{do.call.plotfun}
\alias{do.istat}
\alias{doMultiStraussHard}
\alias{dppmFixAlgorithm}
\alias{dppmFixIntensity}
\alias{emptywindow}
\alias{envelopeEngine}
\alias{envelopeProgressData}
\alias{envelopeTest}
\alias{envelope.hasenvelope}
\alias{envelope.matrix}
\alias{equalpairs}          
\alias{equalpairs.quad}          
\alias{equals.quad}          
\alias{equalsfun.quad}          
\alias{erodemask}
\alias{erode.owin}
\alias{evalCovar}
\alias{evalCovar.ppm}
\alias{evalCovar.lppm}
\alias{evalCovariate}
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{evalPairPotential}
\alias{evalSparse3Dentrywise}
\alias{evaluate2Dkernel}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{existsSpatstatVariable}
\alias{expandSpecialLists}
\alias{expandwinPerfect}
\alias{ExpSmoothLog}
\alias{extractAIC.slrm}
\alias{extractAtomicQtests}
\alias{fakemaintitle}
\alias{family.vblogit} 
\alias{f3engine}
\alias{f3Cengine}
\alias{fasp}
\alias{FDMKERNEL}
\alias{fft2D}
\alias{fftwAvailable}
\alias{fill.coefs}
\alias{findbestlegendpos}
\alias{findCovariate}
\alias{findcbind}
\alias{fii}
\alias{fillNA}
\alias{flatfname}
\alias{flipxypolygon}
\alias{forbid.logi}
\alias{format.numberwithunit}
\alias{FormatFaspFormulae}
\alias{framebottomleft}
\alias{fvexprmap}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{gammabreaks}
\alias{g3engine}
\alias{g3Cengine}
\alias{getdataname}
\alias{getfields}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getlambda.lpp}
\alias{getlastshift}
\alias{getppmdatasubset}
\alias{getppmOriginalCovariates}
\alias{getRandomFieldsModelGen}
\alias{getSpatstatVariable}
\alias{getSumFun}
\alias{geyercounts}
\alias{geyerdelta2}
\alias{GLMpredict}
\alias{good.correction.K}
%\alias{gridadjacencymatrix} %DoNotExport
\alias{gridindex}            
\alias{grid1index}
\alias{grokIndexVector}
\alias{grow.mask}
\alias{hackglmmPQL}
\alias{hasenvelope}
\alias{hasglmfit}
\alias{HermiteCoefs}
\alias{handle.r.b.args}
\alias{handle.rshift.args}
\alias{head.hyperframe}
\alias{hierarchicalordering}
\alias{hiermat}
\alias{ho.engine}
\alias{hsvNA}
\alias{IdenticalRows}
\alias{idorempty}
\alias{illegal.iformula}
\alias{implemented.for.K}
\alias{impliedpresence}
\alias{impliedcoefficients}
\alias{influence.ppmInfluence}
\alias{inpoint}
\alias{instantiate.interact}
\alias{interactionfamilyname}
\alias{intermaker}
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}
\alias{invokeColourmapRule}
\alias{is.atomicQtest}
\alias{is.cadlag}
\alias{is.col.argname}
\alias{is.data}
\alias{is.expandable}
\alias{is.expandable.ppm}
\alias{is.expandable.rmhmodel}
\alias{is.fv}
\alias{is.hyperframe}
\alias{is.imlist}
\alias{is.infline}
\alias{is.interact}
\alias{is.marked.default}    
\alias{is.marked.msr}    
\alias{is.marked.psp}    
\alias{is.marked.quad}
\alias{is.mppm}
\alias{is.multitype.msr}    
\alias{is.multitype.quad}    
\alias{is.multitype.default}    
\alias{is.poisson.mppm}
\alias{is.ppplist}
\alias{is.pp3}
\alias{is.ppx}
\alias{is.psp}
\alias{is.quad}
\alias{is.scov}
\alias{is.sob}
\alias{is.solist}
\alias{is.tess}
\alias{is.vanilla}
\alias{k3engine}
\alias{Kborder.engine}
\alias{Knone.engine}
\alias{Krect.engine}
\alias{Kount}
\alias{Kwtsum}               
\alias{Kpcf.kppm}               
\alias{killinteraction}
\alias{km.rs.opt}
\alias{kppmComLik}
\alias{kppmMinCon}
\alias{kppmPalmLik}
\alias{kraever}
\alias{kraeverRandomFields}
\alias{labels.ppm}
\alias{ldtEngine}
\alias{levels.im}
\alias{levels<-.im}
\alias{levelsAsFactor}
\alias{leverage.ppmInfluence}
\alias{linearKengine}
\alias{linearKmulti}
\alias{linearKmulti.inhom}
\alias{linearKmultiEngine}
\alias{linearpcfengine}
\alias{linearpcfmulti}
\alias{linearpcfmulti.inhom}
\alias{linearPCFmultiEngine}
\alias{listof}
\alias{localKengine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{local2lpp}
\alias{logicalIndex}
\alias{logi.dummy}
\alias{logi.engine}
\alias{logLik.vblogit}
\alias{makeLinnetTolerance}
\alias{maskLaslett}
\alias{match2DkernelName}
\alias{parbreak}
\alias{plan.legend.layout}
\alias{PairPotentialType}
\alias{PDEdensityLPP}
\alias{PoisSaddle}
\alias{PoisSaddleArea}
\alias{PoisSaddleGeyer}
\alias{PoisSaddlePairwise}
\alias{polyLaslett}
\alias{polytileareaEngine}
\alias{positiveIndex}
\alias{PPMmodelmatrix}
\alias{putSpatstatVariable}
\alias{lookup.im}
\alias{lookup2DkernelInfo}
\alias{majorminorversion}
\alias{make.even.breaks}
\alias{makefvlabel}
\alias{makeunitname}
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
\alias{\%mapp\%} %DoNotExport
%NAMESPACE export("%mapp%")
\alias{markappendop}
\alias{marksubset}
\alias{markreplicateop}
\alias{\%mrep\%} %DoNotExport
%NAMESPACE export("%mrep%")
\alias{marksubsetop}
\alias{\%msub\%} %DoNotExport
%NAMESPACE export("%msub%")
\alias{mask2df}
\alias{match.kernel}
\alias{maxflow}
\alias{mctestSigtraceEngine}
\alias{meanlistfv}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{model.se.image}
\alias{modelFrameGam}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{mpl.usable}
\alias{MultiPair.checkmatrix}
\alias{multiply.only.finite.entries}
\alias{multiplicityNumeric}
\alias{na.handle.im}
\alias{names.hyperframe}
\alias{names<-.fv}
\alias{names<-.hyperframe}
\alias{nearest.pixel}
\alias{nearest.valid.pixel}
\alias{newformula}
\alias{newstyle.coeff.handling}
\alias{nncleanEngine}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{n.quad}
\alias{numberwithunit}
\alias{numeric.columns}
\alias{objsurfEngine}
\alias{onecolumn}
\alias{optimStatus}
\alias{outdated.interact}
\alias{oversize.quad}    
\alias{owinpolycheck}
\alias{owinpoly2mask}
\alias{owin2polypath}
\alias{pairs.listof}
\alias{pairs.solist}
\alias{param.quad}
\alias{partialModelMatrix}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{pickoption}
\alias{plotEachLayer}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot3Dpoints}
\alias{plotPolygonBdry}
\alias{plot.addvar}
\alias{plot.barplotdata}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.localpcfmatrix}
\alias{plot.lurk}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.plotpairsim}
\alias{plot.pppmatching}
\alias{plot.profilepl}
\alias{plot.qqppm}
\alias{plot.spatialcdf}
\alias{plot.studpermutest}
\alias{ppllengine}
\alias{ppm.default}
\alias{ppmCovariates}
\alias{ppmDerivatives}
\alias{ppmInfluenceEngine}
\alias{pppdist.mat}
\alias{pppdist.prohorov}
\alias{ppsubset}
\alias{predict.vblogit}
\alias{prefixfv}
\alias{prepareTitle}
\alias{printStatus}
\alias{print.addvar}    
\alias{print.anylist}    
\alias{print.autoexec}    
\alias{print.bt.frame}
\alias{print.bugtable}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.colourmap}
\alias{print.diagppm}
\alias{print.distfun}
\alias{print.detpointprocfamily}
\alias{print.detpointprocfamilyfun}
\alias{print.envelope}
\alias{print.ewcdf}
\alias{print.fasp}       
\alias{print.fv}       
\alias{print.fvfun}       
\alias{print.funxy}       
\alias{print.hasenvelope}       
\alias{print.hierarchicalordering}
\alias{print.hyperframe}
\alias{print.indicfun}       
\alias{print.influence.ppm}       
\alias{print.interact}       
\alias{print.intermaker}       
\alias{print.isf}
\alias{print.laslett}
\alias{print.layered}
\alias{print.leverage.ppm}
\alias{print.lintess}
\alias{print.localpcfmatrix}
\alias{print.lurk}
\alias{print.lut}
\alias{print.minconfit}
\alias{print.mppm}
\alias{print.msr}
\alias{print.nnfun}
\alias{print.numberwithunit}
\alias{print.onearrow}
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
\alias{print.rppm}
\alias{print.splitppp}
\alias{print.simplepanel}
\alias{print.Smoothfun}       
\alias{print.solist}
\alias{print.splitppx}
\alias{print.summary.distfun}
\alias{print.summary.funxy}
\alias{print.summary.hyperframe}
\alias{print.summary.listof}
\alias{print.summary.linim}
\alias{print.summary.linnet}
\alias{print.summary.lintess}
\alias{print.summary.logiquad}
\alias{print.summary.lut}
\alias{print.summary.mppm}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{print.summary.rmhexpand}
\alias{print.summary.ssf}
\alias{print.summary.solist}
\alias{print.summary.splitppp}
\alias{print.summary.splitppx}
\alias{print.summary.unitname}
\alias{print.symbolmap}       
\alias{print.textstring}
\alias{print.texturemap}
\alias{print.tess}
\alias{print.timed}
\alias{print.vblogit}
\alias{print.yardstick}
\alias{project3Dhom}
\alias{putlastshift}
\alias{quad}
\alias{quad.mppm}
\alias{quadBlockSizes}
\alias{RandomFieldsSafe}
\alias{rasterfilter}
\alias{ratfv}
\alias{recognise.spatstat.type}
\alias{rectquadrat.breaks}
\alias{rectquadrat.countEngine}
\alias{reduceformula}
\alias{reheat}
\alias{RelevantDeviation}
\alias{repair.image.xycoords}
\alias{replacementIndex}
\alias{representativeRows}
\alias{resolveEinfo}
\alias{resolve.vargamma.shape}
\alias{rgbNA}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{rMaternInhibition}
\alias{RmhExpandRule}
\alias{rmhsnoop}
\alias{rocData}
\alias{rocModel}
\alias{roseContinuous}
\alias{ruletextline}
\alias{quadrat.testEngine}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{pointgrid}
\alias{rastersample}
\alias{rasterx.mask}
\alias{rastery.mask}
\alias{rasterxy.mask}
\alias{rasterx.im}
\alias{rastery.im}
\alias{rasterxy.im}
\alias{rebadge.fv}
\alias{rebadge.as.crossfun}
\alias{rebadge.as.dotfun}
\alias{rebound}
\alias{rebound.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{reconcile.fv}
\alias{rename.fv}
\alias{repair.old.factor.image}
\alias{reincarnate.interact}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.2D.kernel}
\alias{restrict.mask}
\alias{reversePolyclipArgs}
\alias{rmax.Rigid}
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
\alias{safedeldir}
\alias{safelookup}
\alias{scalardilate.breakpts}
\alias{scalardilate.diagramobj}
\alias{scalardilate.msr}
\alias{scanmeasure}
\alias{scanmeasure.ppp}
\alias{scanmeasure.im}
\alias{scanBinomLRTS}
\alias{scanPoisLRTS}
\alias{second.moment.calc}
\alias{second.moment.engine}
\alias{sewpcf}
\alias{sewsmod}
\alias{shift.diagramobj}              
\alias{shift.influence.ppm}              
\alias{shift.leverage.ppm}              
\alias{shift.msr}              
\alias{shift.quadratcount}              
\alias{shift.quadrattest}              
\alias{shiftxy}              
\alias{shiftxypolygon}              
\alias{signalStatus}
\alias{simulate.profilepl}
\alias{simulrecipe}              
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{Smooth.solist}
\alias{smoothcrossEngine}              
\alias{smoothpointsEngine}
\alias{smudge}
\alias{sort.im}
\alias{sortalongsegment}
\alias{spatstat.xy.coords}
\alias{spatstatClusterModelInfo}
\alias{spatstatDPPModelInfo}
\alias{spatstatRmhInfo}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{spatialCDFtestCalc}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{sphere.volume}
\alias{store.versionstring.spatstat}
\alias{str.hyperframe}
\alias{strausscounts}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summarise.trend}
\alias{summary.envelope}
\alias{summary.hyperframe}
\alias{summary.lintess}
\alias{summary.logiquad}
\alias{summary.lut}
\alias{summary.mppm}
\alias{summary.msr}
\alias{summary.profilepl}
\alias{summary.pppmatching}
\alias{summary.ppx}
\alias{summary.splitppx}
\alias{summary.rmhexpand}
\alias{summary.vblogit}
\alias{sumsymouter}
\alias{superimposeMarks}
\alias{symbolmapdomain}
\alias{symbolmaptype}
\alias{tail.hyperframe}
\alias{tensor1x1}
\alias{thinjump}
\alias{tilecentroids}        
\alias{trianglediameters}
\alias{trim.mask}        
\alias{tweak.fv.entry}
\alias{tweak.ratfv.entry}
\alias{tweak.coefs}
\alias{twostage.test}
\alias{\%unit\%} %DoNotExport
%NAMESPACE export("%unit%")
\alias{unitname.default}
\alias{unitname<-.default}
\alias{unstackFilter}
\alias{update.im}
\alias{update.ippm}
\alias{update.rmhstart}
\alias{validradius}
\alias{validate2Dkernel}
\alias{validate.angles}        
\alias{validate.lpp.coords}
\alias{validate.mask}        
\alias{validate.quad}        
\alias{validate.weights}        
\alias{vanilla.fv}
\alias{varcountEngine}
%\alias{vblogit} %DoNotExport
%\alias{vblogit.fmla} %DoNotExport
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{versionstring.spatstat}
\alias{verifyclass}
\alias{vnnFind}
\alias{Window.lintess}
\alias{Window<-.linnet}
\alias{Window<-.lpp}
\alias{warn.once}
\alias{waxlyrical}
\alias{weightedclosepairs}
\alias{windows.mppm}
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xy.grid}
\alias{X2testEngine}
\alias{xtfrm.im}
\alias{xypolygon2psp}
\alias{xypolyselfint}
%%%% sparse 3D arrays
\alias{sparse3Darray}
\alias{as.sparse3Darray}
\alias{dim.sparse3Darray}
\alias{dim<-.sparse3Darray}
\alias{dimnames.sparse3Darray}
\alias{dimnames<-.sparse3Darray}
\alias{print.sparse3Darray}
\alias{aperm.sparse3Darray}
\alias{as.array.sparse3Darray}
\alias{[.sparse3Darray}
\alias{[<-.sparse3Darray}
\alias{anyNA.sparse3Darray}
\alias{RelevantZero}
\alias{RelevantEmpty}
\alias{isRelevantZero}
\alias{unionOfSparseIndices}
\alias{Math.sparse3Darray}
\alias{Ops.sparse3Darray}
\alias{Summary.sparse3Darray}
\alias{inside3Darray}
\alias{SparseEntries}
\alias{SparseIndices}
\alias{EntriesToSparse}
\alias{mapSparseEntries}
\alias{applySparseEntries}
\alias{sumsymouterSparse}
\alias{tenseur}
\alias{marginSums}
\alias{rbindCompatibleDataFrames}
\alias{bind.sparse3Darray}
%%%%
\alias{spatstatDiagnostic}
%%
\alias{as.ppplist}
\alias{as.imlist}
\alias{pointsAlongNetwork}
\alias{expandSparse}
\alias{allElementsIdentical}
\alias{resampleNetworkDataFrame}
\alias{sparseVectorCumul}
%%
\alias{as.ppm.lppm}
\alias{as.ppm.rppm}
\alias{predict.profilepl}
%%%%%%%
\description{
  Internal spatstat functions.
}
\usage{
\method{[}{splitppx}(x, \dots)
\method{[}{splitppx}(x, \dots) <- value
\method{[}{diagramobj}(x, \dots)
\method{[}{rat}(x, \dots)
acedist.show(X, Y, n, d, timelag)
acedist.noshow(X, Y, n, d)
active.interactions(object)
adaptcoef(new.coef, fitcoef, drop)
adjust.ratfv(f, columns, numfactor, denfactor)
affinexy(X, mat, vec, invert)
affinexypolygon(p, mat, vec, detmat)
ang2rad(ang, unit, start, clockwise)
anycrossing.psp(A,B)
ApplyConnected(X, Engine, r, \dots, rule, auxdata)
applytolayers(L, FUN, \dots)
applyPolyclipArgs(x, p)
areadelta2(X, r, \dots, sparseOK)
areaGain.diri(u, X, r, \dots, W, verbose)
areaGain.grid(u, X, r, \dots, W, ngrid)
areaLoss.diri(X, r, \dots, W, subset)
areaLoss.grid(X, r, \dots, W, subset,
                         method = c("count", "distmap"),
                         ngrid = spatstat.options("ngrid.disc"),
                         exact = FALSE)
areaLoss.poly(X, r, \dots, W, subset, splitem)
assemble.plot.objects(xlim, ylim, \dots, lines, polygon)
AsymmDistance.psp(X, Y, metric, method)
as.breakpts(\dots)
\method{as.character}{unitname}(x, \dots)
\method{as.data.frame}{bw.optim}(x, \dots)
\method{as.data.frame}{fv}(x, \dots)
\method{as.double}{im}(x, \dots)
\method{as.linfun}{linfun}(X, \dots)
\method{as.list}{hyperframe}(x, \dots)
as.listof(x)
\method{as.owin}{lintess}(W, \dots)
as.unitname(s)
augment.msr(x, \dots, sigma, recompute)
avenndist(X)
BartCalc(fY, fK)
bbEngine(\dots)
bermantestCalc(fram, which, alternative, \dots)
bermantestEngine(model, covariate, which, alternative, \dots,
                 modelname, covname, dataname)
bdry.mask(W)
bind.ratfv(x, numerator, denominator, labl, desc, preferred,
           ratio, quotient)
blankcoefnames(x)
bounding.box3(\dots)
break.holes(x, splitby, depth, maxdepth)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
bt.frame(Q, trend, interaction, \dots, covariates,
         correction, rbord, use.gam, allcovar)
bw.optim(cv, h, iopt, \dots, cvname, hname, criterion, unitname)
calc.DR(COV, z, Dim)
calc.NNIR(COV, z, pos, Dim)
calc.SAVE(COV, z, Dim)
calc.SIR(COV, z)
calc.TSE(COV, z, pos, Dim1, Dim2)
cannot.update(\dots)
cartesian(pp, markset, fac = TRUE)
cellmiddles(W, nx, ny, npix, distances)
censtimeCDFest(o, cc, d, breaks, \dots,
     KM, RS, HAN, RAW, han.denom, tt, pmax)
change.default.expand(x, newdefault)
checkbigmatrix(n, m, fatal, silent)
checkfields(X,L)
checksolve(M, action, descrip, target)
check.finespacing(r, eps, win, rmaxdefault, context, action, rname)
check.hist.lengths(hist,breaks)
check.mat.mul(A, B, Acols, Brows, fatal)
check.testfun(f, f1, X)
circticks(R, at, unit, start, clockwise, labels)
clarkevansCalc(X, correction, clipregion, working)
clip.psp(x, window, check, fragments)
cliprect.psp(x, window, fragments)
clippoly.psp(s, window, fragments)
closethresh(X,R,S,twice,\dots)
\method{coef}{summary.kppm}(object, \dots)
\method{coef}{summary.ppm}(object, \dots)
\method{coef}{vblogit}(object, \dots) 
coerce.marks.numeric(X, warn)
\method{compatible}{rat}(A, B, \dots) 
compileCDF(D, B, r, \dots, han.denom, check)
conform.ratfv(x)
crosspairquad(Q,rmax,what)
cobble.xy(x, y, f, fatal, \dots)
codetime(x, hms, what)
col.args.to.grey(x, \dots)
commonPolyclipArgs(\dots, p)
conform.imagelist(X, Zlist)
countingweights(id, areas, check = TRUE)
CressieReadStatistic(OBS,EXP,lambda)
CressieReadSymbol(lambda)
CressieReadName(lambda)
CVforPCF(bw, stuff)
damaged.ppm(object)
data.mppm(x)
datagen.runifpointOnLines(n, L)
datagen.runifpoisppOnLines(lambda, L)
datagen.rpoisppOnLines(lambda, L, lmax, \dots, check)
default.clipwindow(object, epsilon)
default.linnet.tolerance(L)
default.n.tiling(X, nd, ntile, npix, eps, random, quasi, verbose)
default.ntile(X)
deltasuffstat(model, \dots,
              restrict, dataonly, sparseOK, quadsub,
              force, warn.forced, verbose, use.special)
Deviation(x, ref, leaveout, n, xi)
\method{dfbetas}{ppmInfluence}(model, \dots)
dflt.redraw(button, name, env)
densitycrossEngine(Xdata, Xquery, sigma, \dots,
                    weights, edge, varcov,
                    diggle, sorted)
densitypointsEngine(x, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted, spill, cutoff)
diagnose.ppm.engine(object, \dots, type, typename, opt,
                         sigma, rbord, compute.sd, compute.cts,
                         envelope, nsim, nrank,
                         rv, oldstyle, splineargs, verbose)
diagramobj(X, \dots)
digestCovariates(\dots, W)
digital.volume(range, nval, vside)
dilate.owin(\dots)
\method{dim}{fasp}(x)
\method{dim}{hyperframe}(x)
\method{dim}{im}(x)
\method{dim}{msr}(x)
\method{dim}{owin}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
\method{dimnames}{msr}(x)
distributecbind(x)
dist2dpath(dist, method="C")
do.as.im(x, action, \dots, W, eps, dimyx, xy, na.replace)
do.call.plotfun(fun, arglist, \dots)
do.istat(panel)
doMultiStraussHard(iradii, hradii, types)
dppmFixIntensity(DPP, lambda, po)
dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
emptywindow(w)
envelopeEngine(X, fun, simul,
           nsim, nrank, \dots, funargs, funYargs,
           verbose, clipdata, 
           transform, global, ginterval, use.theory,
           alternative, scale, clamp,
           savefuns, savepatterns, saveresultof,
           weights,
           nsim2, VARIANCE, nSD,
           Yname, maxnerr, internal, cl,
           envir.user, expected.arg, do.pwrong,
           foreignclass, collectrubbish)
envelopeProgressData(X, fun, \dots, exponent,
                     alternative, leaveout, scale, clamp,
                     normalize, deflate, rmin, 
                     save.envelope, savefuns, savepatterns)
envelopeTest(X, \dots, exponent, alternative,
            rinterval, leaveout, scale, clamp, tie.rule,
            interpolate, save.interpolant,
            save.envelope, savefuns, savepatterns,
            Xname, verbose)
\method{envelope}{hasenvelope}(Y, \dots, Yname)
\method{envelope}{matrix}(Y, \dots, rvals, observed, theory, funX,
  nsim, nsim2, jsim, jsim.mean,
  type, alternative, scale, clamp, csr, use.theory, nrank, ginterval, nSD,
  savefuns, check, Yname, do.pwrong, weights, precomputed)
equalpairs(U, X, marked=FALSE)
equalpairs.quad(Q)
equals.quad(Q)          
equalsfun.quad(Q)
erodemask(w,r,strict)
erode.owin(\dots)
evalCovar(model, covariate, \dots)
\method{evalCovar}{ppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, modelname, covname, dataname,
          subset) 
\method{evalCovar}{lppm}(model, covariate, \dots, lambdatype, 
          eps, nd, interpolate, jitter, modelname, covname, dataname,
          subset)
evalCovariate(covariate, locations)
evalInteraction(X,P,E,interaction,correction,splitInf,\dots,
                precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,splitInf,\dots,
                Reach,precomputed,savecomputed)
evalPairPotential(X,P,E,pairpot,potpars,R)
evalSparse3Dentrywise(expr, envir)
evaluate2Dkernel(kernel, x, y, sigma, varcov, \dots, scalekernel)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(w)
existsSpatstatVariable(name)
expandSpecialLists(x, special)
expandwinPerfect(W, expand, amount)
ExpSmoothLog(X, \dots, at, weights)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, \dots)
extractAtomicQtests(x)
fakemaintitle(bb, main, \dots)
\method{family}{vblogit}(object, \dots)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
FDMKERNEL(lppobj, sigma, dtt, weights, iterMax, sparse, dtx)
fft2D(z, inverse, west)
fftwAvailable()
fill.coefs(coefs, required)
findbestlegendpos(\dots)
findCovariate(covname, scope, scopename=NULL)
findcbind(root, depth, maxdepth)
fii(interaction, coefs, Vnames, IsOffset)
fillNA(x, value)
flatfname(x)
flipxypolygon(p)
forbid.logi(object)
\method{format}{numberwithunit}(x, \dots, collapse, modifier)
FormatFaspFormulae(f, argname)
framebottomleft(w)
fvexprmap(x)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
gammabreaks(ra, n, gamma)
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
getdataname(defaultvalue, \dots, dataname)
getfields(X, L, fatal = TRUE)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getlambda.lpp(lambda, X, subset, \dots,
              update, leaveoneout, loo.given, lambdaname)
getlastshift(X)
getppmdatasubset(object)
getppmOriginalCovariates(object)
getRandomFieldsModelGen(model)
getSpatstatVariable(name)
getSumFun(abbreviation, classname, ismarked, fatal)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
geyerdelta2(X,r,sat,\dots,sparseOK, correction)
GLMpredict(fit, data, coefs, changecoef, type)
good.correction.K(X)
%gridadjacencymatrix(dims)
gridindex(x, y, xrange, yrange, nx, ny)            
grid1index(x, xrange, nx)
grokIndexVector(ind, len, nama)
grow.mask(M, xmargin=0, ymargin=xmargin)
hackglmmPQL(fixed, random, family, data, correlation, weights,
            control, niter, verbose, subset, \dots, reltol)
hasenvelope(X, E)
hasglmfit(object)
HermiteCoefs(order)
handle.r.b.args(r = NULL, breaks = NULL, window, pixeps = NULL, rmaxdefault)
handle.rshift.args(W, \dots, radius, width, height, edge, clip,
edgedefault)
\method{head}{hyperframe}(x,n,\dots)
hierarchicalordering(i, s)
hiermat(x, h)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
hsvNA(h, s, v, alpha)
IdenticalRows(i,j,a,b)
idorempty(w, r, caller)
illegal.iformula(ifmla, itags, dfvarnames)
implemented.for.K(correction, windowtype, explicit)
impliedpresence(tags, formula, df, extranames=character(0))
impliedcoefficients(object, tag)
\method{influence}{ppmInfluence}(model, \dots)
inpoint(W)
instantiate.interact(x, par)
interactionfamilyname(x)
intermaker(f, blank)
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
invokeColourmapRule(colfun, x, \dots, zlim, colargs)
is.atomicQtest(x)
is.cadlag(s)
is.col.argname(x)
is.data(Q)
is.expandable(x)
\method{is.expandable}{ppm}(x)  
\method{is.expandable}{rmhmodel}(x) 
is.fv(x)
is.hyperframe(x)
is.imlist(x)
is.infline(x)
is.interact(x)
\method{is.marked}{default}(\dots)  
\method{is.marked}{msr}(X, \dots)
\method{is.marked}{psp}(X, \dots)
\method{is.marked}{quad}(X, na.action="warn", \dots)
is.mppm(x)
\method{is.multitype}{default}(X, \dots)  
\method{is.multitype}{msr}(X, \dots)
\method{is.multitype}{quad}(X, na.action="warn", \dots)
\method{is.poisson}{mppm}(x)
is.ppplist(x)
is.pp3(x)
is.ppx(x)
is.psp(x)
is.quad(x)
is.scov(x)
is.solist(x)
is.sob(x)
is.tess(x)
is.vanilla(u)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Krect.engine(X, rmax, nr, correction, weights, ratio, fname)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks)
Kpcf.kppm(model, what)
killinteraction(model)
km.rs.opt(o, cc, d, breaks, KM, RS)
kppmComLik(X, Xname, po, clusters, control, weightfun, rmax,
           algorithm, DPP, \dots)
kppmMinCon(X, Xname, po, clusters, control, statistic, statargs,
           algorithm, DPP, \dots)
kppmPalmLik(X, Xname, po, clusters, control, weightfun, rmax, 
           algorithm, DPP, \dots)
kraever(package, fatal)
kraeverRandomFields()
\method{labels}{ppm}(object, \dots)
ldtEngine(nv, ns, from, to, seglen, huge,
          coUXord, vnndist, vnnwhich, vnnlab)
\method{levels}{im}(x)
\method{levels}{im}(x) <- value
levelsAsFactor(x)
\method{leverage}{ppmInfluence}(model, \dots)
linearKengine(X, \dots, r, reweight, denom, correction,
              ratio, showworking)
linearKmulti(X, I, J, r, \dots, correction)
linearKmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots, correction,
             normalise)
linearpcfengine(X, \dots, r, reweight, denom, correction, ratio)
linearpcfmulti(X, I, J, r, \dots, correction)
linearpcfmulti.inhom(X, I, J, lambdaI, lambdaJ, r, \dots,
                     correction, normalise)
linearKmultiEngine(X, I, J, \dots,
                   r, reweight, denom, correction, showworking)
linearPCFmultiEngine(X, I, J, \dots,
                   r, reweight, denom, correction, showworking)
listof(\dots)
localKengine(X, \dots, wantL, lambda, correction, verbose, rvalue)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
local2lpp(L, seg, tp, X, df.only)
logicalIndex(i, nama, len)
logi.dummy(X, dummytype, nd, mark.repeat, \dots)
logi.engine(Q, trend, interaction, \dots, 
            covariates, subsetexpr, clipwin,
            correction, rbord, covfunargs, allcovar, 
            vnamebase, vnameprefix, justQ, savecomputed, precomputed,
            VB)
\method{logLik}{vblogit}(object, \dots) 
makeLinnetTolerance
maskLaslett(X, \dots, eps, dimyx, xy, oldX, verbose, plotit)
match2DkernelName(kernel)
parbreak(terse)
plan.legend.layout(B, \dots, side, sep, size, sep.frac, size.frac,
                   started, map)
PairPotentialType(pairpot)
PDEdensityLPP(x, sigma, \dots, weights, dx, dt,
              fun, finespacing, finedata)
PoisSaddle(beta, fi)
PoisSaddleArea(beta, fi)
PoisSaddleGeyer(beta, fi)
PoisSaddlePairwise(beta, fi)
polyLaslett(X, \dots, oldX, verbose, plotit)
polytileareaEngine(P, xrange, yrange, nx, ny, DivideByPixelArea)
positiveIndex(i, nama, len)
PPMmodelmatrix(object, data, \dots, subset, Q, keepNA, irregular,
               splitInf)
\method{print}{localpcfmatrix}(x, \dots)
\method{plot}{localpcfmatrix}(x, \dots)
putSpatstatVariable(name, value)
\method{[}{localpcfmatrix}(x, i, \dots)
\method{[}{pp3}(x, i, drop, \dots)
lookup.im(Z, x, y, naok, strict)
lookup2DkernelInfo(kernel)
majorminorversion(v)
make.even.breaks(bmax, npos, bstep)
makefvlabel(op, accent, fname, sub, argname)
makeunitname(sing, plur, mul)
markappend(\dots)
markcbind(\dots)
markformat(x)
\method{markformat}{ppp}(x) 
\method{markformat}{ppx}(x) 
\method{markformat}{psp}(x) 
\method{markformat}{default}(x) 
mark.scale.default(marx, w, \dots, markscale, maxsize, meansize, characters)
markspace.integral(X)
\method{marks}{default}(x, \dots)
\method{marks}{quad}(x, dfok=FALSE, \dots)
markappendop(x, y)
x \%mapp\% y
marksubset(x, index, format)
marksubsetop(x, i)
x \%msub\% i
markreplicateop(x, n)
x \%mrep\% n
mask2df(w)
match.kernel(kernel)
maxflow(costm)
mctestSigtraceEngine(R, devdata, devsim, \dots,
     interpolate, confint, alpha, exponent, unitname)
meanlistfv(z, \dots)
meanX.owin(w)            
meanY.owin(w)
model.se.image(fit, W, \dots, what)
modelFrameGam(formula, \dots)
mpl.engine(Q, trend, interaction, \dots,
         covariates, subsetexpr, clipwin, covfunargs, correction,
	 rbord, use.gam, gcontrol,
         GLM, GLMfamily, GLMcontrol, famille,
         forcefit, nd, eps, allcovar, callstring,
         precomputed, savecomputed, preponly,
         rename.intercept, justQ, weightfactor)
mpl.get.covariates(covariates, locations, type, covfunargs, need.deriv)
mpl.prepare(Q, X, P, trend, interaction, covariates, 
            want.trend, want.inter, correction, rbord, Pname,
            callstring, \dots,
            subsetexpr,
            covfunargs, allcovar, precomputed, savecomputed,
            vnamebase, vnameprefix, warn.illegal, warn.unidentifiable,
            weightfactor, skip.border, clip.interaction, splitInf)
mpl.usable(x)
MultiPair.checkmatrix(mat, n, matname, naok, zerook, asymmok)
multiplicityNumeric(x)
multiply.only.finite.entries(x, a)
na.handle.im(X, na.replace)
\method{names}{fv}(x) <- value
\method{names}{hyperframe}(x)
\method{names}{hyperframe}(x) <- value
nearest.pixel(x, y, Z)
nearest.valid.pixel(x, y, Z)
newformula(old, change, eold, enew)
newstyle.coeff.handling(object)
nncleanEngine(kthNND, k, d, \dots, tol, maxit,
              plothist, lineargs, verbose, Xname)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
n.quad(Q)
numberwithunit(x, u)
numeric.columns(M, logical, others)
objsurfEngine(objfun, optpar, objargs,
              \dots, dotargs, objname,
              ngrid, ratio, verbose)
onecolumn(m)
optimStatus(x, call)
printStatus(x, errors.only)
signalStatus(x, errors.only)
outdated.interact(object)
oversize.quad(Q, \dots, nU, nX, p)
owinpolycheck(W, verbose=TRUE)
owinpoly2mask(w, rasta, check=TRUE)
owin2polypath(w)
\method{pairs}{listof}(\dots, plot=TRUE)
\method{pairs}{solist}(\dots, plot=TRUE)
param.quad(Q)
partialModelMatrix(X,D,model,callstring,\dots)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, \dots,
               r = NULL, breaks = NULL,
               kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
               correction = c("translate", "Ripley"),
               sigma = NULL, varcov = NULL,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J")
pickoption(what="option", key, keymap, \dots,
           exact=FALSE, list.on.err=TRUE, die=TRUE, multi=FALSE,
           allow.all=TRUE)
plotEachLayer(x, \dots, main, plotargs, add, show.all, do.plot)
ploterodewin(W1, W2, col.edge, col.inside, do.plot, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ, do.plot)
plot3Dpoints(xyz, eye, org,
             \dots,
             type, xlim, ylim, zlim,
             add, box, main, cex, box.back, box.front)
plotPolygonBdry(x, \dots)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{barplotdata}(x, \dots)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{lurk}(x, \dots, shade)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{pppmatching}(x, addmatch = NULL, main = NULL, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{profilepl}(x, \dots, add, main, tag, coeff, xvariable,
                         col, lty, lwd, col.opt, lty.opt, lwd.opt)
\method{plot}{qqppm}(x, \dots, limits=TRUE,
           monochrome=spatstat.options('monochrome'),
           limcol=if(monochrome) "black" else "red")
\method{plot}{spatialcdf}(x, \dots, xlab, ylab)
\method{plot}{studpermutest}(x, fmla, \dots,
        lty = NULL, col = NULL, lwd = NULL,
        lty.theo = NULL, col.theo = NULL, lwd.theo = NULL,
        lwd.mean = if (meanonly) 1 else NULL, lty.mean = lty, col.mean = col, 
        separately = FALSE, meanonly = FALSE,
        main = if (meanonly) "group means" else NULL,
        xlim = NULL, ylim = NULL, ylab = NULL, legend = !add,
        legendpos = "topleft", lbox = FALSE, add = FALSE)
ppllengine(X, Y, action="project", check=FALSE)
\method{ppm}{default}(Q, trend, interaction,
       \dots, covariates, data, covfunargs, subset, clipwin,
       correction, rbord, use.gam, method, forcefit, emend, project,
       prior.mean, prior.var,
       nd, eps, gcontrol, nsim, nrmh, start, control,
       verb, callstring)
ppmCovariates(model)
ppmDerivatives(fit, what, Dcovfun, loc, covfunargs)
ppmInfluenceEngine(fit, what, \dots, iScore, iHessian, iArgs,
              drop, method, fine, precomputed, sparseOK,
              fitname, multitypeOK, entrywise, matrix.action,
              dimyx, eps,
              geomsmooth)
pppdist.mat(X, Y, cutoff = 1, q = 1, matching = TRUE,
            precision = 9, approximation = 10)
pppdist.prohorov(X, Y, n, dfix, type, cutoff, matching,
            ccode, auction, precision, approximation) 
ppsubset(X, I, Iname, fatal)
\method{predict}{vblogit}(object, newdata, type, se.fit, dispersion,
                          terms, na.action, \dots)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
prepareTitle(main)
\method{print}{addvar}(x, \dots)
\method{print}{anylist}(x, \dots)
\method{print}{autoexec}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{bugtable}(x, \dots)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{colourmap}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{distfun}(x, \dots)
\method{print}{detpointprocfamily}(x, \dots)
\method{print}{detpointprocfamilyfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{ewcdf}(x, digits, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{funxy}(x, \dots)
\method{print}{fv}(x, \dots, tight)
\method{print}{fvfun}(x, \dots)
\method{print}{hasenvelope}(x, \dots)
\method{print}{hierarchicalordering}(x, \dots)
\method{print}{hyperframe}(x, \dots)
\method{print}{indicfun}(x, \dots)
\method{print}{influence.ppm}(x, \dots)
\method{print}{interact}(x, \dots, family, brief, banner)       
\method{print}{intermaker}(x, \dots)
\method{print}{isf}(x, \dots)
\method{print}{laslett}(x, \dots)
\method{print}{layered}(x, \dots)
\method{print}{leverage.ppm}(x, \dots)
\method{print}{lintess}(x, \dots)
\method{print}{lurk}(x, \dots)
\method{print}{lut}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{mppm}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{nnfun}(x, \dots)
\method{print}{numberwithunit}(x, \dots)
\method{print}{onearrow}(x, \dots)
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
\method{print}{rppm}(x, \dots)
\method{print}{simplepanel}(x, \dots)
\method{print}{Smoothfun}(x, \dots)
\method{print}{solist}(x, \dots)
\method{print}{splitppp}(x, \dots)
\method{print}{splitppx}(x, \dots)
\method{print}{summary.distfun}(x, \dots)
\method{print}{summary.funxy}(x, \dots)
\method{print}{summary.hyperframe}(x, \dots)
\method{print}{summary.linim}(x, \dots)
\method{print}{summary.linnet}(x, \dots)
\method{print}{summary.lintess}(x, \dots)
\method{print}{summary.listof}(x, \dots)
\method{print}{summary.logiquad}(x, \dots, dp=3)
\method{print}{summary.lut}(x, \dots)
\method{print}{summary.mppm}(x, \dots, brief)
\method{print}{summary.owin}(x, \dots)
\method{print}{summary.ppp}(x, \dots, dp)
\method{print}{summary.psp}(x, \dots)
\method{print}{summary.rmhexpand}(x, \dots)
\method{print}{summary.ssf}(x, \dots)
\method{print}{summary.splitppp}(x, \dots)
\method{print}{summary.solist}(x, \dots)
\method{print}{summary.splitppx}(x, \dots)
\method{print}{summary.unitname}(x, \dots)
\method{print}{symbolmap}(x, \dots)
\method{print}{textstring}(x, \dots)
\method{print}{texturemap}(x, \dots)
\method{print}{tess}(x, \dots, brief=FALSE)
\method{print}{timed}(x, \dots)
\method{print}{vblogit}(x, \dots)
\method{print}{yardstick}(x, \dots)
project3Dhom(xyz, eye, org, vert)
putlastshift(X, vec)
quad(data, dummy, w, param)
quad.mppm(x)
quadBlockSizes(nX, nD, p, nMAX, announce)
RandomFieldsSafe()
rasterfilter(X, f)
ratfv(df, numer, denom, \dots, ratio)
recognise.spatstat.type(x)
rectquadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
rectquadrat.countEngine(x, y, xbreaks, ybreaks, weights)
reduceformula(fmla, deletevar, verbose)
reheat(model, invtemp)
RelevantDeviation(x, alternative, clamp, scaling)
repair.image.xycoords(x)
replacementIndex(ii, stuff)
representativeRows(x)
resolveEinfo(x, what, fallback, warn, atomic)
resolve.vargamma.shape(\dots, nu.ker, nu.pcf, default = FALSE)
rgbNA(red, green, blue, alpha, maxColorValue)
rhohatEngine(model, covariate, reference, volume, \dots,
               subset, weights, method, horvitz, smoother, resolution, 
               n, bw, adjust, from, to, 
               bwref, covname, covunits, confidence, modelcall, callstring)
rhohatCalc(ZX, Zvalues, lambda, denom, \dots,
           weights, lambdaX,
           method, horvitz, smoother,
           n, bw, adjust, from, to, 
           bwref, covname, confidence,
           covunits, modelcall, callstring, savestuff)
rMaternInhibition(type, kappa, r, win, stationary, \dots, nsim, drop)
RmhExpandRule(nama)
rocData(covariate, nullmodel, \dots, high)
rocModel(lambda, nullmodel, \dots, high)
rmhsnoop(\dots, Wsim, Wclip, R, xcoords, ycoords, mlevels, mcodes, irep, itype, 
     proptype, proplocn, propmark, propindx, numerator, denominator)
roseContinuous(ang, rad, unit, \dots,
               start, clockwise, main, labels, at, do.plot)
ruletextline(ch, n, terse)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, CR, 
     \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, Xname, fitname)
quadscheme.replicated(data, dummy, method, \dots)
quadscheme.spatial(data, dummy, method, \dots)
pointgrid(W, ngrid)
rastersample(X, Y)
rasterx.mask(w, drop)
rastery.mask(w, drop)
rasterxy.mask(w, drop)
rasterx.im(x)
rastery.im(x)
rasterxy.im(x, drop)
rebadge.fv(x, new.ylab, new.fname, tags, new.desc, new.labl, new.yexp,
           new.dotnames, new.preferred, new.formula, new.tags)
rebadge.as.crossfun(x, main, sub, i, j)
rebadge.as.dotfun(x, main, sub, i)
rebound(x, rect)
\method{rebound}{im}(x, rect)  
\method{rebound}{ppp}(x, rect) 
\method{rebound}{psp}(x, rect) 
\method{rebound}{owin}(x, rect)
reconcile.fv(\dots)
rename.fv(x, fname, ylab, yexp)
repair.old.factor.image(x)
reincarnate.interact(object)
resid4plot(RES, plot.neg, plot.smooth,
           spacing, outer, srange, monochrome, main,
           xlab, ylab, rlab, col.neg, col.smooth, \dots)
resid1plot(RES, opt, plot.neg, plot.smooth,
              srange, monochrome, main,
              add, show.all, do.plot, col.neg, col.smooth, \dots)
resid1panel(observedX, observedV,
            theoreticalX, theoreticalV, theoreticalSD,
            xlab,ylab, \dots, do.plot)
resolve.2D.kernel(\dots,
            sigma, varcov, x, mindist, adjust, bwfun, allow.zero)
restrict.mask(M, W)
reversePolyclipArgs(x, p)
rmax.Rigid(X, g)
rmax.rule(fun, W, lambda)
rotxy(X, angle = pi/2)
rotxypolygon(p, angle = pi/2)
rmhResolveControl(control, model)
rmhResolveExpansion(win, control, imagelist, itype)
rmhResolveTypes(model, start, control)
rmhSnoopEnv(Xinit, Wclip, R)
\method{rmhcontrol}{rmhcontrol}(\dots) 
\method{rmhcontrol}{list}(\dots) 
rmhEngine(InfoList, \dots, verbose, kitchensink, preponly, snoop,
                           overrideXstart, overrideclip) 
\method{rmhmodel}{rmhmodel}(model, \dots) 
\method{rmhstart}{rmhstart}(start, \dots) 
\method{rmhstart}{list}(start, \dots) 
rmpoint.I.allim(n, f, types)
\method{row.names}{hyperframe}(x)
\method{row.names}{hyperframe}(x) <- value
rpoint.multi(n, f, fmax, marks, win, giveup, verbose, warn, nsim, drop)
runifpoispp(lambda, win, \dots, nsim, drop)
runifpoisppOnLines(lambda, L, nsim)
runifrect(n, win, nsim, drop)
safedeldir(X)
safelookup(Z, x, factor, warn)
\method{scalardilate}{breakpts}(X, f, \dots)
\method{scalardilate}{diagramobj}(X, f, \dots)
\method{scalardilate}{msr}(X, f, \dots)
scanmeasure(X, \dots)
\method{scanmeasure}{ppp}(X, r, \dots, method) 
\method{scanmeasure}{im}(X, r, \dots)
scanPoisLRTS(nZ, nG, muZ, muG, alternative)
scanBinomLRTS(nZ, nG, muZ, muG, alternative)
second.moment.calc(x, sigma, edge, what, \dots,
                   varcov=NULL, expand=FALSE, debug=FALSE)
second.moment.engine(x, sigma, edge, what, \dots,
      kernel, obswin, varcov, npts, debug)
sewpcf(d, w, denargs, lambda2area, divisor)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
\method{shift}{diagramobj}(X, \dots)
\method{shift}{influence.ppm}(X, \dots)
\method{shift}{leverage.ppm}(X, \dots)
\method{shift}{msr}(X, \dots)
\method{shift}{quadratcount}(X, \dots)
\method{shift}{quadrattest}(X, \dots)
shiftxy(X, vec = c(0, 0))
shiftxypolygon(p, vec = c(0, 0))
\method{simulate}{profilepl}(object, \dots)
simulrecipe(type, expr, envir, csr, pois, constraints)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
\method{Smooth}{solist}(X, \dots)
smoothcrossEngine(Xdata, Xquery, values, sigma, \dots,
                    weights, varcov, sorted)
smoothpointsEngine(x, values, sigma, \dots,
                    weights, varcov, leaveoneout, sorted, cutoff)
smudge(X)
\method{sort}{im}(x, \dots)
sortalongsegment(df)
spatstat.xy.coords(x, y)
spatstatClusterModelInfo(name, onlyPCP)
spatstatDPPModelInfo(model)
spatstatRmhInfo(cifname)
spatialCDFframe(model, covariate, \dots)
spatialCDFtest(model, covariate, test, \dots,
         dimyx, eps, interpolate, jitter,
         nsim, verbose, modelname, covname, dataname)
spatialCDFtestCalc(fra, test, \dots, details)
sphere.volume(range, nval = 10)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
store.versionstring.spatstat()
\method{str}{hyperframe}(object, \dots)
strausscounts(U,X,r,EqualPairs)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
summarise.trend(trend, w, a)
\method{summary}{envelope}(object,\dots)
\method{summary}{hyperframe}(object, \dots, brief=FALSE)
\method{summary}{lintess}(object, \dots)
\method{summary}{logiquad}(object, \dots, checkdup=FALSE)
\method{summary}{lut}(object, \dots)
\method{summary}{mppm}(object, \dots, brief=FALSE)
\method{summary}{msr}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{pppmatching}(object, \dots)
\method{summary}{ppx}(object, \dots)
\method{summary}{rmhexpand}(object, \dots)
\method{summary}{splitppx}(object, \dots)
\method{summary}{vblogit}(object, \dots)
sumsymouter(x, w)
superimposeMarks(arglist, nobj)
symbolmapdomain(x)
symbolmaptype(x)
\method{tail}{hyperframe}(x,n,\dots)
tensor1x1(A,B)
thinjump(n, p)
tilecentroids(W, nx, ny)
trianglediameters(iedge, jedge, edgelength, \dots,
                  nvert, dmax, check)
trim.mask(M, R, tolerant)
tweak.fv.entry(x, current.tag, new.labl, new.desc, new.tag)
tweak.ratfv.entry(x, \dots)
tweak.coefs(model, new.coef)
twostage.test(X, \dots, exponent, nsim, nsimsub,
              alternative, reuse, leaveout, interpolate,
              savefuns, savepatterns, verbose, testblurb)
x \%unit\% u
\method{unitname}{default}(x) 
\method{unitname}{default}(x) <- value 
unstackFilter(x)
\method{update}{im}(object, \dots) 
\method{update}{ippm}(object, \dots, envir) 
\method{update}{rmhstart}(object, \dots)
validradius(r, caller)
validate2Dkernel(kernel, fatal)
validate.angles(angles, unit, guess)
validate.lpp.coords(X, fatal, context)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
validate.weights(x, recip, how, allowzero, allowinf)
vanilla.fv(x)
varcountEngine(g, B, lambdaB, f)
%vblogit(y, X, offset, eps, m0, S0, S0i, xi0, verb, maxiter, \dots)
%vblogit.fmla(formula, offset, data, subset, weights, verbose, epsilon, \dots)
versionstring.interact(object)
versionstring.ppm(object)
versionstring.spatstat()
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax)
\method{Window}{lintess}(X, \dots)
\method{Window}{linnet}(X, \dots, check=TRUE) <- value
\method{Window}{lpp}(X, \dots, check=TRUE) <- value
warn.once(key, \dots)
waxlyrical(type, terse)
weightedclosepairs(X, r, correction, what)
windows.mppm(x)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
xy.grid(xr, yr, nx, ny, dx, dy)
X2testEngine(OBS, EXP, \dots, method, CR, df, nsim,
     conditional, alternative, testname, dataname)
\method{xtfrm}{im}(x)
xypolyselfint(p, eps, proper, yesorno, checkinternal)
xypolygon2psp(p, w, check)
%%% sparse 3D arrays
sparse3Darray(i,j,k,x,dims,dimnames,strict,nonzero)
as.sparse3Darray(x, \dots)
\method{dim}{sparse3Darray}(x)
\method{dim}{sparse3Darray}(x) <- value
\method{dimnames}{sparse3Darray}(x)
\method{dimnames}{sparse3Darray}(x) <- value
\method{print}{sparse3Darray}(x, \dots)
\method{aperm}{sparse3Darray}(a, perm, resize, \dots)
\method{as.array}{sparse3Darray}(x, \dots)
\method{[}{sparse3Darray}(x, i, j, k, drop, \dots)
\method{[}{sparse3Darray}(x, i, j, k, \dots) <- value
\method{anyNA}{sparse3Darray}(x, recursive)
RelevantZero(x)
RelevantEmpty(x)
isRelevantZero(x)
unionOfSparseIndices(A,B)
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Summary(\dots, na.rm = FALSE)}
%NAMESPACE S3method("Math", "sparse3Darray")
%NAMESPACE S3method("Ops", "sparse3Darray")
%NAMESPACE S3method("Summary", "sparse3Darray")
inside3Darray(d, i, j, k)
SparseEntries(x)
SparseIndices(x)
EntriesToSparse(df, dims)
mapSparseEntries(x, margin, values, conform, across)
applySparseEntries(x, f, \dots)
sumsymouterSparse(x, w, dbg)
tenseur(A, B, alongA, alongB)
marginSums(X, MARGIN)
rbindCompatibleDataFrames(x)
bind.sparse3Darray(A, B, along)
%%
spatstatDiagnostic(msg)
%%
as.ppplist(x, check)
as.imlist(x, check)
pointsAlongNetwork(L, delta)
expandSparse(x, n, across)
allElementsIdentical(x, entry)
resampleNetworkDataFrame(df, template)
sparseVectorCumul(x, i, length)
%%
\method{as.ppm}{lppm}(object)
\method{as.ppm}{rppm}(object)
\method{predict}{profilepl}(object, \dots)
%%%%%%%
}
\details{
  These internal \pkg{spatstat} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat} to the next.
}
\keyword{internal}
