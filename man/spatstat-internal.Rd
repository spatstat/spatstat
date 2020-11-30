\name{spatstat-internal} 
\title{Internal spatstat functions}
\alias{[.diagramobj}
\alias{[.pp3}
\alias{[.rat}
\alias{[.splitppx}
\alias{[<-.splitppx}
\alias{acedist.show}
\alias{acedist.noshow}
\alias{adjust.ratfv}
\alias{affinexy}
\alias{affinexypolygon}
\alias{allElementsIdentical}
\alias{anycrossing.psp}
\alias{applytolayers}
\alias{applyPolyclipArgs}
\alias{assemble.plot.objects}       
\alias{AsymmDistance.psp}
\alias{as.breakpts}
\alias{as.character.unitname}
\alias{as.data.frame.fv}
\alias{as.double.im}
\alias{as.imlist}
\alias{as.linimlist}
\alias{as.list.hyperframe}
\alias{as.listof}
\alias{as.ppplist}
\alias{as.unitname}
\alias{avenndist}
\alias{bbEngine}
\alias{bdry.mask}
\alias{bind.ratfv}
\alias{bounding.box3}
\alias{break.holes}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{cartesian}
\alias{cellmiddles}
\alias{checkbigmatrix}          
\alias{checkfields}          
\alias{check.arc}
\alias{check.finespacing}
\alias{check.hist.lengths}
\alias{circunion}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closethresh}
\alias{cocoEngine}
\alias{coerce.marks.numeric}
\alias{compatible.rat}
\alias{conform.ratfv}
\alias{crosspairquad}
\alias{cobble.xy}
\alias{codetime}
\alias{col.args.to.grey}
\alias{commonPolyclipArgs}
\alias{conform.imagelist}
\alias{countingweights}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{dimnames.hyperframe}
\alias{dimnames<-.hyperframe}
\alias{dflt.redraw}
\alias{diagramobj}
\alias{dim.fasp}               
\alias{dim.hyperframe}               
\alias{dim.im}               
\alias{dim.owin}
\alias{dimnames.fasp}               
\alias{dimnames<-.fasp}
\alias{distributecbind}
\alias{dist2dpath}
\alias{do.as.im}
\alias{do.call.plotfun}
\alias{drawSignedPoly}
\alias{emptywindow}
\alias{equalpairs.quad}          
\alias{equals.quad}          
\alias{equalsfun.quad}          
\alias{erodemask}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{existsSpatstatVariable}
\alias{expandSpecialLists}
\alias{fakemaintitle}
\alias{fasp}
\alias{fft2D}
\alias{fftwAvailable}
\alias{findbestlegendpos}
\alias{findcbind}
\alias{fillNA}
\alias{flatfname}
\alias{flipxypolygon}
\alias{format.numberwithunit}
\alias{FormatFaspFormulae}
\alias{framebottomleft}
\alias{fvexprmap}
\alias{fvlabels}
\alias{fvlabels<-}
\alias{fvlabelmap}
\alias{fvlegend}
\alias{gammabreaks}
\alias{genericNNdistBy}
\alias{getfields}
\alias{getlastshift}
\alias{getSpatstatVariable}
\alias{gridindex}            
\alias{grid1index}
\alias{grow.mask}
\alias{handle.r.b.args}
\alias{head.hyperframe}
\alias{hsvNA}
\alias{IdenticalRowPair}
\alias{IdenticalRows}
\alias{idorempty}
\alias{imageOp}
\alias{imagelistOp}
\alias{inpoint}
\alias{inside.arc}
\alias{interpretAsOrigin}
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}
\alias{invokeColourmapRule}
\alias{is.col.argname}
\alias{is.data}
\alias{is.fv}
\alias{is.hyperframe}
\alias{is.imlist}
\alias{is.infline}
\alias{is.marked.default}    
\alias{is.marked.psp}    
\alias{is.marked.quad}
\alias{is.multitype.quad}    
\alias{is.multitype.default}    
\alias{is.ppplist}
\alias{is.pp3}
\alias{is.ppx}
\alias{is.psp}
\alias{is.quad}
\alias{is.sob}
\alias{is.solist}
\alias{is.tess}
\alias{is.vanilla}
\alias{levels.im}
\alias{levels<-.im}
\alias{levelsAsFactor}
\alias{listof}
\alias{logi.dummy}
\alias{markvaluetype}
\alias{maskLaslett}
\alias{packupNNdata}
\alias{parbreak}
\alias{plan.legend.layout}
\alias{polyLaslett}
\alias{polytileareaEngine}
\alias{putSpatstatVariable}
\alias{lookup.im}
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
\alias{maxflow}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{MinimalTess}
\alias{multiplicityNumeric}
\alias{multiply.only.finite.entries}
\alias{na.handle.im}
\alias{names.hyperframe}
\alias{names<-.fv}
\alias{names<-.hyperframe}
\alias{nearest.pixel}
\alias{nearest.valid.pixel}
\alias{n.quad}
\alias{numberwithunit}
\alias{numeric.columns}
\alias{onecolumn}
\alias{owinpolycheck}
\alias{owinpoly2mask}
\alias{owin2polypath}
\alias{param.quad}
\alias{pickoption}
\alias{plotEachLayer}
\alias{plot3Dpoints}
\alias{plotPolygonBdry}
\alias{plotWidthMap}
\alias{plot.barplotdata}
\alias{plot.indicfun}
\alias{ppllengine}
\alias{pppdist.mat}
\alias{pppdist.prohorov}
\alias{ppsubset}
\alias{prefixfv}
\alias{prepareTitle}
\alias{print.anylist}    
\alias{print.colourmap}
\alias{print.distfun}
\alias{print.ewcdf}
\alias{print.fasp}       
\alias{print.funxy}       
\alias{print.fv}       
\alias{print.fvfun}       
\alias{print.hyperframe}
\alias{print.indicfun}       
\alias{print.laslett}
\alias{print.layered}
\alias{print.lut}
\alias{print.nnfun}
\alias{print.numberwithunit}
\alias{print.onearrow}
\alias{print.pppmatching}
\alias{print.rat}
\alias{print.simplepanel}
\alias{print.solist}
\alias{print.splitppp}
\alias{print.splitppx}
\alias{print.summary.distfun}
\alias{print.summary.funxy}
\alias{print.summary.hyperframe}
\alias{print.summary.listof}
\alias{print.summary.logiquad}
\alias{print.summary.lut}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{print.summary.solist}
\alias{print.summary.splitppp}
\alias{print.summary.splitppx}
\alias{print.summary.unitname}
\alias{print.symbolmap}       
\alias{print.textstring}
\alias{print.texturemap}
\alias{print.tess}
\alias{print.timed}
\alias{print.yardstick}
\alias{project3Dhom}
\alias{putlastshift}
\alias{qtPrepareCoordinate}
\alias{quad}
\alias{rasterfilter}
\alias{ratfv}
\alias{recognise.spatstat.type}
\alias{rectquadrat.breaks}
\alias{rectquadrat.countEngine}
\alias{repair.image.xycoords}
\alias{resolve.stringsAsFactors}
\alias{rgbNA}
\alias{ruletextline}
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
\alias{rebadgeLabels}
\alias{rebound}
\alias{rebound.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{reconcile.fv}
\alias{rename.fv}
\alias{repair.old.factor.image}
\alias{restrict.mask}
\alias{reversePolyclipArgs}
\alias{rotxy}
\alias{rotxypolygon}
\alias{row.names.hyperframe}
\alias{row.names<-.hyperframe}
\alias{runifrect}
\alias{safedeldir}
\alias{scalardilate.breakpts}
\alias{scalardilate.diagramobj}
\alias{shift.diagramobj}              
\alias{shift.quadratcount}              
\alias{shiftxy}              
\alias{shiftxypolygon}              
\alias{simulationresult}
\alias{sort.im}
\alias{spatstatDiagnostic}
\alias{spatstat.deldir.setopt}
\alias{spatstat.xy.coords}
\alias{store.versionstring.spatstat}
\alias{str.hyperframe}
\alias{summary.hyperframe}
\alias{summary.logiquad}
\alias{summary.lut}
\alias{summary.pppmatching}
\alias{summary.ppx}
\alias{summary.splitppx}
\alias{superimposeMarks}
\alias{symbolmapdomain}
\alias{symbolmaptype}
\alias{tail.hyperframe}
\alias{thickSegments}
\alias{tilecentroids}        
\alias{trianglediameters}
\alias{trim.mask}
\alias{tweak.closepairs}
\alias{tweak.fv.entry}
\alias{tweak.ratfv.entry}
\alias{\%unit\%} %DoNotExport
%NAMESPACE export("%unit%")
\alias{unitname.default}
\alias{unitname<-.default}
\alias{unstackFilter}
\alias{update.im}
\alias{validradius}
\alias{validate.mask}        
\alias{validate.quad}        
\alias{vanilla.fv}
\alias{versioncurrency.spatstat}
\alias{versionstring.spatstat}
\alias{verifyclass}
\alias{veryunique}
\alias{warn.once}
\alias{waxlyrical}
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xy.grid}
\alias{xtfrm.im}
\alias{xypolygon2psp}
\alias{xypolyselfint}
\alias{ZeroValue}
\alias{ZeroValue.im}
\alias{[.localpcfmatrix}
\alias{accumulateStatus}
\alias{active.interactions}
\alias{adaptcoef}
\alias{affine.msr}
\alias{ang2rad}
\alias{areadelta2}
\alias{areaGain.diri}       
\alias{areaGain.grid}       
\alias{areaLoss.diri}       
\alias{areaLoss.grid}
\alias{areaLoss.poly}
\alias{as.data.frame.bw.optim}
\alias{augment.msr}
\alias{bandwidth.is.infinite}
\alias{BartCalc}
\alias{bermantestCalc}
\alias{bermantestEngine}
\alias{bigvaluerule}
\alias{blankcoefnames}
\alias{bt.frame}
\alias{bw.optim}
\alias{calc.DR}
\alias{calc.NNIR}
\alias{calc.SAVE}
\alias{calc.SIR}
\alias{calc.TSE}
\alias{cannot.update}
\alias{censtimeCDFest}
\alias{change.default.expand}          
\alias{check.separable}
\alias{check.testfun}
\alias{circticks}
\alias{clarkevansCalc}
\alias{coef.summary.kppm}
\alias{coef.summary.ppm}
\alias{coef.vblogit} 
\alias{compileCDF}
\alias{CressieReadStatistic}
\alias{CressieReadSymbol}
\alias{CressieReadName}
\alias{cutoff2Dkernel}
\alias{CVforPCF}
\alias{damaged.ppm}
\alias{data.mppm}
\alias{datagen.runifpointOnLines}
\alias{datagen.runifpoisppOnLines}
\alias{datagen.rpoisppOnLines}
\alias{default.clipwindow}
\alias{deltasuffstat}
\alias{Deviation}
\alias{dfbetas.ppmInfluence}
\alias{densitycrossEngine}
\alias{densitypointsEngine}
\alias{diagnose.ppm.engine}
\alias{digestCovariates}
\alias{digital.volume}
\alias{dim.msr}
\alias{dimnames.msr}
\alias{doMultiStraussHard}
\alias{dppmFixAlgorithm}
\alias{dppmFixIntensity}
\alias{envelopeEngine}
\alias{envelopeProgressData}
\alias{envelopeTest}
\alias{envelope.hasenvelope}
\alias{envelope.matrix}
\alias{equalpairs}          
\alias{evalCovar}
\alias{evalCovar.ppm}
\alias{evalCovariate}
\alias{evalInteraction}
\alias{evalInterEngine}
\alias{evalPairPotential}
\alias{evaluate2Dkernel}
\alias{expandwinPerfect}
\alias{ExpSmoothLog}
\alias{extractAIC.slrm}
\alias{extractAtomicQtests}
\alias{fakeNeyScot}
\alias{family.vblogit} 
\alias{f3engine}
\alias{f3Cengine}
\alias{fill.coefs}
\alias{findCovariate}
\alias{fii}
\alias{flipxy.msr}
\alias{forbid.logi}
\alias{g3engine}
\alias{g3Cengine}
\alias{getdataname}
\alias{getglmdata}
\alias{getglmfit}
\alias{getglmsubset}
\alias{getppmdatasubset}
\alias{getppmOriginalCovariates}
\alias{getRandomFieldsModelGen}
\alias{getSumFun}
\alias{geyercounts}
\alias{geyerdelta2}
\alias{GLMpredict}
\alias{good.correction.K}
%\alias{gridadjacencymatrix} %DoNotExport
\alias{hackglmmPQL}
\alias{hasenvelope}
\alias{hasglmfit}
\alias{HermiteCoefs}
\alias{handle.rshift.args}
\alias{hierarchicalordering}
\alias{hiermat}
\alias{ho.engine}
\alias{illegal.iformula}
\alias{implemented.for.K}
\alias{impliedpresence}
\alias{impliedcoefficients}
\alias{influence.ppmInfluence}
\alias{instantiate.interact}
\alias{interactionfamilyname}
\alias{intermaker}
\alias{is.atomicQtest}
\alias{is.cadlag}
\alias{is.expandable}
\alias{is.expandable.ppm}
\alias{is.expandable.rmhmodel}
\alias{is.interact}
\alias{is.marked.mppm}    
\alias{is.marked.msr}    
\alias{is.mppm}
\alias{is.multitype.mppm}    
\alias{is.multitype.msr}    
\alias{is.poisson.mppm}
\alias{is.scov}
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
\alias{leverage.ppmInfluence}
\alias{localKengine}
\alias{localKmultiEngine}
\alias{localpcfengine}
\alias{localpcfmatrix}
\alias{logi.engine}
\alias{logLik.vblogit}
\alias{lookup2DkernelInfo}
\alias{match2DkernelName}
\alias{match.kernel}
\alias{mctestSigtraceEngine}
\alias{meanlistfv}
\alias{model.se.image}
\alias{modelFrameGam}
\alias{mpl.engine}
\alias{mpl.get.covariates}
\alias{mpl.prepare}
\alias{mpl.usable}
\alias{MultiPair.checkmatrix}
\alias{newformula}
\alias{newstyle.coeff.handling}
\alias{nncleanEngine}
\alias{nndcumfun}
\alias{no.trend.ppm}
\alias{objsurfEngine}
\alias{optimConverged}
\alias{optimStatus}
\alias{outdated.interact}
\alias{oversize.quad}    
\alias{pairs.listof}
\alias{pairs.solist}
\alias{PairPotentialType}
\alias{partialModelMatrix}
\alias{pcf3engine}
\alias{pcfmulti.inhom}
\alias{ploterodewin}
\alias{ploterodeimage}
\alias{plot.addvar}
\alias{plot.bw.frac}
\alias{plot.bw.optim}
\alias{plot.localpcfmatrix}
\alias{plot.lurk}
\alias{plot.minconfit}
\alias{plot.parres}
\alias{plot.plotpairsim}
\alias{plot.qqppm}
\alias{plot.spatialcdf}
\alias{pointweights}
\alias{PoisSaddle}
\alias{PoisSaddleArea}
\alias{PoisSaddleGeyer}
\alias{PoisSaddlePairwise}
\alias{PPMmodelmatrix}
\alias{ppm.default}
\alias{ppmCovariates}
\alias{ppmDerivatives}
\alias{ppmInfluenceEngine}
\alias{predict.profilepl}
\alias{predict.vblogit}
\alias{printStatus}
\alias{printStatusList}
\alias{print.addvar}    
\alias{print.autoexec}    
\alias{print.bt.frame}
\alias{print.bugtable}
\alias{print.bw.frac}
\alias{print.bw.optim}
\alias{print.diagppm}
\alias{print.densityfun}
\alias{print.detpointprocfamily}
\alias{print.detpointprocfamilyfun}
\alias{print.envelope}
\alias{print.hasenvelope}       
\alias{print.hierarchicalordering}
\alias{print.influence.ppm}       
\alias{print.interact}       
\alias{print.intermaker}       
\alias{print.isf}
\alias{print.leverage.ppm}
\alias{print.localpcfmatrix}
\alias{print.lurk}
\alias{print.minconfit}
\alias{print.mppm}
\alias{print.msr}
\alias{print.parres}
\alias{print.plotpairsim}
\alias{print.plotppm}
\alias{print.profilepl}
\alias{print.quadrattest}
\alias{print.qqppm}
\alias{print.rmhcontrol}
\alias{print.rmhexpand}
\alias{print.rmhmodel}
\alias{print.rmhstart}
\alias{print.rmhInfoList}
\alias{print.rppm}
\alias{print.Smoothfun}       
\alias{print.summary.mppm}
\alias{print.summary.rmhexpand}
\alias{print.summary.ssf}
\alias{print.vblogit}
\alias{quad.mppm}
\alias{quadBlockSizes}
\alias{quadrat.testEngine}
\alias{RandomFieldsSafe}
\alias{reduceformula}
\alias{reheat}
\alias{reincarnate.interact}
\alias{RelevantDeviation}
\alias{rescale.msr}
\alias{resid4plot}
\alias{resid1plot}
\alias{resid1panel}
\alias{resolve.2D.kernel}
\alias{resolveEinfo}
\alias{resolve.foxall.window}
\alias{resolve.lambda}
\alias{resolve.lambda.cross}
\alias{resolve.vargamma.shape}
\alias{rhohatEngine}
\alias{rhohatCalc}
\alias{rMaternInhibition}
\alias{rmax.Rigid}
\alias{rmax.rule}
\alias{RmhExpandRule}
\alias{rmhsnoop}
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
\alias{rocData}
\alias{rocModel}
\alias{roseContinuous}
\alias{rotate.msr}
\alias{rpoint.multi}
\alias{runifpoispp}          
\alias{runifpoisppOnLines}          
\alias{safelookup}
\alias{safevalue}
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
\alias{shift.influence.ppm}              
\alias{shift.leverage.ppm}              
\alias{shift.quadrattest}              
\alias{shift.msr}              
\alias{signalStatus}
\alias{simulate.profilepl}
\alias{simulrecipe}              
\alias{slr.prepare}
\alias{slrAssemblePixelData}
\alias{Smooth.solist}
\alias{smoothcrossEngine}              
\alias{smoothpointsEngine}
\alias{smudge}
\alias{spatstatClusterModelInfo}
\alias{spatstatDPPModelInfo}
\alias{spatstatRmhInfo}
\alias{spatialCDFframe}
\alias{spatialCDFtest}
\alias{spatialCDFtestCalc}
\alias{sphere.volume}
\alias{splitHybridInteraction}
\alias{sp.foundclass}
\alias{sp.foundclasses}
\alias{strausscounts}
\alias{suffloc}
\alias{suffstat.generic}
\alias{suffstat.poisson}
\alias{summarise.trend}
\alias{summary.envelope}
\alias{summary.mppm}
\alias{summary.msr}
\alias{summary.profilepl}
\alias{summary.rmhexpand}
\alias{summary.vblogit}
\alias{thinjump}
\alias{tweak.coefs}
\alias{twostage.test}
\alias{twostage.envelope}
\alias{unitname.msr}
\alias{unitname<-.msr}
\alias{update.ippm}
\alias{update.msr}
\alias{update.rmhstart}
\alias{validate2Dkernel}
\alias{validate.angles}        
\alias{validate.weights}        
\alias{varcountEngine}
%\alias{vblogit} %DoNotExport
%\alias{vblogit.fmla} %DoNotExport
\alias{versionstring.interact}
\alias{versionstring.ppm}
\alias{weightedclosepairs}
\alias{windows.mppm}
\alias{X2testEngine}

%% Linear networks code is corralled here
\alias{ApplyConnected}
\alias{DoCountCrossEnds}
\alias{DoCountEnds}
\alias{FDMKERNEL}
\alias{as.linfun.linfun}
\alias{as.owin.lintess}
\alias{default.linnet.tolerance}
\alias{makeLinnetTolerance}
\alias{print.lintess}
\alias{print.summary.linim}
\alias{print.summary.linnet}
\alias{print.summary.lintess}
\alias{resolve.heat.steps}
\alias{summary.lintess}
\alias{Window.lintess}
\alias{Window<-.linnet}
\alias{Window<-.lpp}
\alias{densitypointsLPP}
\alias{evalCovar.lppm}
\alias{flatdensityfunlpp}
\alias{flatdensityatpointslpp}
\alias{getlambda.lpp}
\alias{local2lpp}
\alias{looHeatLPP}
\alias{looVoronoiLPP}
\alias{validate.lpp.coords}
\alias{as.ppm.lppm}
\alias{as.ppm.rppm}
\alias{pointsAlongNetwork}
\alias{linearKengine}
\alias{linearKmulti}
\alias{linearKmulti.inhom}
\alias{linearKmultiEngine}
\alias{linearpcfengine}
\alias{linearpcfmulti}
\alias{linearpcfmulti.inhom}
\alias{linearPCFmultiEngine}
\alias{resampleNetworkDataFrame}
\alias{sortalongsegment}
\alias{vnnFind}
\alias{ldtEngine}
\alias{qkdeEngine}
\alias{Math.linimlist}
\alias{Ops.linimlist}
\alias{Summary.linimlist}
\alias{Complex.linimlist}
\alias{LinimOp}
\alias{LinimListOp}

%%%%%%%
\description{
  Internal spatstat functions.
}
\usage{
\method{[}{diagramobj}(x, \dots)
\method{[}{pp3}(x, i, drop, \dots)
\method{[}{rat}(x, \dots)
\method{[}{splitppx}(x, \dots)
\method{[}{splitppx}(x, \dots) <- value
acedist.show(X, Y, n, d, timelag)
acedist.noshow(X, Y, n, d)
adjust.ratfv(f, columns, numfactor, denfactor)
affinexy(X, mat, vec, invert)
affinexypolygon(p, mat, vec, detmat)
allElementsIdentical(x, entry)
anycrossing.psp(A,B)
applytolayers(L, FUN, \dots)
applyPolyclipArgs(x, p)
assemble.plot.objects(xlim, ylim, \dots, lines, polygon)
AsymmDistance.psp(X, Y, metric, method)
as.breakpts(\dots)
\method{as.character}{unitname}(x, \dots)
\method{as.data.frame}{fv}(x, \dots)
\method{as.double}{im}(x, \dots)
as.imlist(x, check)
as.linimlist(x, check)
\method{as.list}{hyperframe}(x, \dots)
as.listof(x)
as.ppplist(x, check)
as.unitname(s)
avenndist(X)
bbEngine(\dots)
bdry.mask(W)
bind.ratfv(x, numerator, denominator, labl, desc, preferred,
           ratio, quotient)
bounding.box3(\dots)
break.holes(x, splitby, depth, maxdepth)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
cartesian(pp, markset, fac = TRUE)
cellmiddles(W, nx, ny, npix, distances)
checkbigmatrix(n, m, fatal, silent)
checkfields(X,L)
check.arc(arc, fatal)
check.finespacing(r, eps, win, rmaxdefault, context, action, rname)
check.hist.lengths(hist,breaks)
circunion(arcs)
clip.psp(x, window, check, fragments)
cliprect.psp(x, window, fragments)
clippoly.psp(s, window, fragments)
closethresh(X,R,S,twice,\dots)
cocoEngine(nv, ie, je, algoname)
coerce.marks.numeric(X, warn)
\method{compatible}{rat}(A, B, \dots) 
conform.ratfv(x)
crosspairquad(Q,rmax,what)
cobble.xy(x, y, f, fatal, \dots)
codetime(x, hms, what)
col.args.to.grey(x, \dots)
commonPolyclipArgs(\dots, p)
conform.imagelist(X, Zlist)
countingweights(id, areas, check = TRUE)
default.n.tiling(X, nd, ntile, npix, eps, random, quasi, verbose)
default.ntile(X)
\method{dimnames}{hyperframe}(x)
\method{dimnames}{hyperframe}(x) <- value
dflt.redraw(button, name, env)
diagramobj(X, \dots)
\method{dim}{fasp}(x)
\method{dim}{hyperframe}(x)
\method{dim}{im}(x)
\method{dim}{owin}(x)
\method{dimnames}{fasp}(x)
\method{dimnames}{fasp}(x) <- value
distributecbind(x)
dist2dpath(dist, method="C")
do.as.im(x, action, \dots, W, eps, dimyx, xy, na.replace)
do.call.plotfun(fun, arglist, \dots, envir)
drawSignedPoly(x,y,pars,sgn)
emptywindow(w)
equalpairs.quad(Q)
equals.quad(Q)          
equalsfun.quad(Q)
erodemask(w,r,strict)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(w)
existsSpatstatVariable(name)
expandSpecialLists(x, special)
fakemaintitle(bb, main, \dots)
fasp(fns, which, formulae, dataname, title, rowNames, colNames, checkfv)
fft2D(z, inverse, west)
fftwAvailable()
findbestlegendpos(\dots)
findcbind(root, depth, maxdepth)
fillNA(x, value)
flatfname(x)
flipxypolygon(p)
\method{format}{numberwithunit}(x, \dots, collapse, modifier)
FormatFaspFormulae(f, argname)
framebottomleft(w)
fvexprmap(x)
fvlabels(x, expand=FALSE)
fvlabels(x) <- value
fvlabelmap(x, dot=TRUE)
fvlegend(object, elang)
gammabreaks(ra, n, gamma)
genericNNdistBy(X, by, k)
getfields(X, L, fatal = TRUE)
getlastshift(X)
getSpatstatVariable(name)
gridindex(x, y, xrange, yrange, nx, ny)            
grid1index(x, xrange, nx)
grow.mask(M, xmargin=0, ymargin=xmargin)
handle.r.b.args(r = NULL, breaks = NULL, window, pixeps = NULL, rmaxdefault)
\method{head}{hyperframe}(x,n,\dots)
hsvNA(h, s, v, alpha)
IdenticalRowPair(i,j,a,b)
IdenticalRows(i,j,a,b)
idorempty(w, r, caller)
imageOp(e1, e2, op)
imagelistOp(e1, e2, op)
inpoint(W)
inside.arc(theta, arc)
interpretAsOrigin(x, W)
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
invokeColourmapRule(colfun, x, \dots, zlim, colargs)
is.col.argname(x)
is.data(Q)
is.fv(x)
is.hyperframe(x)
is.imlist(x)
is.infline(x)
\method{is.marked}{default}(\dots)  
\method{is.marked}{psp}(X, \dots)
\method{is.marked}{quad}(X, na.action="warn", \dots)
\method{is.multitype}{default}(X, \dots)  
\method{is.multitype}{quad}(X, na.action="warn", \dots)
is.ppplist(x)
is.pp3(x)
is.ppx(x)
is.psp(x)
is.quad(x)
is.solist(x)
is.sob(x)
is.tess(x)
is.vanilla(u)
\method{levels}{im}(x)
\method{levels}{im}(x) <- value
levelsAsFactor(x)
listof(\dots)
logi.dummy(X, dummytype, nd, mark.repeat, \dots)
maskLaslett(X, \dots, eps, dimyx, xy, oldX, verbose, plotit)
markvaluetype(x)
packupNNdata(NND, NNW, what, k)
parbreak(terse)
plan.legend.layout(B, \dots, side, sep, size, sep.frac, size.frac,
                   started, map)
polyLaslett(X, \dots, oldX, verbose, plotit)
polytileareaEngine(P, xrange, yrange, nx, ny, DivideByPixelArea)
putSpatstatVariable(name, value)
lookup.im(Z, x, y, naok, strict)
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
maxflow(costm)
meanX.owin(w)            
meanY.owin(w)
MinimalTess(W, \dots)
multiplicityNumeric(x)
multiply.only.finite.entries(x, a)
na.handle.im(X, na.replace)
\method{names}{fv}(x) <- value
\method{names}{hyperframe}(x)
\method{names}{hyperframe}(x) <- value
nearest.pixel(x, y, Z)
nearest.valid.pixel(x, y, Z, method, nsearch)
n.quad(Q)
numberwithunit(x, u)
numeric.columns(M, logical, others)
onecolumn(m)
owinpolycheck(W, verbose=TRUE)
owinpoly2mask(w, rasta, check=TRUE)
owin2polypath(w)
param.quad(Q)
pickoption(what="option", key, keymap, \dots,
           exact=FALSE, list.on.err=TRUE, die=TRUE, multi=FALSE,
           allow.all=TRUE)
plotEachLayer(x, \dots, main, plotargs, add, show.all, do.plot)
plot3Dpoints(xyz, eye, org,
             \dots,
             type, xlim, ylim, zlim,
             add, box, main, cex, box.back, box.front)
plotPolygonBdry(x, \dots)
plotWidthMap(bb.leg, zlim, phys.scale, leg.scale, leg.side, leg.args, grafpar)
\method{plot}{barplotdata}(x, \dots)
\method{plot}{indicfun}(x, W, \dots, main)
ppllengine(X, Y, action="project", check=FALSE)
pppdist.mat(X, Y, cutoff = 1, q = 1, matching = TRUE,
            precision = 9, approximation = 10)
pppdist.prohorov(X, Y, n, dfix, type, cutoff, matching,
            ccode, auction, precision, approximation) 
ppsubset(X, I, Iname, fatal)
prefixfv(x, tagprefix, descprefix, lablprefix, whichtags)
prepareTitle(main)
\method{print}{anylist}(x, \dots)
\method{print}{colourmap}(x, \dots)
\method{print}{distfun}(x, \dots)
\method{print}{ewcdf}(x, digits, \dots)
\method{print}{fasp}(x, \dots)
\method{print}{funxy}(x, \dots)
\method{print}{fv}(x, \dots, tight)
\method{print}{fvfun}(x, \dots)
\method{print}{hyperframe}(x, \dots)
\method{print}{indicfun}(x, \dots)
\method{print}{laslett}(x, \dots)
\method{print}{layered}(x, \dots)
\method{print}{lut}(x, \dots)
\method{print}{nnfun}(x, \dots)
\method{print}{numberwithunit}(x, \dots)
\method{print}{onearrow}(x, \dots)
\method{print}{pppmatching}(x, \dots)
\method{print}{rat}(x, \dots)
\method{print}{simplepanel}(x, \dots)
\method{print}{solist}(x, \dots)
\method{print}{splitppp}(x, \dots)
\method{print}{splitppx}(x, \dots)
\method{print}{summary.distfun}(x, \dots)
\method{print}{summary.funxy}(x, \dots)
\method{print}{summary.hyperframe}(x, \dots)
\method{print}{summary.listof}(x, \dots)
\method{print}{summary.logiquad}(x, \dots, dp=3)
\method{print}{summary.lut}(x, \dots)
\method{print}{summary.owin}(x, \dots)
\method{print}{summary.ppp}(x, \dots, dp)
\method{print}{summary.psp}(x, \dots)
\method{print}{summary.splitppp}(x, \dots)
\method{print}{summary.solist}(x, \dots)
\method{print}{summary.splitppx}(x, \dots)
\method{print}{summary.unitname}(x, \dots)
\method{print}{symbolmap}(x, \dots)
\method{print}{textstring}(x, \dots)
\method{print}{texturemap}(x, \dots)
\method{print}{tess}(x, \dots, brief=FALSE)
\method{print}{timed}(x, \dots)
\method{print}{yardstick}(x, \dots)
project3Dhom(xyz, eye, org, vert)
putlastshift(X, vec)
qtPrepareCoordinate(covname, W, origin)
quad(data, dummy, w, param)
rasterfilter(X, f)
ratfv(df, numer, denom, \dots, ratio)
recognise.spatstat.type(x)
rectquadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
rectquadrat.countEngine(x, y, xbreaks, ybreaks, weights)
repair.image.xycoords(x)
resolve.stringsAsFactors(stringsAsFactors)
rgbNA(red, green, blue, alpha, maxColorValue)
ruletextline(ch, n, terse)
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
rebadgeLabels(x, new.fname)
rebound(x, rect)
\method{rebound}{im}(x, rect)  
\method{rebound}{ppp}(x, rect) 
\method{rebound}{psp}(x, rect) 
\method{rebound}{owin}(x, rect)
reconcile.fv(\dots)
rename.fv(x, fname, ylab, yexp)
repair.old.factor.image(x)
restrict.mask(M, W)
reversePolyclipArgs(x, p)
rotxy(X, angle = pi/2)
rotxypolygon(p, angle = pi/2)
\method{row.names}{hyperframe}(x)
\method{row.names}{hyperframe}(x) <- value
runifrect(n, win, nsim, drop)
safedeldir(X)
\method{scalardilate}{breakpts}(X, f, \dots)
\method{scalardilate}{diagramobj}(X, f, \dots)
\method{shift}{diagramobj}(X, \dots)
\method{shift}{quadratcount}(X, \dots)
shiftxy(X, vec = c(0, 0))
shiftxypolygon(p, vec = c(0, 0))
simulationresult(resultlist, nsim, drop, NameBase)
\method{sort}{im}(x, \dots)
spatstatDiagnostic(msg)
spatstat.deldir.setopt(use.trigrafS, use.trigraf, debug.delaunay)
spatstat.xy.coords(x, y)
store.versionstring.spatstat()
\method{str}{hyperframe}(object, \dots)
\method{summary}{hyperframe}(object, \dots, brief=FALSE)
\method{summary}{logiquad}(object, \dots, checkdup=FALSE)
\method{summary}{lut}(object, \dots)
\method{summary}{pppmatching}(object, \dots)
\method{summary}{ppx}(object, \dots)
\method{summary}{splitppx}(object, \dots)
superimposeMarks(arglist, nobj)
symbolmapdomain(x)
symbolmaptype(x)
\method{tail}{hyperframe}(x,n,\dots)
thickSegments(x, widths, \dots,
             add, main, do.plot, show.all, show.window,
             scale, adjust, negative.args,
             legend, leg.side, leg.sep, leg.wid, leg.args, leg.scale,
             zlim, box)
tilecentroids(W, nx, ny)
trianglediameters(iedge, jedge, edgelength, \dots,
                  nvert, dmax, check)
trim.mask(M, R, tolerant)
tweak.closepairs(cl, rmax, i, deltax, deltay, deltaz)
tweak.fv.entry(x, current.tag, new.labl, new.desc, new.tag)
tweak.ratfv.entry(x, \dots)
x \%unit\% u
\method{unitname}{default}(x) 
\method{unitname}{default}(x) <- value 
unstackFilter(x)
\method{update}{im}(object, \dots) 
validradius(r, caller)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
vanilla.fv(x)
versioncurrency.spatstat(today, checkR)
versionstring.spatstat()
veryunique(z)
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
warn.once(key, \dots)
waxlyrical(type, terse)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
xy.grid(xr, yr, nx, ny, dx, dy)
\method{xtfrm}{im}(x)
xypolyselfint(p, eps, proper, yesorno, checkinternal)
xypolygon2psp(p, w, check)
ZeroValue(x)
\method{ZeroValue}{im}
\method{[}{localpcfmatrix}(x, i, \dots)
accumulateStatus(x, stats)
active.interactions(object)
adaptcoef(new.coef, fitcoef, drop)
\method{affine}{msr}(X, mat, vec, \dots)
ang2rad(ang, unit, start, clockwise)
areadelta2(X, r, \dots, sparseOK)
areaGain.diri(u, X, r, \dots, W, verbose)
areaGain.grid(u, X, r, \dots, W, ngrid)
areaLoss.diri(X, r, \dots, W, subset)
areaLoss.grid(X, r, \dots, W, subset,
                         method = c("count", "distmap"),
                         ngrid = spatstat.options("ngrid.disc"),
                         exact = FALSE)
areaLoss.poly(X, r, \dots, W, subset, splitem)
\method{as.data.frame}{bw.optim}(x, \dots)
augment.msr(x, \dots, sigma, recompute)
bandwidth.is.infinite(sigma)
BartCalc(fY, fK)
bermantestCalc(fram, which, alternative, \dots)
bermantestEngine(model, covariate, which, alternative, \dots,
                 modelname, covname, dataname)
bigvaluerule(objfun, objargs, startpar, \dots)
blankcoefnames(x)
bt.frame(Q, trend, interaction, \dots, covariates,
         correction, rbord, use.gam, allcovar)
bw.optim(cv, h, iopt, \dots, cvname, hname, criterion,
         optimum, warnextreme, hargnames, unitname)
calc.DR(COV, z, Dim)
calc.NNIR(COV, z, pos, Dim)
calc.SAVE(COV, z, Dim)
calc.SIR(COV, z)
calc.TSE(COV, z, pos, Dim1, Dim2)
cannot.update(\dots)
censtimeCDFest(o, cc, d, breaks, \dots,
     KM, RS, HAN, RAW, han.denom, tt, pmax)
change.default.expand(x, newdefault)
check.separable(dmat, covname, isconstant, fatal)
check.testfun(f, f1, X)
circticks(R, at, unit, start, clockwise, labels)
clarkevansCalc(X, correction, clipregion, working)
\method{coef}{summary.kppm}(object, \dots)
\method{coef}{summary.ppm}(object, \dots)
\method{coef}{vblogit}(object, \dots) 
compileCDF(D, B, r, \dots, han.denom, check)
CressieReadStatistic(OBS,EXP,lambda,normalise,named)
CressieReadSymbol(lambda)
CressieReadName(lambda)
cutoff2Dkernel(kernel, sigma, varcov, \dots, scalekernel, cutoff, fatal)
CVforPCF(bw, stuff)
damaged.ppm(object)
data.mppm(x)
datagen.runifpointOnLines(n, L)
datagen.runifpoisppOnLines(lambda, L)
datagen.rpoisppOnLines(lambda, L, lmax, \dots, check)
default.clipwindow(object, epsilon)
deltasuffstat(model, \dots,
              restrict, dataonly, sparseOK, quadsub,
              force, warn.forced, verbose, use.special)
Deviation(x, ref, leaveout, n, xi)
\method{dfbetas}{ppmInfluence}(model, \dots)
densitycrossEngine(Xdata, Xquery, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    diggle, sorted, cutoff)
densitypointsEngine(x, sigma, \dots,
                    kernel, scalekernel,
                    weights, edge, varcov,
                    leaveoneout, diggle, sorted, spill, cutoff,
                    debug)
diagnose.ppm.engine(object, \dots, type, typename, opt,
                         sigma, rbord, compute.sd, compute.cts,
                         envelope, nsim, nrank,
                         rv, oldstyle, splineargs, verbose)
digestCovariates(\dots, W)
digital.volume(range, nval, vside)
\method{dim}{msr}(x)
\method{dimnames}{msr}(x)
doMultiStraussHard(iradii, hradii, types)
dppmFixIntensity(DPP, lambda, po)
dppmFixAlgorithm(algorithm, changealgorithm, clusters, startpar)
envelopeEngine(X, fun, simul,
           nsim, nrank, \dots, funargs, funYargs,
           verbose, clipdata, 
           transform, global, ginterval, use.theory,
           alternative, scale, clamp,
           savefuns, savepatterns, saveresultof,
           weights,
           nsim2, VARIANCE, nSD,
           Yname,
           maxnerr, rejectNA, silent, maxerr.action,
           internal, cl,
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
            Xname, badXfatal, verbose)
\method{envelope}{hasenvelope}(Y, \dots, Yname)
\method{envelope}{matrix}(Y, \dots, rvals, observed, theory, funX,
  nsim, nsim2, jsim, jsim.mean,
  type, alternative, scale, clamp, csr, use.theory, nrank, ginterval, nSD,
  savefuns, check, Yname, do.pwrong, weights, precomputed, gaveup)
equalpairs(U, X, marked=FALSE)
evalCovar(model, covariate, \dots)
\method{evalCovar}{ppm}(model, covariate, \dots, lambdatype,
          dimyx, eps, interpolate, jitter, modelname, covname, dataname,
          subset) 
evalCovariate(covariate, locations)
evalInteraction(X,P,E,interaction,correction,splitInf,\dots,
                precomputed,savecomputed)
evalInterEngine(X,P,E,interaction,correction,splitInf,\dots,
                Reach,precomputed,savecomputed)
evalPairPotential(X,P,E,pairpot,potpars,R)
evaluate2Dkernel(kernel, x, y, sigma, varcov, \dots, scalekernel)
expandwinPerfect(W, expand, amount)
ExpSmoothLog(X, \dots, at, weights)
\method{extractAIC}{slrm}(fit, scale = 0, k = 2, \dots)
extractAtomicQtests(x)
fakeNeyScot(Y, lambda, win, saveLambda, saveparents)
\method{family}{vblogit}(object, \dots)
f3engine(x, y, z, box, vside, range, nval, correction)
f3Cengine(x, y, z, box, vside, rmax, nrval)
fill.coefs(coefs, required)
findCovariate(covname, scope, scopename=NULL)
fii(interaction, coefs, Vnames, IsOffset, vnameprefix)
\method{flipxy}{msr}(X)
forbid.logi(object)
g3engine(x, y, z, box, rmax, nrval, correction)
g3Cengine(x, y, z, box, rmax, nrval)
getdataname(defaultvalue, \dots, dataname)
getglmdata(object, drop=FALSE)
getglmfit(object)
getglmsubset(object)
getppmdatasubset(object)
getppmOriginalCovariates(object)
getRandomFieldsModelGen(model)
getSumFun(abbreviation, classname, ismarked, fatal)
geyercounts(U,X,r,sat,Xcounts,EqualPairs)
geyerdelta2(X,r,sat,\dots,sparseOK, correction)
GLMpredict(fit, data, coefs, changecoef, type)
good.correction.K(X)
%gridadjacencymatrix(dims)
hackglmmPQL(fixed, random, family, data, correlation, weights,
            control, niter, verbose, subset, \dots, reltol)
hasenvelope(X, E)
hasglmfit(object)
HermiteCoefs(order)
handle.rshift.args(W, \dots, radius, width, height,
                   edge, clip, edgedefault)
hierarchicalordering(i, s)
hiermat(x, h)
ho.engine(model, \dots, nsim, nrmh, start, control, verb)
illegal.iformula(ifmla, itags, dfvarnames)
implemented.for.K(correction, windowtype, explicit)
impliedpresence(tags, formula, df, extranames=character(0))
impliedcoefficients(object, tag)
\method{influence}{ppmInfluence}(model, \dots)
instantiate.interact(x, par)
interactionfamilyname(x)
intermaker(f, blank)
is.atomicQtest(x)
is.cadlag(s)
is.expandable(x)
\method{is.expandable}{ppm}(x)  
\method{is.expandable}{rmhmodel}(x) 
is.interact(x)
\method{is.marked}{mppm}(X, \dots)
\method{is.marked}{msr}(X, \dots)
is.mppm(x)
\method{is.multitype}{mppm}(X, \dots)
\method{is.multitype}{msr}(X, \dots)
\method{is.poisson}{mppm}(x)
is.scov(x)
k3engine(x, y, z, box, rmax, nrval, correction)
Kborder.engine(X, rmax, nr, correction, weights, ratio)
Knone.engine(X, rmax, nr, weights, ratio)
Krect.engine(X, rmax, nr, correction, weights,
             ratio, fname, use.integers)
Kount(dIJ, bI, b, breaks)
Kwtsum(dIJ, bI, wIJ, b, w, breaks, fatal)
Kpcf.kppm(model, what)
killinteraction(model)
km.rs.opt(o, cc, d, breaks, KM, RS)
kppmComLik(X, Xname, po, clusters, control, weightfun, rmax,
           algorithm, DPP, \dots, pint)
kppmMinCon(X, Xname, po, clusters, control, statistic, statargs,
           algorithm, DPP, \dots)
kppmPalmLik(X, Xname, po, clusters, control, weightfun, rmax, 
           algorithm, DPP, \dots, pint)
kraever(package, fatal)
kraeverRandomFields()
\method{labels}{ppm}(object, \dots)
\method{leverage}{ppmInfluence}(model, \dots)
localKengine(X, \dots, wantL, lambda, rmax, correction, verbose, rvalue)
localKmultiEngine(X, from, to,
           lambdaFrom, lambdaTo, \dots, rmax, wantL,
           correction, verbose, rvalue, sigma, varcov,
           lambdaX, update, leaveoneout,
           Iexplain, Jexplain, Ikey, Jkey, miss.update, miss.leave)
localpcfengine(X, \dots, delta, rmax, nr, stoyan, lambda)
localpcfmatrix(X, i, \dots, lambda, delta, rmax, nr, stoyan)
logi.engine(Q, trend, interaction, \dots, 
            covariates, subsetexpr, clipwin,
            correction, rbord, covfunargs, allcovar, 
            vnamebase, vnameprefix, justQ, savecomputed, precomputed,
            VB)
\method{logLik}{vblogit}(object, \dots)
lookup2DkernelInfo(kernel)
match.kernel(kernel)
match2DkernelName(kernel)
mctestSigtraceEngine(R, devdata, devsim, \dots,
     interpolate, confint, alpha, exponent, unitname)
meanlistfv(z, \dots)
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
newformula(old, change, eold, enew, expandpoly)
newstyle.coeff.handling(object)
nncleanEngine(kthNND, k, d, \dots, tol, maxit,
              plothist, lineargs, verbose, Xname)
nndcumfun(X, \dots, r)
no.trend.ppm(x)
objsurfEngine(objfun, optpar, objargs,
              \dots, dotargs, objname,
              ngrid, ratio, verbose)
optimConverged(x)
optimStatus(x, call)
outdated.interact(object)
oversize.quad(Q, \dots, nU, nX, p)
PairPotentialType(pairpot)
\method{pairs}{listof}(\dots, plot=TRUE)
\method{pairs}{solist}(\dots, plot=TRUE)
\method{plot}{localpcfmatrix}(x, \dots)
PoisSaddle(beta, fi)
PoisSaddleArea(beta, fi)
PoisSaddleGeyer(beta, fi)
PoisSaddlePairwise(beta, fi)
\method{print}{localpcfmatrix}(x, \dots)
PPMmodelmatrix(object, data, \dots, subset, Q, keepNA, irregular,
               splitInf)
printStatus(x, errors.only)
printStatusList(stats)
partialModelMatrix(X,D,model,callstring,\dots)
pcf3engine(x, y, z, box, rmax, nrval, correction, delta)
pcfmulti.inhom(X, I, J, lambdaI = NULL, lambdaJ = NULL, \dots,
               r = NULL, breaks = NULL,
               kernel = "epanechnikov", bw = NULL, stoyan = 0.15,
               correction = c("translate", "Ripley"),
               sigma = NULL, varcov = NULL,
               Iname = "points satisfying condition I",
               Jname = "points satisfying condition J")
ploterodewin(W1, W2, col.edge, col.inside, do.plot, \dots)
ploterodeimage(W, Z, \dots, Wcol, rangeZ, colsZ, do.plot)
\method{plot}{addvar}(x, \dots, do.points=FALSE)
\method{plot}{bw.frac}(x, \dots)
\method{plot}{bw.optim}(x, \dots, showopt, optargs)
\method{plot}{lurk}(x, \dots, shade)
\method{plot}{minconfit}(x, \dots)
\method{plot}{parres}(x, \dots)
\method{plot}{plotpairsim}(x, \dots)
\method{plot}{qqppm}(x, \dots, limits=TRUE,
           monochrome=spatstat.options('monochrome'),
           limcol=if(monochrome) "black" else "red")
\method{plot}{spatialcdf}(x, \dots, xlab, ylab)
pointweights(X, \dots, weights, parent)
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
\method{predict}{vblogit}(object, newdata, type, se.fit, dispersion,
                          terms, na.action, \dots)
\method{predict}{profilepl}(object, \dots)
\method{print}{addvar}(x, \dots)
\method{print}{autoexec}(x, \dots)
\method{print}{bt.frame}(x, \dots)
\method{print}{bugtable}(x, \dots)
\method{print}{bw.frac}(x, \dots)
\method{print}{bw.optim}(x, \dots)
\method{print}{densityfun}(x, \dots)
\method{print}{diagppm}(x, \dots)
\method{print}{detpointprocfamily}(x, \dots)
\method{print}{detpointprocfamilyfun}(x, \dots)
\method{print}{envelope}(x, \dots)
\method{print}{hasenvelope}(x, \dots)
\method{print}{hierarchicalordering}(x, \dots)
\method{print}{influence.ppm}(x, \dots)
\method{print}{interact}(x, \dots, family, brief, banner)       
\method{print}{intermaker}(x, \dots)
\method{print}{isf}(x, \dots)
\method{print}{leverage.ppm}(x, \dots)
\method{print}{lurk}(x, \dots)
\method{print}{minconfit}(x, \dots)
\method{print}{mppm}(x, \dots)
\method{print}{msr}(x, \dots)
\method{print}{parres}(x, \dots)
\method{print}{plotppm}(x, \dots)
\method{print}{plotpairsim}(x, \dots)
\method{print}{profilepl}(x, \dots)
\method{print}{quadrattest}(x, \dots)
\method{print}{qqppm}(x, \dots)
\method{print}{rmhcontrol}(x, \dots)
\method{print}{rmhexpand}(x, \dots, prefix=TRUE)
\method{print}{rmhmodel}(x, \dots)
\method{print}{rmhstart}(x, \dots)
\method{print}{rmhInfoList}(x, \dots)
\method{print}{rppm}(x, \dots)
\method{print}{Smoothfun}(x, \dots)
\method{print}{summary.mppm}(x, \dots, brief)
\method{print}{summary.rmhexpand}(x, \dots)
\method{print}{summary.ssf}(x, \dots)
\method{print}{vblogit}(x, \dots)
quad.mppm(x)
quadBlockSizes(nX, nD, p, nMAX, announce)
RandomFieldsSafe()
reduceformula(fmla, deletevar, verbose)
reheat(model, invtemp)
RelevantDeviation(x, alternative, clamp, scaling)
\method{rescale}{msr}(X, s, unitname)
resolveEinfo(x, what, fallback, warn, atomic)
resolve.foxall.window(X, Y, W, warn.trim)
resolve.lambda(X, lambda, \dots, sigma, varcov, leaveoneout, update)
resolve.lambda.cross(X, I, J, lambdaI, lambdaJ, \dots, lambdaX,
                     sigma, varcov, leaveoneout, update, lambdaIJ,
                     Iexplain, Jexplain, miss.update, miss.leave,
                     caller)
resolve.vargamma.shape(\dots, nu.ker, nu.pcf, default)
roseContinuous(ang, rad, unit, \dots,
               start, clockwise, main, labels, at, do.plot)
\method{rotate}{msr}(X, angle, \dots, centre)
quadrat.testEngine(X, nx, ny, alternative, method, conditional, CR, 
     \dots, nsim, Xcount, xbreaks, ybreaks, tess, fit, df.est, Xname, fitname)
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
rhohatEngine(model, covariate, reference, volume, \dots,
               subset, weights, method, horvitz, smoother, resolution, 
               n, bw, adjust, from, to, 
               bwref, covname, covunits, confidence, modelcall, callstring)
rhohatCalc(ZX, Zvalues, lambda, denom, \dots,
           weights, lambdaX,
           method, horvitz, smoother,
           n, bw, adjust, from, to, 
           bwref, covname, confidence, positiveCI, markovCI,
           covunits, modelcall, callstring, savestuff)
rMaternInhibition(type, kappa, r, win, stationary, \dots, nsim, drop)
rmax.Rigid(X, g)
rmax.rule(fun, W, lambda)
\method{rmhcontrol}{rmhcontrol}(\dots) 
\method{rmhcontrol}{list}(\dots) 
rmhEngine(InfoList, \dots, verbose, kitchensink, preponly, snoop,
                           overrideXstart, overrideclip) 
RmhExpandRule(nama)
rmhResolveControl(control, model)
rmhResolveExpansion(win, control, imagelist, itype)
rmhResolveTypes(model, start, control)
rmhsnoop(\dots, Wsim, Wclip, R, xcoords, ycoords, mlevels, mcodes, irep, itype, 
     proptype, proplocn, propmark, propindx, numerator, denominator,
     panel.only)
rmhSnoopEnv(Xinit, Wclip, R)
\method{rmhmodel}{rmhmodel}(model, \dots) 
\method{rmhstart}{rmhstart}(start, \dots) 
\method{rmhstart}{list}(start, \dots) 
rmpoint.I.allim(n, f, types)
rocData(covariate, nullmodel, \dots, high)
rocModel(lambda, nullmodel, \dots, high)
rpoint.multi(n, f, fmax, marks, win, giveup, verbose, warn, nsim, drop)
runifpoispp(lambda, win, \dots, nsim, drop)
runifpoisppOnLines(lambda, L, nsim, drop)
safelookup(Z, x, factor, warn)
safevalue(x, default)
\method{scalardilate}{msr}(X, f, \dots)
scanmeasure(X, \dots)
\method{scanmeasure}{ppp}(X, r, \dots, method) 
\method{scanmeasure}{im}(X, r, \dots)
scanPoisLRTS(nZ, nG, muZ, muG, alternative)
scanBinomLRTS(nZ, nG, muZ, muG, alternative)
second.moment.calc(x, sigma, edge, what, \dots,
                   varcov, expand, obswin, npts, debug)
second.moment.engine(x, sigma, edge, what, \dots,
      kernel, scalekernel, obswin, varcov, npts, debug, fastgauss)
sewpcf(d, w, denargs, lambda2area, divisor)
sewsmod(d, ff, wt, Ef, rvals, method="smrep", \dots, nwtsteps=500)
\method{shift}{influence.ppm}(X, \dots)
\method{shift}{leverage.ppm}(X, \dots)
\method{shift}{msr}(X, \dots)
\method{shift}{quadrattest}(X, \dots)
signalStatus(x, errors.only)
\method{simulate}{profilepl}(object, \dots)
simulrecipe(type, expr, envir, csr, pois, constraints)
slr.prepare(CallInfo, envir, data, dataAtPoints, splitby, clip)
slrAssemblePixelData(Y, Yname, W, covimages, dataAtPoints, pixelarea)
\method{Smooth}{solist}(X, \dots)
smoothcrossEngine(Xdata, Xquery, values, sigma, \dots,
                    weights, varcov,
                    kernel, scalekernel, sorted, cutoff)
smoothpointsEngine(x, values, sigma, \dots,
                   kernel, scalekernel,
                   weights, varcov, leaveoneout, sorted, cutoff, debug)
smudge(X)
spatstatClusterModelInfo(name, onlyPCP)
spatstatDPPModelInfo(model)
spatstatRmhInfo(cifname)
spatialCDFframe(model, covariate, \dots, jitter)
spatialCDFtest(model, covariate, test, \dots,
         dimyx, eps, interpolate, jitter,
         nsim, verbose, modelname, covname, dataname)
spatialCDFtestCalc(fra, test, \dots, details)
sphere.volume(range, nval = 10)
splitHybridInteraction(coeffs, inte)
sp.foundclass(cname, inlist, formalname, argsgiven)             
sp.foundclasses(cnames, inlist, formalname, argsgiven)
strausscounts(U,X,r,EqualPairs)
suffloc(object)
suffstat.generic(model, X, callstring)
suffstat.poisson(model, X, callstring)
summarise.trend(trend, w, a)
\method{summary}{envelope}(object,\dots)
\method{summary}{mppm}(object, \dots, brief=FALSE)
\method{summary}{msr}(object, \dots)
\method{summary}{profilepl}(object, \dots)
\method{summary}{rmhexpand}(object, \dots)
\method{summary}{vblogit}(object, \dots)
thinjump(n, p)
tweak.coefs(model, new.coef)
twostage.test(X, \dots, exponent, nsim, nsimsub,
              alternative, reuse, leaveout, interpolate,
              savefuns, savepatterns, verbose, badXfatal, testblurb)
twostage.envelope(X, \dots, nsim, nsimsub, nrank,
                  alternative, reuse, leaveout, interpolate,
                  savefuns, savepatterns, verbose, badXfatal, testlabel)
\method{unitname}{msr}(x) 
\method{unitname}{msr}(x) <- value 
\method{update}{ippm}(object, \dots, envir) 
\method{update}{msr}(object, \dots) 
\method{update}{rmhstart}(object, \dots)
validate2Dkernel(kernel, fatal)
validate.angles(angles, unit, guess)
validate.weights(x, recip, how, allowzero, allowinf)
varcountEngine(g, B, lambdaB, f)
%vblogit(y, X, offset, eps, m0, S0, S0i, xi0, verb, maxiter, \dots)
%vblogit.fmla(formula, offset, data, subset, weights, verbose, epsilon, \dots)
versionstring.interact(object)
versionstring.ppm(object)
weightedclosepairs(X, r, correction, what)
windows.mppm(x)
X2testEngine(OBS, EXP, \dots, method, CR, df, nsim,
     conditional, alternative, testname, dataname)

%% Linear networks code is corralled here
ApplyConnected(X, Engine, r, \dots, rule, auxdata)
DoCountEnds(X, D, toler)
DoCountCrossEnds(X, I, J, DIJ, toler)
FDMKERNEL(lppobj, dtt, dtx, M, nsave, weights,
          stepnames, setuponly, verbose)
\method{as.linfun}{linfun}(X, \dots)
\method{as.owin}{lintess}(W, \dots)
default.linnet.tolerance(L)
makeLinnetTolerance(toler)
\method{print}{lintess}(x, \dots)
\method{print}{summary.linim}(x, \dots)
\method{print}{summary.linnet}(x, \dots)
\method{print}{summary.lintess}(x, \dots)
\method{summary}{lintess}(object, \dots)
\method{Window}{lintess}(X, \dots)
\method{Window}{linnet}(X, \dots, check=TRUE) <- value
\method{Window}{lpp}(X, \dots, check=TRUE) <- value
densitypointsLPP(x, sigma, \dots,
                 weights, nsigma, leaveoneout, fast,
                 fastmethod, floored,
                 dx, dt, iterMax, verbose, debug)
\method{evalCovar}{lppm}(model, covariate, \dots, lambdatype, 
          eps, nd, interpolate, jitter, modelname, covname, dataname,
          subset)
flatdensityfunlpp(X, \dots, disconnect, weights, what)
flatdensityatpointslpp(X, \dots, leaveoneout, disconnect, weights, what)
getlambda.lpp(lambda, X, subset, \dots,
              update, leaveoneout, loo.given, lambdaname)
local2lpp(L, seg, tp, X, df.only)
looHeatLPP(U0, Amatrix, npts, niter, nsave,
           lixelweight, lixelmap, verbose) 
looVoronoiLPP(X)
validate.lpp.coords(X, fatal, context)
\method{as.ppm}{lppm}(object)
\method{as.ppm}{rppm}(object)
pointsAlongNetwork(L, delta)
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
resampleNetworkDataFrame(df, template)
sortalongsegment(df)
vnnFind(seg, tp, ns, nv, from, to, seglen, huge, tol, kmax)
ldtEngine(nv, ns, from, to, seglen, huge,
          coUXord, vnndist, vnnwhich, vnnlab)
resolve.heat.steps(sigma, \dots, dx, dt,
                   niter, iterMax, nsave,
                   seglengths, maxdegree, AMbound, L,
                   finespacing, fineNsplit, fineNlixels,
                   W, eps, dimyx, xy, 
                   allow.adjust, warn.adjust,
                   verbose, stepnames)
qkdeEngine(X, sigma, \dots, at, what,
           leaveoneout, diggle, raw, edge2D, edge,
           weights, varcov, positive, shortcut,
           precomputed, savecomputed)
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm = FALSE)}
%NAMESPACE S3method("Math", "linimlist")
%NAMESPACE S3method("Ops", "linimlist")
%NAMESPACE S3method("Complex", "linimlist")
%NAMESPACE S3method("Summary", "linimlist")
LinimOp(e1, e2, op)
LinimListOp(e1, e2, op)
}
\details{
  These internal \pkg{spatstat} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat} to the next.
}
\keyword{internal}
