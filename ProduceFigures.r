require(synapseClient)
require(devtools)
require(Biobase)
require(survival)
require(BCC)
require(ggplot2)

install_github(repo="DreamBox7", username="weiyi-bitw", ref="master")
library(DreamBox7)

data(map)
data(attractome.minimalist)

#=============================================
#	Load METABRIC and OSLO data
#============================================

synapseLogin()

#=====================================================================================
#
#  I didn't change this part. It seems I don't have access to syn1710250?
#
#====================================================================================

metabric = loadEntity('syn1465025')$objects
#
# Since I don't have read access to the full metabric now, 
# I used the following to load a local copy
#
#load("~/workspace/data/dream7/fullrenorm/metabric2000.rdata")
clnc = lazyImputeDFClncOslo(metabric$Complete_METABRIC_Clinical_Features_Data)
clinical <- expandClncOslo(clnc)

ge = exprs(metabric$Complete_METABRIC_Expression_Data)
o = CreateMetageneSpace(ge, attractome=attractome.minimalist, map=map)
meta = o$metaSpace

# keep mesenchymal and lymph-specific metagene for conditioning
ls = meta["ls",]
mt = meta["mt",]

# median-center each metagene
meta = t(apply(meta, 1, function(x){ x - median(x)}))
surv = metabric$Complete_METABRIC_Clinical_Survival_Data_DSS

# load OSLOVAL object

intClncDat = loadEntity("syn1449480")$objects$xIntClinDat
intCnvEset = loadEntity("syn1449473")$objects$xCnvDat
intSurvObject = loadEntity("syn1449477")$objects$xIntSurvObj
intExprEset = loadEntity("syn1449475")$objects$xExprDat

oslo = list(expr = intExprEset, cnv = intCnvEset, clnc = intClncDat, surv = intSurvObject)

#================================================================================================

clnc.oslo = lazyImputeDFClncOslo(oslo$clnc)
clinical.oslo <- expandClncOslo(clnc.oslo)

ge.oslo = exprs(oslo$expr)
o = CreateMetageneSpace(ge.oslo, attractome=attractome.minimalist, map=map)
meta.oslo = o$metaSpace

# keep mesenchymal and lymph-specific metagene for conditioning
ls.oslo = meta["ls",]
mt.oslo = meta["mt",]

# median-center each metagene
meta.oslo = t(apply(meta.oslo, 1, function(x){ x - median(x)}))
surv.oslo = oslo$surv

#===============================================
#	survival curves
#===============================================

# for quartz device...
quartzFonts() # lists current mappings
     
# Add a mapping for Myriad:
# pass in a vector of font filenames in the order: Plain, Bold, Italic, Bold-Italic
quartzFonts(
  Myriad = quartzFont(
    c("MyriadPro-Regular", "MyriadPro-Semibold", "MyriadPro-It", "MyriadPro-SemiboldIt")
  )
)
     
dir.create("figs")

#-------------------
#  figure 1: mitotic
#-------------------

#Plot
tiff(file = "figs/fig1.tiff", width = 7.3, height = 3.5, units = "in", res = 300, compression="lzw", pointsize=8)
par(mar = c(2.5,2.5,2.5,1),       #plot margin
mgp = c(1.2, 0.3, 0),
mfrow = c(1, 2))     #axis and label margin

#
#
#  fig1A: metabric mitotic
#
#

#Mitotc
metafeature=meta["mitotic",]
metafeature = metafeature - median(metafeature)
time=surv[,1]
status=surv[,2]

# truncate to 15-year survival
ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.cin <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.cin<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

plot(
fit.cin,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "A",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")


legend(2200, 0.4, c("High", "Low"),title="CIN", lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.5, expression(paste(italic("P"), " value < ", 2 %*% 10^{-16} )),cex=1, family = "Myriad")

#
#
#  fig1B: mitotic oslo
#
#

metafeature=meta.oslo["mitotic",]
metafeature = metafeature - median(metafeature)
time=surv.oslo[,1]
status=surv.oslo[,2]

# truncate to 15-year survival
ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.cin.oslo <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.cin.oslo<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

plot(
fit.cin.oslo,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "B",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")

legend(3500, 0.9, c("High", "Low"),title="CIN", lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.3,substitute(paste( italic("P"), " value = ", k), list(k=round(pval.cin.oslo, 4))), cex=1, family="Myriad")

dev.off()       #Write


#----------------------
#  figure 2: LYM in ER-
#----------------------



#Plot
tiff(file = "figs/fig2.tiff", width = 3.5, height = 10.5, units = "in", res = 300, compression="lzw", pointsize=12)
par(
mar = c(2.5,3,2.5,1),       #plot margin
mgp = c(1.5, 0.3, 0),     #axis and label margin
mfrow = c(3, 1)
)         #2 columns

#
#
# fig2A: LYM | ER- metabric
#
#
##ER- Samples
#idx.erN = (meta["er",] < 0 & meta["erbb2", ] < 0)
idx.erN = clnc$ER.Expr=="-" 
metafeature=meta["ls",idx.erN]
metafeature = metafeature - median(metafeature)
time=surv[idx.erN,1]
status=surv[idx.erN,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.erN <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.erN<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#ER- Samples
plot(
fit.erN,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "A",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1, 
family = "Myriad")

legend(2300, 0.3, c("High", "Low"),title="LYM", lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.4,substitute(paste( italic("P"), " value = ",k ), list(k=round(pval.erN,4))),cex=1, family = "Myriad")

#
#
#  fig2B: LYM | ER- oslo
#
#

idx.erN = clnc.oslo$ER.Expr=="-" 
metafeature=meta.oslo["ls",idx.erN]
metafeature = metafeature - median(metafeature)
time=surv.oslo[idx.erN,1]
status=surv.oslo[idx.erN,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.erN.oslo <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.erN.oslo<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#ER- Samples
plot(
fit.erN.oslo,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "B",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1, 
family = "Myriad")

legend(3000, 0.9, c("High", "Low"),title="LYM", lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.2,substitute(paste( italic("P"), " value = ", k ),list(k=round(pval.erN.oslo, 4))),cex=1, family = "Myriad")

#-------------------
# figure 2C: LYM in ER+ with positive lymph nodes
#-------------------

# NOTE: OSLO does not have enough samples with ER+ / lymNum+ 
# to perform a meaningful analysis

##lymphNode >3 Samples
idx.erP = (clinical[,"lymph_nodes_positive"] > 4 & clnc$ER.Expr=="+")
metafeature=meta["ls",idx.erP]
metafeature = metafeature - median(metafeature)
time=surv[idx.erP,1]
status=surv[idx.erP,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.erP <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.erP<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#Positive lymph node number >3 Samples
plot(
fit.erP,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "METABRIC ER+ | lymph node number > 4",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")

legend(3000, 0.95, c("High", "Low"),title="LYM", lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.1,substitute(paste( italic("P"), " value = ", k ),list(k=round(pval.erP, 4))),cex=1, family = "Myriad" )

dev.off()       #Write

#------------
# figure 3: SUSD3-FGD3 metagene
#------------


tiff(file = "figs/fig3.tiff", width = 7.3, height = 7, units = "in", res = 300, compression="lzw", pointsize=11)
par(
mar = c(2.5,3,2.5,1),       #plot margin
mgp = c(1.5, 0.3, 0),    #axis and label margin
mfrow=c(2, 2)
)


#------------
# figure 3A: SUSD3, FGD3, ESR1
#------------

#
#
# SUSD3 vs FGD3
#
#
plot(
ge["ILMN_1772686",],
ge["ILMN_1785570",],
lwd = 1:1,
pch=20,         #dots
col="blue",
cex=0.5,
#main = "(A) FGD3 vs. SUSD3 expression level",
ylab = expression(italic("SUSD3")),
xlab = expression(italic("FGD3")),
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")


#
#
# fig3B: SUSD3-FGD3 vs ESR1
#
#

plot(
meta["susd3",],
ge["ILMN_1678535",],
lwd = 1:1,
pch=20,         #dots
col="blue",
cex=0.5,
#main = "(B) FGD3-SUSD3 metagene vs. ESR1 expression level",
xlab = expression(paste( italic("FGD3"),"-", italic("SUSD3"), " metagene")),
ylab = expression(italic("ESR1")),
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")


#
#
# fig3C: SUSD3-FGD3 metabric
#
#

#SUSD3-FGD3

metafeature=meta["susd3",]
metafeature = metafeature - median(metafeature)
time=surv[,1]
status=surv[,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.susd3 <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.susd3<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

plot(
fit.susd3,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "(A) METABRIC",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")

legend(900, 0.4, c("High", "Low"),title=expression(paste( italic("FGD3"),"-", italic("SUSD3"), " metagene")), lwd=1:1, col=c("18","20"),cex=1)
text(2500,0.5,expression(paste( italic("P") , " value < ", 2 %*% 10^{-16} )),cex=1, family = "Myriad")

#
#
# fig3D: SUSD3-FGD3 OSLO
#
#

metafeature=meta.oslo["susd3",]
metafeature = metafeature - median(metafeature)
time=surv.oslo[,1]
status=surv.oslo[,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.susd3.oslo <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.susd3.oslo<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

plot(
fit.susd3.oslo,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "(B) OSLOVAL",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")

legend(1500, 0.3, c("High", "Low"),title=expression(paste( italic("FGD3"),"-", italic("SUSD3"), " metagene")), lwd=1:1, col=c("18","20"),cex=1)
text(3000,0.35,substitute(paste( italic("P"), " value = ", k ),list(k=round(pval.susd3.oslo, 4))),cex=1, , family = "Myriad" )
dev.off()       #Write
#---------
#
# fig 4
#
#--------

# load the predictions made by the model
load("removingFeaturePredictions.rda")

metafeature = pmat[1,]
time=surv.oslo[,1]
status=surv.oslo[,2]

ii = time >= 15*365
time[ii] = 15*365
status[ii] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.m <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.m<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]


tiff(file = "figs/fig4.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw", pointsize=8)
par(
mar = c(3,3,1,1),       #plot margin
mgp = c(1.5, 0.3, 0)     #axis and label margin
)
#Positive lymph node number >3 Samples
plot(
fit.m,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
#main = "METABRIC ER+ | lymph node number > 4",
xlab = "Days",
ylab = "Survival (%)",
yscale = 100,
cex.main = 1,
cex.axis = 0.9,
cex.lab = 1,
family = "Myriad")

legend(500, 0.3, c("Poor", "Good"),title="Predicted survival", lwd=1:1, col=c("18","20"),cex=1)
text(3500,0.1,expression(paste( italic("P") , " value < ", 2 %*% 10^{-16} ) ),cex=1, family = "Myriad" )

dev.off()       #Write
