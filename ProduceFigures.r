library(DreamBox7)
library(affy)

setwd("~/Dropbox/CU/PhD/STM/")

data(map)
data(attractome.minimalist)


#=============================================
#	Create working space
#============================================

load("~/workspace/data/dream7/fullrenorm/metabric2000.rdata")
clnc = lazyImputeDFClncOslo(metabric$Complete_METABRIC_Clinical_Features_Data)
clinical <- expandClncOslo(clnc)

ge = exprs(metabric$Complete_METABRIC_Expression_Data)
meta = CreateMetageneSpace(ge, attractome=attractome.minimalist, map=map)
meta = meta$metaSpace

ls = meta["ls",]
meta = t(apply(meta, 1, function(x){ x - median(x)}))

idx = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
mes.lymphneg = meta["mt",] * idx
mes.lymphneg[idx] = mes.lymphneg[idx] - median(mes.lymphneg[idx])
meta = rbind(meta, mes.lymphneg)

idx = (clinical[,"lymph_nodes_positive"] > 3)
ls.lymphpos = ls * idx
ls.lymphpos[idx] = ls.lymphpos[idx] - median(ls.lymphpos[idx])
meta = rbind(meta, ls.lymphpos)

idx = clnc$ER.Expr == "-"
ls.erneg = ls * idx
ls.erneg[idx] = ls.erneg[idx] - median(ls.erneg[idx])
meta = rbind(meta, ls.erneg)

surv.os = metabric$Complete_METABRIC_Clinical_Survival_Data_OS
surv.ds = metabric$Complete_METABRIC_Clinical_Survival_Data_DSS

save(clinical, clnc, meta, surv.os, surv.ds, file="~/workspace/data/dream7/fullrenorm/metabric2000.workspace.rda")

#===============================================
#	directly load working space
#===============================================

load("~/workspace/data/dream7/fullrenorm/metabric2000.workspace.rda")

#===============================================
#	survival curves
#===============================================

#Mitotc
metafeature=meta["mitotic",]
time=surv.ds[,1]
status=surv.ds[,2]

idx = time >= 15*365
time[idx] = 15*365
status[idx] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.cin <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.cin<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#Plot
tiff(file = "figs/mitotic.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0))     #axis and label margin

plot(
fit.cin,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "15-year survival based on mitotic CIN attractor expresssion level",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(2000, 0.4, c("High", "Low"),title="Mitotic CIN attractor expression level", lwd=1:1, col=c("18","20"),cex=0.4)
text(2000,0.6,paste("P-value < 2E-16 "),cex=0.5)
dev.off()       #Write

tiff(file = "figs/lsxlymNum.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0))     #axis and label margin

plot(
fit.cin,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "15-year survival based on LYM*(10-lymNum)",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(2000, 0.4, c("High", "Low"),title="LYM*(10-lymNum)", lwd=1:1, col=c("18","20"),cex=0.4)
text(2000,0.6,paste("P-value < 2E-16 "),cex=0.5)
dev.off()       #Write
# MES

##Early Stage Samples
idx.early = (clinical[,"lymph_nodes_positive"]<1 & clinical[,"size"] < 30)
metafeature=meta["mes.lymphneg",idx.early]
time=surv.ds[idx.early,1]
status=surv.ds[idx.early,2]

idx = time >= 10*365
time[idx] = 10*365
status[idx] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.lymphneg <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.lymphneg<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#Plot
tiff(file = "figs/meslymneg30.metabric.ds.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)     #axis and label margin
)         #2 columns

plot(
fit.lymphneg,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "10-year survival based on MES metagene restricted to early-stage samples",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(1000, 0.4, c("High", "Low"),title="Mesenchymal transition attractor expression level", lwd=1:1, col=c("18","20"),cex=0.4)
text(2500,0.75,paste("P-value = ",round(pval.lymphneg,4)),cex=0.5)
dev.off()       #Write

##ER- Samples
#idx.erN = (meta["er",] < 0 & meta["erbb2", ] < 0)
idx.erN = (clnc$ER.Expr=="-" & clnc$Her2.Expr=="-" & clnc$PR.Expr=="-")
metafeature=meta["ls",idx.erN]
time=surv.ds[idx.erN,1]
status=surv.ds[idx.erN,2]

idx = time >= 15*365
time[idx] = 15*365
status[idx] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.erN <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.erN<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

##lymphNode >3 Samples
idx.lymph3 = (clinical[,"lymph_nodes_positive"] > 4 & clnc$ER.Expr=="+")
metafeature=meta["ls",idx.lymph3]
time=surv.ds[idx.lymph3,1]
status=surv.ds[idx.lymph3,2]

idx = time >= 15*365
time[idx] = 15*365
status[idx] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.lymph3 <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.lymph3<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

#Plot
tiff(file = "figs/lym.TN.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)     #axis and label margin
)         #2 columns

#ER- Samples
plot(
fit.erN,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "15-year survival based on LYM metagene restricted to ER-/HER2- samples",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(2000, 0.3, c("High", "Low"),title="Lymphocyte attractor expression level", lwd=1:1, col=c("18","20"),cex=0.4)
text(3000,0.4,paste("P-value = ",round(pval.erN,4)),cex=0.5)

dev.off()

tiff(file = "figs/lym.erN.lymNum4.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)     #axis and label margin
)
#Positive lymph node number >3 Samples
plot(
fit.lymph3,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "15-year survival based on LYM metagene restricted to ER+ / lymph node number > 4",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(2000, 0.9, c("High", "Low"),title="Lymphocyte attractor expression level", lwd=1:1, col=c("18","20"),cex=0.4)
text(3000,0.1,paste("P-value = ",round(pval.lymph3,4)),cex=0.5)

dev.off()       #Write

#SUSD3-FGD3

metafeature=meta["susd3",]
time=surv.ds[,1]
status=surv.ds[,2]

idx = time >= 15*365
time[idx] = 15*365
status[idx] = 0

X<-cbind(time,status,as.numeric( metafeature<median(metafeature)))
colnames(X)=c("time","status", "x")
fit.susd3 <- survfit(Surv(time, status) ~ x, data = data.frame(X))
pval.susd3<-summary(coxph(Surv(time, status) ~ metafeature, ))$logtest[3]

tiff(file = "figs/susd3.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)     #axis and label margin
)

plot(
fit.susd3,
col = c("18","20"),     #Red, Blue
lwd = 1:1,
mark.time = FALSE,
main = "15-year survival based on SUSD3-FGD3 metagene",
xlab = "Days",
ylab = "% Survived",
yscale = 100,
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

legend(1500, 0.4, c("High", "Low"),title="SUSD3-FGD3 metagene expression level", lwd=1:1, col=c("18","20"),cex=0.4)
text(2000,0.6,paste("P-value < 2E-16 "),cex=0.5)

dev.off()       #Write
#============================
#	Scatter plots
#===========================

load("~/workspace/data/dream7/fullrenorm/ge.49576x1981.rda")

tiff(file = "figs/ER_CIN.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)      #axis and label margin
)

#Mitotic CIN metagene vs. Estrogen Receptor Metagene

plot(
meta["mitotic",],
ge["ILMN_1678535",], # ESR1
lwd = 1:1,
pch=20,         #dots
col="blue",
cex=0.5,
main = "Mitotic CIN metagene expression level vs. ESR1 expression level",
xlab = "Mitotic CIN metagene",
ylab = "Estrogen Receptor",
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)
dev.off()       #Write

#SUSD3 vs ESR1 (single gene)

tiff(file = "figs/SUSD3_ER.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)      #axis and label margin
)

#SUSD3 vs. ESR1
plot(
meta["susd3",],
meta["er",],
lwd = 1:1,
pch=20,         #dots
col="blue",
cex=0.5,
main = "SUSD3-FGD3 metagene expression level vs. ESR1 expression level",
xlab = "SUSD3-FGD3 metagene",
ylab = "ER metagene ",
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)
dev.off()       #Write

#SUSD3 vs. FGD3

tiff(file = "figs/SUSD3_FGD3.metabric.tiff", width = 3.5, height = 3.5, units = "in", res = 300, compression="lzw")
par(
mar = c(2,2,2,1),       #plot margin
mgp = c(1, 0.4, 0)      #axis and label margin
)

#SUSD3 vs. ESR1
plot(
ge["ILMN_1785570",],
ge["ILMN_1772686",],
lwd = 1:1,
pch=20,         #dots
col="blue",
cex=0.5,
main = "Association between expression levels of SUSD3 and FGD3",
xlab = "SUSD3",
ylab = "FGD3",
cex.main = 0.5,
cex.axis = 0.5,
cex.lab = 0.5)

lines(x=c(7, 7), y=c(5, 7), col="red", lty="1")

dev.off()       #Write
