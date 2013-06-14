require(synapseClient)
require(devtools)
require(Biobase)
require(survival)
require(BCC)
require(devtools)

# install DreamBox7 package from github
install_github(repo="DreamBox7", username="weiyi-bitw", ref="master")
library(DreamBox7)

synapseLogin()

# Load metabric data
clnc = loadEntity("syn1710260")$objects[[1]]
cna = loadEntity("syn1710262")$objects[[1]]
survdss = loadEntity("syn1730400")$objects[[1]]
expr = loadEntity("syn1710275")$objects[[1]]
surv = loadEntity("syn1710277")$objects[[1]]

metabric = list(expr = expr, cnv = cna, clnc = clnc, surv = surv, survdss = survdss)


# Load OSLOVAL data
intClncDat = loadEntity("syn1710251")$objects[[1]]
intCnvEset = loadEntity("syn1710253")$objects[[1]]
intSurvObject = loadEntity("syn1710257")$objects[[1]]
intExprEset = loadEntity("syn1710255")$objects[[1]]

oslo = list(expr = intExprEset, cnv = intCnvEset, clnc = intClncDat, surv = intSurvObject)

# source model source file
source_url("https://raw.github.com/weiyi-bitw/BCCModels/master/syn1417992_OSDS.R")

# create new model instance
gdModel = GoldiloxModel$new()

# train model
gdModel$customTrain(
metabric$expr,
metabric$cnv, 
metabric$clnc, 
metabric$surv, 
metabric$survdss
)

# get training score

trainPredictions <- gdModel$customPredict(
metabric$expr, 
metabric$cnv,
metabric$clnc
)

tr.ci = getCCDIdx(trainPredictions, metabric$surv);tr.ci

# Create prediction on OSLOVAL dataset

osloPredictions <- gdModel$customPredict(
oslo$expr, 
oslo$cnv,
oslo$clnc
)

# note that since the weights of the subclassifiers contain random factor, 
# the CI will be slightly different each time
oslo.ci = getCCDIdx(osloPredictions, oslo$surv);oslo.ci



