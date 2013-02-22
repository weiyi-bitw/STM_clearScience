library(synapseClient)
synapseLogin()

if(!require(cafr)){
	syn = loadEntity("syn1446295")
	install.packages(paste(syn$cacheDir, "cafr", sep="/"), repos=NULL)
	library(cafr)
}
if(!require(impute)){
	source("http://bioconductor.org/biocLite.R")
	biocLite(impute)
	library(impute)
}

panCancerAttractorSearch = function(synList, seed, a=5, verbose=TRUE){
	out = list(attractors=NULL, mis=NULL, annotations=NULL, topGenes=NULL)
	for(synID in synList){
	#===== Loading data =====
		cat("Processing " , synID, "...\n", sep="");flush.console()
		syn = loadEntity(synID)
		annot = c(syn$properties$name, syn$properties$tissueType, syn$properties$platform, syn$properties$numSamples)
		out$annotations = rbind(out$annotations, annot)
	#===== Impute missing values =====
		cat("\t", syn$properties$name, "\n");flush.console()
		cat("Imputing missing expression values... \n");flush.console()
		ge = as.matrix(read.delim(file.path(syn$cacheDir, syn$files[[1]]), row.names=1, check.names=F ))
		ge[ge=="null"] = NA
		ge = impute.knn(ge)$data

	#===== Extract only tumor samples =====
		sampleIDs = colnames(ge)
		sampleTypes = sapply(sampleIDs, function(x){as.numeric(substr(strsplit(x, split="-")[[1]][4], 1, 2))})
		idx = sampleTypes < 10
		ge = ge[,idx]

	#===== find attractor and store only top 20 genes =====
		cat("Converge to attractor... \n");flush.console()
		atr = CAFrun(ge, ge[seed,], a = a, verbose=verbose, sorting=T)
		out$attractors = cbind(out$attractors, names(atr)[1:50])
		out$mis = cbind(out$mis, as.numeric(atr[1:50]))
	}
	rownames(out$annotations) = synList
	colnames(out$annotations) = c("Name", "TissueType", "Platform", "NumSamples")
	out$annotations = data.frame(out$annotations)
	colnames(out$attractors) = out$annotations$TissueType
	colnames(out$mis) = out$annotations$TissueType
	out$mis = data.frame(out$mis)

	allGenes = unique( as.vector(out$attractors) )
	mat = matrix(0, nrow=length(allGenes), ncol=ncol(out$attractors))
	rownames(mat) = allGenes
	for(i in 1:ncol(out$attractors)){
		mat[out$attractors[,i],i] = out$mis[,i]
	}
	scores = apply(mat, 1, mean)
	scores = sort(scores, decreasing=T)[1:20]

	t = table(as.vector(out$attractors))

	out$topGenes = cbind( t[names(scores)], round(scores, 3) )
	colnames(out$topGenes) = c("Frequency", "Avg MI")
	cat("Done!\n");flush.console()
	return (out)
}

synList = c("syn363348", "syn363356", "syn363394", "syn363404")

# find the Chromosomal instability attractor using seed gene CENPA
out = panCancerAttractorSearch(synList, "CENPA")
cin = out
save(cin, file="~/Dropbox/workspace/tcgaPanCancer/cin.rdata")

# find the mesenchymal transition attractor using seed gene COL5A2
out = panCancerAttractorSearch(synList, "COL5A2")
mes = out
save(mes, file="~/Dropbox/workspace/tcgaPanCancer/mes.rdata")

# find the lymphocyte-specific attractor using seed gene LAPTM5
out = panCancerAttractorSearch(synList, "LAPTM5")
lym = out
save(lym, file="~/Dropbox/workspace/tcgaPanCancer/lym.rdata")



