genome <- initializeGenomeObject(file = "genome.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))
model <- initializeModelObject(parameter = parameter, model = "ROC")
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
runMCMC(mcmc = mcmc, genome = genome, model = model)

trace <- parameter$getTraceObject()

plot(x = trace, what = "Mutation", mixture = 1)

pdf()
plot(x = trace, what = "Mutation", mixture = 1)
dev.off()

pdf()
AnaCoDa:::plot.Rcpp_Trace(x = trace, what = "Mutation", geneIndex = 1, mixture = 1, log.10.scale = F)
dev.off()

AnaCoDa:::plotCodonSpecificParameters(trace = trace, mixture = 1, type = "Mutation", main = "Mutation Parameter Traces")

pdf(title = "plotCSPtest")
AnaCoDa:::plotCodonSpecificParameters(trace = trace, mixture = 1, type = "Mutation", main = "Mutation Parameter Traces")
dev.off()

#plot Rcpp trace without having AnaCoDa loaded
plot.Rcpp_Trace()

#create alternate library save path
myPaths <- .libPaths()
myPaths <- c(myPaths, "~/Desktop/R/lib-dev")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths) #add new path

library(devtools)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("ebarne24/RibModelFramework", dependencies = TRUE)

usethis::use_build_ignore(c("dev"))

