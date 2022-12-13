#Changes by Elizabeth Barnes, 12/9/2022
##Including only plotting function for ROC parameters

#' Plot Parameter 
#' 
#' @param x A parameter object
#' 
#' @param what Which aspect of the parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#'
#' @param mixture.name a vector with names/descriptions of the mixture distributions in the parameter object
#'
#' @param with.ci Plot with or without confidence intervals. Default value
#' is TRUE
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC or FONSE
#' parameter object for each mixture element.
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
plot.Rcpp_ROCParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, aa.exclude = c("M", "W"), aa.names = NULL, ...)
{
  plotParameterObject(x, what = what, samples= samples, mixture.name=mixture.name, with.ci=with.ci, aa.exclude = aa.exclude, aa.names = aa.names, ...)
}


### NOT EXPOSED
plotParameterObject <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, aa.exclude = c("M", "W"), aa.names = NULL, ...){
  numMixtures <- x$numMixtures
  means <- data.frame(matrix(0,ncol=numMixtures,nrow=40))
  sd.values <- data.frame(matrix(0,ncol=numMixtures*2,nrow=40))
  
  ## Check to ensure aa.names passed are valid
  if(is.null(aa.names)) {
    aa.names <- aminoAcids() }
  else {
    if(length(aa.names) == 0) {
      stop("aa.names not found, length of zero")
    }
    

  ## remove aa we don't want to or can't plot
  aa.names <- aa.names[aa.names %in% aa.exclude]
  paramType <- ifelse(what == "Mutation", 0, 1)
  #cat("ParamType: ", paramType, "\n")
  
  #warning message for any values include in aa.exclude
  if(length(aa.exclude) > 0) {
    warning("Members of aa.names included in aa.exclude will be excluded.")
  }
  
  
  ##NOTE by Elizabeth Barnes 12/8/2022
  ### aa M, W, and X are skipped because n_syn = 1 or 3
  
  for (mixture in 1:numMixtures) {
    # get codon specific parameter
    count <- 1
    for (aa in aa.names) {
      codons <- AAToCodon(aa, T)
      for (i in 1:length(codons))
      {
        means[count,mixture] <- x$getCodonSpecificPosteriorMean(mixture, samples,codons[i], paramType, TRUE, log_scale=FALSE)
        tmp <- x$getCodonSpecificQuantile(mixture,samples, codons[i], paramType, c(0.025, 0.975), TRUE, log_scale=FALSE)
        
        ## This approach to storing the quantiles may seem unconventional, but I actually found it to be the most straight forward approach
        ## for plotting later.
        sd.values[count,mixture] <- tmp[1]
        sd.values[count,mixture+numMixtures] <- tmp[2]
        count <- count + 1
      }
    }
  }}
  ## Begin graphing
  mat <- matrix(rep(0,numMixtures*numMixtures),
                nrow = numMixtures, ncol = numMixtures, byrow = TRUE)
  count <- 1
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i<=j){
        mat[i,j] <-count
        count <- count + 1
      }
    }
  }
  nf <- layout(mat,widths=c(rep(5,numMixtures)),heights=c(rep(5,numMixtures)),respect=FALSE)
  par(mar=c(1,1,1,1))
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i==j)
      {
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="",xaxt='n',yaxt='n',ann=FALSE)
        if(is.null(mixture.name)){
          text(x = 0.5, y = 0.5, paste0("Mixture\nElement",i), 
               cex = 1.6, col = "black")
        }else{
          text(x = 0.5, y = 0.5, mixture.name[i], 
               cex = 1.6, col = "black")
        }
      }
      else if (i < j){
        if(with.ci){
          plot(means[,j],means[,i],ann=FALSE,xlim=range(cbind(sd.values[,j],sd.values[,j+numMixtures])),ylim=range(cbind(sd.values[,i],sd.values[,i+numMixtures])))
          upper.panel.plot(means[,j],means[,i],sd.x=cbind(sd.values[,j],sd.values[,j+numMixtures]),sd.y=cbind(sd.values[,i],sd.values[,i+numMixtures]))
          #confidenceInterval.plot(x = means[,j],y = mean[,i], sd.x=sd.values[,j],sd.y=sd.values[,i])
        } else{
          plot(means[,j],means[,i],ann=FALSE,xlim=range(means[,j]),ylim=range(means[,i]))
          upper.panel.plot(means[,j],means[,i])
        }
      }
    }
  }
}

##Test plots made after modifying code

genome <- initializeGenomeObject(file = "orf_coding.fasta")

parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))

model <- initializeModelObject(parameter = parameter, model = "ROC")

mcmc <- initializeMCMCObject(samples = 500, thinning = 10, adaptive.width = 50)

runMCMC(mcmc = mcmc, genome = genome, model = model)

original_parameter_plot <- plotParameterObject(x = parameter, what = "Mutation", samples = 2000, mixture.name = NULL, with.ci = TRUE)

