#current code, with aa.names = NULL
## add an argument for codon list = NULL
### If codon list = NULL, code will run like normal. Otherwise, will map to codon list 

library(AnaCoDa)

plot.Rcpp_Trace <- function(x, what=c("Mutation", "Selection", "MixtureProbability" ,"Sphi", "Mphi", "Aphi", "Sepsilon", "ExpectedPhi", "Expression","NSEProb","NSERate","InitiationCost","PartitionFunction"), 
                            geneIndex=1, mixture = 1, log.10.scale=F, aa.names = NULL, codon.table = NULL)
{
  if(what[1] == "Mutation")
  {
    plotCodonSpecificParameters(x, mixture, "Mutation", main="Mutation Parameter Traces", aa.names = aa.names)
  }
  if(what[1] == "Selection")
  {
    plotCodonSpecificParameters(x, mixture, "Selection", main="Selection Parameter Traces", aa.names = aa.names)
  }
  if(what[1] == "Alpha")
  {
    plotCodonSpecificParameters(x, mixture, "Alpha", main="Alpha Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names)
  }
  if(what[1] == "Lambda")
  {
    plotCodonSpecificParameters(x, mixture, "Lambda", main="Lambda Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names)
  } 
  
  if(what[1] == "MeanWaitingTime")
  {
    plotCodonSpecificParameters(x, mixture, "MeanWaitingTime", main="Mean Waiting Time Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names)
  }  
  if(what[1] == "VarWaitingTime")
  {
    plotCodonSpecificParameters(x, mixture, "VarWaitingTime", main="Variance Waiting Time Parameter Traces", ROC.or.FONSE=FALSE, aa.names = aa.names)
  }  
  if(what[1] == "NSEProb")
  {
    plotCodonSpecificParameters(x, mixture, "NSEProb", main="Nonsense Error Probability Parameter Traces", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names)
  }  
  if(what[1] == "MixtureProbability")
  {
    plotMixtureProbability(x)
  }
  if(what[1] == "Sphi")
  {
    plotHyperParameterTrace(x, what = what[1]) 
  }
  if(what[1] == "Mphi") 
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "Aphi")
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "InitiationCost")
  {
    plotFONSEHyperParameterTrace(x,what=what[1])
  }
  if(what[1] == "PartitionFunction")
  {
    plotPANSEHyperParameterTrace(x,what=what[1])
  }
  if(what[1] == "Sepsilon") 
  {
    plotHyperParameterTrace(x, what = what[1])
  }
  if(what[1] == "ExpectedPhi")
  {
    plotExpectedPhiTrace(x)
  }
  if(what[1] == "Expression")
  {
    plotExpressionTrace(x, geneIndex)
  }
  if(what[1] == "AcceptanceRatio")
  {
    plotAcceptanceRatios(x)
  }
  if(what[1] == "NSERate")
  {
    plotCodonSpecificParameters(x, mixture, "NSERate", main="NSERate", ROC.or.FONSE=FALSE, log.10.scale=log.10.scale, aa.names = aa.names)
  }  
}
  

#current code for making sure aa.names passed are valid
##currently, aa.names <- aminoAcids()
  #map to codons from codon table
if( is.null(aa.names) ) {
  aa.names <- aminoAcids()
} else {
  aa.match <- (aa.names %in% aminoAcids())
  ## test to ensure there's no aa being called that don't exist in trace
  aa.mismatch <- aa.names[!aa.match]
  if(length(aa.mismatch) > 0){
    warning("Members ", aa.mismatch, "of aa.names argument absent from trace object and will be excluded.",
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
  aa.names <- aa.names[aa.match]
}
with.ref.codon <- ifelse(ROC.or.FONSE, TRUE, FALSE)


#codon table to make sorting easier
for(aa in aa.names)
{
  aa.names <- aminoAcids()
  codon.table <- (lapply(aa.names, function(x) {codons = AAToCodon(x);n_syn = length(codons);  ct_ending = grepl("[CT]$", codons); return(data.frame(aa = x, n_syn, codons, ct_ending))}))
  codon.table <- do.call(rbind, codon.table)
}

##filtering codon table by either AG or CT ending
for(aa in aa.names)
{
  aa.by.table <- codon.table$aa
  aaEndingCT <- aa.by.table[which(codon.table$ct_ending=="TRUE")]
  aaEndingAG <- aa.by.table[which(codon.table$ct_ending=="FALSE")]
}

#argument for passing to codon list
##if codon.list = NULL, use aa.names
if( is.null(codon.table) ) {
  aa.names <- aminoAcids()
  warning("codon.table null, using default aminoAcids() list")
} else { aa.names <- aa.by.table
}


#loading example genome, parameter, and trace to explore arguments within them for obtaining amino acids
#NOTE: noticed typo on help page for "GetCodonCounts" <- "fiven" instead of "given" in description
genome <- initializeGenomeObject(file = "orf_coding.fasta")
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome)))

##some possibilities
genome$getCodonCountsPerGene()
###inside parentheses get notice "expecting a single string value: [type=S4; extent=1]
###perhaps need to use "getGenes" first?

genome$getGenes()
###does not provide information about what to include within parentheses to get to run

trace <- parameter$getTraceObject()
###within trace object, there are arguments for "getCodonSpecific" acceptance rate trace, parameter trace but no information in help pages
trace$getCodonSpecificAcceptanceRateTraceForAA()
trace$getCodonSpecificParameterTrace()

codon_counts <- getCodonCounts(genome)

getCodonCounts

