#Load external libraries
library (AnaCoDa)
library(ggplot2)

##Library help information
anacodaHelp <- ??AnaCoDa

#Define key parameters
roundInitial <- 1
roundMax <- 4
roundMy <- roundInitial

sampleGenome <- TRUE
sampleSize <- 500
sampleSeed <- 06272001

initialSphi <- 1
num.mixtures <- 1
whichModel <- "ROC"
with.phi <- FALSE

makeFitPlots <- TRUE
compareToPrevious <- TRUE
makeComparisonPlots <- TRUE
samples <- 500

##Set up objects
inputDirectory <- "~/Documents/Research/gilchrist-lab-eb/input"
outputDirectory <- "~/Documents/Research/gilchrist-lab-eb/output"

inputRestartFile <- paste0(outputDirectory, "Restart/rstart.round.", roundMy - 1, ".rst")
outputRestartFile <- paste0(outputDirectory, "Restart/rstart.round.", roundMy, ".rst")

fasta.file <- "orf_coding.fasta"

###Create genome object
genome <- initializeGenomeObject(file = fasta.file)

#Begin loop 
if(makeFitPlots) {
  roundMy <- roundInitial
  if(sampleGenome) {
    genomeSampled <- sampleGenomeFunc(genome = genome, sampleSize = sampleSize, sampleSeed = sampleSeed)
    sampleSeed <- sampleSeed + 1
    length(genomeSampled)
  }
  
  #initialize parameter and model objects
  parameter <- initializeParameterObject(genome = genomeSampled, model = whichModel, sphi = initialSphi, num.mixtures = num.mixtures, gene.assignment = rep(1, length(genomeSampled))) 
  
  
  model <- initializeModelObject(parameter = parameter, model = whichModel, with.phi = with.phi)
  
  if(roundInitial = roundMy) {
    outputFile <-   outputFile <- paste0("Graphs/csp.traces.mutation.round-", roundMy, ".pdf")
    
    trace <- parameter$getTraceObject()
    
    pdf(file = outputFile)
    plot(x = trace, what = "Mutation", samples = samples)
    dev.off()
    
    outputFile <-   outputFile <- paste0("Graphs/csp.traces.selection.round-", roundMy, ".pdf")
    
    pdf(file = outputFile)
    plot(x = trace, what = "Selection", samples = samples)
    dev.off()
    
  } else {
    outputFile <-   outputFile <- paste0("Graphs/csp.traces.mutation.round-", roundMy, ".pdf")
    
    trace <- parameter$getTraceObject()
    
    pdf(file = outputFile)
    plot(x = trace, what = "Mutation", samples = samples)
    dev.off()
    
    outputFile <-   outputFile <- paste0("Graphs/csp.traces.selection.round-", roundMy, ".pdf")
    
    pdf(file = outputFile)
    plot(x = trace, what = "Selection", samples = samples)
    dev.off()
    
  }
  roundMy <- roundMy + 1 
}