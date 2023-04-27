#create codon table for mapping split amino acids
codons <- c("AAA","AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC","ATG", "ATT", "CAA", "CAC", "CAG", "CAT","CCA","CCC","CCG","CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT","TCA", "TCC", "TCG", "TCT","TGC", "TGG", "TGT", "TTA", "TTC","TTG", "TTT" )

std <- c("K", "N", "K", "N", "T", "T", "T", "T", "R", "Z", "R", "Z", "I", "I",
         "M", "I", "Q", "H", "Q", "H", "P", "P", "P", "P", "R", "R", "R", "R",
         "L", "L", "L", "L", "E", "D", "E", "D", "A", "A", "A", "A", "G",
         "G", "G", "G","V", "V", "V", "V", "Y", "Y", "S", "S", "S", "S", "C",
         "W", "C", "L", "F", "L", "F")

aa <- c("KR", "NY", "KR", "NY", "TR", "TY", "TR", "TY", "RR", "ZY", "RR", "ZY", "IR", "IY",
        "MR", "IY", "QR", "HY", "QR", "HY", "PR", "PY", "PR", "PY", "RR", "RY", "RR", "RY",
        "LR", "LY", "LR", "LY", "ER", "DY", "ER", "DY", "AR", "AY", "AR", "AY", "GR",
        "GY", "GR", "GY", "VR", "VY", "VR", "VY", "YY", "YY", "SR", "SY", "SR", "SY",
        "CY", "WR", "CY", "LR", "FY", "LR", "FY")

aa.to.codon.table <- data.frame(std, aa, codons)

amino.acids <- aa.to.codon.table$aa #the purine/pyrimidine portions of amino acids

if(is.null(aa.to.codon.table)) {
  ## focal flag really exclude.focal, not include.focal
  localAAToCodon <- function(aa, focal = FALSE) {
    AAToCodon(aa, focal=focal)
  } 
} else {
  localAAToCodon <- function(aa, focal = FALSE) {
    if(length(aa) > 1) {
      warning("Error: only one aa may be inputted to AAToCodon", 
              call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
              domain = NULL)
    }
    aa.match <- aa.to.codon.table$aa==aa
    codons <- aa.to.codon.table[aa.match, "codons"]
    return(codons)
    }
  }


localAAToCodon(aa = "AR")
aa.names <- aminoAcids()

plotCodonSpecificParameters.backup(trace = trace, mixture = 1, type = "Mutation", main = "Mutation Parameter Traces", ROC.or.FONSE = TRUE, aa.names = aa.names)


pdf("CSPtest.pdf")
plotCodonSpecificParameters.backup(trace = trace, mixture = 1, type = "Mutation", main = "Mutation Parameter Traces", ROC.or.FONSE = TRUE)
dev.off()



# Called from Plot Trace Object (plot for trace)
# NOT EXPOSED
# 
#' Plot Codon Specific Parameter
#' @param trace An Rcpp trace object initialized with \code{initializeTraceObject}.
#'
#' @param mixture The mixture for which to plot values.
#'
#' @param type A string containing one of the following to graph: \code{Mutation, Selection, Alpha, LambdaPrime, MeanWaitingTime, VarWaitingTime}. 
#'
#' @param main The title of the plot.
#'
#' @param ROC.or.FONSE A logical value determining if the Parameter was ROC/FONSE or not.
#'
#' @param log.10.scale A logical value determining if figures should be plotted on the log.10.scale (default=F). Should not be applied to mutation and selection parameters estimated by ROC/FONSE.
#'
#' @return This function has no return value.
#' 
#' @description Plots a codon-specific set of traces, specified with the \code{type} parameter.
#'
plotCodonSpecificParameters.backup <- plotCodonSpecificParameters <- function(trace, mixture, type="Mutation", main="Mutation Parameter Traces", ROC.or.FONSE=TRUE, log.10.scale=FALSE, aa.names = aminoAcids())
{
  opar <- par(no.readonly = TRUE)
  ### Trace plot.
  if (ROC.or.FONSE)
  {
    nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
                 rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  }else
  {    
    nf <- layout(matrix(c(rep(1, 4), 2:25), nrow = 7, ncol = 4, byrow = TRUE),
                 rep(1, 4), c(2, 8, 8, 8, 8, 8, 8), respect = FALSE) 
  }
  ### Plot title.
  if (ROC.or.FONSE){
    par(mar = c(0, 0, 0, 0))
  }else{
    par(mar = c(1,1,1,1))
  }
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  ### TODO change to groupList -> checks for ROC like model is not necessary!
  
  ## Check to ensure aa.names passed are valid
  ## To have this work with non-standard AA codes, we need to add a genetic.code
  ## argument when calling the function.
  aa.match <- (aa.names %in% aminoAcids())
  ## test to ensure there's no aa being called that don't exist in trace
  aa.mismatch <- aa.names[!aa.match]
  if(length(aa.mismatch) > 0){
    warning("Members ", aa.mismatch, "of aa.names not found in aminoAcids() object and will be excluded.",
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
  aa.names <- aa.names[aa.match]
  
  with.ref.codon <- ifelse(ROC.or.FONSE, TRUE, FALSE)
  for(aa in aa.names)
  { 
    codons <- AAToCodon(aa, with.ref.codon)
    if(length(codons) == 0) next
    if (!ROC.or.FONSE){
      if(aa == "X") next
    }
    if (ROC.or.FONSE){
      if(aa == "X" || aa == "M" || aa == "W") next
    }
    cur.trace <- vector("list", length(codons))
    paramType <- 0
    if(type == "Mutation"){
      ylab <- expression(Delta~"M")
      paramType <- 0
      special <- FALSE
    }else if (type == "Selection"){
      ylab <- expression(Delta~eta)
      paramType <- 1
      special <- FALSE
    }else if (type == "Alpha"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*alpha)
      } else{
        ylab <- expression(alpha)
      }
      paramType <- 0
      special <- FALSE
    }else if (type == "Lambda"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*lambda)
      } else{
        ylab <- expression(lambda)
      }
      paramType <- 1
      special <- FALSE
    }else if (type == "MeanWaitingTime"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*alpha/lambda)
      }else{
        ylab <- expression(alpha/lambda)
      }
      special <- TRUE
    }else if (type == "VarWaitingTime"){
      ylab <- expression(alpha/lambda^"2")
      special <- TRUE
    }else if (type == "NSEProb"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*"Pr(NSE)")
      } else{
        ylab <- expression("E[Pr(NSE)]")
      }
      special <- TRUE
    }else if (type == "VarNSEProb"){
      ylab <- expression("Var[Pr(NSE)]")
      special <- TRUE
    }else if (type == "NSERate"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*"NSERate")
      } else{
        ylab <- expression("NSERate")
      }
      paramType <- 2
      special <- FALSE
    }else{
      stop("Parameter 'type' not recognized! Must be one of: 'Mutation', 'Selection', 'Alpha', 'Lambda', 'MeanWaitingTime', 'VarWaitingTime', 'NSEProb', 'NSERate'.")
    }
    
    for(i in 1:length(codons)){
      if(special){
        tmpAlpha <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0, with.ref.codon)
        tmpLambdaPrime <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1, with.ref.codon)
        
        if (type == "MeanWaitingTime"){          
          cur.trace[[i]] <- tmpAlpha / tmpLambdaPrime
        }else if (type == "VarWaitingTime"){
          cur.trace[[i]] <- tmpAlpha / (tmpLambdaPrime * tmpLambdaPrime)
        } else if (type == "NSEProb" || type == "VarNSEProb"){
          tmpNSERate <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 2, with.ref.codon)
          if (type == "NSEProb")
          {
            cur.trace[[i]] <- tmpNSERate*(tmpAlpha/tmpLambdaPrime)
          } else {
            cur.trace[[i]] <- tmpNSERate*tmpNSERate*(tmpAlpha/(tmpLambdaPrime * tmpLambdaPrime))
          }
        }
        if (log.10.scale)
        {
          cur.trace[[i]] <- log10(cur.trace[[i]])
        }
      }
      else{
        cur.trace[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], paramType, with.ref.codon)
        if (log.10.scale)
        {
          cur.trace[[i]] <- log10(cur.trace[[i]])
        }
      }
    }
    
    cur.trace <- do.call("cbind", cur.trace)
    if(length(cur.trace) == 0) next
    x <- 1:dim(cur.trace)[1]
    xlim <- range(x)
    ylim <- range(cur.trace, na.rm=T)
    
    main.aa <- aa #TODO map to three leter code
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(cur.trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = cur.trace[, i.codon], col = .codonColors[[codons[i.codon]]])
    }
    colors <- unlist(.codonColors[codons])
    legend("topleft", legend = codons, col = colors, 
           lty = rep(1, length(codons)), bty = "n", cex = 0.75)
  }
  par(opar)
} 


# Called from Plot Trace Object (plot for trace)
# NOT EXPOSED
# 
#' Plot Codon Specific Parameter
#' @param trace An Rcpp trace object initialized with \code{initializeTraceObject}.
#'
#' @param mixture The mixture for which to plot values.
#'
#' @param type A string containing one of the following to graph: \code{Mutation, Selection, Alpha, LambdaPrime, MeanWaitingTime, VarWaitingTime}. 
#'
#' @param main The title of the plot.
#'
#' @param ROC.or.FONSE A logical value determining if the Parameter was ROC/FONSE or not.
#'
#' @param log.10.scale A logical value determining if figures should be plotted on the log.10.scale (default=F). Should not be applied to mutation and selection parameters estimated by ROC/FONSE.
#'
#' @return This function has no return value.
#' 
#' @description Plots a codon-specific set of traces, specified with the \code{type} parameter.
#'
plotCodonSpecificParameters <- function(trace, mixture, type="Mutation", main="Mutation Parameter Traces", ROC.or.FONSE=TRUE, log.10.scale=F, aa.names = NULL)
{
  opar <- par(no.readonly = T) 
  ### Trace plot.
  if (ROC.or.FONSE)
  {
    nf <- layout(matrix(c(rep(1, 4), 2:21), nrow = 6, ncol = 4, byrow = TRUE),
                 rep(1, 4), c(2, 8, 8, 8, 8, 8), respect = FALSE)  
  }else
  {    
    nf <- layout(matrix(c(rep(1, 4), 2:25), nrow = 7, ncol = 4, byrow = TRUE),
                 rep(1, 4), c(2, 8, 8, 8, 8, 8, 8), respect = FALSE) 
  }
  ### Plot title.
  if (ROC.or.FONSE){
    par(mar = c(0, 0, 0, 0))
  }else{
    par(mar = c(1,1,1,1))
  }
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  ### TODO change to groupList -> checks for ROC like model is not necessary!
  
  ## Check to ensure aa.names passed are valid
  if( is.null(aa.names) ) {
    aa.names <- aminoAcids()
  } else {
    aa.match <- (aa.names %in% amino.acids)
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
  
  
  for(aa in aa.names)
  { 
    codons <- localAAToCodon(aa, with.ref.codon)
    if(length(codons) == 0) next
    if (!ROC.or.FONSE){
      if(aa == "X") next
    }
    if (ROC.or.FONSE){
      if(aa == "X" || aa == "M" || aa == "W") next
    }
    cur.trace <- vector("list", length(codons))
    paramType <- 0
    if(type == "Mutation"){
      ylab <- expression(Delta~"M")
      paramType <- 0
      special <- FALSE
    }else if (type == "Selection"){
      ylab <- expression(Delta~eta)
      paramType <- 1
      special <- FALSE
    }else if (type == "Alpha"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*alpha)
      } else{
        ylab <- expression(alpha)
      }
      paramType <- 0
      special <- FALSE
    }else if (type == "Lambda"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*lambda)
      } else{
        ylab <- expression(lambda)
      }
      paramType <- 1
      special <- FALSE
    }else if (type == "MeanWaitingTime"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*alpha/lambda)
      }else{
        ylab <- expression(alpha/lambda)
      }
      special <- TRUE
    }else if (type == "VarWaitingTime"){
      ylab <- expression(alpha/lambda^"2")
      special <- TRUE
    }else if (type == "NSEProb"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*"Pr(NSE)")
      } else{
        ylab <- expression("E[Pr(NSE)]")
      }
      special <- TRUE
    }else if (type == "VarNSEProb"){
      ylab <- expression("Var[Pr(NSE)]")
      special <- TRUE
    }else if (type == "NSERate"){
      if (log.10.scale)
      {
        ylab <- expression("log"[10]*"NSERate")
      } else{
        ylab <- expression("NSERate")
      }
      paramType <- 2
      special <- FALSE
    }else{
      stop("Parameter 'type' not recognized! Must be one of: 'Mutation', 'Selection', 'Alpha', 'Lambda', 'MeanWaitingTime', 'VarWaitingTime', 'NSEProb', 'NSERate'.")
    }
    
    for(i in 1:length(codons)){
      if(special){
        tmpAlpha <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0, with.ref.codon)
        tmpLambdaPrime <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1, with.ref.codon)
        
        if (type == "MeanWaitingTime"){          
          cur.trace[[i]] <- tmpAlpha / tmpLambdaPrime
        }else if (type == "VarWaitingTime"){
          cur.trace[[i]] <- tmpAlpha / (tmpLambdaPrime * tmpLambdaPrime)
        } else if (type == "NSEProb" || type == "VarNSEProb"){
          tmpNSERate <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 2, with.ref.codon)
          if (type == "NSEProb")
          {
            cur.trace[[i]] <- tmpNSERate*(tmpAlpha/tmpLambdaPrime)
          } else {
            cur.trace[[i]] <- tmpNSERate*tmpNSERate*(tmpAlpha/(tmpLambdaPrime * tmpLambdaPrime))
          }
        }
        if (log.10.scale)
        {
          cur.trace[[i]] <- log10(cur.trace[[i]])
        }
      } # end if(special)
      else{
        cur.trace[[i]] <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], paramType, with.ref.codon)
        if (log.10.scale)
        {
          cur.trace[[i]] <- log10(cur.trace[[i]])
        }
      }
    }
    
    cur.trace <- do.call("cbind", cur.trace)
    if(length(cur.trace) == 0) next
    x <- 1:dim(cur.trace)[1]
    xlim <- range(x)
    ylim <- range(cur.trace, na.rm=T)
    
    main.aa <- aa #TODO map to three leter code
    plot(NULL, NULL, xlim = xlim, ylim = ylim,
         xlab = "Samples", ylab = ylab, main = main.aa)
    plot.order <- order(apply(cur.trace, 2, sd), decreasing = TRUE)
    for(i.codon in plot.order){
      lines(x = x, y = cur.trace[, i.codon], col = .codonColors[[codons[i.codon]]])
    }
    colors <- unlist(.codonColors[codons])
    legend("topleft", legend = codons, col = colors, 
           lty = rep(1, length(codons)), bty = "n", cex = 0.75)
  }
  par(opar)
} 


#######################



