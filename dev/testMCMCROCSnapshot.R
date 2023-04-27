
#loading libraries
library(AnaCoDa)
library(testthat)

rm(list=ls(all.names=TRUE))
context("MCMC with ROC")

fileName = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
expressionFile = file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR_phi_withPhiSet.csv")
selectionMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_1.csv")
selectionHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "selection_2.csv")
mutationMainFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_1.csv")
mutationHtFile = file.path("UnitTestingData", "testMCMCROCFiles", "mutation_2.csv")
mcmcSaveFile <- file.path("UnitTestingOut", "testMCMCROCobject.Rda")

sphi_init <- c(1,1)
numMixtures <- 2
mixDef <- "allUnique"

samples <- 10
thinning <- 10
adaptiveWidth <- 10
divergence.iteration <- 0

set.seed(446141)

#creating data frame containing possibilities for estimating or fixing phi, dM, and dE
##true or false factor is created simply to allow me to replicate the list more easily
T_or_F_factor <- c(TRUE,FALSE)
T_or_F <- data.frame(est_phi = rep(T_or_F_factor, each = 4),
                     est_dM = c(rep(T, 2), rep(F,2)),
                     est_dE = c(T,F))

##objects containing true and false possibilities from columns of T_or_F data frame
est_phi <- T_or_F$est_phi
est_dM <- T_or_F$est_dM
est_dE <- T_or_F$est_dE

#for loop begins

for(est_phi) {
  for(est_dM) {
    for(est_dE) {
        
      genome <- initializeGenomeObject(file = fileName, observed.expression.file = expressionFile, match.expression.by.id=FALSE)
      
      geneAssignment <- sample(c(1,2), size = length(genome), replace = TRUE, prob = c(0.3, 0.7)) #c(rep(1,500), rep(2,500))
      
      parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
      parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
      parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)
      
      model <- initializeModelObject(parameter, "ROC", with.phi = est_phi) 
      
      outFile <- file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")
      
      sink(outFile)
      runMCMC(mcmc, genome, model, 1, divergence.iteration)
      sink()
      
      #get log likelihood trace when 'with phi' is true for tenth sample
      actual_loglik <- mcmc$getLogLikelihoodTrace()[samples])
      
###testing begins for log likelihood
      test_that("get expected log likelihood value when 'with phi' is true, not expected log likelihood when 'with phi' is false", {
        parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, split.serine = TRUE, mixture.definition = mixDef)
        parameter$initSelectionCategories(c(selectionMainFile, selectionHtFile), 2,F)
        parameter$initMutationCategories(c(mutationMainFile, mutationHtFile), 2,F)
        
        model <- initializeModelObject(parameter, "ROC", with.phi = TRUE) 
        
        outFile <- file.path("UnitTestingOut", "testMCMCROCLogPhi.txt")
        
        sink(outFile)
        runMCMC(mcmc, genome, model, 1, divergence.iteration)
        sink()
        
        expect_loglik <- mcmc$getLogLikelihoodTrace([samples])
        expect_value(expect_loglik, actual_loglik )
      }
      
      
      
    }
  }
}







             
