#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   btm_rphg_tests.R
# Description:   rPHG tests for the Wheat PHG workshop (2020)
# Author:        Brandon Monier
# Created:       2020-02-27 at 17:49:29
# Last Modified: 2020-02-27 at 18:25:07
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to test out new ideas and
#    workflows for the wheat PHG workshop (2020).
#
#    This script makes use of the test small seq database in maize.
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages ----
library(rPHG)


## Make configuration file ----
rPHG:::configFileMaker(
    dbName = "phgSmallSeq.db",
    dbType = "sqlite",
    exportPath = paste0(getwd(), "/data")
)


## Get config file name
configFilePath <- paste0("data/configFilePHG.txt")


## Various graph building methods ----
method1 <- "GATK_PIPELINE"
method2 <- "CONSENSUS"



# === Tests - get PHG object and functions ==========================

## Start logger ----
rPHG::startLogger(fullPath = paste0(getwd(), "/output"))


## Get methods ----
phgMethods <- rPHG::showPHGMethods(configFile = configFilePath)


## Build the PHG (see rPHG graphBuilder code) ----
phgObject <- rPHG::graphBuilder(
    configFile = configFilePath,
    methods = method2 # CONSENSUS
)


## Get path matrix ----
phgPathMat <- rPHG::pathMatrix(
    configFile = configFilePath,
    pathDir = "/home/bm646/Temporary/test_rphg/phgSmallSeq/hap_count_best_path/"
)


## Get number of haplotypes in range ----
numHaplo <- rPHG::numHaploPerRange(phgObject = phgObject)



# === New tests (Wheat Hackathon - 2020-02-20) ======================

## Get path matrix ----
mat_test <- rPHG::pathsForMethod(
    configFile = configFilePath,
    pathMethod = "PATH_METHOD"
)


## Get read mapping data frame ----
rmDF <- rPHG::readMappingsForLineName(
    configFile = configFilePath,
    lineName = "RefA1_gbs",
    readMappingMethodName = "HAP_COUNT_METHOD",
    haplotypeMethodName = "CONSENSUS"
)



# === Visualization - 2020-02-25 ====================================

## Number of haplotypes per ref range
phgObject %>%
    rPHG::numHaploPerRange(start = 4e3, end = 5e4) %>%
    rPHG::plotNumHaplo()










