# Title:  rPHG flow tests - DO NOT RUN
# Author: Brandon Monier

# === Preamble (change depending on user) ===========================

## Load packages ----
library(rPHG)


## Get config path ----
configFilePath <- system.file(
    "extdata",
    "configSQLite.txt",
    package = "rPHG"
)

configFilePath <- "/home/bm646/Temporary/test_rphg/wheat_cap_2019_dbs/congif_sql_btm.txt"


## Various graph building methods ----
method1 <- "GATK_PIPELINE"
method2 <- "CONSENSUS"



# === Tests - get PHG object and functions ==========================

## Start logger ----
rPHG::startLogger(fullPath = "/home/bm646/Temporary/test_rphg/")


## Get methods ----
phgMethods <- rPHG::showPHGMethods(configFile = configFilePath)


## Build the PHG (see rPHG graphBuilder code) ----
phgObject <- rPHG::graphBuilder(
    configFile = configFilePath,
    methods = phgMethods$method_name[5], # CONSENSUS
    chrom = NULL,
    includeSequence = FALSE,
    includeVariant = FALSE
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
    rPHG::numHaploPerRange(start = 153e6, end = 154e6) %>%
    rPHG::plotNumHaplo()










