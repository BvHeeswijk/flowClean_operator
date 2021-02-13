library(tercen)
library(dplyr)
library(tibble)
library(tidyverse)
library(flowCore)
library(flowClean)

# Minimum amount of cells is 30.000
# QC cutoff value = 10.000


matrix2flowset <- function(a_matrix){ 
  
  minRange <- matrixStats::colMins(a_matrix)
  maxRange <- matrixStats::colMaxs(a_matrix)
  rnge <- maxRange - minRange
  
  df_params <- data.frame(
    name = colnames(a_matrix),
    desc = colnames(a_matrix),
    range = rnge,
    minRange = minRange,
    maxRange = maxRange
  )
  
  params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(params) <- df_params
  Biobase::varMetadata(params) <- data.frame(
    labelDescription = c("Name of Parameter",
                         "Description of Parameter",
                         "Range of Parameter",
                         "Minimum Parameter Value after Transformation",
                         "Maximum Parameter Value after Transformation")
  )
  flowset <- flowCore::flowFrame(a_matrix, params)
  
  return(flowset)
}

synPerturbedcsv <- "f5b451f5-2953-43ed-adbc-56ff1ec4b190"
testcsv <- "6a29c574-2528-4891-9128-8b94760d087a"



ctx <- tercenCtx(workflowId = "affe1383f9318d3b8141a4bfa700a359",
                 stepId = testcsv)

time <- ctx$cselect(ctx$cnames[[1]])
data <- ctx$as.matrix() %>% t()
data <- as.matrix(cbind(data, time))
# Testing if the amount of input cells are >30000
if (nrow(data)<30000) {stop("flowClean requires >30000 cells to be measured.")}

fc_frame <- matrix2flowset(data)

cutoffvar = ifelse((ctx$op.value('cutoff') == 10), "median", as.double(ctx$op.value('cutoff')))

qc_list <- 
  flowClean::clean(
    fc_frame,
    #filePrefixWithDir = "QC_output",
    #vectMarkers = "markercolumns",
    #ext = "fcs",
    binSize = 0.01,
    nCellCutoff = 500,
    #announce = TRUE,
    cutoff = cutoffvar,
    #diagnostic = TRUE,
    fcMax = 1.3,
    returnVector = TRUE,
    nstable = 5
  )

flag <- ifelse(qc_list >= 10000, "fail", "pass")
qc_df <- data.frame(flag, .ci = (0:(length(flag)-1)))
qc_result <- ctx$addNamespace(qc_df)
ctx$save(qc_result)
