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

library(tercen)
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(magrittr)
library(tibble)
library(flowCore)
library(flowClean)

# Minimum amount of cells is 30.000
# QC cutoff value = 10.000


matrix2flowFrame <- function(a_matrix){ 
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
  flowFrame <- flowCore::flowFrame(a_matrix, params)
  
  return(flowFrame)
}

flowclean_QC <- function(flowframe, input.pars){
  qc_list <- flowClean::clean(
    flowframe,
    #filePrefixWithDir = "QC_output",
    vectMarkers = vectormarkers,
    #ext = "fcs",
    binSize = input.pars$binSize,
    nCellCutoff = input.pars$nCellCutoff,
    #announce = TRUE,
    cutoff = input.pars$cutoff,
    #diagnostic = TRUE,
    fcMax = input.pars$fcMax,
    returnVector = TRUE,
    nstable = input.pars$nstable)
  return(ifelse(qc_list >= 10000, "fail", "pass"))
}

ctx <- tercenCtx(workflowId = "36321b1d18264d80a1106708ce034a11",
                 stepId = "6bc01961-4601-420a-9067-eba8a427dddb")

if(ctx$cnames[1] == "filename") {filename <- TRUE
if(ctx$cnames[2] != "Time") stop("Time not detected in the second column.")
}else{filename <- FALSE
if(ctx$cnames[1] != "Time") stop("filename or Time not detected in the top column.")
}

celldf <- ctx %>% dplyr::select(.ri, .ci) 
if(nrow(celldf) != length(table(celldf)))stop("There are multiple values in one of the cells.")

input.pars <- list(
  nstable = ifelse(is.null(ctx$op.value('nstable')), 5, as.double(ctx$op.value('nstable'))),
  fcMax = ifelse(is.null(ctx$op.value('fcMax')),  1.3, as.double(ctx$op.value('fcMax'))),
  nCellCutoff = ifelse(is.null(ctx$op.value('nCellCutoff')),  500, as.double(ctx$op.value('nCellCutoff'))),
  binSize = ifelse(is.null(ctx$op.value('binSize')),  0.01, as.double(ctx$op.value('binSize'))),
  cutoff = ifelse(is.null(ctx$op.value('cutoff')),  10, as.double(ctx$op.value('cutoff')))
)


channels <-ctx$rselect(ctx$rnames[1]) %>% t()
vectorindex <- c(1:length(channels))
fsccolumns <- grepl("FSC", channels)
ssccolumns <- grepl("SSC", channels)
vectormarkers <- vectorindex[ifelse((fsccolumns+ssccolumns) == 0, TRUE, FALSE)]

if(filename == TRUE){
  data <- ctx$as.matrix() %>% t() %>% cbind((ctx$cselect(ctx$cnames[[2]]))) %>% cbind((ctx$cselect(ctx$cnames[[1]])))
}
if(filename == FALSE){
  data <- ctx$as.matrix() %>% t() %>% cbind((ctx$cselect(ctx$cnames[[1]])))
  data$filename <- "singlefile"
}
filenames <- unique(data$filename)
qc_df <- data.frame(matrix(ncol=0, nrow=nrow(data)))

QC_allfiles <- lapply(filenames, function(x) {
  singlefiledata <- data[data$filename == x,]
  if(nrow(singlefiledata)<30000){
    warning(sprintf("%s does not have the required >30.000 events.", x))
    return(rep('failed_QC', nrow(singlefiledata)))}
  else{
    singlefileflowframe <- singlefiledata[1:(ncol(singlefiledata)-1)] %>% as.matrix() %>% matrix2flowFrame()
    singlefileQC_vector <- flowclean_QC(singlefileflowframe, input.pars)
  }
})

qc_df$QC_flag <- do.call(c, QC_allfiles)
flowClean_QC <- cbind(qc_df, .ci = (0:(nrow(qc_df)-1)))
ctx$addNamespace(flowClean_QC) %>% ctx$save()
ctx <- tercenCtx(workflowId = "36321b1d18264d80a1106708ce034a11",
                 stepId = "16291844-2e39-40ba-8891-b7111ebe88ad")

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
