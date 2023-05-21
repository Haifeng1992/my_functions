my_zsore <- function(data_column) {
  data_column = (data_column-mean(data_column))/sd(data_column)
  return(data_column)
}

## preselect gene expression explaining 50% variations over all samples
low_expressed_gene_filter = function(log_data, percentage){
  var.x1 <- apply(t(log_data), 2, var)
  var.sort <- sort(var.x1, decreasing=TRUE)
  sum.a <- cumsum(var.sort)
  half.a <- var.sort[which(sum.a>(sum.a[length(var.x1)]/percentage))[1]]
  x1.cov <- t(log_data)[,var.x1>=half.a]
  x1.cov = log2(x1.cov+1)
  return(x1.cov)
}
try_rfe_if_work <- function(x.train, y.train){
  tryCatch(
    # This is what I want to do...
    {
      tmp1 = caret::rfe(x.train, y.train, 
                        rfeControl = rfeControl(functions = caretFuncs, 
                                                method = "cv",
                                                numbers = 2),
                        metric = "AUC",
                        #sizes = c(1:20), method = "svmRadial")
                        sizes = c(2, 5, 10), method = "svmRadial")
      return(tmp1)
    },
    # ... but if an error occurs, tell me what happened: 
    error=function(error_message) {
      message("This is my custom message.")
      message("And below is the error message from R:")
      message(error_message)
      return(NA)
    }
  )
}

initialization <- function() {
  library("glmnet")
  #library("pcr")
  library("RColorBrewer")
  library("ComplexHeatmap")
  library("circlize")
  library("tidyverse")
  #library("hrbrthemes")
  library("viridis")
  library("stats")
  library("caret")
  library("pROC")
  #library("DMwR")
  #library("PerfMeas")
  #library("qpgraph")
  library("PRROC")
  library("e1071")
  library("ROCR")
  library("unikn")
  library("animation")
  library("MASS")
  library("ggplot2")
  library("reshape2")
  library(caret)
  library(ggfortify)
  library(gridExtra)
  library("dplyr")
  library("randomForest")
  library("edgeR")
  library("readxl")
  library("ggpubr")
}

label_to_numeric <- function(test_label, predicted_results){
  test_label = as.character(test_label)
  test_label[test_label == "IBC"] = 0
  test_label[test_label == "DCIS"] = 1
  test_label = as.numeric(test_label)
  #test_label = as.factor(test_label)
  
  predicted_results = as.character(predicted_results)
  predicted_results[predicted_results == "IBC"] = 0
  predicted_results[predicted_results == "DCIS"] = 1
  predicted_results = as.numeric(predicted_results)
  #predicted_results = as.factor(predicted_results)
  
  tmp = list()
  tmp[[1]] = test_label
  tmp[[2]] = predicted_results
  
  return(tmp)
}


map_between_probe_and_gene <- function(probe_list, to_gene_or_to_probe) {
  #probe_list = GeneList
  load("/Users/xuhaifeng/Documents/phd1/data/tcga_entrez_gene_pos.Rdata")
  
  if(to_gene_or_to_probe == "to_gene"){
    tmp = tcga_entrez_gene_pos
    tmp = tmp[is.element(tmp$EntrezID, probe_list),]
    tmp = tmp$Gene
    return(tmp)
  }else if(to_gene_or_to_probe == "to_probe"){
    tmp = tcga_entrez_gene_pos
    tmp = tmp[is.element(tmp$Gene, probe_list),]
    tmp = tmp$EntrezID
    return(tmp)
  }else{
    print("wrong input! to_gene_or_to_probe must be set to 'to_gene' or 'to_probe' ")
  }
  
}