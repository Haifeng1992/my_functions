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


svm_models_cll <- function(data, label, genes, rfe_or_not, folds) {
  set.seed(4399)
  existed_genes = intersect(genes, rownames(data))
  if(length(genes) == length(existed_genes)){
    all_genes_existed = "all genes existed"
  }else{
    all_genes_existed = genes[is.element(genes, rownames(data)) == FALSE]
  }

  #pFLOW_data = data[existed_genes,]
  pFLOW_data = data[existed_genes,]
  pFLOW_data_t = as.data.frame(t(pFLOW_data))
  #pFLOW_data_t = as.data.frame(t(pFLOW_data))
  label = as.character(label)
  label = as.factor(label)
  
  Response = label
  pflow_baseline = cbind(Response, pFLOW_data_t)
  pflow_baseline$Response = as.character(pflow_baseline$Response)
  pflow_baseline$Response = factor(pflow_baseline$Response, ordered = TRUE)
  test_folds = createFolds(y = pflow_baseline$Response, k = folds, list = TRUE, returnTrain = FALSE)
  
  GeneList = list()
  auc_list = list()
  accuracy_list = list()
  ba_list = list()
  all_results_list = list()
  
  if(rfe_or_not == TRUE){
    for (i in 1:folds) {
      test_data <- pflow_baseline[as.numeric(test_folds[[i]]), ]
      train_data = pflow_baseline[-as.numeric(test_folds[[i]]), ]
      
      x.train = train_data[,2:ncol(train_data)]
      y.train = train_data$Response
      y.train = as.character(y.train)
      y.train = factor(y.train, ordered = TRUE)
      
      x.test = test_data[,2:ncol(test_data)]
      y.test = test_data$Response
      y.test = as.character(y.test)
      y.test = factor(y.test, ordered = TRUE)
      
      ################## RFE feature selection  #######################################
      tmp1 = try(caret::rfe(x = x.train, y = y.train, 
                            rfeControl = rfeControl(functions = caretFuncs,
                                                    number = folds-1),
                            #sizes = c(1:20), method = "svmRadial")
                            sizes = c(2, 5, 10, 20), method = "svmRadial")
      )
      if(is.list(tmp1) == TRUE){
        a = tmp1[["optVariables"]]
        genes = a
        GeneList[[i]] = genes
        ################## Select Certain Genes #######################################
        # 100 tests common genes
        
        
        x.train = x.train[,genes]
        x.test = x.test[,genes]
        
        model = svm(x.train, y.train)
        pred_binary = predict(model,x.test)
        pred_binary = factor(pred_binary, ordered = TRUE)
        
        CM = confusionMatrix(pred_binary, y.test)
        accuracy_list[[i]] = CM$overall[1]
        ba_list[[i]] = CM$byClass[11]
        
        model = svm(x.train, y.train, probability = TRUE)
        pred_binary = predict(model,x.test, probability = TRUE)
        y.test = as.character(y.test)
        y.test[y.test == "PR"] = 1
        y.test[y.test == "SD"] = 0
        pred_binary = attr(pred_binary,"probabilities")
        pred_binary = as.data.frame(pred_binary)
        pred_binary = pred_binary$PR
        
        tt = roc(response = y.test, predictor =pred_binary, direction = "<")
        #plot(tt)
        auc_list[[i]] = auc(tt)
      }else{
        model = svm(x.train, y.train)
        pred_binary = predict(model,x.test)
        pred_binary = factor(pred_binary, ordered = TRUE)
        
        CM = confusionMatrix(pred_binary, y.test)
        accuracy_list[[i]] = CM$overall[1]
        ba_list[[i]] = CM$byClass[11]
        
        model = svm(x.train, y.train, probability = TRUE)
        pred_binary = predict(model,x.test, probability = TRUE)
        y.test = as.character(y.test)
        y.test[y.test == "PR"] = 1
        y.test[y.test == "SD"] = 0
        pred_binary = attr(pred_binary,"probabilities")
        pred_binary = as.data.frame(pred_binary)
        pred_binary = pred_binary$PR
        
        tt = roc(response = y.test, predictor =pred_binary, direction = "<")
        plot(tt)
        auc_list[[i]] = auc(tt)
        
        GeneList[[i]] = "no gene selected"
      }
    }
  }else{
    for (i in 1:folds) {
      test_data <- pflow_baseline[as.numeric(test_folds[[i]]), ]
      train_data = pflow_baseline[-as.numeric(test_folds[[i]]), ]
      
      x.train = train_data[,2:ncol(train_data)]
      y.train = train_data$Response
      y.train = as.character(y.train)
      y.train = factor(y.train, ordered = TRUE)
      
      x.test = test_data[,2:ncol(test_data)]
      y.test = test_data$Response
      y.test = as.character(y.test)
      y.test = factor(y.test, ordered = TRUE)
      
      model = svm(x.train, y.train)
      pred_binary = predict(model,x.test)
      pred_binary = factor(pred_binary, ordered = TRUE)
      
      CM = confusionMatrix(pred_binary, y.test)
      accuracy_list[[i]] = CM$overall[1]
      ba_list[[i]] = CM$byClass[11]
      
      model = svm(x.train, y.train, probability = TRUE)
      pred_binary = predict(model,x.test, probability = TRUE)
      y.test = as.character(y.test)
      y.test[y.test == "PR"] = 1
      y.test[y.test == "SD"] = 0
      pred_binary = attr(pred_binary,"probabilities")
      pred_binary = as.data.frame(pred_binary)
      pred_binary = pred_binary$PR
      
      y.test = as.numeric(y.test)
      
      tt = roc(response = y.test, predictor =pred_binary, direction = "<")
      plot(tt)
      auc_list[[i]] = auc(tt)
      
      GeneList[[i]] = "no gene selected"
    }
  }

  kkk = unlist(accuracy_list)
  all_results_list[[1]] = mean(kkk)
  
  this_ba = unlist(ba_list)
  all_results_list[[2]] = mean(this_ba)
  
  ccc = unlist(auc_list)
  all_results_list[[3]] = mean(ccc)
  
  these_genes = unlist(GeneList)
  all_results_list[[4]] = these_genes
  
  all_results_list[[5]] = all_genes_existed
  names(all_results_list)[1:5] = c("accuracy", "balanced_accuracy", "AUC", "RFE_genes", "missing_genes")
  return(all_results_list)
}

svm_models_cll_sub <- function(data, label, genes, folds, tuning) {
  set.seed(23432)
  existed_genes = intersect(genes, colnames(data))
  if(length(genes) == length(existed_genes)){
    all_genes_existed = "all genes existed"
  }else{
    all_genes_existed = genes[is.element(genes, rownames(data)) == FALSE]
  }
  
  #pFLOW_data = data[existed_genes,]
  pFLOW_data = data[,existed_genes]
  pFLOW_data_t = pFLOW_data
  #pFLOW_data_t = as.data.frame(t(pFLOW_data))
  label = as.character(label)
  label = as.factor(label)
  
  Response = label
  pflow_baseline = cbind(Response, pFLOW_data_t)
  pflow_baseline$Response = as.character(pflow_baseline$Response)
  pflow_baseline$Response = factor(pflow_baseline$Response, ordered = TRUE)
  test_folds = createFolds(y = pflow_baseline$Response, k = folds, list = TRUE, returnTrain = FALSE)
  
  GeneList = list()
  auc_list = list()
  accuracy_list = list()
  ba_list = list()
  all_results_list = list()
  for (i in 1:folds) {
    test_data <- pflow_baseline[as.numeric(test_folds[[i]]), ]
    train_data = pflow_baseline[-as.numeric(test_folds[[i]]), ]
    
    x.train = train_data[,2:ncol(train_data)]
    y.train = train_data$Response
    y.train = as.character(y.train)
    y.train = factor(y.train, ordered = TRUE)
    
    x.test = test_data[,2:ncol(test_data)]
    y.test = test_data$Response
    y.test = as.character(y.test)
    y.test = factor(y.test, ordered = TRUE)
    
    if(tuning == TRUE){
      obj = e1071::tune(svm, train.x = x.train, train.y = y.train,
                        ranges = list(gamma = 2^(-1:1), cost = 2^(2:6)),
                        tunecontrol = tune.control(nrepeat = 5, cross = 2)
      )
      model = svm(x.train, y.train, probability = TRUE, cost = obj$best.parameters$cost, gamma = obj$best.parameters$gamma)
    }else{
      model = svm(x.train, y.train, probability = TRUE)
    }
    
    #model = svm(x.train, y.train)
    pred_binary = predict(model,x.test)
    pred_binary = factor(pred_binary, ordered = TRUE)
    
    CM = confusionMatrix(pred_binary, y.test, positive = "PR")
    accuracy_list[[i]] = CM$overall[1]
    ba_list[[i]] = CM$byClass[11]
    
    #model = svm(x.train, y.train, probability = TRUE)
    pred_binary = predict(model,x.test, probability = TRUE)
    y.test = as.character(y.test)
    y.test[y.test == "PR"] = 1
    y.test[y.test == "SD"] = 0
    pred_binary = attr(pred_binary,"probabilities")
    pred_binary = as.data.frame(pred_binary)
    pred_binary = pred_binary$PR
    
    y.test = as.numeric(y.test)
    
    tt = roc(response = y.test, predictor =pred_binary, direction = "<")
    plot(tt)
    auc_list[[i]] = auc(tt)
    
    GeneList[[i]] = "no gene selected"
  }
  
  kkk = unlist(accuracy_list)
  all_results_list[[1]] = mean(kkk)
  
  this_ba = unlist(ba_list)
  all_results_list[[2]] = mean(this_ba)
  
  ccc = unlist(auc_list)
  all_results_list[[3]] = mean(ccc)
  
  these_genes = unlist(GeneList)
  all_results_list[[4]] = these_genes
  
  all_results_list[[5]] = all_genes_existed
  names(all_results_list)[1:5] = c("accuracy", "balanced_accuracy", "AUC", "RFE_genes", "missing_genes")
  return(all_results_list)
}

lasso_models_cll <- function(data, label, genes, folds) {
  
  existed_genes = intersect(genes, rownames(data))
  if(length(genes) == length(existed_genes)){
    all_genes_existed = "all genes existed"
  }else{
    all_genes_existed = genes[is.element(genes, rownames(data)) == FALSE]
  }
  
  pFLOW_data = data[existed_genes,]
  pFLOW_data_t = as.data.frame(t(pFLOW_data))
  Response = label
  pflow_baseline = cbind(Response, pFLOW_data_t)
  pflow_baseline$Response = as.character(pflow_baseline$Response)
  pflow_baseline$Response = factor(pflow_baseline$Response, ordered = TRUE)
  test_folds = createFolds(y = pflow_baseline$Response, k = folds, list = TRUE, returnTrain = FALSE)
  
  GeneList = list()
  auc_list = list()
  accuracy_list = list()
  ba_list = list()
  all_results_list = list()
  

  for (i in 1:folds) {
    test_data <- pflow_baseline[as.numeric(test_folds[[i]]), ]
    train_data = pflow_baseline[-as.numeric(test_folds[[i]]), ]
    
    x.train = train_data[,2:ncol(train_data)]
    y.train = train_data$Response
    y.train = as.character(y.train)
    y.train = factor(y.train, ordered = TRUE)
    
    x.test = test_data[,2:ncol(test_data)]
    y.test = test_data$Response
    y.test = as.character(y.test)
    y.test = factor(y.test, ordered = TRUE)
    
    ################## RFE feature selection  #######################################
    cv_output <- cv.glmnet(x = as.matrix(x.train), y = y.train, 
                           alpha = 1, family = "binomial", type.measure="auc", nfolds = 8)
    #plot(cv_output)
    model= glmnet(x = as.matrix(x.train), y = as.matrix(y.train), 
                  alpha = 1, family = "binomial", type.measure="auc")
    
    #get the common genes
    betaGene = model$beta[,model$lambda==cv_output$lambda.min]
    betaGene = betaGene[betaGene != 0]
    betaGene = as.data.frame(betaGene)
    
    
    
    result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
    label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
    result = as.numeric(result)
    roc = roc(response = y.test, predictor = result, quiet = TRUE)
    #, print.auc = TRUE, plot=TRUE, legacy.axes = TRUE)
    #plot(roc)
    auc = auc(roc)
    CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
    
    auc_list[[i]] = auc
    accuracy_list[[i]] = CM$overall[1]
    ba_list[[i]] = CM$byClass[11]
    GeneList[[i]] = rownames(betaGene)
  }
  
  kkk = unlist(accuracy_list)
  all_results_list[[1]] = mean(kkk)
  
  this_ba = unlist(ba_list)
  all_results_list[[2]] = mean(this_ba)
  
  ccc = unlist(auc_list)
  all_results_list[[3]] = mean(ccc)
  
  these_genes = unlist(GeneList)
  all_results_list[[4]] = these_genes
  
  all_results_list[[5]] = all_genes_existed
  names(all_results_list)[1:5] = c("accuracy", "balanced_accuracy", "AUC", "RFE_genes", "missing_genes")
  return(all_results_list)
}

lasso_models_cll_with_celltypes_only <- function(data, label, genes, folds) {
  
  existed_genes = intersect(genes, rownames(data))
  if(length(genes) == length(existed_genes)){
    all_genes_existed = "all genes existed"
  }else{
    all_genes_existed = genes[is.element(genes, rownames(data)) == FALSE]
  }
  
  pFLOW_data = data[existed_genes,]
  pFLOW_data_t = as.data.frame(t(pFLOW_data))
  Response = label
  pflow_baseline = cbind(Response, pFLOW_data_t)
  pflow_baseline$Response = as.character(pflow_baseline$Response)
  pflow_baseline$Response = factor(pflow_baseline$Response, ordered = TRUE)
  test_folds = createFolds(y = pflow_baseline$Response, k = folds, list = TRUE, returnTrain = FALSE)
  
  GeneList = list()
  auc_list = list()
  accuracy_list = list()
  ba_list = list()
  all_results_list = list()
  
  
  for (i in 1:folds) {
    test_data <- pflow_baseline[as.numeric(test_folds[[i]]), ]
    train_data = pflow_baseline[-as.numeric(test_folds[[i]]), ]
    
    x.train = train_data[,2:ncol(train_data)]
    y.train = train_data$Response
    y.train = as.character(y.train)
    y.train = factor(y.train, ordered = TRUE)
    
    x.test = test_data[,2:ncol(test_data)]
    y.test = test_data$Response
    y.test = as.character(y.test)
    y.test = factor(y.test, ordered = TRUE)
    
    ################## RFE feature selection  #######################################
    cv_output <- cv.glmnet(x = as.matrix(x.train), y = y.train, 
                           alpha = 1, family = "binomial", type.measure="auc", nfolds = 8)
    #plot(cv_output)
    model= glmnet(x = as.matrix(x.train), y = as.matrix(y.train), 
                  alpha = 1, family = "binomial", type.measure="auc")
    
    #get the common genes
    betaGene = model$beta[,model$lambda==cv_output$lambda.min]
    betaGene = betaGene[betaGene != 0]
    betaGene = as.data.frame(betaGene)
    
    
    
    result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
    label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
    result = as.numeric(result)
    roc = roc(response = y.test, predictor = result, quiet = TRUE)
    #, print.auc = TRUE, plot=TRUE, legacy.axes = TRUE)
    #plot(roc)
    auc = auc(roc)
    CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
    
    auc_list[[i]] = auc
    accuracy_list[[i]] = CM$overall[1]
    ba_list[[i]] = CM$byClass[11]
    GeneList[[i]] = rownames(betaGene)
  }
  
  kkk = unlist(accuracy_list)
  all_results_list[[1]] = mean(kkk)
  
  this_ba = unlist(ba_list)
  all_results_list[[2]] = mean(this_ba)
  
  ccc = unlist(auc_list)
  all_results_list[[3]] = mean(ccc)
  
  these_genes = unlist(GeneList)
  all_results_list[[4]] = these_genes
  
  all_results_list[[5]] = all_genes_existed
  names(all_results_list)[1:5] = c("accuracy", "balanced_accuracy", "AUC", "RFE_genes", "missing_genes")
  return(all_results_list)
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

ten_fold_cv_greedy_feature <- function(tmp_data){
  LumB = tmp_data[tmp_data$subtype == "LumB",]
  LumA = tmp_data[tmp_data$subtype == "LumA",]
  Her2 = tmp_data[tmp_data$subtype == "Her2",]
  Basal = tmp_data[tmp_data$subtype == "Basal",]
  Normal = tmp_data[tmp_data$subtype == "Normal",]
  
  test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)
  
  auc_results = matrix(0,1,10)
  
  for (i in 1:10) {
    test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                       LumA[as.numeric(test_index_LumA[[i]]), ],
                       Her2[as.numeric(test_index_Her2[[i]]), ],
                       Basal[as.numeric(test_index_Basal[[i]]), ],
                       Normal[as.numeric(test_index_Normal[[i]]), ]
    )
    train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                        LumA[as.numeric(-test_index_LumA[[i]]), ],
                        Her2[as.numeric(-test_index_Her2[[i]]), ],
                        Basal[as.numeric(-test_index_Basal[[i]]), ],
                        Normal[as.numeric(-test_index_Normal[[i]]), ]
    ) 
    model = lda(x = as.matrix(train_data[,3]), grouping = train_data$label)
    results = predict(model, as.matrix(test_data[,3]), type = "class",)
    
    labels = label_to_numeric(test_data$label, results$class)
    
    roc = roc(response = labels[[1]], predictor = labels[[2]], quiet = TRUE)
    auc = auc(roc)
    auc_results[1,i] = auc
  }
  auc_mean = mean(auc_results)
  
  return(auc_mean)
}

ten_fold_cv_greedy_feature_prauc <- function(tmp_data){
  LumB = tmp_data[tmp_data$subtype == "LumB",]
  LumA = tmp_data[tmp_data$subtype == "LumA",]
  Her2 = tmp_data[tmp_data$subtype == "Her2",]
  Basal = tmp_data[tmp_data$subtype == "Basal",]
  Normal = tmp_data[tmp_data$subtype == "Normal",]
  
  test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)
  
  auc_results = matrix(0,1,10)
  
  for (i in 1:10) {
    test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                       LumA[as.numeric(test_index_LumA[[i]]), ],
                       Her2[as.numeric(test_index_Her2[[i]]), ],
                       Basal[as.numeric(test_index_Basal[[i]]), ],
                       Normal[as.numeric(test_index_Normal[[i]]), ]
    )
    train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                        LumA[as.numeric(-test_index_LumA[[i]]), ],
                        Her2[as.numeric(-test_index_Her2[[i]]), ],
                        Basal[as.numeric(-test_index_Basal[[i]]), ],
                        Normal[as.numeric(-test_index_Normal[[i]]), ]
    ) 
    model = lda(x = as.matrix(train_data[,3]), grouping = train_data$label)
    results = predict(model, as.matrix(test_data[,3]), type = "class",)
    
    y.test = test_data$label
    pr_label = vector(mode = "logical", length = length(y.test))
    for (pr_i in 1:length(y.test)) {
      if(y.test[pr_i] == as.factor("DCIS")){
        pr_label[pr_i] = TRUE
      }
    }
    
    wpr<-pr.curve(results$posterior[,1], weights.class0 = pr_label, curve = TRUE,
                  max.compute = T, min.compute = T, rand.compute = T)
    au_prc = wpr$auc.integral
    auc_results[1,i] = au_prc
  }
  auc_mean = mean(auc_results)
  
  return(auc_mean)
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

nested_cv <- function(expr_data, type, purpose) {
  #probe_list = GeneList
  #expr_data = mlbased
  #expr_data = traditional
  
  LumB = expr_data[expr_data$subtype == "LumB",]
  LumA = expr_data[expr_data$subtype == "LumA",]
  Her2 = expr_data[expr_data$subtype == "Her2",]
  Basal = expr_data[expr_data$subtype == "Basal",]
  Normal = expr_data[expr_data$subtype == "Normal",]
  
  test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)
  
  # initialization of the nested cross-validation
  results = matrix(0,5,10)
  mse_results = matrix(0,1,10)
  GeneList = list()
  i = 1
  par(mfrow=c(3,4))
  #par(mfrow=c(1,1))
  tmp = matrix(as.factor(0), 1, 3)
  #my_pal <- seecol(pal_unikn_light)
  
  sens_new = list()
  spec_new = list()
  precision_new = list()
  recall_new = list()
  
  #### nested cross-validation ####
  for (i in 1:10) {
    rando = sample.int(100,1)
    set.seed(rando)
    test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                       LumA[as.numeric(test_index_LumA[[i]]), ],
                       Her2[as.numeric(test_index_Her2[[i]]), ],
                       Basal[as.numeric(test_index_Basal[[i]]), ],
                       Normal[as.numeric(test_index_Normal[[i]]), ]
    )
    train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                        LumA[as.numeric(-test_index_LumA[[i]]), ],
                        Her2[as.numeric(-test_index_Her2[[i]]), ],
                        Basal[as.numeric(-test_index_Basal[[i]]), ],
                        Normal[as.numeric(-test_index_Normal[[i]]), ]
    )   
    #test_data <- expr_data[as.numeric(test_index[[i]]), ]
    #train_data <- expr_data[-as.numeric(test_index[[i]]), ]
    x.test = test_data[, 3:ncol(expr_data)]
    y.test = test_data[, 1]
    x.train = train_data[, 3:ncol(expr_data)]
    y.train = train_data[, 1]
    
    y.train = as.character(y.train)
    y.train[y.train == "IBC"] = 0
    y.train[y.train == "IDC"] = 0
    y.train[y.train == "DCIS"] = 1
    y.train = as.numeric(y.train)
    y.train = as.factor(y.train)
    y.test = as.character(y.test)
    y.test[y.test == "IBC"] = 0
    y.test[y.test == "IDC"] = 0
    y.test[y.test == "DCIS"] = 1
    y.test = as.numeric(y.test)
    
    
    
    
    if(type == "normal"){
      pr_label = vector(mode = "logical", length = length(y.test))
      for (pr_i in 1:length(y.test)) {
        if(y.test[pr_i] == 1){
          pr_label[pr_i] = TRUE
        }
      }
      y.test = as.factor(y.test)
      cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), 
                             alpha = 1, family = "binomial", type.measure="auc")
      #plot(cv_output)
      model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                    alpha = 1, family = "binomial", type.measure="auc")
      
      #get the common genes
      betaGene = model$beta[,model$lambda==cv_output$lambda.min]
      betaGene = betaGene[betaGene != 0]
      betaGene = as.data.frame(betaGene)
      GeneList[[i]] = rownames(betaGene)
      
      result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
      label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
      result = as.numeric(result)
    }else if(type == "one_feature"){
      pr_label = vector(mode = "logical", length = length(y.test))
      for (pr_i in 1:length(y.test)) {
        if(y.test[pr_i] == 1){
          pr_label[pr_i] = TRUE
        }
      }
      data =cbind(y.train, x.train)
      colnames(data) = c("label", "ROR")
      data =as.data.frame(data)
      
      model = stats::glm(label~.,family=gaussian, data=data)
      
      test_data =cbind(y.test, x.test)
      test_data = as.data.frame(test_data)
      colnames(test_data) = c("label", "ROR")
      result = predict(model, test_data, type= "response")
      label_result = result
      label_result[label_result>=0.5] = 1
      label_result[label_result<0.5] = 0
    }else{
      print("the type should be 'normal' or 'one_feature'")
      return()
    }
    
    #calculate average mse for this fold
    mse = cbind(as.data.frame(result), as.data.frame(as.numeric(as.character(y.test))))
    temp_c = as.data.frame(matrix(0,nrow(mse),1))
    mse = cbind(mse, temp_c)
    for (j in 1:nrow(mse)) {
      mse[j,3] = (mse[j,1] - mse[j,2])^2
    }
    mse_results[1,i] = mean(mse[j,3])
    
    roc = roc(response = y.test, predictor = result, quiet = TRUE)
    
    sens_new[[i]] = roc$sensitivities
    spec_new[[i]] = roc$specificities
    
    
    #, print.auc = TRUE, plot=TRUE, legacy.axes = TRUE)
    #plot(roc)
    auc = auc(roc)
    
    #cutoff table with ROCR package
    #pred <- prediction(result, y.test)
    #perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    #plot(perf, col=rainbow(10))
    #str(perf)
    #cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
    #                      tpr=perf@y.values[[1]])
    
    #PR curve
    #if(type == "normal"){
      wpr<-pr.curve(result, weights.class0 = pr_label, curve = TRUE,
                  max.compute = T, min.compute = T, rand.compute = T)
      au_prc = wpr$auc.integral
      #if (i == 1){
      #  plot(wpr, add = FALSE, color = my_pal[i])
      #}else{
      #  plot(wpr, add = FALSE, color = my_pal[i])
      #}
      #no_skill = length(y.test[y.test==1]) / length(y.test)
      #lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
      #text(0.5, 0.07, round(no_skill,4))
    #}
      precision_new[[i]] = wpr$curve[,1]
      recall_new[[i]] = wpr$curve[,2]
    #balanced accuracy at cutoff 0.5
    CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
    
    results[1,i] = auc
    results[2,i] = CM$byClass[11]
    results[3,i] = au_prc
    results[4,i] = CM$byClass[1]
    results[5,i] = CM$byClass[2]
  }
  #tmp = tmp[2:nrow(tmp),]
  
  rownames(results) = c("AUC", "balanced_accuracy","auprc","spec", "sens(DCIS)")
  
  #### calculate overall metrics
  mean_values = matrix(0,nrow(results),1)
  sd_values = matrix(0,nrow(results),1)
  
  mean_values = apply(results, 1, mean)
  sd_values = apply(results, 1, sd)
  results = cbind(results, mean_values, sd_values)
  results = results[,11:12]
  
  #### below for recording predicted probability of each fold
  
  results_new = list()
  
  if(purpose == "draw_mean_curves"){
    print(results)
    results_new = list()
    results_new[[1]] = sens_new
    results_new[[2]] = spec_new
    results_new[[3]] = precision_new
    results_new[[4]] = recall_new
  }else if(purpose == "normal"){
    results_new[[1]] = results
    results_new[[2]] = GeneList
  }else{
    print("purpose should be normal or draw_mean_curves")
  }
  
  return(results_new)
  #return(results)
  
}

nested_cv_lung <- function(expr_data, type, purpose) {
  
  test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)
  
  # initialization of the nested cross-validation
  results = matrix(0,5,10)
  mse_results = matrix(0,1,10)
  GeneList = list()
  i = 1
  par(mfrow=c(3,4))
  #par(mfrow=c(1,1))
  tmp = matrix(as.factor(0), 1, 3)
  #my_pal <- seecol(pal_unikn_light)
  
  sens_new = list()
  spec_new = list()
  precision_new = list()
  recall_new = list()
  
  #### nested cross-validation ####
  for (i in 1:10) {
    rando = sample.int(100,1)
    set.seed(rando)
    test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                       LumA[as.numeric(test_index_LumA[[i]]), ],
                       Her2[as.numeric(test_index_Her2[[i]]), ],
                       Basal[as.numeric(test_index_Basal[[i]]), ],
                       Normal[as.numeric(test_index_Normal[[i]]), ]
    )
    train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                        LumA[as.numeric(-test_index_LumA[[i]]), ],
                        Her2[as.numeric(-test_index_Her2[[i]]), ],
                        Basal[as.numeric(-test_index_Basal[[i]]), ],
                        Normal[as.numeric(-test_index_Normal[[i]]), ]
    )   
    #test_data <- expr_data[as.numeric(test_index[[i]]), ]
    #train_data <- expr_data[-as.numeric(test_index[[i]]), ]
    x.test = test_data[, 3:ncol(expr_data)]
    y.test = test_data[, 1]
    x.train = train_data[, 3:ncol(expr_data)]
    y.train = train_data[, 1]
    
    y.train = as.character(y.train)
    y.train[y.train == "IBC"] = 0
    y.train[y.train == "IDC"] = 0
    y.train[y.train == "DCIS"] = 1
    y.train = as.numeric(y.train)
    y.train = as.factor(y.train)
    y.test = as.character(y.test)
    y.test[y.test == "IBC"] = 0
    y.test[y.test == "IDC"] = 0
    y.test[y.test == "DCIS"] = 1
    y.test = as.numeric(y.test)
    
    
    
    
    if(type == "normal"){
      pr_label = vector(mode = "logical", length = length(y.test))
      for (pr_i in 1:length(y.test)) {
        if(y.test[pr_i] == 1){
          pr_label[pr_i] = TRUE
        }
      }
      y.test = as.factor(y.test)
      cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), 
                             alpha = 1, family = "binomial", type.measure="auc")
      #plot(cv_output)
      model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                    alpha = 1, family = "binomial", type.measure="auc")
      
      #get the common genes
      betaGene = model$beta[,model$lambda==cv_output$lambda.min]
      betaGene = betaGene[betaGene != 0]
      betaGene = as.data.frame(betaGene)
      GeneList[[i]] = rownames(betaGene)
      
      result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
      label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
      result = as.numeric(result)
    }else if(type == "one_feature"){
      pr_label = vector(mode = "logical", length = length(y.test))
      for (pr_i in 1:length(y.test)) {
        if(y.test[pr_i] == 1){
          pr_label[pr_i] = TRUE
        }
      }
      data =cbind(y.train, x.train)
      colnames(data) = c("label", "ROR")
      data =as.data.frame(data)
      
      model = stats::glm(label~.,family=gaussian, data=data)
      
      test_data =cbind(y.test, x.test)
      test_data = as.data.frame(test_data)
      colnames(test_data) = c("label", "ROR")
      result = predict(model, test_data, type= "response")
      label_result = result
      label_result[label_result>=0.5] = 1
      label_result[label_result<0.5] = 0
    }else{
      print("the type should be 'normal' or 'one_feature'")
      return()
    }
    
    #calculate average mse for this fold
    mse = cbind(as.data.frame(result), as.data.frame(as.numeric(as.character(y.test))))
    temp_c = as.data.frame(matrix(0,nrow(mse),1))
    mse = cbind(mse, temp_c)
    for (j in 1:nrow(mse)) {
      mse[j,3] = (mse[j,1] - mse[j,2])^2
    }
    mse_results[1,i] = mean(mse[j,3])
    
    roc = roc(response = y.test, predictor = result, quiet = TRUE)
    
    sens_new[[i]] = roc$sensitivities
    spec_new[[i]] = roc$specificities
    
    
    #, print.auc = TRUE, plot=TRUE, legacy.axes = TRUE)
    #plot(roc)
    auc = auc(roc)
    
    #cutoff table with ROCR package
    #pred <- prediction(result, y.test)
    #perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    #plot(perf, col=rainbow(10))
    #str(perf)
    #cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
    #                      tpr=perf@y.values[[1]])
    
    #PR curve
    #if(type == "normal"){
    wpr<-pr.curve(result, weights.class0 = pr_label, curve = TRUE,
                  max.compute = T, min.compute = T, rand.compute = T)
    au_prc = wpr$auc.integral
    #if (i == 1){
    #  plot(wpr, add = FALSE, color = my_pal[i])
    #}else{
    #  plot(wpr, add = FALSE, color = my_pal[i])
    #}
    #no_skill = length(y.test[y.test==1]) / length(y.test)
    #lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
    #text(0.5, 0.07, round(no_skill,4))
    #}
    precision_new[[i]] = wpr$curve[,1]
    recall_new[[i]] = wpr$curve[,2]
    #balanced accuracy at cutoff 0.5
    CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
    
    results[1,i] = auc
    results[2,i] = CM$byClass[11]
    results[3,i] = au_prc
    results[4,i] = CM$byClass[1]
    results[5,i] = CM$byClass[2]
  }
  #tmp = tmp[2:nrow(tmp),]
  
  rownames(results) = c("AUC", "balanced_accuracy","auprc","spec", "sens(DCIS)")
  
  #### calculate overall metrics
  mean_values = matrix(0,nrow(results),1)
  sd_values = matrix(0,nrow(results),1)
  
  mean_values = apply(results, 1, mean)
  sd_values = apply(results, 1, sd)
  results = cbind(results, mean_values, sd_values)
  results = results[,11:12]
  
  #### below for recording predicted probability of each fold
  
  results_new = list()
  
  if(purpose == "draw_mean_curves"){
    print(results)
    results_new = list()
    results_new[[1]] = sens_new
    results_new[[2]] = spec_new
    results_new[[3]] = precision_new
    results_new[[4]] = recall_new
  }else if(purpose == "normal"){
    results_new[[1]] = results
    results_new[[2]] = GeneList
  }else{
    print("purpose should be normal or draw_mean_curves")
  }
  
  return(results_new)
  #return(results)
  
}


nested_cv_plot <- function(expr_data, type) {
  #probe_list = GeneList
  #expr_data = mlbased
  #expr_data = traditional
  
  LumB = expr_data[expr_data$subtype == "LumB",]
  LumA = expr_data[expr_data$subtype == "LumA",]
  Her2 = expr_data[expr_data$subtype == "Her2",]
  Basal = expr_data[expr_data$subtype == "Basal",]
  Normal = expr_data[expr_data$subtype == "Normal",]
  
  test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
  test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)
  
  # initialization of the nested cross-validation
  results = matrix(0,5,10)
  mse_results = matrix(0,1,10)
  GeneList = list()
  i = 1
  par(mfrow=c(3,4))
  #par(mfrow=c(1,1))
  tmp = matrix(as.factor(0), 1, 3)
  #my_pal <- seecol(pal_unikn_light)
  
  sens_new = list()
  spec_new = list()
  
  #### nested cross-validation ####
  for (i in 1:10) {
    rando = sample.int(100,1)
    set.seed(rando)
    test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                       LumA[as.numeric(test_index_LumA[[i]]), ],
                       Her2[as.numeric(test_index_Her2[[i]]), ],
                       Basal[as.numeric(test_index_Basal[[i]]), ],
                       Normal[as.numeric(test_index_Normal[[i]]), ]
    )
    train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                        LumA[as.numeric(-test_index_LumA[[i]]), ],
                        Her2[as.numeric(-test_index_Her2[[i]]), ],
                        Basal[as.numeric(-test_index_Basal[[i]]), ],
                        Normal[as.numeric(-test_index_Normal[[i]]), ]
    )   
    #test_data <- expr_data[as.numeric(test_index[[i]]), ]
    #train_data <- expr_data[-as.numeric(test_index[[i]]), ]
    x.test = test_data[, 3:ncol(expr_data)]
    y.test = test_data[, 1]
    x.train = train_data[, 3:ncol(expr_data)]
    y.train = train_data[, 1]
    
    y.train = as.character(y.train)
    y.train[y.train == "IBC"] = 0
    y.train[y.train == "IDC"] = 0
    y.train[y.train == "DCIS"] = 1
    y.train = as.numeric(y.train)
    y.train = as.factor(y.train)
    y.test = as.character(y.test)
    y.test[y.test == "IBC"] = 0
    y.test[y.test == "IDC"] = 0
    y.test[y.test == "DCIS"] = 1
    y.test = as.numeric(y.test)
    
    
      pr_label = vector(mode = "logical", length = length(y.test))
      for (pr_i in 1:length(y.test)) {
        if(y.test[pr_i] == 1){
          pr_label[pr_i] = TRUE
        }
      }
      y.test = as.factor(y.test)
      cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), 
                             alpha = 1, family = "binomial", type.measure="auc")
      #plot(cv_output)
      model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                    alpha = 1, family = "binomial", type.measure="auc")
      
      #get the common genes
      betaGene = model$beta[,model$lambda==cv_output$lambda.min]
      betaGene = betaGene[betaGene != 0]
      betaGene = as.data.frame(betaGene)
      GeneList[[i]] = rownames(betaGene)
      
      result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
      label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
      result = as.numeric(result)
    
    #calculate average mse for this fold
    mse = cbind(as.data.frame(result), as.data.frame(as.numeric(as.character(y.test))))
    temp_c = as.data.frame(matrix(0,nrow(mse),1))
    mse = cbind(mse, temp_c)
    for (j in 1:nrow(mse)) {
      mse[j,3] = (mse[j,1] - mse[j,2])^2
    }
    mse_results[1,i] = mean(mse[j,3])
    
    roc = roc(response = y.test, predictor = result, quiet = TRUE)
    auc = auc(roc)
    
    sens_new[[i]] = roc$sensitivities
    spec_new[[i]] = roc$specificities
    
    #, print.auc = TRUE, plot=TRUE, legacy.axes = TRUE)
    #plot(roc, yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
    
    #auc = auc(roc)
    
    #cutoff table with ROCR package
    #pred <- prediction(result, y.test)
    #perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    #plot(perf, col=rainbow(10))
    #str(perf)
    #cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
    #                      tpr=perf@y.values[[1]])
    
    #PR curve
      wpr<-pr.curve(result, weights.class0 = pr_label, curve = TRUE,
                    max.compute = T, min.compute = T, rand.compute = T)
      au_prc = wpr$auc.integral
      
      if(type == "auroc"){
        plot(roc, yaxt = "n",cex.axis=1.5, cex.lab =1.5, main = "")
        
        axis(side = 2, labels = FALSE)
        
        ## Draw the y-axis.
        axis(side = 2,
             ## Rotate the labels.
             las = 2,
             ## Adjust the label position.
             mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab =1.5)
      }else if(type == "auprc"){
        if (i == 1){
          wpr$auc.integral = ""
          plot(wpr, add = FALSE, color = "black", yaxt = "n",cex.axis=1.5, cex.lab =1.5, main = "")
        }else{
          plot(wpr, add = FALSE, color = "black", yaxt = "n",cex.axis=1.5, cex.lab =1.5, main = "")
        }
        title(main = "")
        no_skill = length(y.test[y.test==1]) / length(y.test)
        lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
        text(0.5, 0.07, round(no_skill,4), cex =1.5)
        axis(side = 2, labels = FALSE)
        
        ## Draw the y-axis.
        axis(side = 2,
             ## Rotate the labels.
             las = 2,
             ## Adjust the label position.
             mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab =1.5)
      }else{
        print("auroc or auprc?")
      }
    #balanced accuracy at cutoff 0.5
    CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
    
    results[1,i] = auc
    results[2,i] = CM$byClass[11]
    results[3,i] = au_prc
    results[4,i] = CM$byClass[1]
    results[5,i] = CM$byClass[2]
  }
  #tmp = tmp[2:nrow(tmp),]
  
  rownames(results) = c("AUC", "balanced_accuracy","auprc","spec", "sens(DCIS)")
  
  #### calculate overall metrics
  mean_values = matrix(0,nrow(results),1)
  sd_values = matrix(0,nrow(results),1)
  
  mean_values = apply(results, 1, mean)
  sd_values = apply(results, 1, sd)
  results = cbind(results, mean_values, sd_values)
  results = results[,11:12]
  
  #### below for recording predicted probability of each fold
  
  results_new = list()
  results_new[[1]] = results
  results_new[[2]] = GeneList
  
  #results[[1]] = sens_new
  #results[[2]] = spec_new
  return(results)
  #return(results_new)
  
}

roc_and_prc <- function(iter, y.test, result, sens_new, spec_new, auc_new, precision_new, recall_new, prauc_new, pr_label) {
  new_list = list()
  out <- try(roc(y.test, result, quiet = TRUE), silent = TRUE)
  #temp = class(out) == "try-error"

  if(class(out) == "try-error"){
  }else{
    k = roc(y.test, result, quiet = TRUE)
    sens_new[[iter]] = k$sensitivities
    spec_new[[iter]] = k$specificities
    
    #title("result of random selected 10 features in training set")
    auc_new[[iter]] = auc(k)
    
    wpr = pr.curve(result, weights.class0 = pr_label, curve = TRUE,
                   max.compute = T, min.compute = T, rand.compute = T)
    plot(wpr)
    prauc_new[[iter]]  = wpr$auc.integral
    
    precision_new[[iter]] = wpr$curve[,2]
    recall_new[[iter]] = wpr$curve[,1]
    
    new_list[[1]] = sens_new
    new_list[[2]] = spec_new
    new_list[[3]] = auc_new
    new_list[[4]] = precision_new
    new_list[[5]] = recall_new
    new_list[[6]] = prauc_new
    
    return(new_list)
  }
}