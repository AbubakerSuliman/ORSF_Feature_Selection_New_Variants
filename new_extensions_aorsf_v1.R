#New Extensions of LCIV for Accelerated Oblique Random Survival Forests#
#Abubaker Suliman 06 March 2025, UAEU#
library(survival)
library(mRMRe)
library(glmnet)
library(carSurv)

f_cvglment_se <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  
  #data <- as.data.frame(cbind(y_node, x_node))
  #print(ncol(x_node))
  if(nrow(x_node) <= 25) #| ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  #future::plan("multisession", workers = n_wrk)
  #doMC::registerDoMC(cores = 3)
  suppressWarnings(
    fit <- try(
      glmnet::cv.glmnet(as.matrix(x_node), as.matrix(y_node), family = "cox", nfolds = 3, parallel = TRUE),
      silent = TRUE
    )
  )
  #foreach::registerDoSEQ()
  
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), nrow=ncol(x_node), ncol=1))
  }
  
  #return(matrix(fit$glmnet.fit$beta[, which(fit$lambda==fit$lambda.min)], ncol = 1))
  return(matrix(fit$glmnet.fit$beta[, which(fit$lambda==fit$lambda.1se)], ncol = 1))
  
}

f_cvglment_min <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  
  #data <- as.data.frame(cbind(y_node, x_node))
  #print(ncol(x_node))
  if(nrow(x_node) <= 25) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  #future::plan("multisession", workers = n_wrk)
  suppressWarnings(
    fit <- try(
      glmnet::cv.glmnet(as.matrix(x_node), as.matrix(y_node), family = "cox", nfolds = 3),
      silent = TRUE
    )
  )
  
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), nrow=ncol(x_node), ncol=1))
  }
  
  return(matrix(fit$glmnet.fit$beta[, which(fit$lambda==fit$lambda.min)], ncol = 1))
  #return(matrix(fit$glmnet.fit$beta[, which(fit$lambda==fit$lambda.1se)], ncol = 1))
  
}

f_mRMR_clsc_2q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        mRMR_data <- mRMRe::mRMR.data(data=as.data.frame(df))

        mRMR_classic <- mRMRe::mRMR.classic(data=mRMR_data, target_indices=ncol(df), feature_count=ncol(df)-1)

        vec = setNames(unlist(mRMRe::scores(mRMR_classic)), mRMR_classic@feature_names[unlist(mRMRe::solutions(mRMR_classic))])
        vec = vec[is.finite(vec)]
        
        #print(vec)
        #if(any(str_detect(names(vec),"target"))) 
        #vec = vec[sort(names(vec))] #Problem x1 x10 x11 x12 x2 x3 etc
        sel_var = names(vec[vec > quantile(vec, 0.5)]) #0.75 / 0.25
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(sel_var, collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[sel_var] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  if(inherits(fit, "try-error")){
    #print("Error")
    
    return(matrix(runif(ncol(x_node)), nrow=ncol(x_node), ncol=1))
  }
  
  return(matrix(vec, ncol = 1))
}

f_mRMR_clsc_1q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        mRMR_data <- mRMRe::mRMR.data(data=as.data.frame(df))
        
        mRMR_classic <- mRMRe::mRMR.classic(data=mRMR_data, target_indices=ncol(df), feature_count=ncol(df)-1)
        
        vec = setNames(unlist(mRMRe::scores(mRMR_classic)), mRMR_classic@feature_names[unlist(mRMRe::solutions(mRMR_classic))])
        vec = vec[is.finite(vec)]
        
        #print(vec)
        #if(any(str_detect(names(vec),"target"))) 
        #vec = vec[sort(names(vec))] #Problem x1 x10 x11 x12 x2 x3 etc
        sel_var = names(vec[vec > quantile(vec, 0.25)]) #0.75 / 0.5
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(sel_var, collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[sel_var] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  if(inherits(fit, "try-error")){
    #print("Error")
    
    return(matrix(runif(ncol(x_node)), nrow=ncol(x_node), ncol=1))
  }
  
  return(matrix(vec, ncol = 1))
}

f_mRMR_clsc_3q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        mRMR_data <- mRMRe::mRMR.data(data=as.data.frame(df))
        
        mRMR_classic <- mRMRe::mRMR.classic(data=mRMR_data, target_indices=ncol(df), feature_count=ncol(df)-1)
        
        vec = setNames(unlist(mRMRe::scores(mRMR_classic)), mRMR_classic@feature_names[unlist(mRMRe::solutions(mRMR_classic))])
        vec = vec[is.finite(vec)]
        
        #print(vec)
        #if(any(str_detect(names(vec),"target"))) 
        #vec = vec[sort(names(vec))] #Problem x1 x10 x11 x12 x2 x3 etc
        sel_var = names(vec[vec > quantile(vec, 0.75)]) #0.75 / 0.25
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(sel_var, collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[sel_var] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  if(inherits(fit, "try-error")){
    #print("Error")
    
    return(matrix(runif(ncol(x_node)), nrow=ncol(x_node), ncol=1))
  }
  
  return(matrix(vec, ncol = 1))
}

f_csurv_fndr <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        
        # Estimate CAR scores
        carScores <- carSurv::carSurvScore(obsTime=y_node[ ,'time'], obsEvent=y_node[ ,'status'], X=as.matrix(x_node))
        selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "fndr", plotDiag = FALSE, threshold = 0.05)
        #selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "threshold", plotDiag = FALSE, threshold = 0.5)
        #vec = vec[is.finite(vec)]
        
        if(length(selectedVar)==0) return(matrix(runif(ncol(x_node)), ncol = 1))
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(colnames(x_node)[selectedVar], collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[colnames(x_node)[selectedVar]] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  #browser()
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), ncol = 1))
  }
  return(matrix(vec, ncol = 1))
}

f_csurv_2q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))

  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        
        # Estimate CAR scores
        carScores <- carSurv::carSurvScore(obsTime=y_node[ ,'time'], obsEvent=y_node[ ,'status'], X=as.matrix(x_node))
        #selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "fndr", plotDiag = FALSE, threshold = 0.05)
        selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "threshold", plotDiag = FALSE, threshold = 0.5)
        #vec = vec[is.finite(vec)]
        
        if(length(selectedVar)==0) return(matrix(runif(ncol(x_node)), ncol = 1))
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(colnames(x_node)[selectedVar], collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[colnames(x_node)[selectedVar]] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  #browser()
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), ncol = 1))
  }
  return(matrix(vec, ncol = 1))
}

f_csurv_1q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        
        # Estimate CAR scores
        carScores <- carSurv::carSurvScore(obsTime=y_node[ ,'time'], obsEvent=y_node[ ,'status'], X=as.matrix(x_node))
        #selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "fndr", plotDiag = FALSE, threshold = 0.05)
        selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "threshold", plotDiag = FALSE, threshold = 0.25)
        #vec = vec[is.finite(vec)]
        
        if(length(selectedVar)==0) return(matrix(runif(ncol(x_node)), ncol = 1))
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(colnames(x_node)[selectedVar], collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[colnames(x_node)[selectedVar]] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  #browser()
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), ncol = 1))
  }
  return(matrix(vec, ncol = 1))
}

f_csurv_3q <- function(x_node, y_node, w_node) {
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  if(nrow(x_node) <= 15) # |ncol(x_node)<5
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      {
        df = data.frame(x_node)
        df$target <- survival::Surv(y_node[ ,'time'], y_node[ ,'status'])
        
        # Estimate CAR scores
        carScores <- carSurv::carSurvScore(obsTime=y_node[ ,'time'], obsEvent=y_node[ ,'status'], X=as.matrix(x_node))
        #selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "fndr", plotDiag = FALSE, threshold = 0.05)
        selectedVar <- carSurv::carVarSelect(carSurvScores=carScores, method = "threshold", plotDiag = FALSE, threshold = 0.75)
        #vec = vec[is.finite(vec)]
        
        if(length(selectedVar)==0) return(matrix(runif(ncol(x_node)), ncol = 1))
        
        formula_r = as.formula(sprintf("target ~ %s", paste0(colnames(x_node)[selectedVar], collapse = " + ")))
        adj_model <- eval(bquote(coxph(.(formula_r), data = df)))
        
        vec = setNames(object = vector(mode = "numeric", length = length(colnames(x_node))), nm = colnames(x_node))
        vec[colnames(x_node)[selectedVar]] = coef(adj_model)
        if(any(is.infinite(vec) | is.na(vec))){
          #browser()
          vec[is.na(vec) | is.infinite(vec)] = 0
        }
      },
      silent = TRUE
    )
  )
  #browser()
  if(inherits(fit, "try-error")){
    #print("Error")
    return(matrix(runif(ncol(x_node)), ncol = 1))
  }
  return(matrix(vec, ncol = 1))
}


penalized_fitter <- function(x_node, y_node, w_node, alpha, df_target=3){
  
  #browser()
  colnames(y_node) <- c('time', 'status')
  colnames(x_node) <- paste("x", seq(ncol(x_node)), sep = '')
  
  data <- as.data.frame(cbind(y_node, x_node))
  
  if(nrow(data) <= 10)
    return(matrix(runif(ncol(x_node)), ncol = 1))
  
  suppressWarnings(
    fit <- try(
      glmnet::glmnet(x = x_node,
                     y = survival::Surv(data$time, data$status),
                     #weights = w_node,
                     alpha = 0.5,
                     family = "cox"),
      silent = TRUE
    )
  )
  
  # suppressWarnings(
  #   fit <- try(
  #     glmnet::glmnet(x = x_node,
  #                    y = y_node,
  #                    weights = w_node,
  #                    alpha = alpha,
  #                    family = family),
  #     silent = TRUE
  #   )
  # )
  
  if(aorsf:::is_error(fit)){
    return(matrix(0, nrow=ncol(x_node), ncol=1))
  }
  
  for(i in seq_along(fit$df)){
    if(fit$df[i] >= df_target || i == length(fit$df)){
      return(matrix(fit$beta[, i, drop=TRUE], ncol = 1))
    }
  }
  
}

#f_mRMR_clsc_2q(x_node = smarto1 %>% select(-time, -event) %>% select_if(is.numeric) %>% as.data.frame(), y_node = smarto1 %>% select(time, event) %>% as.data.frame())


