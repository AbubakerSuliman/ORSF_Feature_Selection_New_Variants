define_learners <- function(seed_r, tune_folds=3, n_evals=1, mtry_lwr = 3, mtry_upr = 7){
  
  source("LearnerSurvAorsfNew.R")
  source("new_extensions_aorsf_v1.R")
  
  #browser()
  set.seed(seed_r)
  rsmp_tune_cv = rsmp("holdout")
  trm_tune = trm("evals", n_evals = n_evals)
  
  fallback = lrn("surv.kaplan")
  
  po_cat = po("imputelearner", lrn("classif.rpart"), param_vals =list(affect_columns = selector_grep("fac_")), id="class") %>>% 
    po("encode", method = "treatment", affect_columns = selector_type("factor"), id = "binary_enc") %>>%
    po("removeconstants")
  
  po_con = po("imputelearner", lrn("regr.rpart"), param_vals =list(affect_columns = selector_grep("num_")), id="reg") %>>% 
    po("scale", param_vals =list(affect_columns = selector_grep("num_")))
  
  preproc = po_cat %>>% 
    po_con 

  learners_lst = list()
  
  ###################################
  ### A. Baselearners ###############
  ###################################
  
  learners_lst$kaplan = as_learner(preproc %>>% lrn("surv.kaplan", id = "class_nonpar_kaplan"))
  learners_lst$kaplan$encapsulate("evaluate", fallback)

  learners_lst$cox = as_learner(preproc %>>% lrn("surv.coxph", id = "class_semipar_coxph"))
  learners_lst$cox$encapsulate("evaluate", fallback)
 
  learners_lst$aft = as_learner(preproc %>>% lrn("surv.parametric", id = "class_par_aft"))
  learners_lst$aft$encapsulate("evaluate", fallback)
  
  learners_lst$cv_glmnet = auto_tuner(
    learner = preproc %>>% lrn("surv.cv_glmnet",
                  id = "class_semipar_cvglmnet",
                  alpha = to_tune(p_dbl(0, 1))),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
  ,
    store_models = FALSE
  )

  learners_lst$cv_glmnet$encapsulate("evaluate", fallback)

  
  learners_lst$penalized = auto_tuner(
    learner = preproc %>>% lrn("surv.penalized",
                                    id = "class_semipar_penalized",
                                    lambda1 = to_tune(p_int(0, 10)),
                                    lambda2 = to_tune(p_int(0, 10))),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  learners_lst$penalized$encapsulate("evaluate", fallback)
  
  ###################################
  ### B. RandomForestlearners #######
  ###################################
  ###############################
  #Random Survival Forests #####
  ##############################
  
  learners_lst$at_lrn_rfsrc = auto_tuner(
    learner = preproc %>>% lrn("surv.rfsrc",
                  id = "ml_ranfor_rfsrc_brier",
                  splitrule = "bs.gradient",
                  ntree = to_tune(p_int(250, 750)),
                  mtry = to_tune(p_int(2, 7)),
                  nodesize = to_tune(p_int(15, 25))),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search"),
    store_models = FALSE
  )


  learners_lst$at_lrn_rfsrc$encapsulate("evaluate", fallback)
  
  
  ###################################################
  #Accelerated Oblique Random Survival Forests  #####
  ###################################################

  learners_lst$at_lrn_aorsff = auto_tuner(
    learner = preproc %>>% lrn("surv.aorsf",
                  id = "ml_ranfor_aorsf_fast",
                  control_type = "fast", #cph, net
                  split_rule = 'cstat', #'logrank'
                  n_thread = 1,
                  control_net_df_target = 3,
                  n_tree = to_tune(p_int(lower = 250, upper = 750)),
                  leaf_min_obs = to_tune(p_int(lower = 15, upper = 25)),
                  mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
    ),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search"),
    store_models = FALSE
  )

  learners_lst$at_lrn_aorsff$encapsulate("evaluate", fallback)
  
  learners_lst$at_lrn_aorsfc = auto_tuner(
    learner = preproc %>>% lrn("surv.aorsf",
                  id = "ml_ranfor_aorsf_cph",
                  control_type = "cph", #cph, net
                  split_rule = 'cstat', #'logrank'
                  n_thread = 1,
                  control_net_df_target = 3,
                  n_tree = to_tune(p_int(lower = 250, upper = 750)),
                  leaf_min_obs = to_tune(p_int(lower = 15, upper = 25)),
                  mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
    ),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search"),
    store_models = FALSE
  )
  learners_lst$at_lrn_aorsfc$encapsulate("evaluate", fallback)

  
  learners_lst$at_lrn_aorsfn = auto_tuner(
    learner = preproc %>>% lrn("surv.aorsf",
                  id = "ml_ranfor_aorsf_net",
                  control_type = "net", #cph, net
                  split_rule = 'cstat', #'logrank'
                  n_thread = 1,
                  control_net_df_target = 3,
                  n_tree = to_tune(p_int(lower = 250, upper = 750)),
                  leaf_min_obs = to_tune(p_int(lower = 15, upper = 25)),
                  mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
    ),
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  learners_lst$at_lrn_aorsfn$encapsulate("evaluate", fallback)
  
  
  learner_min = LearnerSurvAorsfNew$new()
  learner_min$id = "surv.aorsfnew_min"
  learner_min$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_min$param_set$values$control_type = "custom"
  learner_min$param_set$values$control_custom_fun = f_cvglment_min
  learner_min$param_set$values$split_rule = 'cstat'
  learner_min$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_min$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  learners_lst$at_lrn_aorsf_min = auto_tuner(
    learner = preproc %>>% learner_min,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  learners_lst$at_lrn_aorsf_min$encapsulate("evaluate", fallback)
  
  
  learner_se = LearnerSurvAorsfNew$new()
  learner_se$id = "surv.aorsfnew_se"
  learner_se$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_se$param_set$values$control_type = "custom"
  learner_se$param_set$values$control_custom_fun = f_cvglment_se
  learner_se$param_set$values$split_rule = 'cstat'
  learner_se$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_se$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_se = auto_tuner(
    learner = preproc %>>% learner_se,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_se$encapsulate("evaluate", fallback)
  
  learner_mr_2q = LearnerSurvAorsfNew$new()
  learner_mr_2q$id = "surv.aorsfnew_mr_2q"
  learner_mr_2q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_mr_2q$param_set$values$control_type = "custom"
  learner_mr_2q$param_set$values$control_custom_fun = f_mRMR_clsc_2q
  learner_mr_2q$param_set$values$split_rule = 'cstat'
  learner_mr_2q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_mr_2q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_mr_2q = auto_tuner(
    learner = preproc %>>% learner_mr_2q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_mr_2q$encapsulate("evaluate", fallback)
  
  learner_mr_1q = LearnerSurvAorsfNew$new()
  learner_mr_1q$id = "surv.aorsfnew_mr_1q"
  learner_mr_1q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_mr_1q$param_set$values$control_type = "custom"
  learner_mr_1q$param_set$values$control_custom_fun = f_mRMR_clsc_1q
  learner_mr_1q$param_set$values$split_rule = 'cstat'
  learner_mr_1q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_mr_1q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_mr_1q = auto_tuner(
    learner = preproc %>>% learner_mr_1q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_mr_1q$encapsulate("evaluate", fallback)
  
  
  learner_mr_3q = LearnerSurvAorsfNew$new()
  learner_mr_3q$id = "surv.aorsfnew_mr_3q"
  learner_mr_3q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_mr_3q$param_set$values$control_type = "custom"
  learner_mr_3q$param_set$values$control_custom_fun = f_mRMR_clsc_3q
  learner_mr_3q$param_set$values$split_rule = 'cstat'
  learner_mr_3q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_mr_3q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_mr_3q = auto_tuner(
    learner = preproc %>>% learner_mr_3q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_mr_3q$encapsulate("evaluate", fallback)
  
  #########################################################################################################
  #########################################################################################################
  #########################################################################################################
  
  
  learner_cv_fndr = LearnerSurvAorsfNew$new()
  learner_cv_fndr$id = "surv.aorsfnew_cv_fndr"
  learner_cv_fndr$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_cv_fndr$param_set$values$control_type = "custom"
  learner_cv_fndr$param_set$values$control_custom_fun = f_csurv_fndr
  learner_cv_fndr$param_set$values$split_rule = 'cstat'
  learner_cv_fndr$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_cv_fndr$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_cv_fndr = auto_tuner(
    learner = preproc %>>% learner_cv_fndr,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_cv_fndr$encapsulate("evaluate", fallback)
  
  
  learner_cv_1q = LearnerSurvAorsfNew$new()
  learner_cv_1q$id = "surv.aorsfnew_cv_1q"
  learner_cv_1q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_cv_1q$param_set$values$control_type = "custom"
  learner_cv_1q$param_set$values$control_custom_fun = f_csurv_1q
  learner_cv_1q$param_set$values$split_rule = 'cstat'
  learner_cv_1q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_cv_1q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_cv_1q = auto_tuner(
    learner = preproc %>>% learner_cv_1q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_cv_1q$encapsulate("evaluate", fallback)
  
  
  learner_cv_2q = LearnerSurvAorsfNew$new()
  learner_cv_2q$id = "surv.aorsfnew_cv_2q"
  learner_cv_2q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_cv_2q$param_set$values$control_type = "custom"
  learner_cv_2q$param_set$values$control_custom_fun = f_csurv_2q
  learner_cv_2q$param_set$values$split_rule = 'cstat'
  learner_cv_2q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_cv_2q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_cv_2q = auto_tuner(
    learner = preproc %>>% learner_cv_2q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_cv_2q$encapsulate("evaluate", fallback)
  
  
  learner_cv_3q = LearnerSurvAorsfNew$new()
  learner_cv_3q$id = "surv.aorsfnew_cv_3q"
  learner_cv_3q$param_set$values$mtry = to_tune(p_int(lower = mtry_lwr, upper = mtry_upr))
  learner_cv_3q$param_set$values$control_type = "custom"
  learner_cv_3q$param_set$values$control_custom_fun = f_csurv_3q
  learner_cv_3q$param_set$values$split_rule = 'cstat'
  learner_cv_3q$param_set$values$n_tree = to_tune(p_int(lower = 250, upper = 750))
  learner_cv_3q$param_set$values$leaf_min_obs = to_tune(p_int(lower = 15, upper = 25))
  
  learners_lst$at_lrn_aorsf_cv_3q = auto_tuner(
    learner = preproc %>>% learner_cv_3q,
    resampling = rsmp_tune_cv,
    measure = msr("surv.cindex"),
    terminator = trm_tune,
    tuner = tnr("random_search")
    ,
    store_models = FALSE
  )
  
  learners_lst$at_lrn_aorsf_cv_3q$encapsulate("evaluate", fallback)


  return(learners_lst)
  
}
