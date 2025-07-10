########################################
#A. Load Libraries       ###############
########################################
library(readr)
library(tidyverse)
library(mlr3pipelines)
# install learners
library(mlr3learners)
library(mlr3proba)
library(mlr3extralearners)
library(paradox) #ParamSet
library(mlr3tuning) #AutoTuner
library(mlr3pipelines)
library(reticulate)
library(survex)
library(prodlim)
library(epiR)
# control Pyhon warnings through reticulate:
warnings <- import("warnings")
# ignore / suppress:
warnings$simplefilter("ignore")

#Load Functions
source("define_learners_novel_v1.R")

#PATH
main_dir = getwd()
date_r = "03May2025"

cc <- commandArgs(trailingOnly  = TRUE)

sim_scenario = as.character(cc[1]) #"overall"
sim_id = as.character(cc[2]) #123
n_workers = as.integer(cc[3]) #10
n_evals = as.integer(cc[4]) #1

# sim_scenario = "grace"
# sim_id = 123
# n_workers = 20
# n_evals = 30
########################################
#B. Read Dataset         ###############
########################################
smarto <- read_csv("1_data/1_original/smarto.csv", na = "missing")

grace <- read_csv("1_data/1_original/grace.csv")

cost <- read_csv("1_data/1_original/cost.csv")

#The TRACE study group.
# Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L., F. E. Nielsen, Melchior, T., Joen, T. and
# P. K. Andersen (1997), Does in-hospital ventricular fibrillation affect prognosis after myocardial
# infarction?, European Heart Journal 18, 919924
#TRACE The TRACE study group of myocardial infarction
#status a numeric vector code. Survival status. 9: dead from myocardial infarction, 0: alive, 7,8: dead from other causes.

trace <- read_csv("1_data/1_original/TRACE.csv")

whas500 <- read_csv("1_data/1_original/whas500.csv")

leukemia <- read_csv("1_data/1_original/LeukSurv.csv")

prostate <- read_csv("1_data/1_original/prostate.csv", na = "missing")

colon <- read_csv("1_data/1_original/colon.csv", na = "missing")

breast <- read_csv("1_data/1_original/GBSG2.csv")

#survival::rotterdam
rotterdam <- read_csv("1_data/1_original/rotterdam.csv")

#####################################
#C. Data Preprocessiong #############
#####################################
missing_cols = colnames(smarto)[colSums(is.na(smarto))/nrow(smarto)>0.15]
smarto1 = smarto %>% select(-all_of(c('pid', missing_cols))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

#Time in days and Censoring indicator. 1 = Death. 0 = Censored
grace1 = grace %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

cost1 = cost %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

#time a numeric vector. Survival time in years.
#Medain Flw (IQR) 7.0 years (6.5 - 7.5) 
trace1 = trace %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor)

whas5001 = whas500 %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

leukemia1 = leukemia %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(fac_district = fct_collapse(fac_district, 
                                      A = c("1", "2", "3", "4", "5", "6", "7", "8"),
                                      B = c("9", "10", "11", "12", "13", "14", "15", "16"),
                                      C = c("17", "18", "19", "20", "21", "22", "23", "24")),
         time = time/365.25)

prostate1 = prostate %>% select(-all_of(c('pid'))) %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  #Months of Follow-up 
  mutate(time = ifelse(time==0, 0.5, time)/12)

#Death table(survival::colon$etype, survival::colon$status)
colon1 = colon %>% select(-all_of(c('pid')))  %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

# time recurrence free survival time (in days); cens censoring indicator (0- censored, 1- event).
breast1 = breast %>% select(-all_of(c('pid')))  %>% mutate_if(str_detect(names(.), "fac"), factor) %>%
  mutate(time = time/365.25)

rotterdam1 = rotterdam %>% select(-all_of(c('pid', 'year', 'death', 'dtime'))) %>%
  rename(time=rtime, event=recur) %>%
  mutate_at(vars(meno, size, grade, hormon, chemo), factor) %>%
  rename_at(vars(meno, size, grade, hormon, chemo), ~paste0("fac_", .)) %>%
  rename_at(vars(age, nodes, pgr, er), ~paste0("num_", .)) %>%
  mutate(time = time/365.25)
  

#########################################################################
#D. Benchmarking (Statistical Models - ML Models )           ############
#########################################################################

########################
#D.1 Tasks     #########
########################

if(sim_scenario == "smarto"){
  tsk_train = as_task_surv(smarto1, time = "time", event = "event", type = "right", id = "smarto")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}


if(sim_scenario == "grace"){
  tsk_train = as_task_surv(grace1, time = "time", event = "event", type = "right", id = "grace")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 5
}


if(sim_scenario == "cost"){
  tsk_train = as_task_surv(cost1, time = "time", event = "event", type = "right", id = "cost")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}


if(sim_scenario == "trace"){
  tsk_train = as_task_surv(trace1, time = "time", event = "event", type = "right", id = "trace")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 6
}

if(sim_scenario == "whas"){
  tsk_train = as_task_surv(whas5001, time = "time", event = "event", type = "right", id = "whas")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

if(sim_scenario == "breast"){
  tsk_train = as_task_surv(breast1, time = "time", event = "event", type = "right", id = "breast")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

if(sim_scenario == "colon"){
  tsk_train = as_task_surv(colon1, time = "time", event = "event", type = "right", id = "colon")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

if(sim_scenario == "prostate"){
  tsk_train = as_task_surv(prostate1, time = "time", event = "event", type = "right", id = "prostate")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

if(sim_scenario == "leukemia"){
  tsk_train = as_task_surv(leukemia1, time = "time", event = "event", type = "right", id = "leukemia1")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

if(sim_scenario == "rotterdam"){
  tsk_train = as_task_surv(rotterdam1, time = "time", event = "event", type = "right", id = "rotterdam")
  tsk_train$set_col_roles("event", c("target", "stratum"))
  
  mtry_lwr_v = 3
  mtry_upr_v = 7
}

########################
#D.2 Learners     ######
########################

learners_lst = define_learners(seed_r=sim_id, n_evals = n_evals, mtry_lwr = mtry_lwr_v, mtry_upr = mtry_upr_v)
#learners_lst_red = learners_lst[c("at_lrn_aorsff", "at_lrn_aorsfc", "at_lrn_aorsfn")]
#learners_lst_red = learners_lst[c("at_lrn_aorsf_mr_1q", "at_lrn_aorsf_mr_3q", "at_lrn_aorsf_cv_1q", "at_lrn_aorsf_cv_2q", 
 #                                 "at_lrn_aorsf_cv_3q")]# "at_lrn_aorsf_mr_2q", "at_lrn_aorsf_mr_1q", "at_lrn_aorsf_mr_3q", "at_lrn_aorsf_cv_fndr", "at_lrn_aorsf_cv_1q", "at_lrn_aorsf_cv_2q", "at_lrn_aorsf_cv_3q"


########################
#D.3 Benchmark     #####
########################
set.seed(sim_id)
rsmp_outer_cv = rsmp("cv", folds = 3)

# set logger to maximum information - save to scratch
file_path = sprintf("%s/3_results/%s/%s", main_dir, date_r, sim_scenario) 
dir.create(file_path, showWarnings = FALSE, recursive = TRUE)
options("mlr3.debug" = TRUE)
logfile = sprintf("%s/log.txt", file_path)
logger = lgr::get_logger("mlr3")
logger$set_threshold("debug")
logger$add_appender(lgr::AppenderJson$new(logfile), name = "json")

future::plan("multisession", workers = n_workers)
message("Number of parallel workers: ", future::nbrOfWorkers())

start <- Sys.time()
set.seed(sim_id)
bmr = benchmark(benchmark_grid(tsk_train, learners_lst, rsmp_outer_cv), store_models = FALSE)

end <- Sys.time()
difftime(end, start, units="hours")

#bmr$resample_results$resample_result[[1]]$errors

# msr_txt = c("surv.rcll", "surv.cindex", "surv.dcalib")
# measures = msrs(msr_txt)
msrs =  list(
  msr("surv.dcalib", id = "D-Calibration"),
  msr("surv.graf", id = "Brier score"),
  msr("surv.cindex", id = "C_Harrell"),
  msr("surv.cindex", weight_meth = "G2", id = "C_Uno"),
  msr("time_train", id = "time_train"),
  msr("time_predict", id = "time_predict")
)
#bmr$aggregate(msrs)
Score = bmr$score(msrs)
Score = Score %>% select(-c(task, learner, resampling, prediction_test))
saveRDS(Score, sprintf("%s/Score.rds", file_path))

future::plan("sequential")



