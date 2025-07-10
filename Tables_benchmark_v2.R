library(tidyverse)
library(ggpubr)
library(rstatix)
setwd("~/SimulationMar1/3_results/03May2025")
# Score_grace <- readRDS("~/SimulationMar1/3_results/grace_19APR25/Score.rds")
# Score_smarto <- readRDS("~/SimulationMar1/3_results/smarto_18APR25/Score.rds")
# Score_cost <- readRDS("~/SimulationMar1/3_results/cost_15APR25/Score.rds")
# Score_grace <- readRDS("grace/Score.rds")
# Score_smarto <- readRDS("smarto/Score.rds")
# Score_cost <- readRDS("cost/Score.rds")
# Score_trace <- readRDS("trace/Score.rds")
# df = bind_rows(list(grace = Score_grace, smarto=Score_smarto, cost=Score_cost, trace=Score_trace), .id = "Dataset")

lst_fls = grep("benchmark.*", list.files(), value = TRUE, invert = TRUE)

df = map_dfr(lst_fls, ~{
  if(any(grep("Score", list.files(.x))))
  df_scr = readRDS(sprintf("%s/Score.rds", .x))
})


df = df %>% mutate(learner_id1 = str_replace_all(learner_id, pattern = "class\\.binary_enc\\.removeconstants\\.reg\\.scale\\.|surv\\.|class_|nonpar_|semipar_|par_|ml_|ann_|ranfor_|\\.tuned|new", ""),
                   learner_id1 = str_replace_all(learner_id1, c("_cv_"="_cars_", "_mr_"="_mrmr_", "aorsf_min"="aorsf_lasso_min", "aorsf_se"="aorsf_lasso_net", "rfsrc_brier"="RSF")),
                   learner_id1 = str_replace_all(learner_id1, c("_"="-", "aorsf"="ORSF", "penalized" = "Penalized Cox PH")),
                   learner_id2 = factor(learner_id1, levels=c("coxph", "Penalized Cox PH", "RSF", "ORSF-fast", "ORSF-cph", "ORSF-net", "ORSF-lasso-min", "ORSF-lasso-net", "ORSF-mrmr-1q", "ORSF-mrmr-2q", "ORSF-mrmr-3q", "ORSF-cars-fndr", "ORSF-cars-1q", "ORSF-cars-2q", "ORSF-cars-3q"))) %>%
  filter(!str_detect(learner_id1, "kaplan|cvglmnet|aft"))

df = df %>% 
  mutate(task_id = factor(task_id, levels=c("whas", "cost", "grace", "trace", "smarto", "prostate", "breast", "colon", "leukemia1", "rotterdam"), 
                          labels=c("whas (N=461)", "COST (N=518)", "GRACE (N=1000)", "TRACE (N=1878)", "SMARTO (N=3873)", "prostate (N=502)", "GBSG2 (N=686)", "colon (N=929)", "leukemia (N=1043)", "rotterdam (N=2982)")))

df_res = df %>% group_by(task_id, learner_id2) %>% 
  summarise_at(vars(`D-Calibration` : time_predict), mean) %>% 
  pivot_wider(id_cols = learner_id2, 
              names_from = task_id, 
              names_glue = "{task_id}_{.value}", 
              values_from = `D-Calibration`:time_predict)

write.csv(df, sprintf("benchmark results long %s.csv", Sys.Date()))
write.csv(df_res, sprintf("benchmark results wide %s.csv", Sys.Date()))

df_plot = df %>% 
  filter(str_detect(learner_id2,"ORSF-net|ORSF-lasso|ORSF-mrmr|ORSF-cars" )) %>%
  mutate(learner_cat = case_when(str_detect(learner_id2, "ORSF-net") ~ "ORSF-net",
                                 str_detect(learner_id2, "ORSF-lasso") ~ "ORSF-lasso",
                                 str_detect(learner_id2, "ORSF-mrmr") ~ "ORSF-mrmr",
                                 str_detect(learner_id2, "ORSF-cars") ~ "ORSF-cars"),
         learner_cat = factor(learner_cat, levels=c("ORSF-net", "ORSF-lasso", "ORSF-mrmr", "ORSF-cars")),
         Disease = ifelse(str_detect(tolower(task_id), "cost|whas|grace|trace|smarto"), "CVD", "Cancer"))

p_cvd = ggplot(data = df_plot %>% filter(Disease=="CVD"), 
       aes(x = learner_id2, y = time_train/60/60, color=learner_cat)) + 
  geom_boxplot() + 
  facet_wrap(task_id~., nrow = 5) +
  theme_bw() + 
  labs(x="Methods", y="Trainin time (in Hrs)", color="Method category") +
  theme(axis.text.x = element_text(angle = 15),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.position = "bottom") + 
  scale_y_continuous(limits = c(0, 13), expand = c(0.025, 0.025))

p_onc = ggplot(data = df_plot %>% filter(Disease=="Cancer"), 
               aes(x = learner_id2, y = time_train/60/60, color=learner_cat)) + 
  geom_boxplot() + 
  facet_wrap(task_id~., nrow = 5) +
  theme_bw() + 
  labs(x="Methods", y="Trainin time (in Hrs)", color="Method category") +
  theme(axis.text.x = element_text(angle = 15),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.position = "bottom") + 
  scale_y_continuous(limits = c(0, 13), expand = c(0.025, 0.025))

annotate_figure(ggarrange(p_cvd, p_onc, ncol = 2, common.legend = TRUE, 
                          #widths = c(4, 4), 
                          legend = "bottom", labels = c("A)", "B)"), font.label = list(size = 12, color = "black", face = "plain", family = NULL), hjust = -0.32, vjust = 1.5))

