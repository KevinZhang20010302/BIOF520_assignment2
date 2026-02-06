Final_BIOF520_Assignment2
================
Kevin Zhang
2026-02-04

### read data

``` r
knowles <- readRDS("~/Documents/BIOF520/Assignment2/knowles_matched_TaLG_final.rds")
UROMOL <- readRDS("~/Documents/BIOF520/Assignment2/UROMOL_TaLG.teachingcohort.rds")

knowles_clinical <- knowles %>% select(knowles_ID,
                                       Age,
                                       Sex,
                                       Concomitant.CIS,
                                       UROMOL2021.classification,
                                       BCG,
                                       FUtime_days.,
                                       Recurrence,
                                       RFS_time
                                       ) 

UROMOL_clinical <- UROMOL %>% select(UROMOL.ID,
                                       Age,
                                       Sex,
                                       Concomitant.CIS,
                                       UROMOL2021.classification,
                                       BCG,
                                       FUtime_days.,
                                       Recurrence,
                                       RFS_time
                                       )

knowles_gene <- knowles$exprs
UROMOL_gene <- UROMOL$exprs
```

### remove batch effect (normalization for cross-platform (Microarray vs RNA-seq))

``` r
knowles_gene_t <- t(knowles_gene) 
UROMOL_gene_t <- t(UROMOL_gene) #convert gene name to rkw

both_gene_name <- intersect(rownames(knowles_gene_t), rownames(UROMOL_gene_t))

knowles_UROMOL_gene_t <- cbind(knowles_gene_t[both_gene_name, ],
                               UROMOL_gene_t[both_gene_name, ])

knowles_UROMOL_batch <- c(rep("knowles", ncol(knowles_gene_t)), rep("UROMOL", ncol(UROMOL_gene_t)))

normalized_knowles_UROMOL_gene_t <- removeBatchEffect(knowles_UROMOL_gene_t, batch = knowles_UROMOL_batch)
```

    ## design matrix of interest not specified. Assuming a one-group experiment.

``` r
normalized_knowles_gene_t <- normalized_knowles_UROMOL_gene_t[, knowles_UROMOL_batch == "knowles"]
normalized_UROMOL_gene_t <- normalized_knowles_UROMOL_gene_t[, knowles_UROMOL_batch == "UROMOL"]

normalized_knowles_gene <- t(normalized_knowles_gene_t)
normalized_UROMOL_gene <- t(normalized_UROMOL_gene_t)
```

### clean the clinical data based on survival data

``` r
### first convert RFS_time by FUtime_days/30 if recurrence is not null but RFS_time is null
knowles_index <- !is.na(knowles_clinical$Recurrence) & is.na(knowles_clinical$RFS_time)
knowles_clinical_final <- knowles_clinical
knowles_clinical_final$RFS_time[knowles_index] <- knowles_clinical_final$FUtime_days.[knowles_index] / 30

UROMOL_index <- !is.na(UROMOL_clinical$Recurrence) & is.na(UROMOL_clinical$RFS_time)
UROMOL_clinical_final <- UROMOL_clinical
UROMOL_clinical_final$RFS_time[UROMOL_index] <- UROMOL_clinical_final$FUtime_days.[UROMOL_index] / 30


### remove NA in both RFS and recurrence
remove_k_index <- is.na(knowles_clinical$Recurrence) & is.na(knowles_clinical$RFS_time)
knowles_clinical_final <- knowles_clinical_final[!remove_k_index, ]
keep_k_rows <- rownames(knowles_clinical_final)
knowles_gene_final <- normalized_knowles_gene[keep_k_rows, ]

remove_u_index <- is.na(UROMOL_clinical$Recurrence) & is.na(UROMOL_clinical$RFS_time)
UROMOL_clinical_final <- UROMOL_clinical_final[!remove_u_index, ]
keep_u_rows <- rownames(UROMOL_clinical_final)
UROMOL_gene_final <- normalized_UROMOL_gene[keep_u_rows, ]
```

### further clean clinical data

``` r
knowles_clinical_final_final <- knowles_clinical_final %>%
  mutate(Sex = factor(Sex),
         Concomitant.CIS = factor(Concomitant.CIS, levels = c("No", "Yes")),
         UROMOL2021.classification = factor(UROMOL2021.classification, levels = c("Class_1", "Class_2a", "Class_2b", "Class_3"),
                                            labels = c("Class 1", "Class 2a", "Class 2b", "Class 3")),
         BCG = factor(BCG, levels = c("0", "1")),
         risk_classification = case_when(RFS_time > 24 ~ "low",
                                         RFS_time <= 24 & Recurrence == 1 ~ "high",
                                         RFS_time <= 24 & Recurrence == 0 ~ "unknown")
         )

UROMOL_clinical_final_final <- UROMOL_clinical_final %>%
  mutate(Sex = factor(Sex),
         Concomitant.CIS = factor(Concomitant.CIS),
         UROMOL2021.classification = factor(UROMOL2021.classification),
         BCG = factor(BCG),
         risk_classification = case_when(RFS_time > 24 ~ "low",
                                         RFS_time <= 24 & Recurrence == 1 ~ "high",
                                         RFS_time <= 24 & Recurrence == 0 ~ "unknown")
         )
```

### make the x variables (meanwhile remove NAs)

``` r
X_clinical_k <- model.matrix(~ Age + Sex + Concomitant.CIS + UROMOL2021.classification + BCG, data = knowles_clinical_final_final)[,-1]

X_clinical_u <- model.matrix(~ Age + Sex + Concomitant.CIS + UROMOL2021.classification + BCG, data = UROMOL_clinical_final_final)[,-1]


keep_k_rows <- rownames(X_clinical_k)
knowles_clinical_remove_NA <- knowles_clinical_final_final[keep_k_rows,]
knowles_gene_keep <- knowles_gene_final[keep_k_rows, ]
knowles_full <- cbind(X_clinical_k, knowles_gene_keep)

keep_u_rows <- rownames(X_clinical_u)
UROMOL_clinical_remove_NA <- UROMOL_clinical_final_final[keep_u_rows,]
UROMOL_gene_keep <- UROMOL_gene_final[keep_u_rows, ]
UROMOL_full <- cbind(X_clinical_u, UROMOL_gene_keep)
```

### make the y variables

``` r
survival_data_k <- with(
  knowles_clinical_remove_NA,
  Surv(time = RFS_time, event = Recurrence)
)


survival_data_u <- with(
  UROMOL_clinical_remove_NA,
  Surv(time = RFS_time, event = Recurrence)
)
```

### train data and test data spilt

``` r
set.seed(1234)
u_length <- nrow(X_clinical_u)

train_size <- floor(0.8 * u_length)

train_index <- sample(seq_len(u_length), size = train_size)


UROMOL_train <- UROMOL_full[train_index, ]
UROMOL_test <- UROMOL_full[-train_index, ]

UROMOL_survival_train <- survival_data_u[train_index]
UROMOL_survival_test <- survival_data_u[-train_index]
```

### cox model fit

``` r
set.seed(1234)
cox_cv_fit <- cv.glmnet(
    x = UROMOL_train,
    y = UROMOL_survival_train,
    family = "cox",
    alpha = 1,
    nfold = 10
  )


cox_penalty <- glmnet(
  x = UROMOL_train,
  y = UROMOL_survival_train,
  family = "cox",
  alpha = 1,
  lambda = cox_cv_fit$lambda.min
)

U_train_risk_score<- predict(cox_penalty, newx = UROMOL_train, type = "link")

no_recurrence_model <- coxph(UROMOL_survival_train ~ U_train_risk_score)

coefficient_matrix <- as.matrix(coef(cox_penalty))
non_zero_coefficient <- rownames(coefficient_matrix)[coefficient_matrix!=0]
non_zero_coefficient 
```

    ## [1] "SLC44A1"  "CDC37L1"  "GLE1"     "AGAP2"    "PLAA"     "CDKN2B"   "MLLT3"   
    ## [8] "AATK"     "ARHGEF25"

\###test on test data

``` r
### test risk score prediction
U_test_risk_score<- predict(cox_penalty, newx = UROMOL_test, type = "link")
### classification

new_data_df <- data.frame(U_train_risk_score = as.numeric(U_test_risk_score))

U_test_prediction_curves <- survfit(no_recurrence_model, newdata = new_data_df)

U_test_24_month_probability <- summary(U_test_prediction_curves, times = 24)$surv

UROMOL_test_result_df <- as.data.frame(UROMOL_test)

UROMOL_test_result_df$predict_classification <- ifelse(as.numeric(U_test_24_month_probability) > 0.5, "low", "high")
```

### KM plot

``` r
u_predict_classification <- UROMOL_test_result_df$predict_classification

km_curve <- survfit(UROMOL_survival_test ~ u_predict_classification)

km_plot <- ggsurvplot(km_curve,
                      data = data.frame(UROMOL_survival_test, u_predict_classification),
                      risk.table = TRUE,
                      pval = TRUE,
                      conf.int = TRUE,
                      palette = c("red", "blue"),
                      title = "Kaplan-Meier plot for internal test",
                      legend.labs = c("High Risk", "Low Risk"),
                      xlab = "Time (Months)",
                      ylab = "Probabilty Of No Recurrent",
                      font.title = 10,
                      font.subtitle = 10,
                      font.caption = 10,
                      font.x = 10,
                      font.y = 10,
                      font.tickslab = 10,
                      font.legend = 10,
                      risk.table.fontsize = 3.5)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the ggpubr package.
    ##   Please report the issue at <https://github.com/kassambara/ggpubr/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
km_plot
```

    ## Ignoring unknown labels:
    ## • colour : "Strata"

![](Untitled_files/figure-gfm/unnamed-chunk-10-1.png)<!-- --> \###
c-index calculation for internal test

``` r
flip_u_risk_score <- -1 * U_test_risk_score

u_test_c_index <- concordance(UROMOL_survival_test~ flip_u_risk_score)$concordance

u_test_c_index
```

    ## [1] 0.5585586

### km of EAU risk classification

``` r
test_rowname <- rownames(UROMOL_test)

EAU_test<- UROMOL[test_rowname,]$EAU.risk

EAU_no_na_index <- is.na(EAU_test)

EAU_test_remove_NA <- EAU_test[!EAU_no_na_index]

EAU_binary <- ifelse(EAU_test_remove_NA == "Low", "Low Risk", "High/Intermediate Risk")

UROMOL_survival_test_EAU <- UROMOL_survival_test[!EAU_no_na_index]

EAU_km_curve <-  survfit(UROMOL_survival_test_EAU  ~ EAU_binary)

EAU_km_plot <- ggsurvplot(EAU_km_curve,
                      data = data.frame(UROMOL_survival_test_EAU, EAU_binary),
                      risk.table = TRUE,
                      pval = TRUE,
                      conf.int = TRUE,
                      palette = c("red", "blue"),
                      title = "Kaplan-Meier plot for internal test EAU",
                      legend.labs = c("High/Intermediate Risk", "Low Risk"),
                      xlab = "Time (Months)",
                      ylab = "Probabilty Of No Recurrent",
                      font.title = 10,
                      font.subtitle = 10,
                      font.caption = 10,
                      font.x = 10,
                      font.y = 10,
                      font.tickslab = 10,
                      font.legend = 10,
                      risk.table.fontsize = 3.5)

EAU_km_plot
```

    ## Ignoring unknown labels:
    ## • colour : "Strata"

![](Untitled_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### AUC for my classification

``` r
u_test_id <- rownames(UROMOL_test_result_df)

UROMOL_ground_truth <- UROMOL_clinical_final_final[u_test_id, ]$risk_classification

u_test_roc_data <- data.frame(
  ground_truth_risk = UROMOL_ground_truth,
  predicted_risk = UROMOL_test_result_df$predict_classification
)

u_test_roc_data_nounkown <- u_test_roc_data[u_test_roc_data$ground_truth_risk != "unknown", ]

u_test_roc_data_nounkown <- u_test_roc_data_nounkown %>% mutate(
  ground_truth_risk = ifelse(ground_truth_risk == "high", 1, 0),
  predicted_risk = ifelse(predicted_risk == "high", 1, 0),
)

u_test_roc <- roc(u_test_roc_data_nounkown$ground_truth_risk, u_test_roc_data_nounkown$predicted_risk)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
roc_plot<-ggroc(u_test_roc, color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey")+
  theme_minimal() +
  labs(title = "ROC plot for new classification model",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  annotate("text", x=0.3, y = 0.2, label = paste("AUC =", round(auc(u_test_roc), 2)), size = 6, color = "red") +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

roc_plot
```

![](Untitled_files/figure-gfm/unnamed-chunk-13-1.png)<!-- --> \### AUC
for old classification

``` r
UROMOL_ground_truth_EAU <- UROMOL_ground_truth[!EAU_no_na_index]


eau_u_test_roc_data <- data.frame(
  ground_truth_risk = UROMOL_ground_truth_EAU,
  predicted_risk = EAU_binary
)

eau_u_test_roc_data_nounkown <- eau_u_test_roc_data[eau_u_test_roc_data$ground_truth_risk != "unknown", ]

eau_u_test_roc_data_nounkown <- eau_u_test_roc_data_nounkown %>% mutate(
  ground_truth_risk = ifelse(ground_truth_risk == "high", 1, 0),
  predicted_risk = ifelse(predicted_risk == "High/Intermediate Risk", 1, 0),
)

eau_u_test_roc <- roc(eau_u_test_roc_data_nounkown$ground_truth_risk, eau_u_test_roc_data_nounkown$predicted_risk)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
eau_roc_plot<-ggroc(eau_u_test_roc, color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey")+
  theme_minimal() +
  labs(title = "ROC plot for EAU classification model",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  annotate("text", x=0.3, y = 0.2, label = paste("AUC =", round(auc(eau_u_test_roc), 2)), size = 6, color = "red") +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

eau_roc_plot
```

![](Untitled_files/figure-gfm/unnamed-chunk-14-1.png)<!-- --> \###
external validation

``` r
k_test_risk_score<- predict(cox_penalty, newx = knowles_full, type = "link")
### classification

new_data_df <- data.frame(U_train_risk_score = as.numeric(k_test_risk_score))

K_test_prediction_curves <- survfit(no_recurrence_model, newdata = new_data_df)

K_test_24_month_probability <- summary(K_test_prediction_curves, times = 24)$surv

knowles_test_result_df <- as.data.frame(knowles_full)

knowles_test_result_df$predict_classification <- ifelse(as.numeric(K_test_24_month_probability) > 0.5, "low", "high")
```

### km plot for external

``` r
k_predict_classification <- knowles_test_result_df$predict_classification

k_km_curve <- survfit(survival_data_k ~ k_predict_classification)

k_km_plot <- ggsurvplot(k_km_curve,
                      data = data.frame(survival_data_k, k_predict_classification),
                      risk.table = TRUE,
                      pval = TRUE,
                      conf.int = TRUE,
                      palette = c("red", "blue"),
                      title = "Kaplan-Meier plot for external test",
                      legend.labs = c("High Risk", "Low Risk"),
                      xlab = "Time (Months)",
                      ylab = "Probabilty Of No Recurrent",
                      font.title = 10,
                      font.subtitle = 10,
                      font.caption = 10,
                      font.x = 10,
                      font.y = 10,
                      font.tickslab = 10,
                      font.legend = 10,
                      risk.table.fontsize = 3.5)

k_km_plot
```

    ## Ignoring unknown labels:
    ## • colour : "Strata"

![](Untitled_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### c-index for external

``` r
flip_k_risk_score <- -1 * k_test_risk_score

k_test_c_index <- concordance(survival_data_k~ flip_k_risk_score)$concordance

k_test_c_index
```

    ## [1] 0.5414938

### ROC plot for external

``` r
k_id <- rownames(knowles_test_result_df)

knowles_ground_truth <- knowles_clinical_final_final[k_id, ]$risk_classification

k_test_roc_data <- data.frame(
  ground_truth_risk = knowles_ground_truth,
  predicted_risk = knowles_test_result_df$predict_classification
)

k_test_roc_data_nounkown <- k_test_roc_data[k_test_roc_data$ground_truth_risk != "unknown", ]

k_test_roc_data_nounkown <- k_test_roc_data_nounkown %>% mutate(
  ground_truth_risk = ifelse(ground_truth_risk == "high", 1, 0),
  predicted_risk = ifelse(predicted_risk == "high", 1, 0),
)

k_test_roc <- roc(k_test_roc_data_nounkown$ground_truth_risk, k_test_roc_data_nounkown$predicted_risk)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
k_roc_plot<-ggroc(k_test_roc, color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey")+
  theme_minimal() +
  labs(title = "ROC plot for new classification model",
       x = "False Positive Rate",
       y = "True Positive Rate") +
  annotate("text", x=0.3, y = 0.2, label = paste("AUC =", round(auc(k_test_roc), 2)), size = 6, color = "red") +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))

k_roc_plot
```

![](Untitled_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
### package a final function
risk_classifer <- function(x, t, p) {
  risk_score <- predict(cox_penalty, newx = x, type = "link")
  new_data_df <- data.frame(U_train_risk_score = as.numeric(risk_score))
  prediction_curves <- survfit(no_recurrence_model, newdata = new_data_df)
  time_month_probability <- summary(prediction_curves, times = t)$surv
  result_df <- as.data.frame(x)
  result_df$predict_probability <- as.numeric(time_month_probability)
  result_df$predict_classification <- ifelse(as.numeric(time_month_probability) > p, "low", "high")
  return(result_df)
}
```

### test final packed model correct or not

``` r
test_final_u <- risk_classifer(UROMOL_test, 24, 0.5)
test_final_k <- risk_classifer(knowles_full, 24, 0.5)
identical(test_final_u$predict_classification, UROMOL_test_result_df$predict_classification)
```

    ## [1] TRUE

``` r
identical(test_final_k$predict_classification, knowles_test_result_df$predict_classification)
```

    ## [1] TRUE

``` r
test_final_u$predict_probability
```

    ##  [1] 0.4969424874 0.8342755501 0.7899691469 0.4816987431 0.3380995233
    ##  [6] 0.2984850599 0.6921437339 0.4869617286 0.2031952117 0.5875558826
    ## [11] 0.5152606010 0.6853978433 0.3166444752 0.4238734040 0.3321353835
    ## [16] 0.7920531641 0.4730714409 0.4544371089 0.5488801260 0.3229843558
    ## [21] 0.6338207986 0.4237195967 0.3969447053 0.4094875181 0.5072108285
    ## [26] 0.4305053109 0.3113408948 0.4348153050 0.3171814405 0.5178015655
    ## [31] 0.7126199832 0.5654587474 0.5020090663 0.4273239153 0.4375161918
    ## [36] 0.6277062696 0.5452856783 0.6216707703 0.1664813351 0.5782210283
    ## [41] 0.9436741475 0.5899967097 0.1340072769 0.4759869617 0.6022242842
    ## [46] 0.6204262683 0.0003852435 0.2169915957 0.8644861375 0.2243567829
    ## [51] 0.7130960330 0.4732772551 0.4388523372 0.3110031510 0.4306962930

``` r
test_final_k$predict_probability
```

    ##  [1] 0.4842936 0.4967047 0.5050541 0.5365009 0.4298365 0.4368755 0.4936332
    ##  [8] 0.5497534 0.3270710 0.5148844 0.5834642 0.4801497 0.3458208 0.5566127
    ## [15] 0.5677590 0.3859472 0.5239131 0.5393914 0.4221532 0.5119697 0.4887504
    ## [22] 0.4192202 0.4856726 0.3370702 0.5027674 0.5932516 0.5806820 0.3313785
    ## [29] 0.4486333 0.4935558 0.5294791 0.3149634 0.3335036 0.4723144 0.5105396
    ## [36] 0.5456937 0.4544004 0.4855861 0.5058032 0.4693813 0.5731007 0.5384838
    ## [43] 0.5262153 0.5465407 0.4042575 0.5234577 0.6177719 0.5345024 0.3758389
    ## [50] 0.5301829 0.5347118 0.5209364 0.4523243 0.4967796 0.4970146 0.5374184
    ## [57] 0.4431884 0.4858449 0.5352796 0.5123137 0.3409271 0.5840917 0.4976846
    ## [64] 0.5317745 0.3129040 0.4588091 0.3957418 0.4923618 0.3273526 0.5466448
    ## [71] 0.5321372 0.4238651 0.3449979
