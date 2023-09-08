rm(list = ls())

# Packages ----
library(dplyr)
library(Boruta)
library(ranger)
# library(randomForest)
library(caret)

library(ggplot2)
theme_set(theme_light())

COL = c(black = "black"
        ,red = rgb(100, 38, 33, maxColorValue = 100)
        ,green = rgb(38, 77, 19, maxColorValue = 100)
        ,blue = rgb(28, 24, 61, maxColorValue = 100)
        ,purple = rgb(76, 32, 72, maxColorValue = 100)
        ,cyan = rgb(21, 75, 87, maxColorValue = 100)
        ,dark_cyan = rgb(0, 47, 59, maxColorValue = 100)
        ,yellow = rgb(99, 90, 13, maxColorValue = 100)
)
# par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0))
# barplot(rep(1, length(COL)), col = COL, border = NA, axes = FALSE, space = c(0, 0))


# Load data ----
if (TRUE) {
    tmp <- base::sort(list.files(pattern = "datalists_",
                                 path = "./data",
                                 full.names = TRUE),
                      decreasing = TRUE)
    # print(tmp[1])
    load(tmp[1])
    tmp <- base::sort(list.files(pattern = "datacompile_",
                                 path = "./data",
                                 full.names = TRUE),
                      decreasing = TRUE)
    # print(tmp[1])
    DATA <- readRDS(tmp[1])


    ## About missing values ----
    countNA <- apply(DATA, 2, function(x) sum(is.na(x)))
    sort(countNA[countNA > 0], decreasing = TRUE)

    ## Log-transform precipitation ----
    tmp <- grep("precip", names(DATA), ignore.case = TRUE, value = TRUE)
    DATA[, tmp] <- log10(DATA[, tmp] + 1)

    ## Convert all character columns to factors ----
    DATA <- as.data.frame(unclass(DATA), stringsAsFactors = TRUE)
    summary(DATA[, predictors0])

    ## Remove predictors that are not varying ----
    # (if any)
    tmp <- apply(DATA, 2, function(x) all(x == x[1]))
    (pred2remove <- names(tmp)[tmp & !is.na(tmp)])
    predictors <- base::setdiff(predictors, pred2remove)

    DATAnoNA <- na.omit(DATA[, c(RESPONSE, predictors)])
}


# Cross-validation ----
if (TRUE) {
    metrics <- function(obs, pred, w = NULL) {
        # obs = observed data
        # pred = predicted data
        # w = weights for weighted calculations (e.g., geographic weights)

        # If weights are not given, calculate with equal weights
        if (is.null(w)) {
            w <- rep.int(1L, length(obs))
        }
        wsum <- sum(w)
        obsmean <- sum(obs * w) / wsum
        e <- obs - pred
        MAE <- sum(abs(e) * w) / wsum
        MSE <- sum(e^2 * w) / wsum
        RMSE <- sqrt(MSE)
        R2 <- 1 - sum(e^2 * w) / sum((obs - obsmean)^2 * w)
        c(MAE = MAE, RMSE = RMSE, R2 = R2)
    }

    # Split the data into folds for k-fold cross-validation
    K = 10
    set.seed(123456789)
    FOLDS <- caret::createFolds(DATAnoNA[, RESPONSE], k = K)
    # Actual size of each fold
    n_fold <- sapply(FOLDS, length)

    # Number of trees in random forests
    NTREE = 500

    ## Option 1: All variables ----
    set.seed(123456789)
    M <- numeric()
    for (k in 1:K) { # k = 1
        rf <- ranger(dependent.variable.name = RESPONSE,
                     data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, predictors)],
                     respect.unordered.factors = 'order',
                     num.trees = NTREE)
        yhat <- predict(rf, data = DATAnoNA[FOLDS[[k]], predictors])$predictions
        M <- rbind(M,
                   metrics(DATAnoNA[FOLDS[[k]], RESPONSE], yhat))
    }
    M1 <- as_tibble(M) %>%
        mutate(p = length(predictors),
               n = n_fold,
               Version = "All")
    saveRDS(M1, paste0("./dataderived/CV_All_", Sys.Date(), ".rds"))


    ## Option 2: Boruta ----
    set.seed(123456789)
    M <- numeric()
    for (k in 1:K) { # k = 1
        B <- Boruta::Boruta(as.formula(paste(RESPONSE, ".", sep = " ~ ")),
                            doTrace = 0, maxRuns = 500,
                            data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, predictors)])
        st <- attStats(B)
        # print(st)
        # v = rownames(st)[st$decision != "Rejected"]
        v <- rownames(st)[st$decision == "Confirmed"]
        rf <- ranger(dependent.variable.name = RESPONSE,
                     data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, v)],
                     respect.unordered.factors = 'order',
                     num.trees = NTREE)
        yhat <- predict(rf, data = DATAnoNA[FOLDS[[k]], v])$predictions
        M <- rbind(M,
                   c(metrics(DATAnoNA[FOLDS[[k]], RESPONSE], yhat), p = length(v)))
    }
    M2 <- as_tibble(M) %>%
        mutate(n = n_fold,
               Version = "Boruta")
    saveRDS(M2, paste0("./dataderived/CV_Boruta_", Sys.Date(), ".rds"))


    ## Option 3: InteractionForest ----
    library(diversityForest)
    set.seed(123456789)
    M <- numeric()
    for (k in 1:K) { # k = 1
        rfinter <- interactionfor(dependent.variable.name = RESPONSE,
                                  data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, predictors)],
                                  num.trees = NTREE)
        v <- predictors
        yhat <- predict(rfinter, data = DATAnoNA[FOLDS[[k]], v])$predictions
        M <- rbind(M,
                   c(metrics(DATAnoNA[FOLDS[[k]], RESPONSE], yhat), p = length(v)))
    }
    M3 <- as_tibble(M) %>%
        mutate(n = n_fold,
               Version = "InteractionForest")
    saveRDS(M3, paste0("./dataderived/CV_InteractionForest_", Sys.Date(), ".rds"))


    ## Option 4: InteractionForest Recursive ----
    library(diversityForest)
    set.seed(123456789)
    M <- numeric()
    for (k in 1:K) { # k = 1
        print(k)
        rfinter <- interactionfor(dependent.variable.name = RESPONSE,
                                  data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, predictors)],
                                  num.trees = NTREE)
        while (any(rfinter$eim.univ.sorted < 0)) {
            v_univ <- names(rfinter$eim.univ.sorted)[rfinter$eim.univ.sorted > 0]
            v <- c(v_univ
                   # , v_quant, v_qual
            ) %>% unique() %>% sort()
            rfinter <- interactionfor(dependent.variable.name = RESPONSE,
                                      data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, v)],
                                      num.trees = NTREE)
        }
        yhat <- predict(rfinter, data = DATAnoNA[FOLDS[[k]], v])$predictions
        M <- rbind(M,
                   c(metrics(DATAnoNA[FOLDS[[k]], RESPONSE], yhat), p = length(v)))
    }
    M4 <- as_tibble(M) %>%
        mutate(n = n_fold,
               Version = "InteractionForestRec")
    saveRDS(M4, paste0("./dataderived/CV_InteractionForestRec_", Sys.Date(), ".rds"))


    ## Option 5: Boruta Recursive ----
    set.seed(123456789)
    M <- numeric()
    for (k in 1:K) { # k = 1
        B <- Boruta::Boruta(as.formula(paste(RESPONSE, ".", sep = " ~ ")),
                            doTrace = 0, maxRuns = 500,
                            data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, predictors)])
        st <- attStats(B)
        # print(st)
        v <- rownames(st)[st$decision != "Rejected"]
        # v <- rownames(st)[st$decision == "Confirmed"]
        while (any(st$decision == "Rejected")) {
            B <- Boruta::Boruta(as.formula(paste(RESPONSE, ".", sep = " ~ ")),
                                doTrace = 0, maxRuns = 500,
                                data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, v)])
            st <- attStats(B)
            v <- rownames(st)[st$decision != "Rejected"]
        }
        rf <- ranger(dependent.variable.name = RESPONSE,
                     data = DATAnoNA[-FOLDS[[k]], c(RESPONSE, v)],
                     respect.unordered.factors = 'order',
                     num.trees = NTREE)
        yhat <- predict(rf, data = DATAnoNA[FOLDS[[k]], v])$predictions
        M <- rbind(M,
                   c(metrics(DATAnoNA[FOLDS[[k]], RESPONSE], yhat), p = length(v)))
    }
    M5 <- as_tibble(M) %>%
        mutate(n = n_fold,
               Version = "BorutaRec")
    saveRDS(M5, paste0("./dataderived/CV_BorutaRec_", Sys.Date(), ".rds"))
}

# See the updated code for plots in "hatchery_RF.qmd"
M <- dplyr::bind_rows(M1, M2, M3, M4, M5)
M_long <- reshape2::melt(M,
                         variable.name = "Metric",
                         id.vars = c("n", "p", "Version"))

ggplot(M_long, aes(x = Version, y = value, fill = Version)) +
    geom_boxplot() +
    facet_wrap(vars(Metric), ncol = 3, scales = "free_y") +
    xlab("") + ylab("Value") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Dark2")

ggplot(M_long, aes(x = Version, y = p, fill = Version)) +
    geom_boxplot() +
    xlab("") + ylab("p") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Dark2")

tapply(M_long$p, M_long$Version, mean)


am = "BH"
attach(M)
pairwise.t.test(MAE, Version,
                paired = TRUE,
                p.adjust.method = am)
pairwise.t.test(RMSE, Version,
                paired = TRUE,
                p.adjust.method = am)
pairwise.t.test(R2, Version,
                paired = TRUE,
                p.adjust.method = am)
detach(M)


# https://www.datanovia.com/en/lessons/pairwise-t-test/
library(rstatix)
library(ggpubr)
pwc <- M %>%
    pairwise_t_test(MAE ~ Version,
                    paired = TRUE,
                    p.adjust.method = am)
# pwc
pwc <- pwc %>% add_xy_position(x = "Version", step.increase = 1)
ggboxplot(M, x = "Version", y = "MAE") +
    stat_pvalue_manual(pwc, hide.ns = TRUE)


# Variable selection ----

# Continue with one strategy (Boruta recursive)
set.seed(123456789)

B <- Boruta::Boruta(as.formula(paste(RESPONSE, ".", sep = " ~ ")),
                    doTrace = 0, maxRuns = 500,
                    data = DATAnoNA[, c(RESPONSE, predictors)])
st <- attStats(B)
# print(st)
v <- rownames(st)[st$decision != "Rejected"]
# v <- rownames(st)[st$decision == "Confirmed"]
while (any(st$decision == "Rejected")) {
    B <- Boruta::Boruta(as.formula(paste(RESPONSE, ".", sep = " ~ ")),
                        doTrace = 0, maxRuns = 500,
                        data = DATAnoNA[, c(RESPONSE, v)])
    st <- attStats(B)
    v <- rownames(st)[st$decision != "Rejected"]
}
saveRDS(B, paste0("./dataderived/VarSel_BorutaRec_B_", Sys.Date(), ".rds"))
saveRDS(v, paste0("./dataderived/VarSel_BorutaRec_v_", Sys.Date(), ".rds"))


# RF tuning and summary ----
# Run "Load data" first
# rm(list = ls())

library(ranger)
library(caret)

pfun <- function(object, newdata) {  # prediction wrapper
    unname(predict(object, data = newdata)$predictions)
}

# Number of trees in random forests
NTREE = 500

# Load most recent variable selection
tmp <- base::sort(list.files(pattern = "BorutaRec_v_",
                             path = "./dataderived",
                             full.names = TRUE),
                  decreasing = TRUE)
predictors <- readRDS(tmp[1])

# Redefine working data subset
DATAnoNA <- na.omit(DATA[, c(RESPONSE, predictors)])

# Re-estimate random forest
set.seed(123456789)
(rf1 <- ranger(Yield ~ .,
               num.trees = NTREE,
               write.forest = TRUE,
               data = DATAnoNA[, c(RESPONSE, predictors)]))
saveRDS(rf1, paste0("./dataderived/RF_rf1_", Sys.Date(), ".rds"))

# Settings for cross-validation
set.seed(123456789)
train.control.reg <- caret::trainControl(method = "repeatedcv",
                                         number = 10,
                                         repeats = 1)
rfGrid.reg <- expand.grid(splitrule = "variance",
                          mtry = c(5, 7, 10, 15, 20, 25, 30, 40, 50),
                          min.node.size = c(1, 3, 5, 10))
(rfCaret <- train(Yield ~ .,
                  data = DATAnoNA[sample(1:nrow(DATAnoNA), replace = FALSE),
                                  c(RESPONSE, predictors)],
                  method = "ranger",
                  trControl = train.control.reg,
                  tuneGrid = rfGrid.reg,
                  num.trees = NTREE,
                  importance = "none",
                  metric = "Rsquared"))
saveRDS(rfCaret, paste0("./dataderived/RF_rfCaret_", Sys.Date(), ".rds"))
(rf2 <- rfCaret$finalModel)


# Shapley ----
# Run "Load data" first
# rm(list = ls())

library(fastshap)
library(shapviz)
library(ranger)

pfun <- function(object, newdata) {  # prediction wrapper
    unname(predict(object, data = newdata)$predictions)
}

# Number of trees in random forests
NTREE = 500

# Load most recent variable selection
# tmp <- base::sort(list.files(pattern = "BorutaRec_B_",
#                              path = "./dataderived",
#                              full.names = TRUE),
#                   decreasing = TRUE)
# B <- readRDS(tmp[1])
tmp <- base::sort(list.files(pattern = "BorutaRec_v_",
                             path = "./dataderived",
                             full.names = TRUE),
                  decreasing = TRUE)
predictors <- readRDS(tmp[1])

# Redefine working data subset
DATAnoNA <- na.omit(DATA[, c(RESPONSE, predictors)])

tmp <- base::sort(list.files(pattern = "RF_rfCaret_",
                             path = "./dataderived",
                             full.names = TRUE),
                  decreasing = TRUE)
rfCaret <- readRDS(tmp[1])
rf2 <- rfCaret$finalModel

(baseline <- mean(pfun(rf2, newdata = DATAnoNA)))
mean(DATAnoNA$Yield)
saveRDS(baseline, paste0("./dataderived/RF_baseline_", Sys.Date(), ".rds"))


library(doParallel)
cl <- parallel::makeCluster(parallel::detectCores())
registerDoParallel(cl)
# Load ranger to use it in parallel
# https://stackoverflow.com/questions/34022083/parallel-computation-loading-packages-in-each-thread-only-once
# https://tshafer.com/blog/2020/08/r-packages-s3-methods
clusterCall(cl, function() library(ranger))

# Feature importance scores for the whole dataset
set.seed(123456789)
shap_rf1 <- fastshap::explain(rf2,
                              X = DATAnoNA[, predictors],
                              pred_wrapper = pfun,
                              nsim = 100,
                              parallel = TRUE,
                              adjust = TRUE,
                              shap_only = FALSE)
saveRDS(shap_rf1, paste0("./dataderived/shap_rfall_", Sys.Date(), ".rds"))

# Feature importance scores for well-predicted crashes
# Select well-predicted cases
DATAnoNA <- DATAnoNA %>%
    mutate(Fitted = predict(rf2, data = DATAnoNA)$predictions) %>%
    mutate(AbsError = abs(Yield - Fitted))
iwp <- with(DATAnoNA, which(AbsError < 3 & Yield < 1))
length(iwp)
set.seed(123456789)
shap_rfcrashes <- fastshap::explain(rf2,
                                    X = DATAnoNA[iwp, predictors],
                                    pred_wrapper = pfun,
                                    nsim = 100,
                                    parallel = TRUE,
                                    adjust = TRUE,
                                    shap_only = FALSE)
shap_rfcrashes <- list(shap = shap_rfcrashes,
                       iwp = iwp,
                       DATAnoNA = DATAnoNA)
saveRDS(shap_rfcrashes, paste0("./dataderived/shap_rfcrashes_", Sys.Date(), ".rds"))


# Feature importance scores for well-predicted high yields
# Select well-predicted cases
iwp <- with(DATAnoNA, which(AbsError < 3 & Yield > 20))
length(iwp)
set.seed(123456789)
shap_rfhigh <- fastshap::explain(rf2,
                                 X = DATAnoNA[iwp, predictors],
                                 pred_wrapper = pfun,
                                 nsim = 100,
                                 parallel = TRUE,
                                 adjust = TRUE,
                                 shap_only = FALSE)
shap_rfhigh <- list(shap = shap_rfhigh,
                    iwp = iwp,
                    DATAnoNA = DATAnoNA)
saveRDS(shap_rfhigh, paste0("./dataderived/shap_rfhigh_", Sys.Date(), ".rds"))
