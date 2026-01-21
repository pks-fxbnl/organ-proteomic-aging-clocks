###########################################
## Demo code for proteomic aging models  ##
## - Data generation (train & test)      ##
## - Conventional age model (all proteins)
## - Organ-specific age model (30 proteins)
###########################################

### ----------------------------------------
### 0. Load required packages
### ----------------------------------------
library(glmnet)

### ----------------------------------------
### 1. Generate demo proteomics datasets
### ----------------------------------------

set.seed(123)  # for reproducibility

# Number of proteins
n_proteins <- 1254

# Protein names (fake but realistic)
protein_names <- paste0("Protein_", seq_len(n_proteins))

### 1.1 Training dataset: 98 healthy participants
n_train <- 98

train_demo <- data.frame(
  SampleID = paste0("Healthy_", seq_len(n_train)),
  matrix(
    rnorm(n_train * n_proteins, mean = 0, sd = 1),  # standardized values
    nrow = n_train,
    ncol = n_proteins
  )
)

colnames(train_demo)[-1] <- protein_names

### 1.2 Test dataset: 1850 HIV participants
n_test <- 1850

test_demo <- data.frame(
  SampleID = paste0("HIV_", seq_len(n_test)),
  matrix(
    rnorm(n_test * n_proteins, mean = 0, sd = 1),  # standardized values
    nrow = n_test,
    ncol = n_proteins
  )
)

colnames(test_demo)[-1] <- protein_names

### 1.3 Convert SampleID to rownames and drop column

# Training dataset
rownames(train_demo) <- train_demo$SampleID
train_demo$SampleID <- NULL

# Test dataset
rownames(test_demo) <- test_demo$SampleID
test_demo$SampleID <- NULL

# Quick checks (optional)
# dim(train_demo)  # 98 × 1254
# dim(test_demo)   # 1850 × 1254
# head(train_demo[, 1:5])

### 1.4 Generate demo chronological ages

# Train: 98 healthy participants (age 20–80)
set.seed(123)  # for reproducibility
train_age_demo <- round(runif(n_train, min = 20, max = 80))

# Test: 1850 HIV participants (not used in training, only as background demo)
set.seed(456)
test_age_demo <- round(runif(n_test, min = 20, max = 80))

# Quick checks (optional)
# length(train_age_demo)  # 98
# length(test_age_demo)   # 1850

### ----------------------------------------
### 2. Helper function: lambda-ensemble LASSO
### ----------------------------------------

lambda_ensemble_lasso <- function(x_train, y_train, x_test,
                                  n_bootstrap = 500,
                                  seed = 1234) {
  set.seed(seed)
  
  n_train <- nrow(x_train)
  n_test  <- nrow(x_test)
  
  # Matrices to store predictions from each bootstrap model
  train_pred_mat <- matrix(NA_real_, nrow = n_train, ncol = n_bootstrap)
  test_pred_mat  <- matrix(NA_real_, nrow = n_test,  ncol = n_bootstrap)
  
  # Vector to store lambda from each bootstrap iteration
  lambda_vec <- numeric(n_bootstrap)
  
  # Data frame to store loss metrics per iteration
  loss_table <- data.frame(
    Iteration = integer(n_bootstrap),
    Lambda    = numeric(n_bootstrap),
    CV_Error  = numeric(n_bootstrap),
    R2        = numeric(n_bootstrap)
  )
  
  # List to store coefficients (including intercept) from each bootstrap model
  coef_list <- vector("list", n_bootstrap)
  
  for (b in seq_len(n_bootstrap)) {
    # 1. Bootstrap resampling of the training set
    sample_idx <- sample.int(n_train, size = n_train, replace = TRUE)
    X_boot <- x_train[sample_idx, , drop = FALSE]
    y_boot <- y_train[sample_idx]
    
    # 2. Cross-validated LASSO on bootstrap sample to select optimal lambda
    cv_fit <- cv.glmnet(
      X_boot, y_boot,
      alpha        = 1,
      nfolds       = 5,
      type.measure = "mse"
    )
    
    lambda_opt <- cv_fit$lambda.min
    lambda_vec[b] <- lambda_opt
    
    # 3. Fit LASSO model on bootstrap sample with selected lambda
    fit_boot <- glmnet(
      X_boot, y_boot,
      alpha  = 1,
      lambda = lambda_opt
    )
    
    # 4. Store coefficients
    coef_vec <- as.vector(coef(fit_boot))
    names(coef_vec) <- rownames(coef(fit_boot))
    coef_list[[b]] <- data.frame(
      Feature     = names(coef_vec),
      Coefficient = coef_vec,
      Iteration   = b,
      row.names   = NULL
    )
    
    # 5. Compute R² on bootstrap sample
    y_boot_pred <- as.numeric(predict(fit_boot, X_boot))
    r2_boot <- 1 - sum((y_boot - y_boot_pred)^2) / sum((y_boot - mean(y_boot))^2)
    
    # 6. Predict on full training data using this bootstrap model
    train_pred_mat[, b] <- as.numeric(predict(fit_boot, x_train))
    
    # 7. Refit on full training data using the same lambda
    fit_full <- glmnet(
      x_train, y_train,
      alpha  = 1,
      lambda = lambda_opt
    )
    
    # 8. Predict on test data using full-data model
    test_pred_mat[, b] <- as.numeric(predict(fit_full, x_test))
    
    # 9. Record loss metrics
    loss_table$Iteration[b] <- b
    loss_table$Lambda[b]    <- lambda_opt
    loss_table$CV_Error[b]  <- min(cv_fit$cvm, na.rm = TRUE)
    loss_table$R2[b]        <- r2_boot
  }
  
  # Ensemble predictions (mean over bootstrap models)
  train_pred_mean <- rowMeans(train_pred_mat, na.rm = TRUE)
  test_pred_mean  <- rowMeans(test_pred_mat,  na.rm = TRUE)
  
  list(
    train_pred     = train_pred_mean,
    test_pred      = test_pred_mean,
    lambda_vec     = lambda_vec,
    loss_table     = loss_table,
    coefficients   = coef_list,
    train_pred_mat = train_pred_mat,
    test_pred_mat  = test_pred_mat
  )
}

### ----------------------------------------
### 3. Conventional proteomic age model (all proteins)
### ----------------------------------------

# Convert demo proteomic data to matrices
x_train_full <- as.matrix(train_demo)   # 98 × 1254
x_test_full  <- as.matrix(test_demo)    # 1850 × 1254
y_train      <- as.numeric(train_age_demo)

# Fit conventional age model and predict on test set
conv_res <- lambda_ensemble_lasso(
  x_train = x_train_full,
  y_train = y_train,
  x_test  = x_test_full,
  n_bootstrap = 500,
  seed = 2025
)

# Predicted conventional proteomic age on test set
conventional_age_demo <- conv_res$test_pred

# Optional checks
# summary(conventional_age_demo)
# head(conventional_age_demo)
# head(conv_res$loss_table)

### ----------------------------------------
### 4. Organ-specific age model (random 30 proteins)
### ----------------------------------------

# Randomly select 30 proteins as "organ-specific proteins" for this demo
set.seed(999)  # to make the subset reproducible
n_features <- ncol(x_train_full)

organ_feature_idx <- sample(seq_len(n_features), size = 30, replace = FALSE)
organ_specific_proteins <- colnames(x_train_full)[organ_feature_idx]

# Subset training and test data to these 30 proteins
x_train_organ <- x_train_full[, organ_feature_idx, drop = FALSE]
x_test_organ  <- x_test_full[,  organ_feature_idx, drop = FALSE]

# Fit organ-specific age model and predict on test set
organ_res <- lambda_ensemble_lasso(
  x_train = x_train_organ,
  y_train = y_train,
  x_test  = x_test_organ,
  n_bootstrap = 500,
  seed = 2026
)

# Predicted organ-specific proteomic age on test set
organ_specific_age_demo <- organ_res$test_pred

# Optional checks
# summary(organ_specific_age_demo)
# head(organ_specific_age_demo)
# organ_specific_proteins   # the 30 proteins used in this demo
# head(organ_res$loss_table)

#result:
#organ_specific_age_demo
#conventional_age_demo

###########################################
## End of demo script
###########################################


