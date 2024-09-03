#' Likelihood ratio test for testing equality means against umbrella-ordered alternatives with unknown peak in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @return Peak numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_1 = mu_2 = ... = mu_k vs H_1:mu_1 <=.....<= mu_(h-1)<= mu_h >= mu_(h+1)>=....>= mu_k (at least one strict inequality, with unknown peak), where mu_i represents the population means of the i-th treatment. The input consists of two variables: sample_data and significance_level. The output consists of the Peak of the umbrella, critical value, the UniLRT* test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @importFrom Iso pava
#' @importFrom stats qchisq
#' @importFrom stats quantile rnorm var
#' @author Subha Halder
UniLRTstar <- function(sample_data, significance_level){
  set.seed(456)
  R_MLE <- function(X, n) {
    X1 <- X[-1]
    n1 <- n[-1]
    sorted_indices <- order(X1)
    X1_sorted <- X1[sorted_indices]
    n1_sorted <- n1[sorted_indices]
    A <- numeric(length(X1_sorted))
    for (j in 2:length(X)) {
      A[j-1] <- (n[1] * X[1] + sum(n1_sorted[1:(j - 1)] * X1_sorted[1:(j - 1)])) /
        (n[1] + sum(n1_sorted[1:(j - 1)]))
    }
    if (all(X1 >= X[1])) {
      new_X <- X
    } else if (A[length(X)-2] >= X1_sorted[length(X)-1]) {
      X <- rep(A[length(X)-1], length(X))
      new_X <- X
    } else {
      comparisons <- logical(length(X1_sorted) - 1)
      comparisons1 <- logical(length(X1_sorted) - 1)
      stored_values <- numeric(0)
      for (k in 1:(length(X1_sorted) - 1)) {
        comparisons[k] <- A[k] < X1_sorted[k + 1]
        if(comparisons1[k] <- A[k] < X1_sorted[k + 1]) {
          for (s in 1:k) {
            stored_values[s] <- X1_sorted[s]
          }
          break
        }
      }
      selected_A_values <- A[comparisons]
      X[1] <- selected_A_values[1]
      for (l in 2:length(X)) {
        if (X[l] %in% stored_values) {
          X[l] <- selected_A_values[1]
        }
      }
      new_X <- X
    }
    return(new_X)
  }

  unimodal <- function(X, w, lmode) {
    if(lmode==length(X)){
      new_X <- Iso::pava(X,w)
    } else if (lmode==1){
      new_X <- -(Iso::pava(-X,w))
    } else{
      X1 <- X[1:(lmode-1)]
      X2 <- X[(lmode+1):length(X)]
      w1 <- w[1:(lmode-1)]
      w2 <- w[(lmode+1):length(X)]
      newX1 <- Iso::pava(X1,w1)
      newX2 <- -(Iso::pava(-X2,w2))
      Y <- c(X[lmode],newX1,newX2)
      w_new<- c(w[lmode],w1,w2)
      v <- -(R_MLE(-Y,w_new))
      max_v <- v[-1]
      new_X1 <- max_v[1:length(X1)]
      new_X2 <- max_v[(length(X1)+1):length(max_v)]
      new_X <- c(new_X1,v[1],new_X2)
    }
    return(new_X)
  }
  unimod_peak <- function(sample_data_list) {
    n <- sapply(sample_data_list, length)
    mu0 <- sapply(sample_data_list, mean)
    var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - mu0[i])^2) / n[i])
    w0 <- n / var0
    Uni_star <- list()
    res <- numeric()
    for (p in 1:length(mu0)) {
      repeat {
        new_mu0 <- unimodal(X = sapply(sample_data_list, mean), w = w0, lmode = p)
        new_var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - new_mu0[i])^2) / n[i])
        new_w0 <- n / new_var0
        if (max(abs(new_mu0 - mu0)) <= 0.000001) {
          break  # Exit the loop if the difference is less than epsilon
        }
        w0 <- new_w0
        mu0 <- new_mu0
        var0 <- new_var0
      }
      Uni_star[[p]] <- new_mu0
      res[p] <- sum(new_w0*(sapply(sample_data_list, mean)-Uni_star[[p]])^2)
    }
    return(which.min(res))
  }
  pk <-  unimod_peak(sample_data)
  LRT_H0_new <- function(sample_data_list) {
    means <- sapply(sample_data_list, mean)
    sample_sizes <- sapply(sample_data_list, length)
    S <- unlist(sample_data_list)
    mu1 <- mean(S)
    var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - means[i])^2)) / sample_sizes[i])
    u1 <- sample_sizes / var1

    repeat {
      new_mu1 <- (sum(u1 * means)) / sum(u1)
      new_var1 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu1)^2)) / sample_sizes[i])
      new_u1 <- sample_sizes / new_var1

      if (max(abs(new_mu1 - mu1)) <= 0.0000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      u1 <- new_u1
      mu1 <- new_mu1
      var1 <- new_var1
    }

    return(var1)
  }

  LRT_H1_new <- function(sample_data_list,h) {
    n <- sapply(sample_data_list, length)
    mu0 <- sapply(sample_data_list, mean)
    var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - mu0[i])^2)) / n[i])
    w0 <- n / var0
    repeat {
      new_mu0 <- unimodal(X = sapply(sample_data_list, mean), w = w0, lmode = h)
      new_var0 <- sapply(1:length(sample_data_list), function(i) (sum((sample_data_list[[i]] - new_mu0[i])^2)) / n[i])
      new_w0 <- n / new_var0

      if (max(abs(new_mu0 - mu0)) <= 0.0000001) {
        break  # Exit the loop if the difference is less than epsilon
      }

      w0 <- new_w0
      mu0 <- new_mu0
      var0 <- new_var0
    }

    return(var0)
  }
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 20000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  lambda_values_star <- numeric(num_samples)
  for (i in 1:num_samples) {
    bootstrap_samples <- lapply(sample_data, function(x) rnorm(n = length(x), mean = 0, sd = sqrt(var(x))))
    V_R_star <- LRT_H1_new(bootstrap_samples,h=pk) / LRT_H0_new(bootstrap_samples)
    weights <- sapply(1:num_datasets, function(i) V_R_star[i]^(length(bootstrap_samples[[i]]) / 2))
    lambda_values_star[i] <- prod(weights)
  }
  sort_lambda_star <- sort(lambda_values_star)
  quantile_value <- quantile(sort_lambda_star, probs = significance_level)
  V_R <- LRT_H1_new(sample_data, h=pk) / LRT_H0_new(sample_data)
  weights <- sapply(1:num_datasets, function(i) V_R[i]^(length(sample_data[[i]]) / 2))
  lambda <- prod(weights)
  if (lambda < quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste( "Peak:", pk, "; Critical value:", quantile_value, "; UniLRT Test statistic:", lambda, "; Result:", result))
}






