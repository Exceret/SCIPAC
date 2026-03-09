#' The core function to obtain Lambda using different regression methods
#'
#' @param bulk.dat a dimension-reduced bulk data. Rows stand for samples and columns stand for PCs.
#' @param y class labels or survival time.
#' If family = "binomial", "cumulative", or "gaussian", y represents class labels.
#' If family = "cox", y is a data frame with two columns. Row names of y are the same as bulk.dat's and the two columns are named ("time", "status").
#' This renaming is required for the algorithm used.
#' @param family this argument takes one of the following input:
#' (1). "binomial": logistic regression with elastic net;
#' (2). "cumulative": reversed proportional log odds ratio model with elastic net;
#' (3). "gaussian": linear regression with elastic net; and
#' (4). "cox": cox regression with elastic net.
#' @param K.means.res a list contains (1). k, the number of clusters;
#' (2). ct.assignment, a data frame with one column indicating cluster assignment. Row names of ct.assignment are cell names
#' (3). centers, cluster centroids. Rows are for PCs and columns are for clusters.
#' Each cluster centroid is calculated by taking the average value of all the cells in the cluster.
#' @param ela.net.alpha the parameter alpha used for the elastic net. 0 for ridge and 1 for lasso. The default is 0.4.
#' @param nfold The number of folds used in cross-validation for regression models. The default is \code{nfold = 10}.
#' @export
#' @return a data frame with two columns. One column named ct.assignment is the cluster assignment, and the other one named
#' Lambda is the calculated Lambda for each cell. Cells from the same cluster have the same value of Lambda.

classifier.Lambda.core <- function(
  bulk.dat,
  y,
  family,
  K.means.res,
  ela.net.alpha,
  nfold
) {
  K.means.cent <- K.means.res$centers
  K <- K.means.res$k
  ct.assign <- K.means.res$ct.assignment

  if ((family == "binomial") | (family == "cumulative")) {
    if (is.numeric(y)) {
      y <- factor(y)
    } else if (is.factor(y)) {
      y <- as.numeric(y)
      y <- factor(y, levels = c(1:length(unique(y))))
    }

    class.lab <- unique(y)
    n.class <- length(class.lab)

    new.sample <- c()
    new.y <- c()
    for (i in 1:n.class) {
      lab.i <- class.lab[i]
      lab.idx <- which(y == lab.i)
      lab.row.names <- rownames(bulk.dat[lab.idx, ])
      lab.sample <- sample(lab.row.names, length(lab.idx), replace = TRUE)
      lab.dat <- bulk.dat[lab.sample, ]

      y.bt <- rep(lab.i, length(lab.idx))

      new.sample <- rbind(new.sample, lab.dat)
      new.y <- append(new.y, y.bt)
    }
    new.sample.ave <- colSums(new.sample) / nrow(new.sample)

    if (family == "binomial") {
      # Apply logistic regression
      cv.ela.net <- glmnet::cv.glmnet(
        as.matrix(new.sample),
        new.y,
        family = "binomial",
        alpha = ela.net.alpha,
        nfolds = nfold
      )
      logit.mol.ela.net <- glmnet::glmnet(
        new.sample,
        new.y,
        family = "binomial",
        alpha = ela.net.alpha,
        lambda = cv.ela.net$lambda.min,
        standardize = TRUE
      )
      beta <- c(stats::coef(logit.mol.ela.net)[-1])
    } else if (family == "cumulative") {
      # Apply penalized ordinal regression
      pen.logit.ela.net <- ordinalNet::ordinalNet(
        as.matrix(new.sample),
        new.y,
        family = "cumulative",
        link = "logit",
        alpha = ela.net.alpha,
        parallelTerms = TRUE,
        nonparallelTerms = FALSE,
        standardize = TRUE,
        reverse = TRUE
      )
      beta <- stats::coef(pen.logit.ela.net)[-c(1:(n.class - 1))]
    }

    # Calculate Lambda
    Lambda <- sapply(c(1:K), function(k) {
      chosen.cen <- K.means.cent[, k]
      x <- chosen.cen - new.sample.ave
      logit <- t(beta) %*% x
    })
    ct.idx <- ct.assign$cluster_assignment
    ct.assign$Lambda <- Lambda[ct.idx]
    # Do the standardization with consistent handling
    if (length(ct.assign$Lambda) == 0) {
      stop("Lambda calculation failed - no values produced")
    }

    if (stats::sd(ct.assign$Lambda) == 0) {
      warning(
        "Lambda values have zero variance in binomial/cumulative regression"
      )
      ct.assign$Lambda <- rep(0, length(ct.assign$Lambda))
    } else {
      ct.assign$Lambda <- scale(ct.assign$Lambda, center = TRUE, scale = TRUE)
      # Ensure it's a vector, not a matrix
      ct.assign$Lambda <- as.numeric(ct.assign$Lambda)
    }
    return(ct.assign)
  } else if (family == "gaussian") {
    n.row <- nrow(bulk.dat)
    re.sample <- sample(1:n.row, n.row, replace = TRUE)
    new.sample <- bulk.dat[re.sample, ]

    new.y <- y[re.sample]

    new.sample.ave <- colSums(new.sample) / nrow(new.sample)

    # Apply linear regression with elastic net
    cv.ela.net <- glmnet::cv.glmnet(
      as.matrix(new.sample),
      new.y,
      family = "gaussian",
      alpha = ela.net.alpha,
      nfolds = nfold
    )
    linear.mol.ela.net <- glmnet::glmnet(
      new.sample,
      new.y,
      family = "gaussian",
      alpha = ela.net.alpha,
      lambda = cv.ela.net$lambda.min,
      standardize = TRUE
    )
    beta <- c(stats::coef(linear.mol.ela.net)[-1])

    # Calculate Lambda
    Lambda <- sapply(c(1:K), function(k) {
      chosen.cen <- K.means.cent[, k]
      x <- chosen.cen - new.sample.ave
      logit <- t(beta) %*% x
    })
    ct.idx <- ct.assign$cluster_assignment
    ct.assign$Lambda <- Lambda[ct.idx]
    # Do the standardization with consistent handling
    if (length(ct.assign$Lambda) == 0) {
      stop("Lambda calculation failed - no values produced")
    }

    if (stats::sd(ct.assign$Lambda) == 0) {
      warning(
        "Lambda values have zero variance in binomial/cumulative regression"
      )
      ct.assign$Lambda <- rep(0, length(ct.assign$Lambda))
    } else {
      ct.assign$Lambda <- scale(ct.assign$Lambda, center = TRUE, scale = TRUE)
      # Ensure it's a vector, not a matrix
      ct.assign$Lambda <- as.numeric(ct.assign$Lambda)
    }
    return(ct.assign)
  } else if (family == "cox") {
    # Re-sampling the whole survival data
    all.rownames <- rownames(bulk.dat)
    re.sample <- sample(all.rownames, length(all.rownames), replace = TRUE)
    re.sample.dat <- bulk.dat[re.sample, ]

    # Handle matrix format for y
    if (is.matrix(y)) {
      y <- as.data.frame(y)
    }

    # Use match to ensure proper alignment of row names
    matched_indices <- match(re.sample, rownames(y))
    new.y <- y[matched_indices, ]
    new.sample.ave <- colSums(re.sample.dat) / nrow(re.sample.dat)

    # Adjust nfolds based on number of events to prevent CV errors
    n_events <- sum(new.y$status == 1)
    adjusted_nfold <- min(nfold, n_events)

    # Create Surv object for cox regression
    surv_obj <- survival::Surv(time = new.y$time, event = new.y$status)

    # Apply cox regression with adjusted nfolds
    cv.ela.net <- glmnet::cv.glmnet(
      re.sample.dat,
      surv_obj,
      family = "cox",
      type.measure = "C",
      alpha = ela.net.alpha,
      nfolds = adjusted_nfold
    )
    cox.ela.net <- glmnet::glmnet(
      re.sample.dat,
      surv_obj,
      family = "cox",
      alpha = ela.net.alpha,
      lambda = cv.ela.net$lambda.min,
      standardize = TRUE
    )
    beta <- as.matrix(stats::coef(cox.ela.net))

    # Calculate Lambda
    Lambda <- sapply(c(1:K), function(k) {
      chosen.cen <- K.means.cent[, k]
      x <- chosen.cen - new.sample.ave
      logit <- t(beta) %*% x
    })
    ct.idx <- ct.assign$cluster_assignment
    ct.assign$Lambda <- Lambda[ct.idx]
    # Do the standardization
    ct.assign$Lambda <- scale(ct.assign$Lambda, center = TRUE, scale = TRUE)
    return(ct.assign)
  } else {
    stop(
      "Please choose a valid family argument from 'binomial', 'cumulative', 'gaussian' or 'cox'"
    )
  }
}

#' Apply parallel computing to calculate classifier.Lambda.core
#'
#' @param bulk.dat a dimension-reduced bulk data. Rows stand for samples and columns stand for PCs.
#' @param y class labels or survival time.
#' If family = "binomial", "cumulative", or "gaussian", y represents class labels.
#' If family = "cox", y is a data frame with two columns. Row names of y are the same as bulk.dat's and the two columns are named ("time", "status").
#' This renaming is required for the algorithm used.
#' @param family this argument takes one of the following input:
#' (1). "binomial": logistic regression with elastic net;
#' (2). "cumulative": reversed proportional log odds ratio model with elastic net;
#' (3). "gaussian": linear regression with elastic net; and
#' (4). "cox": cox regression with elastic net.
#' @param K.means.res a list contains (1). k, the number of clusters;
#' (2). ct.assignment, a data frame with one column indicating cluster assignment. Row names of ct.assignment are cell names
#' (3). centers, cluster centroids. Rows are for PCs and columns are for clusters.
#' Each cluster centroid is calculated by taking the average value of all the cells in the cluster.
#' @param ela.net.alpha the parameter alpha used for the elastic net. 0 for ridge and 1 for lasso. The default is 0.4.
#' @param bt.size the number of bootstrap samples. The default is 50.
#' @param nfold The number of folds used in cross-validation for regression models. The default is \code{nfold = 10}.
#' @param numCores the number of cores used for parallel computing. The default is 7.
#' @param debug logical. Whether to print debug information. The default is FALSE.
#' @export
#' @return A data frame whose rows are cells and columns are bootstrap samples.

classifier.Lambda <- function(
  bulk.dat,
  y,
  family,
  K.means.res,
  ela.net.alpha,
  bt.size,
  nfold = nfold,
  numCores,
  debug = FALSE
) {
  fx <- function(seed) {
    set.seed(seed)
    tryCatch(
      {
        if (debug && seed == 1) {
          message("Debug: Running first bootstrap sample...")
          message(
            "  Bulk data dimensions: ",
            nrow(bulk.dat),
            " x ",
            ncol(bulk.dat)
          )
          message("  Y dimensions: ", nrow(y), " x ", ncol(y))
          message("  Number of clusters: ", K.means.res$k)
        }
        result <- classifier.Lambda.core(
          bulk.dat,
          y,
          family,
          K.means.res,
          ela.net.alpha = ela.net.alpha,
          nfold = nfold
        )
        if (debug && seed == 1) {
          message("  Lambda calculation successful")
          message("  Result structure: ", paste(names(result), collapse = ", "))
          if (!is.null(result$Lambda)) {
            message("  Lambda length: ", length(result$Lambda))
          }
        }
        return(result)
      },
      error = function(e) {
        if (debug) {
          message("Error in bootstrap sample ", seed, ": ", e$message)
        }
        return(list(error = TRUE, message = as.character(e), seed = seed))
      }
    )
  }
  K <- K.means.res$k
  ct.assign <- K.means.res$ct.assignment

  seed.ls <- c(1:bt.size)

  if (debug || numCores == 1) {
    # Use lapply for better error reporting in debug mode
    Lambda.tab <- lapply(seed.ls, fx)
  } else {
    Lambda.tab <- parallel::mclapply(seed.ls, fx, mc.cores = numCores)
  }

  # Check for errors in parallel execution
  errors <- sapply(Lambda.tab, function(x) {
    inherits(x, "try-error") || (!is.null(x$error) && x$error)
  })

  if (all(errors)) {
    # If all bootstrap samples failed, check the first error message
    first_error <- Lambda.tab[[1]]
    if (!is.null(first_error$message)) {
      stop("All bootstrap samples failed. First error: ", first_error$message)
    } else {
      stop("All bootstrap samples failed in parallel execution")
    }
  }

  # Get valid results only
  valid_idx <- which(!errors)
  if (length(valid_idx) < bt.size * 0.5) {
    warning(sprintf(
      "Only %d out of %d bootstrap samples succeeded. Results may be unreliable.",
      length(valid_idx),
      bt.size
    ))
  }

  Lambda.res <- matrix(NA, nrow = nrow(ct.assign), ncol = length(valid_idx))
  for (i in seq_along(valid_idx)) {
    idx <- valid_idx[i]
    if (!is.null(Lambda.tab[[idx]]$Lambda)) {
      Lambda.res[, i] <- Lambda.tab[[idx]]$Lambda
    }
  }

  # Delete columns with NAs.
  if (sum(is.na(Lambda.res[1, ])) == 0) {
    Lambda.res <- Lambda.res
  } else {
    na.idx <- which(is.na(Lambda.res[1, ]))
    Lambda.res <- Lambda.res[, -na.idx]
  }

  if (ncol(Lambda.res) == 0) {
    stop("No valid bootstrap samples were produced")
  }

  return(Lambda.res)
}


#' Summarize Lambda values
#'
#' @param Lambda.res Lambda values for each cell with multiple bootstrap samples.
#' @param K.means.res a list contains (1). k, the number of clusters;
#' (2). ct.assignment, a data frame with one column indicating cluster assignment. Row names of ct.assignment are cell names
#' (3). centers, cluster centroids. Rows are for PCs and columns are for clusters.
#' Each cluster centroid is calculated by taking the average value of all the cells in the cluster.
#' @param CI.alpha significance level used to decide significantly positive/negative results. The default is 0.05.
#' @return A data frame whose rows are cells and columns are bootstrap samples.
#' @export
#' @importFrom stats var

obtain.ct.Lambda <- function(Lambda.res, K.means.res, CI.alpha = 0.05) {
  K <- K.means.res$k
  prob <- 1 - CI.alpha / 2

  Lambda.est <- rowSums(Lambda.res) / ncol(Lambda.res)
  Lambda.upper <- Lambda.est +
    stats::qnorm(prob) * sqrt(apply(Lambda.res, 1, var))
  Lambda.lower <- Lambda.est -
    stats::qnorm(prob) * sqrt(apply(Lambda.res, 1, var))

  Lambda.std <- sqrt(apply(Lambda.res, 1, var))
  Lambda.z <- apply(Lambda.res, 1, mean) / Lambda.std
  Lambda.z.sign <- sign(Lambda.z)
  Lambda.pval.nlog <- -log10(
    2 * stats::pnorm(q = abs(Lambda.z), lower.tail = FALSE)
  )
  Lambda.pval <- Lambda.z.sign * Lambda.pval.nlog

  Lambda.sig <- rep(NA, nrow(Lambda.res))

  sign.res <- sign(Lambda.upper) + sign(Lambda.lower)
  Lambda.sig[sign.res == 0] <- "Neutral"
  Lambda.sig[sign.res > 0] <- "Positive"
  Lambda.sig[sign.res < 0] <- "Negative"

  ct.assign <- K.means.res$ct.assignment
  ct.assign$SCIPAC_Lambda_est <- Lambda.est
  ct.assign$SCIPAC_Lambda_upper <- Lambda.upper
  ct.assign$SCIPAC_Lambda_lower <- Lambda.lower
  ct.assign$SCIPAC <- Lambda.sig
  ct.assign$SCIPAC_log_pval <- Lambda.pval

  ct.assign$SCIPAC <- factor(
    ct.assign$SCIPAC,
    levels = c("Positive", "Negative", "Neutral")
  )

  return(ct.assign)
}

#' Obtain phenotype-associated cells using SCIPAC
#'
#' @param bulk.dat a dimension-reduced bulk data. Rows stand for samples and columns stand for PCs.
#' @param y class labels or survival time.
#' \itemize{
#' \item If \code{family = "binomial"} or \code{"cumulative"}, \code{y} represents class labels.
#' \item If \code{family = "gaussian"}, \code{y} represents a continuous variable.
#' \item If \code{family = "cox"}, \code{y} is a data frame with two columns. Row names of \code{y} are the same as \code{bulk.dat}'s and the two columns are named \code{("time", "status")}.
#' This renaming is required for the algorithm used.
#' }
#' @param family this argument takes one of the following input:
#' \itemize{
#' \item \code{family = "binomial"}: logistic regression with elastic net;
#' \item \code{family = "cumulative"}: reversed proportional log odds ratio model with elastic net;
#' \item \code{family = "gaussian"}: linear regression with elastic net; and
#' \item \code{family = "cox"}: cox regression with elastic net.
#' }
#' @param ct.res a list contains
#' \itemize{
#' \item \code{k}, the number of clusters;
#' \item \code{ct.assignment}, a data frame with one column indicating cluster assignment. Row names of \code{ct.assignment} are cell names
#' \item \code{centers}, cluster centroids. A data frame whose rows are for PCs and columns are for clusters. Each cluster centroid is calculated by taking the average value of all the cells in the cluster.
#' }
#' @param ela.net.alpha the parameter alpha used for the elastic net. \code{0} for ridge and \code{1} for lasso. The default is \code{ela.net.alpha = 0.4}.
#' @param bt.size the number of bootstrap samples. The default is \code{bt.size = 50}.
#' @param numCores the number of cores used for parallel computing. The default is \code{numCores = 7}.
#' @param CI.alpha significance level used to decide significantly positive/negative results. The default is \code{CI.alpha = 0.05}.
#' @param nfold The number of folds used in cross-validation for regression models. The default is \code{nfold = 10}.
#' @param debug If TRUE, runs in debug mode with single core and verbose output. Default is FALSE.
#' @return A data frame with six columns. Row names are cells. The six columns are
#' \itemize{
#' \item (1). \code{cluster_assignment}: the cluster assignment of each cell;
#' \item (2). \code{Lambda.est}: the estimated Lambda;
#' \item (3). \code{Lambda.upper}: the upper value of the confidence interval;
#' \item (4). \code{Lambda.lower}: the lower value of the confidence interval;
#' \item (5). \code{sig}: significance identification, includes "Positive", "Negative", and "Neutral";
#' \item (6). \code{log.pval}: log10 p-values
#' }
#' @export

SCIPAC <- function(
  bulk.dat,
  y,
  family,
  ct.res,
  ela.net.alpha = 0.4,
  bt.size = 50,
  numCores = 7,
  CI.alpha = 0.05,
  nfold = 10,
  debug = FALSE
) {
  if (debug) {
    message("Running SCIPAC in debug mode with single core...")
    numCores <- 1
  }

  Lambda.res <- classifier.Lambda(
    bulk.dat,
    y,
    family,
    ct.res,
    ela.net.alpha = ela.net.alpha,
    bt.size = bt.size,
    nfold = nfold,
    numCores = numCores,
    debug = debug
  )
  ct.assign <- obtain.ct.Lambda(Lambda.res, ct.res, CI.alpha = CI.alpha)
  return(ct.assign)
}
