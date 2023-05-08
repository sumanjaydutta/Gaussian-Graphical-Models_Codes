#synthetic experiment for seeing source of errors for Shrnikage compared with GGMs
#we compare four methods here, LW, RBLW, GReedyPrune and HybridMB
library(rlist)
library(glasso)
library(MASS)
library(doParallel)
library(Matrix)
library(slam)
library(glmnet)
library(gurobi)
library(viridis)
library(ggplot2)
library(igraph)
library(nlshrink)
library(CovTools)
library(corpcor)

cl = makeCluster(4) # use 4 processes at once
registerDoParallel(cl)
# L1 norm of matrices as vectors
exp1_norm = function(prec1, prec2) {
  return(sum(abs(prec1 - prec2)))
}

# Unused.
# given the covariance of an attractive model,
# create a walk-summable model by putting random signs
random_sign = function(input_cov) {
  prec = solve(input_cov)
  dims = dim(prec)
  signs = matrix(2 * rbinom(dims[1] * dims[2], 1, 0.5) - 1,
                 nrow = dims[1],
                 ncol = dims[2])
  signs = signs * t(signs) # ensures symmetry and positive diagonal
  prec.new = prec * signs
  return(solve(prec.new))
}

cov_bm = function(n) {
  initial_variance = 0.5
  cov = matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      cov[i, j] = initial_variance + min(i, j) / n
    }
  }
  return(cov)
}

# path laplacian plus a tiny "mass" at 0 to make it invertible
cov_path_free = function(m) {
  adj = matrix(0,nrow=m,ncol=m)
  for (i in 1:(m - 1)) {
    adj[i, i + 1] = 1
    adj[i + 1, i] = 1
  }
  laplacian = diag(colSums(adj)) - adj
  laplacian[1,1] = laplacian[1,1] + 1/m
  cov = solve(laplacian)
  return(standardize(as.matrix(cov)))
}

# Covariance of a path and some cliques (more complicated version)
cov_path_2cliques = function(m, n1, n2, d) {
  cov1 = standardize(cov_bm(m)) # CHANGED
  prec2_basic = diag(d) + matrix(-(0.99 / d), nrow = d, ncol = d)
  prec3_basic = diag(d) + matrix(-(0.5 / d), nrow = d, ncol = d)
  cov2_basic = standardize(solve(prec2_basic))
  cov3_basic = standardize(solve(prec3_basic))
  cov2 = as.matrix(.bdiag(replicate(n1 / d, cov2_basic, simplify = FALSE)))
  cov3 = as.matrix(.bdiag(replicate(n2 / d, cov3_basic, simplify = FALSE)))
  return(as.matrix(bdiag(cov1, cov2, cov3)))
}

# Covariance of a path and some cliques
cov_path_and_cliques = function(m, n, d, rho = 0.95) {
  cov1 = standardize(cov_bm(m)) 
  prec_basic = diag(d) + matrix(-(rho / d), nrow = d, ncol = d)
  cov_basic = standardize(solve(prec_basic)) 
  cov2 = as.matrix(.bdiag(replicate(n / d, cov_basic, simplify = FALSE)))
  return(as.matrix(bdiag(cov1, cov2)))
}


# path of length m and sqrt(n) paths of length sqrt(n)
cov_path_and_paths = function(m, n) {
  cov1 = cov_bm(m)
  l = sqrt(n)
  cov_basic = cov_bm(l)
  cov2 = as.matrix(.bdiag(replicate(n / l, cov_basic, simplify = FALSE)))
  return(as.matrix(bdiag(cov1, cov2)))
}

# renormalize to unit variance. Without this rescaling,
# many of the older algorithms will not work well.
standardize = function(M) {
  #renormalize to diagonal 1
  diags = diag(M)
  M = t(M / sqrt(diags)) / sqrt(diags)
  return(M)
}

# make a rect grid
cov_rect_grid = function(a, b) {
  adj = matrix(0, a * b, a * b)
  for (i in 1:a) {
    for (j in 1:b) {
      row = (i - 1) * b + j
      if (i < a) {
        adj[row, row + b] = 1
      }
      if (j < b) {
        adj[row, row + 1] = 1
      }
    }
  }
  adj = adj + t(adj)
  laplacian = diag(apply(adj, 1, sum)) - adj
  prec = laplacian[1:(a * b - 1), 1:(a * b - 1)]
  cov = solve(prec)
  return(standardize(cov))
}


# Some examples from the "hard examples" appendix.
cov_copy_example = function(n, delta) {
  v_basic = rbind(
    c(3 / 4, -1 / 4, -1 / 4, -1 / 4, delta, 0, 0, 0, 0, 0, 0, 0),
    c(3 / 4, -1 / 4, -1 / 4, -1 / 4, 0, delta, 0, 0, 0, 0, 0, 0),
    c(-1 / 4, 3 / 4, -1 / 4, -1 / 4, 0, 0, delta, 0, 0, 0, 0, 0),
    c(-1 / 4, 3 / 4, -1 / 4, -1 / 4, 0, 0, 0, delta, 0, 0, 0, 0),
    c(-1 / 4, -1 / 4, 3 / 4, -1 / 4, 0, 0, 0, 0, delta, 0, 0, 0),
    c(-1 / 4, -1 / 4, 3 / 4, -1 / 4, 0, 0, 0, 0, 0, delta, 0, 0),
    c(-1 / 4, -1 / 4, -1 / 4, 3 / 4, 0, 0, 0, 0, 0, 0, delta, 0),
    c(-1 / 4, -1 / 4, -1 / 4, 3 / 4, 0, 0, 0, 0, 0, 0, 0, delta)
  )
  cov_basic = v_basic %*% t(v_basic)
  cov_full = as.matrix(.bdiag(replicate(n / 8, cov_basic, simplify = FALSE)))
  return(cov_full)
}

cov_weird_copy = function(n, delta) {
  v_basic = rbind(
    c(3 / 4, -1 / 4, -1 / 4, -1 / 4, delta, 0, 0, 0, 0, 0, 0, 0),
    c(3 / 4, -1 / 4, -1 / 4, -1 / 4, 0, delta, 0, 0, 0, 0, 0, 0),
    c(-1 / 4, 3 / 4, -1 / 4, -1 / 4, 0, 0, delta, 0, 0, 0, 0, 0),
    c(-1 / 4, 3 / 4, -1 / 4, -1 / 4, 0, 0, 0, delta, 0, 0, 0, 0),
    c(-1 / 4, -1 / 4, 3 / 4, -1 / 4, 0, 0, 0, 0, delta, 0, 0, 0),
    c(-1 / 4, -1 / 4, 3 / 4, -1 / 4, 0, 0, 0, 0, 0, delta, 0, 0),
    c(-1 / 4, -1 / 4, -1 / 4, 3 / 4, 0, 0, 0, 0, 0, 0, delta, 0),
    c(-1 / 4, -1 / 4, -1 / 4, 3 / 4, 0, 0, 0, 0, 0, 0, 0, delta)
  )
  cov_basic = v_basic %*% t(v_basic)
  prec1 = solve(cov_basic)
  prec2 = prec1[-1, -1]
  cov2 = solve(prec2)
  
  # rescale to constant-order variances
  cov3 = cov2
  v = cov3[1, 1]
  cov3[1, ] = cov3[1, ] / sqrt(v)
  cov3[, 1] = cov3[, 1] / sqrt(v)
  
  cov_full =  as.matrix(.bdiag(replicate(n / 7, cov3, simplify = FALSE)))
  return(cov_full)
}


# We could use glassopath, however we found glasso was faster in our setup
# and more reliable. Hence the current implementation.
best_glasso <-
  function(emp_cov, true_prec, rho_list, kappa_eval = NULL) {
    # corresponds to result with infinite penalty
    best_prec = matrix(0, nrow = nrow(emp_cov), ncol = ncol(emp_cov))
    best_err = Inf
    best_rho = Inf
    glasso_prec_list = foreach(rho = rho_list, .packages = "glasso") %dopar% {
      glasso(emp_cov, rho)[["wi"]]
    }
    for (idx in 1:length(rho_list)) {
      rho = rho_list[idx]
      glasso_prec = glasso_prec_list[[idx]]
      if (EXPERIMENT == 1) {
        err = exp1_norm(glasso_prec, true_prec)#norm(glasso_prec - true_prec,type="f")^2
      } else {
        err = get_adj_err(true_prec, glasso_prec, kappa = kappa_eval)
      }
      if (err < best_err) {
        best_err = err
        best_prec = glasso_prec
        best_rho = rho
      }
    }
    return(list(rho = best_rho, prec = best_prec, err = best_err))
  }

## CHECKED: returns same thing as greedy.old
# run greedy for k iterations from i
# rows of X are samples
greedy = function(X, k, i) {
  # make sure k is valid
  k = min(k, dim(cov) - 1)
  S = vector()
  X.S = X # projected version
  
  # Current implementation is basically Gram-Schmidt.
  for (step in 1:k) {
    dots = (X.S[, i] %*% X.S)[1, ]
    norms = sqrt(colSums(X.S ^ 2))
    scores = abs(dots) / norms
    scores[i] = NA
    scores[S] = NA
    
    if (all(is.na(scores))) {
      break     # no valid pivots remain, can happen due to numerics
    }
    pivot = which.max(scores)
    #print(scores[pivot]/sqrt(dim(X)[1]))
    dots.pivot = (X.S[, pivot] %*% X.S)[1, ]
    X.S = X.S - outer(X.S[, pivot], dots.pivot / (norms[pivot] ^ 2))
    #X.S[,pivot] = NA # uneeded with new code
    S = c(S, pivot)
  }
  return(S)
}

# implementation using NA: apparently NA makes
# R not use optimized routines so it is bad
greedy.old2 = function(X, k, i) {
  # make sure k is valid
  k = min(k, dim(cov) - 1)
  S = vector()
  X.S = X # projected version
  
  # Current implementation is basically Gram-Schmidt. Think if
  # this is the best way for a small # of steps
  for (step in 1:k) {
    dots = (X.S[, i] %*% X.S)[1, ]
    norms = sqrt(colSums(X.S ^ 2))
    scores = abs(dots) / norms
    scores[i] = NA
    
    if (all(is.na(scores))) {
      break     # no valid pivots remain, can happen due to numerics
    }
    pivot = which.max(scores)
    #print(scores[pivot]/sqrt(dim(X)[1]))
    dots.pivot = (X.S[, pivot] %*% X.S)[1, ]
    X.S = X.S - outer(X.S[, pivot], dots.pivot / (norms[pivot] ^ 2))
    X.S[, pivot] = NA
    S = c(S, pivot)
  }
  return(S)
}

# run greedy for k iterations from i
# still can be useful if we have a lot of samples
greedy.old = function(emp_cov, k, i) {
  cov = emp_cov
  
  # make sure k is valid
  k = min(k, dim(cov) - 1)
  S = vector()
  for (step in 1:k) {
    # adding abs reduced warnings from cols of S
    # (which don't matter)
    scores = abs(cov[i, ]) / sqrt(abs(diag(cov)))
    scores[i] = NA
    scores[S] = NA
    
    if (all(is.na(scores))) {
      break     # no valid pivots remain, can happen due to numerics
    }
    pivot = which.max(scores)
    #print(scores[pivot])
    
    #print(pivot)
    # compute the next covariance matrix,
    # covariance matrix of (X_j - Cov[j,pivot]/Cov[pivot,pivot] * X_pivot)
    # = Cov(j,k) - Cov(j, pivot)Cov(k,pivot)/Cov(pivot,pivot) - Cov(pivot,k)Cov(pivot,j)/Cov(pivot,pivot)
    #            + Cov[j,pivot]Cov[pivot,k]/Cov(pivot,pivot])
    cov = cov - outer(cov[, pivot], cov[, pivot]) / cov[pivot, pivot]
    # NA values make matrix operations slower
    #cov[pivot,] = NA
    #cov[,pivot] = NA
    #image(cov)
    #print(cov[,i])
    S = c(S, pivot)
  }
  return(S)
}

prune = function(X, S, i, threshold) {
  # minimize (e_i - w_S e_S)^T Cov (e_i - w_S e_S)
  # = e_i^T Cov e_i - 2 w_S e_S^T Cov e_i + w_S e_S^T Cov(w_S e_S)
  # solution at 0 = Cov e_i - Cov (w_S e_S)
  S_final = vector()
  M = t(X[, S]) %*% X[, S]
  v = t(X[, S]) %*% X[, i]
  
  w = solve(M, v)
  
  err = norm(X[, i] - X[, S] %*% w, type = "f") ^ 2

  for (idx in 1:length(S)) {
    j = S[idx]
    S_noj = S[-idx]
    
    w_noj = solve(M[-idx, -idx], v[-idx])
    
    err_noj = norm(X[, i] - X[, S_noj] %*% w_noj, type = "f") ^ 2
    if (err_noj > (1 + threshold) * err) {
      S_final = append(S_final, j)
    }
  }
  return(S_final)
}

# REIMPLEMENTED THIS TO BE FASTER IT IS MUCH SLOWER THAN GREEDY
prune.old = function(X, S, i, threshold) {
  # minimize (e_i - w_S e_S)^T Cov (e_i - w_S e_S)
  # = e_i^T Cov e_i - 2 w_S e_S^T Cov e_i + w_S e_S^T Cov(w_S e_S)
  # solution at 0 = Cov e_i - Cov (w_S e_S)
  cov = emp_cov
  S_final = vector()
  w = ginv(X[, S]) %*% X[, i]
  #print(w)
  err = norm(X[, i] - X[, S] %*% w, type = "f") ^ 2
  #print(paste("with all: ", err))M
  
  for (idx in 1:length(S)) {
    j = S[idx]
    S_noj = S[-idx]
    w_noj = ginv(X[, S_noj]) %*% X[, i]
    err_noj = norm(X[, i] - X[, S_noj] %*% w_noj, type = "f") ^ 2
    #print(paste("without ",j,":", err_noj))
    if (err_noj > (1 + threshold) * err) {
      S_final = append(S_final, j)
    }
  }
  return(S_final)
}

# Note: doesn't work if k=2 due to list->number conversions in R
greedy_and_prune = function(X, emp_cov, i, k, threshold) {
  # If there are many rows, writing down the whole covariance
  # matrix will be faster (the greedy.old method)
  # they should always return the same result
  S = if (nrow(X) < ncol(X) * 0.7) {
    greedy(X, k, i)
  } else {
    greedy.old(emp_cov, k, i)
  }
  return(prune(X, S, i, threshold))
}

full_greedy_and_prune = function(X, emp_cov, k, threshold) {
  n = ncol(X)
  S = list()
  for (i in 1:n) {
    S[[i]] = greedy_and_prune(X, emp_cov, i, k, threshold)
  }
  # Make the neighborhood sets consistent by intersecting
  for (i in 1:n) {
    #print(S[[i]])
    S_cleaned = vector()
    for (j in S[[i]]) {
      if (i %in% S[[j]]) {
        S_cleaned = c(S_cleaned, j)
      }
    }
    S[[i]] = S_cleaned
    #print(S[[i]])
  }

  prec_ = matrix(0, nrow = n, ncol = n)
  # estimate coefficients using linear regression
  for (i in 1:n) {
    U = S[[i]]
    if (length(U) == 0) {
      prec_[i, i] = 1 / emp_cov[i, i]
    } else {
      w = ginv(X[, U]) %*% X[, i]
      # estimate for conditional variance
      var_i = norm(X[, i] - X[, U] %*% w, type = "f") ^ 2 / nrow(X)
      prec_[i, i] = 1 / var_i
      for (idx in 1:length(U)) {
        j = U[idx]
        w_j = w[idx]
        # E[X_i | X_{sim i}] = \sum_i
        
        # average the estimates both ways, remembering the correct sign.
        #prec_[i,j] = prec_[i,j] - (1/2) * (w_j * prec_[i,i])
        #prec_[j,i] = prec_[j,i] - (1/2) * (w_j * prec_[i,i])
        
        #TODO: properly test to check this is better
        # instead, take the estimate of min-norm
        new = -(w_j * prec_[i, i])
        if (prec_[i, j] == 0 || abs(prec_[i, j]) > abs(new)) {
          prec_[i, j] = new
        }
        prec_[j, i] = prec_[i, j]
      }
    }
  }
  return(prec_)
}

best_greedy <-
  function(X,
           emp_cov,
           true_prec,
           k_list,
           threshold_list,
           kappa_eval = NULL) {
    best_prec = matrix(0, nrow = nrow(emp_cov), ncol = ncol(emp_cov))
    best_err = Inf
    best_threshold = -1
    best_k = -1
    if (FALSE) {
      for (k in k_list) {
        ## ONLY FOR DEBUGGING
        for (threshold in threshold_list) {
          full_greedy_and_prune(X, emp_cov, k, threshold)
        }
      }
      print("DEBUG PHASE OVER")
    }
    greedy_prec_list = foreach(k = k_list) %:% foreach(
      threshold = threshold_list,
      .export = c(
        "full_greedy_and_prune",
        "greedy_and_prune",
        "greedy",
        "greedy.old",
        "prune"
      ),
      .packages = "MASS"
    ) %dopar% {
      full_greedy_and_prune(X, emp_cov, k, threshold)
    }
    for (idx_k in 1:length(greedy_prec_list)) {
      k = k_list[idx_k]
      for (idx in 1:length(threshold_list)) {
        threshold = threshold_list[idx]
        greedy_prec = greedy_prec_list[[idx_k]][[idx]]#full_greedy_and_prune(X,emp_cov,k,threshold)
        if (EXPERIMENT == 1) {
          err = exp1_norm(greedy_prec, true_prec)#norm(greedy_prec - true_prec,type="f")^2
        } else {
          err = get_adj_err(true_prec, greedy_prec, kappa = kappa_eval)
        }
        if (err < best_err) {
          best_err = err
          best_prec = greedy_prec
          best_threshold = threshold
          best_k = k
        }
      }
    }
    #}
    return(list(
      threshold = best_threshold,
      prec = best_prec,
      err = best_err,
      k = best_k
    ))
  }

best_clime_old <- function(X, lambdas, true_prec) {
  best_prec = matrix(0, nrow = nrow(emp_cov), ncol = ncol(emp_cov))
  best_err = norm(true_prec, type = "f") ^ 2
  best_lambda = Inf
  clime_prec_list = foreach(lambda = lambdas, .packages = "clime") %dopar% {
    clime(X, lambda = lambda)$Omega[[1]]
  }
  for (idx in 1:length(lambdas)) {
    lambda = lambdas[idx]
    clime_prec = clime_prec_list[[idx]]
    #glasso_prec = glasso(emp_cov, rho)[["wi"]]
    err = norm(clime_prec - true_prec, type = "f") ^ 2
    if (err < best_err) {
      best_err = err
      best_prec = clime_prec
      best_lambda = lambda
    }
  }
  return(list(
    lambda = best_lambda,
    prec = best_prec,
    err = best_err
  ))
}

best_clime <- function(emp_cov, lambdas, true_prec, kappa_eval = NULL) {
  best_prec = matrix(0, nrow = nrow(emp_cov), ncol = ncol(emp_cov))
  best_err = Inf
  best_lambda = Inf
  # CLIME is now parallelized itself
  #clime_prec_list = foreach(lambda=lambdas,
  #                          .export = c("my.clime","basis_vector"),
  #                          .packages=c("gurobi","Matrix","doParallel")) %dopar% {
  #  my.clime(emp_cov,lambda=lambda)
  #}
  for (idx in 1:length(lambdas)) {
    lambda = lambdas[idx]
    clime_prec = my.clime(emp_cov, lambda)#clime_prec_list[[idx]]
    if (EXPERIMENT == 1) {
      err = exp1_norm(clime_prec, true_prec)#norm(clime_prec - true_prec,type="f")^2
    } else {
      err = get_adj_err(true_prec, clime_prec, kappa = kappa_eval)
    }
    
    if (err < best_err) {
      best_err = err
      best_prec = clime_prec
      best_lambda = lambda
    }
  }
  return(list(
    lambda = best_lambda,
    prec = best_prec,
    err = best_err
  ))
}


sparsity_threshold = function(Theta, kappa) {
  Adj = Theta # gets the right dimensions
  # return the adjacency matrix using the rule
  # edge <->
  # |\Theta_{ij}|/\sqrt{\Theta_{ii} \Theta_{jj}} > kappa
  for (i in 1:dim(Theta)[1]) {
    Adj[i, ] = abs(Theta[i, ]) / sqrt(Theta[i, i] * diag(Theta)) > kappa
    Adj[i, i] = 0
  }
  
  # Some estimators like ACLIME can put 0 on diagonal
  # which necessitates this cleanup step.
  for (i in 1:dim(Theta)[1]) {
    for (j in 1:dim(Theta)[2]) {
      if (is.na(Adj[i, j])) {
        Adj[i, j] = 0
      }
    }
  }
  return(Adj)
}

get_kappa = function(prec, cleaning = 0.0001) {
  kappa = 1
  
  # due to finite precision arithmetic, there may be tiny entries
  # which should be zero, that would screw up our calculation of kappa.
  # Hence we use cleaning as an "epsilon" to safely delete those.
  Theta = prec
  Theta[abs(prec) < cleaning] = 0
  for (i in 1:dim(Theta)[1]) {
    entries = abs(Theta[i, ]) / sqrt(Theta[i, i] * diag(Theta))
    kappa = min(kappa, min(entries[entries > 0]))
  }
  return(kappa)
}

get_adj_err = function(prec, result_prec, kappa = NULL) {
  if (is.null(kappa)) {
    kappa = get_kappa(prec) / 2
  }
  truth = sparsity_threshold(prec, kappa)
  edge_err = sum(abs(sparsity_threshold(result_prec, kappa) - truth))
  return(edge_err)
}


basis_vector = function(i, n) {
  v = rep(0, n)
  v[i] = 1
  return(v)
}

# the symmetrization procedure recommended in the CLIME paper.
clime.symmetrize = function(omega.tilde) {
  p = dim(omega.tilde)[1]
  for (i in 1:p) {
    for (j in 1:p) {
      if (abs(omega.tilde[i, j]) < abs(omega.tilde[j, i])) {
        omega.tilde[j, i] = omega.tilde[i, j]
      }
    }
  }
  return(omega.tilde)
}

# Better implementation: solve small instead of big LPs.
# (like 3x faster in a simple example, matches result of other implementation)
my.clime = function(Sigma, lambda) {
  p = dim(Sigma)[1]
  
  omega.tilde = foreach (
    j = 1:p,
    .combine = cbind,
    .export = c("basis_vector"),
    .packages = c("gurobi")
  ) %dopar% {
    #print(j)
    model = list()
    
    params = 2 * p
    model$modelsense = 'min'
    model$obj = rep(1, params) # corresponds to l1 norm
    model$sense = c(rep(c("<", ">"), p), "=")
    rhs = rep(0, params + 1)
    rhs[j * 2 - 1] = -1
    rhs[j * 2] = -1
    
    for (i in 1:p) {
      # first case is abs is negative, second case is abs is positive
      rhs[i * 2 - 1] = rhs[i * 2 - 1] + lambda
      rhs[i * 2] = rhs[i * 2] - lambda
    }
    
    rhs[params + 1] = 0
    model$rhs = rhs
    
    A = matrix(0, params + 1, params)
    Sigma.extended = cbind(Sigma, Sigma)
    b_vec = c(rep(1, p), rep(0, p))
    c_vec = c(rep(0, p), rep(1, p))
    
    for (i in 1:p) {
      # Note: use ENTRY-wise multiplication below!
      A[2 * i - 1, ] = -Sigma.extended[i, ] * (b_vec - c_vec)
      A[2 * i, ] = -Sigma.extended[i, ] * (b_vec - c_vec)
    }
    A[params + 1, ] = basis_vector(j + p, params) # force c_jj = 0
    model$A = A
    #gurobi_write(model,"model.rlp")
    
    # Method=1 uses only dual simplex which is fastest in our case.
    params = list(Method = 1, OutputFlag = 0)
    result = gurobi(model, params)
    if (result$status == "INFEASIBLE") {
      # Unspecified case: pick a safe default (as we assume a scaling where variance is 1).
      rslt = rep(0, p)
      rslt[j] = 1
      return(rslt)
    }
    else {
      return(result$x[1:p] - result$x[(p + 1):(2 * p)])
    }
  }
  
  ## Recommended symmetrization step from CLIME paper.
  ## STEP 2.5
  # Omega_{ij} = smaller of omega.tilde_{ij} and omega.tilde_{ji} in norm.
  omega.tilde = clime.symmetrize(omega.tilde)

  return(omega.tilde)
}


# For only this method, we use classical notation to match the ACLIME paper:
# n = number of samples, p = number of parameters
aclime = function(Sigma.star, n) {
  p = dim(Sigma.star)[1]
  Sigma.hat = Sigma.star + diag(p) / n
  sigma = diag(Sigma.star)
  delta = 2 # As recommended in the ACLIME paper.
  lambda.n = delta * sqrt(log(p) / n)
  ## STEP 1: for every j, solve
  # omega.hat_{. j} = argmin |b_j|_1
  # s.t. |Sigma.hat b_j - e_j|_{\infty} <= lambda.n * max(sigma_i,sigma_j) * b_jj
  #      b_jj > 0
  
  # rewrite this as
  # omega.hat_{. j} = argmin \sum_i (b_{ij} + c_{ij})
  # s.t. -lambda.n * max(sigma_i,sigma_j) * b_jj <= (Sigma.hat (b_j - c_j) - e_j)_i
  #                                                 <= lambda.n * max(sigma_i,sigma_j) * b_jj
  # c_jj = 0 and b,c >= 0
  
  # rewrite main constraint as (abs. value is negative case)
  # -lambda.n * max(sigma_i,sigma_j) * b_jj - (Sigma.hat (b_j - c_j))_i <= (-e_j)_i
  # and (abs. value is positive case)
  # lambda.n * max(sigma_i,sigma_j) * b_jj - (Sigma.hat (b_j - c_j))_i >= -(e_j)_i
  
  # This is what is passed on to STEP 2
  #omega.breve = rep(0,p)
  omega.breve = foreach (
    j = 1:p,
    .export = c("basis_vector"),
    .combine = c,
    .packages = c("gurobi")
  ) %dopar% {
    model = list()
    
    params = 2 * p
    model$modelsense = 'min'
    model$obj = rep(1, params) # corresponds to l1 norm
    model$sense = c(rep(c("<", ">"), p), "=")
    rhs = rep(0, params + 1)
    rhs[j * 2 - 1] = -1
    rhs[j * 2] = -1
    rhs[params + 1] = 0
    model$rhs = rhs
    
    A = matrix(0, params + 1, params)
    Sigma.hat.extended = cbind(Sigma.hat, Sigma.hat)
    for (i in 1:p) {
      b_vec = c(rep(1, p), rep(0, p))
      c_vec = c(rep(0, p), rep(1, p))
      # Note: use ENTRY-wise multiplication below!
      A[2 * i - 1, ] = (
        -lambda.n * max(sigma[i], sigma[j]) * basis_vector(j, params)
        - Sigma.hat.extended[i, ] * (b_vec - c_vec)
      )
      A[2 * i, ] = (
        lambda.n * max(sigma[i], sigma[j]) * basis_vector(j, params)
        - Sigma.hat.extended[i, ] * (b_vec - c_vec)
      )
    }
    A[params + 1, ] = basis_vector(j + p, params) # force c_jj = 0
    model$A = A
    #gurobi_write(model,"model.rlp")
    
    #params = list(OutputFlag=0)
    param = list(Method = 1, OutputFlag = 0)
    result = gurobi(model, param)
    
    ## STEP 1.5:
    # omega.breve_jj = if(sigma_j <= sqrt(n/log(p))) { omega.hat_{jj}} else { sqrt(log(p)/n) }
    
    omega.hat = result$x
    #omega.breve[j] = if(sigma[j] <= sqrt(n/log(p))) { omega.hat[j] - omega.hat[j + p] } else { sqrt(log(p)/n) }
    return(if (sigma[j] <= sqrt(n / log(p))) {
      omega.hat[j] - omega.hat[j + p]
    } else {
      sqrt(log(p) / n)
    })
    #print(result$x)
  }
  
  
  ## STEP 2:
  # omega.tilde_{. j} = argmin |b|_1
  # s.t. for all $i$, |(Sigma.hat b - e_j)_i| \le lambda.n * sqrt(sigma_i * omega.breve_jj)
  
  # i.e. same as before but
  # we replace lambda.n * max(sigma_i,sigma_j) * b_jj by lambda.n * sqrt(sigma_i * omega.breve_jj)
  
  omega.tilde = foreach (
    j = 1:p,
    .combine = cbind,
    .export = c("basis_vector"),
    .packages = c("gurobi")
  ) %dopar% {
    model = list()
    
    params = 2 * p
    model$modelsense = 'min'
    model$obj = rep(1, params) # corresponds to l1 norm
    model$sense = c(rep(c("<", ">"), p), "=")
    rhs = rep(0, params + 1)
    rhs[j * 2 - 1] = -1
    rhs[j * 2] = -1
    for (i in 1:p) {
      # first case is abs is negative, second case is abs is positive
      rhs[i * 2 - 1] = rhs[i * 2 - 1] + lambda.n * sqrt(sigma[i] * omega.breve[j])
      rhs[i * 2] = rhs[i * 2] - lambda.n * sqrt(sigma[i] * omega.breve[j])
    }
    rhs[params + 1] = 0
    model$rhs = rhs
    
    A = matrix(0, params + 1, params)
    Sigma.hat.extended = cbind(Sigma.hat, Sigma.hat)
    for (i in 1:p) {
      b_vec = c(rep(1, p), rep(0, p))
      c_vec = c(rep(0, p), rep(1, p))
      # Note: use ENTRY-wise multiplication below!
      A[2 * i - 1, ] = -Sigma.hat.extended[i, ] * (b_vec - c_vec)
      A[2 * i, ] = -Sigma.hat.extended[i, ] * (b_vec - c_vec)
    }
    A[params + 1, ] = basis_vector(j + p, params) # force c_jj = 0
    model$A = A
    #gurobi_write(model,"model.rlp")
    
    #params = list(OutputFlag=0)
    param = list(Method = 1, OutputFlag = 0)
    result = gurobi(model, param)
    #omega.tilde[,j] = result$x[1:p] - result$x[(p + 1):(2 * p)]
    return(result$x[1:p] - result$x[(p + 1):(2 * p)])
  }
  
  ## STEP 2.5
  # Omega_{ij} = smaller of omega.tilde_{ij} and omega.tilde_{ji} in norm.
  for (i in 1:p) {
    for (j in 1:p) {
      if (abs(omega.tilde[i, j]) < abs(omega.tilde[j, i])) {
        omega.tilde[j, i] = omega.tilde[i, j]
      }
    }
  }
  
  return(omega.tilde)
}

MB.singlenode = function(X, i, lambda) {
  rslt = glmnet(X[, -i], X[, i], lambda = lambda, standardize = TRUE)
  
  sigma.hat = norm(X[, i] - X[, -i] %*% rslt$beta, type = "f") / sqrt(dim(X)[1])
  # w = -\Theta_{ij}/\Theta_{ii} = \Theta_{ij} * sigma^2
  # so \Theta_{ij} = -w/sigma^2
  if (i > length(rslt$beta)) {
    #annoying corner case, R's indexing is bad
    theta.hat = c(-rslt$beta[, 1], 1) / sigma.hat ^ 2
  } else {
    theta.hat = c(-rslt$beta[-(i:length(rslt$beta))], 1,-rslt$beta[i:length(rslt$beta)]) / 
                sigma.hat ^2
  }
  return (theta.hat)
}

MB = function(X, lambda) {
  p = dim(X)[2]
  if (FALSE) {
    Theta = matrix(nrow = p, ncol = p)
    for (i in 1:p) {
      #print(i)
      Theta[i, ] = MB.singlenode(X, i, lambda)
    }
  } else {
    Theta = foreach (
      i = 1:p,
      .combine = rbind,
      .export = c("MB.singlenode"),
      .packages = c("glmnet")
    ) %dopar% {
      MB.singlenode(X, i, lambda)
    }
  }
  return(Theta)
}

best_MB <- function(X, true_prec, lambda_list, kappa_eval = NULL) {
  # corresponds to result with infinite penalty
  p = ncol(X)
  best_prec = matrix(0, nrow = p, ncol = p)
  best_err = Inf
  best_lambda = Inf
  # MB itself is parallelized now
  prec_list = foreach(
    lambda = lambda_list,
    .export = c("MB", "MB.singlenode"),
    .packages = c("glmnet")
  ) %do% {
    #glasso(emp_cov, rho,penalize.diagonal=FALSE)[["wi"]]
    MB(X, lambda)
  }
  for (idx in 1:length(lambda_list)) {
    lambda = lambda_list[idx]
    prec = prec_list[[idx]]
    if (EXPERIMENT == 1) {
      err = exp1_norm(prec, true_prec) # (prec - true_prec,type="f")^2
    } else {
      err = get_adj_err(true_prec, prec, kappa = kappa_eval)
    }
    if (err < best_err) {
      best_err = err
      best_prec = prec
      best_lambda = lambda
    }
  }
  return(list(
    lambda = best_lambda,
    prec = best_prec,
    err = best_err
  ))
}

hybrid.MB.singlenode.v2 = function(X, i, lambda) {
  # pick out single node
  j = greedy(X, 1, i)[1]
  dots.j = (X[, j] %*% X)[1, ]
  # project away column j
  X.proj = X - outer(X[, j], dots.j / (dots.j[j]))
  
  X_proj_norms = sqrt(colSums(X.proj ^ 2) / nrow(X))
  print(X_proj_norms)
  
  # set column j back to original
  X_reg = t(t(X) / X_proj_norms)
  X_reg[, j] = X[, j]
  print(X_reg)
  
  penalty.factor = replicate(dim(X_reg)[2], 1)
  penalty.factor[i] = Inf
  penalty.factor[j] = 0
  
  rslt = glmnet(
    X_reg,
    X[, i],
    lambda = lambda,
    penalty.factor = penalty.factor,
    standardize = FALSE
  )#,exclude=c(i))
  
  sigma.hat = norm(X[, i] - X %*% rslt$beta, type = "f") / sqrt(dim(X)[1])
  
  
  # w = -\Theta_{ij}/\Theta_{ii} = \Theta_{ij} * sigma^2
  # so \Theta_{ij} = -w/sigma^2
  theta.hat = -rslt$beta[, 1] / sigma.hat ^ 2
  theta.hat[i] = 1 / sigma.hat ^ 2
  return (theta.hat)
}

# Unused
l1_constrained_lsq = function(Y, X, lambda) {
  eps = 1e-10 # handle numerical issues w.r.t. psdness
  
  # implement l1 constraint by having signed copies of variables
  p = ncol(X)
  #print(p)
  #rslt = lsqlincon(cbind(X,-X), Y, A = matrix(replicate(p * 2,1),nrow=1), b = lambda)
  C = cbind(X, -X)
  # | Y - A X |_2^2 = <X, A^T A X> - 2 <Y, A X> + <Y,Y>
  model = list()
  
  # gurobi does NOT divide the quadratic term by 2, so we mult. linear by 2 equivalently
  model$modelsense = 'min'
  model$Q = t(C) %*% C + eps * diag(1, nrow = dim(C)[2])
  model$obj = -2 * (t(C) %*% Y)[, 1] # corresponds to l1 norm
  model$A = matrix(replicate(p * 2, 1), nrow = 1)
  model$sense = c("<")
  model$rhs = c(lambda)
  model$lb = 0
  params = list(OutputFlag = 0, Threads = 1)
  
  rslt = gurobi(model, params)
  
  #rslt = quadprog(t(C) %*% C + eps * diag(1,nrow=dim(C)[2]), -(t(C) %*% Y)[,1],
  #                A = matrix(replicate(p * 2,1),nrow=1), b = lambda, lb = 0)
  #print(sum(abs(rslt$xmin)))
  unduplicated = rslt$x[1:p] - rslt$x[(p + 1):(p * 2)]
  return(unduplicated)
}



hybrid.MB.singlenode = function(X,
                                i,
                                gamma = 1,
                                mult_tolerance = 1.01,
                                d_init_factor = 2) {
  # pick out single node
  j = greedy(X, 1, i)[1]
  #print(j)
  dots.j = (X[, j] %*% X)[1, ]
  # project away column j
  X.proj = X - outer(X[, j], dots.j / (dots.j[j]))
  # set column j back to original
  X_reg = X.proj
  X_reg[, j] = X[, j]
  
  X_reg_norms = sqrt(colSums(X_reg ^ 2))
  
  exclude = c(i) # TESTING: was c[i,j]
  penalty.factor = replicate(dim(X_reg)[2], 1)
  penalty.factor[j] = 0

  X_prediction = t(t(X_reg) / X_reg_norms)
  #X_prediction2 = t(t(X_reg[, c(-i)]) / X_reg_norms[c(-i)])
  Y_prediction = X_reg[, i]
  inc = 4
  
  # FORMULA for glmnet max lambda: https://stackoverflow.com/questions/25257780/how-does-glmnet-compute-the-maximal-lambda-value
  # times two to be safe
  lambda = 2 * max(abs(colSums(X_prediction* Y_prediction)))/100
  # start with lambda that is too big, then stop when it's small enough
  MAX_ITERS = 16 # shouldn't be important
  for (iter in 1:MAX_ITERS) {
    rslt = glmnet(
      X_prediction,
      Y_prediction,
      lambda = lambda,
      standardize = FALSE,
      penalty.factor = penalty.factor,
      exclude = exclude
    )
    sigma.hat = norm(Y_prediction - X_prediction %*% rslt$beta, type = "f") /
      sqrt(dim(X)[1])
    l = sum(abs(rslt$beta))
    if (gamma * sigma.hat <= l) {
      # l1 norm is big, stop shrinking
      break
    }
    lambda = lambda/ inc
  }
  upper = lambda * inc
  lower = lambda
  while (upper > lower * mult_tolerance) {
    mid = (upper + lower)/2
    rslt = glmnet(
      X_prediction,
      Y_prediction,
      lambda = mid,
      standardize = FALSE,
      penalty.factor = penalty.factor,
      exclude = c(i)
    )
    l = sum(abs(rslt$beta))
    #print(rslt$beta[j])
    #print(l)
    sigma.hat = norm(Y_prediction - X_prediction %*% rslt$beta, type = "f") /
      sqrt(dim(X)[1])
    if (gamma * sigma.hat <= l) {
      # l1 norm is big, use more regularization
      lower = mid
    } else {
      upper = mid
    }
    #print(mid)
  }
  beta = rslt$beta[,1]/ X_reg_norms
  a_j =  beta[j] + dots.j[i] / dots.j[j] - beta %*% dots.j / (dots.j[j])

    beta[i] = -1
  beta[j] = a_j

  theta.hat = -beta / sigma.hat ^ 2
  return (theta.hat)
}





hybrid.MB = function(X, gamma = 2) {
  p = dim(X)[2]
  #Theta = matrix(nrow=p,ncol=p)
  if (FALSE) {
    CHUNK_SIZE = 20
    NUM_CHUNKS = ceil(p / CHUNK_SIZE)
    
    Theta = foreach (
      chunk = 0:NUM_CHUNKS,
      .combine = rbind,
      .export = c("hybrid.MB.singlenode", "greedy", "l1_constrained_lsq"),
      .packages = c("glmnet")
    ) %do% {
      foreach (i = (chunk + 1):(min(p, chunk + CHUNK_SIZE)), .combine = rbind) %dopar% {
        hybrid.MB.singlenode(X, i, gamma = gamma)
      }
    }
  } else {
    Theta = foreach (
      i = 1:p,
      .combine = rbind,
      .export = c("hybrid.MB.singlenode", "greedy", "l1_constrained_lsq"),
      .packages = c("glmnet")
    ) %dopar% {  # CHANGED
      #print(i)
      hybrid.MB.singlenode(X, i, gamma = gamma)
    }
  }
  # naive symmetrization
  Theta = (Theta + t(Theta))/2
  return(Theta)
}

best_hybrid.MB <-
  function(X,
           true_prec,
           gamma_list,
           d_list = c(1),
           kappa_eval = NULL) {
    # corresponds to result with infinite penalty
    p = ncol(X)
    best_prec = matrix(0, nrow = p, ncol = p)
    best_err = Inf
    best_lambda = c(Inf, Inf)
    #moved parallelism inside
    prec_list = foreach(
      d = d_list,
      .export = c(
        "hybrid.MB",
        "hybrid.MB.singlenode",
        "greedy",
        "l1_constrained_lsq"
      ),
      .packages = c("gurobi")
    ) %:% foreach(gamma = gamma_list) %do% {
      #glasso(emp_cov, rho,penalize.diagonal=FALSE)[["wi"]]
      hybrid.MB(X, gamma = gamma)
    }
    for (idx_d in 1:length(d_list)) {
      d = d_list[idx_d]
      for (idx in 1:length(gamma_list)) {
        gamma = gamma_list[idx]
        prec = prec_list[[idx_d]][[idx]]
        if (EXPERIMENT == 1) {
          err = exp1_norm(prec, true_prec) # (prec - true_prec,type="f")^2
        } else {
          err = get_adj_err(true_prec, prec, kappa = kappa_eval)
        }
        if (err < best_err) {
          best_err = err
          best_prec = prec
          best_lambda = c(d, gamma)
        }
      }
    }
    return(list(
      #d=best_lambda[1],
      gamma = best_lambda[2],
      prec = best_prec,
      err = best_err
    ))
  }

# Do an e.g. 5-fold crossvalidation
crossvalidate = function(X, method, params, folds = 5) {
  n = nrow(X)
  holdout = split(sample(1:n), 1:folds)
  #print(holdout)
  best_param = params[[1]]
  best_obj = Inf
  num_nz = Inf
  
  for (param in params) {
    print(param)
    #print(param[0])
    #print(param[1])
    objective = 0
    for (test_idx in 1:folds) {
      X_train = X[-holdout[[test_idx]], ]
      X_test = X[holdout[[test_idx]], ]
      result = method(X_train, param)
      
      test_size = nrow(X_test)
      holdout_cov = (1 / test_size) * t(X_test) %*% X_test

      # variance reduction cross-validation as in paper
      p = dim(holdout_cov)[1]
      err = 0
      for (i in 1:p) {
        predictor = result[i, ]
        predictor[i] = 0
        if (result[i,i] != 0) {
          predictor = -predictor / result[i, i]
        }
        err = err + (basis_vector(i, p) - predictor) %*% holdout_cov %*% (basis_vector(i, p) - predictor)
      }
      print(err)
      objective = objective + err / p
    }
    objective = objective / folds
    print(objective)
    if (!is.na(objective) && objective < best_obj) {
      best_obj = objective
      best_param = param
    }
  }
  return(list(param = best_param, obj = best_obj))
}

# Equally spaced grid with k + 1 points
log_grid = function(a, b, k) {
  return(exp(log(a) + (0:k) * (log(b) - log(a)) / k))
}

# This may return a list of length < k + 1 due to rounding
log_grid_discrete = function(a, b, k) {
  return(unique(round(log_grid(a, b, k))))
}


is.walk_summable = function(M) {
  N = -abs(M)
  diag(N) = diag(M)
  eigs = eigen(N)$values
  print(eigs)
  return(!any(eigs < 0))
}

walk_summable_dist = function(M) {
  N = -abs(M)
  diag(N) = diag(M)
  return(nearPD(N)$normF)
}

walk_summable_rel = function(M) {
  return(walk_summable_dist(M) / norm(M, ))
}
#glasso
glasso_train = function(X_train, lambda) {
    # symmetrization in our CV objective. (doesn't hurt glasso performance)
    v = glasso((1 / nrow(X_train)) * t(X_train) %*% X_train, lambda)$wi
    return((v + t(v))/2)
  }
#greedy
greedy_params = list()
  for (k in log_grid_discrete(3,24,6)) {#log_grid_discrete(3, 26, 6)) {
    for (threshold in log_grid(0.001, 0.1, 7)) {
      greedy_params = list.append(greedy_params, c(k, threshold))
   } 
  }
  greedy_train = function(X_train, p) {
    full_greedy_and_prune(X_train, (1 / nrow(X_train)) * t(X_train) %*% X_train, p[1], p[2])
  }
niftyret = read.csv("niftyreturns.csv")
covnif = cov(niftyret)
precnif = solve(covnif)
munif = colMeans(niftyret)
for (i in 101: 200)
{
set.seed(i)
sample_distribution <- mvrnorm(n = 200,
                               mu = munif, 
                               Sigma = covnif)
sample_mat = as.matrix(sample_distribution)
cov_lw = linshrink_cov(sample_mat)
prec_lw = solve(cov_lw)
cov_rb = cov.shrink(sample_mat)
prec_rb = solve(cov_rb)
res_glasso = crossvalidate(scale(sample_mat), glasso_train, log_grid(5e-4,0.4,14))
prec_gl = glasso_train(scale(sample_mat), res_glasso$param)
res_greedy = crossvalidate(scale(sample_mat), greedy_train, greedy_params)
prec_gr = greedy_train(scale(sample_mat), res_greedy$param)
l2_lw = norm(precnif - prec_lw, type = "F")  
l2_rb = norm(precnif - prec_rb, type = "F") 
l2_gl = norm(precnif - prec_gl, type = "F") 
l2_gr = norm(precnif - prec_gr, type = "F") 
frob_vec = c(l2_lw,l2_rb,l2_gl,l2_gr)
print(frob_vec)
}
