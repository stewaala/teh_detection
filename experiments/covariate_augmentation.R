# ==============================================================================
# Add-on: Extra nuisance covariates for the log-link RCT simulator
# ==============================================================================

#' Append P nuisance covariates to an existing simulated dataset
#' (Designed for output of copula_rct_log(); does not modify A, Y, CRTE, CATE)
#'
#' @param dat data.table with columns at least C1..C4, X1..X3, A, Y, CRTE, CATE
#' @param P integer, number of extra covariates to add (e.g., 20, 50, 100)
#' @param z_prefix character, prefix for new columns (default "Z")
#' @param z_seed optional integer; if NULL uses seed + 1e6 in the wrapper
#' @param mix named numeric vec summing to 1: c(risk=., hte=., transforms=., noise=.)
#' @param block_size integer >=1; if >1, induces block collinearity among proxies
#' @param target_corr list with typical weak correlations for proxies:
#'        list(risk = 0.08, hte = 0.08)
#' @param margins list with type proportions (sum to 1):
#'        list(continuous = 0.7, binary = 0.2, count = 0.1)
#' @param keep_latents logical; if TRUE, returns R_latent and H_latent columns
#' @return data.table with P new columns Z001..Zxxx and an attribute "z_meta"
add_noise_covariates <- function(
    dat,
    P,
    z_prefix    = "Z",
    z_seed      = NULL,
    mix         = c(risk = 0.4, hte = 0.2, transforms = 0.2, noise = 0.2),
    block_size  = 5,
    target_corr = list(risk = 0.08, hte = 0.08),
    margins     = list(continuous = 0.7, binary = 0.2, count = 0.1),
    keep_latents = FALSE
) {
  stopifnot(P >= 0)
  if (P == 0L) return(dat)
  
  # Ensure data.table
  if (!data.table::is.data.table(dat)) dat <- data.table::as.data.table(dat)
  
  # RNG
  if (!is.null(z_seed)) set.seed(as.integer(z_seed))
  
  n <- nrow(dat)
  
  # Latent signals (center + scale to unit variance for stable loadings)
  R_latent <- 0.3 * dat[["C1"]] + 0.4 * sin(dat[["C4"]])
  R_latent <- as.numeric(scale(R_latent, center = TRUE, scale = TRUE))
  
  H_latent <- -dat[["C1"]] - dat[["C2"]] + as.numeric(dat[["C3"]] > 0) + dat[["C4"]]^2
  H_latent <- as.numeric(scale(H_latent, center = TRUE, scale = TRUE))
  
  # Mixture sizes (category of how each Z is constructed)
  mix <- pmax(mix, 0); mix <- mix / sum(mix)
  cat_counts <- as.integer(as.vector(rmultinom(1, P, prob = mix)))
  names(cat_counts) <- names(mix)
  # Guarantee exact sum
  while (sum(cat_counts) < P) { cat_counts[which.max(mix)] <- cat_counts[which.max(mix)] + 1 }
  while (sum(cat_counts) > P) { cat_counts[which.max(cat_counts)] <- cat_counts[which.max(cat_counts)] - 1 }
  
  n_risk  <- cat_counts["risk"];        if (is.na(n_risk))  n_risk <- 0L
  n_hte   <- cat_counts["hte"];         if (is.na(n_hte))   n_hte <- 0L
  n_trf   <- cat_counts["transforms"];  if (is.na(n_trf))   n_trf <- 0L
  n_noise <- cat_counts["noise"];       if (is.na(n_noise)) n_noise <- 0L
  
  # Margins (final observed types)
  margins <- pmax(unlist(margins), 0); margins <- margins / sum(margins)
  type_levels <- c("continuous","binary","count")
  types <- sample(type_levels, size = P, replace = TRUE, prob = margins)
  
  # Helper to compute loading 'a' for desired correlation c:
  # If Z = a*S + e, with Var(S)=Var(e)=1 -> Corr(Z,S) = a / sqrt(a^2 + 1).
  loading_for_corr <- function(c) c / sqrt(1 - c^2)
  
  # Block structure (shared noise to induce collinearity within blocks)
  make_block_ids <- function(k) {
    if (block_size <= 1L || k <= 1L) return(rep(NA_integer_, k))
    rep(seq_len(ceiling(k / block_size)), each = block_size)[seq_len(k)]
  }
  
  # Allocate column names
  z_width <- max(3L, nchar(as.character(P)))
  z_names <- sprintf("%s%0*d", z_prefix, z_width, seq_len(P))
  
  Z_mat <- matrix(NA_real_, nrow = n, ncol = P)
  meta  <- vector("list", P)
  
  col_idx <- 0L
  
  # ---- 1) Risk proxies -------------------------------------------------------
  if (n_risk > 0L) {
    block_ids <- make_block_ids(n_risk)
    # Pre-draw shared block noises (N(0,1) vectors)
    block_U <- lapply(unique(block_ids[!is.na(block_ids)]), function(b) rnorm(n))
    names(block_U) <- as.character(unique(block_ids[!is.na(block_ids)]))
    
    for (j in seq_len(n_risk)) {
      col_idx <- col_idx + 1L
      # target correlation around specified mean, with small jitter
      c_tgt <- pmin(pmax(rnorm(1, mean = target_corr$risk, sd = 0.02), 0.03), 0.18)
      a     <- loading_for_corr(c_tgt)
      
      # block collinearity strength
      gamma <- if (!is.na(block_ids[j])) runif(1, 0.3, 0.6) else 0
      uvec  <- if (gamma > 0) block_U[[as.character(block_ids[j])]] else 0
      
      z_lat <- a * R_latent + sqrt(1 - gamma^2) * rnorm(n) + gamma * uvec
      
      Z_mat[, col_idx] <- z_lat
      meta[[col_idx]] <- list(category = "risk", target_corr = c_tgt, loading = a,
                              block_id = block_ids[j], gamma = gamma, base = "R_latent")
    }
  }
  
  # ---- 2) HTE proxies --------------------------------------------------------
  if (n_hte > 0L) {
    block_ids <- make_block_ids(n_hte)
    block_U <- lapply(unique(block_ids[!is.na(block_ids)]), function(b) rnorm(n))
    names(block_U) <- as.character(unique(block_ids[!is.na(block_ids)]))
    
    for (j in seq_len(n_hte)) {
      col_idx <- col_idx + 1L
      c_tgt <- pmin(pmax(rnorm(1, mean = target_corr$hte, sd = 0.02), 0.03), 0.18)
      a     <- loading_for_corr(c_tgt)
      gamma <- if (!is.na(block_ids[j])) runif(1, 0.3, 0.6) else 0
      uvec  <- if (gamma > 0) block_U[[as.character(block_ids[j])]] else 0
      
      z_lat <- a * H_latent + sqrt(1 - gamma^2) * rnorm(n) + gamma * uvec
      
      Z_mat[, col_idx] <- z_lat
      meta[[col_idx]] <- list(category = "hte", target_corr = c_tgt, loading = a,
                              block_id = block_ids[j], gamma = gamma, base = "H_latent")
    }
  }
  
  # ---- 3) Simple transforms of C's -------------------------------------------
  if (n_trf > 0L) {
    C_names <- c("C1","C2","C3","C4")
    g_funs  <- list(
      id  = function(x) x,
      abs = function(x) abs(x),
      sq  = function(x) x^2,
      log1p_abs = function(x) log1p(abs(x))
    )
    for (j in seq_len(n_trf)) {
      col_idx <- col_idx + 1L
      cj  <- sample(C_names, 1)
      gjn <- sample(names(g_funs), 1)
      base <- dat[[cj]]
      # treat C3 as numeric for transforms
      if (cj == "C3") base <- as.numeric(base)
      
      z_lat <- as.numeric(scale(g_funs[[gjn]](base))) + 0.1 * rnorm(n)
      
      Z_mat[, col_idx] <- z_lat
      meta[[col_idx]] <- list(category = "transforms", C = cj, g = gjn)
    }
  }
  
  # ---- 4) Pure noise ----------------------------------------------------------
  if (n_noise > 0L) {
    for (j in seq_len(n_noise)) {
      col_idx <- col_idx + 1L
      # mix of gaussian and heavy-tail
      if (runif(1) < 0.25) {
        z_lat <- rt(n, df = 3)
      } else {
        z_lat <- rnorm(n)
      }
      Z_mat[, col_idx] <- z_lat
      meta[[col_idx]] <- list(category = "noise")
    }
  }
  
  # ---- Apply margins (types): continuous / binary / count --------------------
  # Work column-wise to keep behaviour explicit
  type_info <- vector("list", P)
  for (j in seq_len(P)) {
    zj <- Z_mat[, j]
    
    if (types[j] == "continuous") {
      # optionally produce some log-normal heavy tails
      if (runif(1) < 0.15) {
        zj <- scale(zj)            # center/scale
        zj <- as.numeric(exp(0.2 * zj))  # log-normal-ish, positive skew
        zj <- as.numeric(scale(zj))
      } else {
        zj <- as.numeric(scale(zj))
      }
      Z_mat[, j] <- zj
      type_info[[j]] <- list(type = "continuous")
      
    } else if (types[j] == "binary") {
      # choose prevalence between 0.1 and 0.9
      p <- runif(1, 0.1, 0.9)
      thr <- stats::quantile(zj, probs = 1 - p, names = FALSE, type = 7)
      zb <- as.integer(zj > thr)
      Z_mat[, j] <- zb
      type_info[[j]] <- list(type = "binary", prevalence = round(mean(zb), 3))
      
    } else { # count
      # small-signal log-link Poisson: lambda = exp(mu + s * z_std)
      zstd <- as.numeric(scale(zj))
      s <- runif(1, 0.05, 0.25)                  # weak dependence on latent
      mu <- log(runif(1, 0.3, 1.5))              # average count ~ 0.3 .. 1.5
      lambda <- pmin(exp(mu + s * zstd), 20)     # cap to avoid extremes
      Z_mat[, j] <- rpois(n, lambda)
      type_info[[j]] <- list(type = "count",
                             mean = round(mean(Z_mat[, j]), 3))
    }
  }
  
  # Bind to data.table
  Z_dt <- data.table::as.data.table(Z_mat)
  data.table::setnames(Z_dt, z_names)
  dat_out <- data.table::copy(dat)
  dat_out <- cbind(dat_out, Z_dt)
  
  # Add latents if requested
  if (keep_latents) {
    dat_out[, R_latent := R_latent]
    dat_out[, H_latent := H_latent]
  }
  
  # Store metadata for reproducibility
  z_meta <- list(
    P = P,
    z_prefix = z_prefix,
    z_seed = z_seed,
    mix = mix,
    block_size = block_size,
    target_corr = target_corr,
    margins = margins,
    categories = vapply(meta, function(m) m$category %||% NA_character_, character(1)),
    types = vapply(type_info, function(m) m$type %||% NA_character_, character(1)),
    details = meta,
    type_details = type_info
  )
  attr(dat_out, "z_meta") <- z_meta
  
  dat_out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Wrapper: simulate log-link RCT data and append P nuisance covariates
#' (Calls your original copula_rct_log(); unchanged.)
#'
#' @param n,rho,seed,p_base passed to copula_rct_log()
#' @param P,z_prefix,z_seed,mix,block_size,target_corr,margins,keep_latents
#'        passed to add_noise_covariates()
#' @return data.table with original columns + Z001..Zxxx; attr "z_meta"
copula_rct_log_plusZ <- function(
    n,
    rho         = 0,
    seed        = 111,
    p_base      = .EEXP_H_CONST * exp(-2),
    P           = 0,
    z_prefix    = "Z",
    z_seed      = NULL,
    mix         = c(risk = 0.4, hte = 0.2, transforms = 0.2, noise = 0.2),
    block_size  = 5,
    target_corr = list(risk = 0.08, hte = 0.08),
    margins     = list(continuous = 0.7, binary = 0.2, count = 0.1),
    keep_latents = FALSE
) {
  dat <- copula_rct_log(n = n, rho = rho, seed = seed, p_base = p_base)
  add_noise_covariates(
    dat = dat,
    P = P,
    z_prefix = z_prefix,
    z_seed = z_seed %||% (seed + 1e6),
    mix = mix,
    block_size = block_size,
    target_corr = target_corr,
    margins = margins,
    keep_latents = keep_latents
  )
}
