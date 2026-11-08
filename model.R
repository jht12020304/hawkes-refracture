# Joint model of recurrent fractures and death with a shared gamma frailty Z.
#   fracture  lambda_i(t) = Z_i     * (mu*gam*t^(gam-1)   + eta*beta   * sum_j exp(-beta  *(t - t_j)))
#   death     h_i(t)      = Z_i^nu  * (rho*kap*t^(kap-1)  + etad*betad * sum_j exp(-betad *(t - t_j)))
# Window [0, C_i], C = fracture-surveillance censor. eta = expected fractures triggered per fracture,
# ln2/beta = half-life of the triggering effect; etad = expected deaths triggered per fracture.
# Models: J0 independent (theta = 0); J1 shared Z, nu = 1 (closed form); J2 nu free (quadrature);
#         J2n = J2 without the fracture kernel (eta = 0): frailty only.
suppressMessages(library(dplyr))

## ---- Data structure ----
mk_d <- function(evlist, Cvec, Dvec){
  n_i <- lengths(evlist); pid <- rep(seq_along(evlist), n_i); t <- unlist(evlist, use.names = FALSE)
  idx <- split(seq_along(t), pid)
  prs <- do.call(rbind, lapply(idx, function(ii){            # all (earlier, later) event pairs per person
    n <- length(ii); if (n < 2) return(NULL)
    j <- rep(2:n, times = (2:n) - 1); l <- sequence((2:n) - 1)
    cbind(tar = ii[j], gap = t[ii[j]] - t[ii[l]])
  }))
  prev <- rep(NA_integer_, length(t))
  for (ii in idx) prev[ii] <- c(NA_integer_, head(ii, -1))
  list(t = t, tg = pmax(t, 1e-6), pid = pid, C = Cvec, Cg = pmax(Cvec, 1e-6), D = Dvec,
       n = n_i, nE = length(t), gapC = Cvec[pid] - t, prev = prev, idx = idx,
       pr_tar = prs[, "tar"], pr_gap = prs[, "gap"])
}

# Analysis input: extreme-trauma bones dropped, same-day bones merged into one event.
load_inputs <- function(){
  ev <- readRDS("data/events.rds")
  list(rows = filter(ev, coalesce(xtrauma, 0) != 1),
       P    = filter(readRDS("data/persons.rds"), !is.na(cens_yrs)))
}
merge_events <- function(rows) rows %>% group_by(cohort, ID, event_id) %>%
  summarise(t = min(t_yrs), .groups = "drop")

# p: persons of one cohort; e: events (ID, t). Deaths after the window are treated as censored.
build_d <- function(p, e){
  p <- arrange(p, ID)
  e <- e %>% filter(ID %in% p$ID) %>% arrange(ID, t)
  D <- as.integer(p$death == 1 & p$fu_yrs <= p$cens_yrs + 0.1)
  d <- mk_d(split(e$t, factor(match(e$ID, p$ID), levels = seq_len(nrow(p)))), p$cens_yrs, D)
  d$p <- p; d
}

## ---- Gauss-Legendre nodes on [0,1] (Golub-Welsch) ----
gl_nodes <- function(M){
  i <- 1:(M - 1); b <- i/sqrt(4*i^2 - 1)
  A <- matrix(0, M, M); A[cbind(i, i + 1)] <- b; A[cbind(i + 1, i)] <- b
  e <- eigen(A, symmetric = TRUE); o <- order(e$values)
  list(x = (e$values[o] + 1)/2, w = e$vectors[1, o]^2)
}
QUAD <- gl_nodes(30); QUAD_HI <- gl_nodes(100)   # 30 for the search, 100 for the final fit

## ---- Likelihood (all parameters on log scale) ----
PAR_NAMES <- list(
  J0  = c("mu","gam","eta","beta","rho","kap","etad","betad"),
  J1  = c("mu","gam","eta","beta","rho","kap","etad","betad","theta"),
  J2  = c("mu","gam","eta","beta","rho","kap","etad","betad","theta","nu"),
  J2n = c("mu","gam","rho","kap","etad","betad","theta","nu"))

unpack <- function(par, mod){
  q <- setNames(as.list(exp(par)), PAR_NAMES[[mod]])
  modifyList(list(eta = 0, beta = 1, theta = 0, nu = 1), q)
}

nll_j <- function(par, d, mod, Q = QUAD){
  q <- unpack(par, mod)
  lam <- q$mu*q$gam*d$tg^(q$gam - 1); Lam <- q$mu*d$C^q$gam
  hd  <- q$rho*q$kap*d$Cg^(q$kap - 1); Hd <- q$rho*d$C^q$kap
  if (d$nE){
    if (length(d$pr_tar)){
      s <- rowsum(exp(-q$beta*d$pr_gap), d$pr_tar)
      ii <- as.integer(rownames(s)); lam[ii] <- lam[ii] + q$eta*q$beta*s
    }
    tr <- rowsum(cbind(1 - exp(-q$beta*d$gapC), exp(-q$betad*d$gapC), 1 - exp(-q$betad*d$gapC)), d$pid)
    ii <- as.integer(rownames(tr))
    Lam[ii] <- Lam[ii] + q$eta*tr[, 1]
    hd[ii]  <- hd[ii]  + q$etad*q$betad*tr[, 2]
    Hd[ii]  <- Hd[ii]  + q$etad*tr[, 3]
  }
  base <- sum(log(lam)) + sum(d$D*log(hd))
  ll <- if (mod == "J0") base - sum(Lam) - sum(Hd)
  else if (mod == "J1"){ K <- 1/q$theta; nD <- d$n + d$D
    base + sum(lgamma(nD + K) - lgamma(K) + nD*log(q$theta) - (nD + K)*log1p(q$theta*(Lam + Hd)))
  } else {                                     # integrate over Z ~ Gamma(1/theta, theta) via its quantiles
    z  <- qgamma(Q$x, shape = 1/q$theta, scale = q$theta)
    lg <- outer(d$n + q$nu*d$D, log(z)) - outer(Lam, z) - outer(Hd, z^q$nu)
    lg <- sweep(lg, 2, log(Q$w), "+")
    mx <- lg[, 1]; for (j in 2:ncol(lg)) mx <- pmax(mx, lg[, j])
    base + sum(mx + log(rowSums(exp(lg - mx))))
  }
  if (!is.finite(ll)) 1e10 else -ll
}

## ---- Optimisation ----
polish <- function(start, d, mod, it = 10000){
  o <- optim(start, nll_j, d = d, mod = mod, method = "Nelder-Mead", control = list(maxit = it))
  o <- optim(o$par, nll_j, d = d, mod = mod, method = "BFGS", control = list(maxit = 2000))
  optim(o$par, nll_j, d = d, mod = mod, Q = QUAD_HI, method = "BFGS", hessian = TRUE, control = list(maxit = 2000))
}
# theta and nu tend to stick at boundaries from a single start -> fixed set of starts, keep the best.
polish_multi <- function(starts, d, mod){
  fits <- Filter(function(f) !is.null(f) && is.finite(f$value),
                 lapply(starts, function(s) tryCatch(polish(s, d, mod), error = function(e) NULL)))
  if (!length(fits)) stop(mod, ": all starts failed")
  fits[[which.min(vapply(fits, `[[`, 1, "value"))]]
}
# J0 and J2n are both nested in J2: a lower J2 log-likelihood means the optimiser failed.
check_nesting <- function(f, label){
  ll <- sapply(f, function(m) -m$value)
  if (ll["J2"] < max(ll[intersect(c("J0","J1","J2n"), names(ll))]) - 1e-6)
    stop(label, ": nesting violated (", paste(names(ll), round(ll, 3), collapse = " "), ")")
}
aic <- function(m) 2*m$value + 2*length(m$par)
ci_tab <- function(m, mod){
  V  <- tryCatch(solve(m$hessian), error = function(e) matrix(NA, length(m$par), length(m$par)))
  se <- sqrt(pmax(diag(V), 0))
  data.frame(par = PAR_NAMES[[mod]], est = exp(m$par), lo = exp(m$par - 1.96*se), hi = exp(m$par + 1.96*se))
}

## ---- Fracture-only model W3 (theta = 0, no death part): used for starts and goodness of fit ----
nll_w3 <- function(par, d){
  mu <- exp(par[1]); gam <- exp(par[2]); eta <- exp(par[3]); beta <- exp(par[4])
  lam <- mu*gam*d$tg^(gam - 1); Lam <- mu*d$C^gam
  if (d$nE){
    if (length(d$pr_tar)){
      s <- rowsum(exp(-beta*d$pr_gap), d$pr_tar)
      ii <- as.integer(rownames(s)); lam[ii] <- lam[ii] + eta*beta*s
    }
    tr <- rowsum(1 - exp(-beta*d$gapC), d$pid)
    ii <- as.integer(rownames(tr)); Lam[ii] <- Lam[ii] + eta*tr
  }
  ll <- sum(log(lam)) - sum(Lam)
  if (!is.finite(ll)) 1e10 else -ll
}
fit_w3 <- function(d){
  cr <- log(d$nE/sum(d$C))
  fits <- lapply(c(.06, 1), function(b){
    o <- optim(c(cr, log(1.3), log(.2), log(b)), nll_w3, d = d, method = "Nelder-Mead", control = list(maxit = 8000))
    optim(o$par, nll_w3, d = d, method = "BFGS", control = list(maxit = 2000))
  })
  fits[[which.min(vapply(fits, `[[`, 1, "value"))]]
}
# Right-truncation-corrected sequential PIT of each next event; KS distance to Uniform(0,1).
gof_D <- function(d, par){
  mu <- exp(par[1]); gam <- exp(par[2]); eta <- exp(par[3]); beta <- exp(par[4])
  L0 <- mu*d$tg^gam
  if (length(d$pr_tar)){
    s <- rowsum(1 - exp(-beta*d$pr_gap), d$pr_tar)
    ii <- as.integer(rownames(s)); L0[ii] <- L0[ii] + eta*s
  }
  vC <- eta*(1 - exp(-beta*d$gapC)); SC <- numeric(d$nE)
  for (ii in d$idx) SC[ii] <- c(0, cumsum(vC[ii][-length(ii)]))
  LC  <- mu*d$C[d$pid]^gam + SC
  L0p <- ifelse(is.na(d$prev), 0, L0[d$prev])
  u <- (1 - exp(-(L0 - L0p)))/(1 - exp(-(LC - L0p)))
  keep <- abs(d$t - d$C[d$pid]) > 1e-9          # drop events that define the censor
  c(D = unname(suppressWarnings(ks.test(u[keep], "punif"))$statistic), umed = median(u[keep]))
}

## ---- Simulation from the joint model (branching for fractures; death = earliest competing arrival) ----
sim_joint <- function(q, A, trunc4 = FALSE){
  np <- length(A)
  Z  <- if (q$theta > 0) rgamma(np, 1/q$theta, scale = q$theta) else rep(1, np)
  n0 <- rpois(np, Z*q$mu*A^q$gam); pid <- rep(seq_len(np), n0)
  t  <- A[pid]*runif(sum(n0))^(1/q$gam)
  allp <- pid; allt <- t
  while (length(pid)){
    nk <- rpois(length(pid), Z[pid]*q$eta)
    rp <- rep(pid, nk); rt <- rep(t, nk) + rexp(sum(nk), q$beta)
    ok <- rt < A[rp]; pid <- rp[ok]; t <- rt[ok]
    allp <- c(allp, pid); allt <- c(allt, t)
    if (length(allp) > 3e5) stop("sim_joint: explosive branching")
  }
  zv <- Z^q$nu
  Td <- (rexp(np)/(zv*q$rho))^(1/q$kap)                  # baseline death time
  if (length(allt)){                                       # each fracture adds a death candidate of mass z^nu*etad
    mj <- zv[allp]*q$etad; Ej <- rexp(length(allt))
    Tj <- suppressWarnings(ifelse(Ej < mj, allt - log(1 - Ej/mj)/q$betad, Inf))
    tm <- tapply(Tj, allp, min)
    ii <- as.integer(names(tm)); Td[ii] <- pmin(Td[ii], tm)
  }
  C <- pmin(Td, A)
  keep <- allt <= C[allp]; allp <- allp[keep]; allt <- allt[keep]
  if (trunc4){                                             # MrOS: surveillance stops at the 4th fracture
    o <- order(allp, allt); allp <- allp[o]; allt <- allt[o]
    r  <- unlist(lapply(split(allt, allp), seq_along), use.names = FALSE)
    t4 <- tapply(allt[r == 4], allp[r == 4], min)
    ii <- as.integer(names(t4)); C[ii] <- pmin(C[ii], t4)
    keep <- allt <= C[allp]; allp <- allp[keep]; allt <- allt[keep]
  }
  o <- order(allp, allt)
  mk_d(split(allt[o], factor(allp[o], levels = seq_len(np))), C, as.integer(Td <= C + 1e-12))
}

# Self-check: recover known J1 parameters from a perturbed start; J2 with nu = 1 must equal J1.
sim_check <- function(){
  set.seed(11)
  tq <- list(mu = .02, gam = 1.3, eta = .4, beta = .15, rho = .02, kap = 2, etad = .15, betad = 1, theta = .5, nu = 1)
  d  <- sim_joint(tq, rep(15, 4000))
  truth <- log(unlist(tq[PAR_NAMES$J1]))
  f   <- polish(truth + rnorm(9, 0, .3), d, "J1")
  d12 <- nll_j(c(f$par, 0), d, "J2", Q = QUAD_HI) - f$value
  cat(sprintf("self-check: %d events, %d deaths | max log error %.2f | J2(nu=1) - J1 = %.4f\n",
              d$nE, sum(d$D), max(abs(f$par - truth)), d12))
  stopifnot(max(abs(f$par - truth)) < 0.7, abs(d12) < 1)
}
