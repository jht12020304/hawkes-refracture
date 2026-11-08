# Main analysis: fit J0 / J1 / J2 / J2n per cohort, compare by AIC, parametric-bootstrap goodness of fit.
# Output: results/joint_fits.rds, results/joint_models.csv, results/joint_params.csv
source("model.R")
dir.create("results", showWarnings = FALSE)

cat("== self-check ==\n"); sim_check()

inp <- load_inputs(); E <- merge_events(inp$rows)
B <- 25                                    # GOF bootstrap replicates
res <- list(); models <- list(); params <- list()

for (co in c("MrOS", "SOF")){
  d  <- build_d(filter(inp$P, cohort == co), filter(E, cohort == co))
  w3 <- fit_w3(d)
  cat(sprintf("\n== %s: %d persons, %d events, %d deaths in window ==\n", co, nrow(d$p), d$nE, sum(d$D)))

  J0  <- polish(c(w3$par, log(.01), log(1.8), log(.1), log(1)), d, "J0")
  J1  <- polish_multi(lapply(c(.05, .3, .8, 1.5), function(th) c(J0$par, log(th))), d, "J1")
  J2  <- polish_multi(c(lapply(c(.2, .5, 1), function(nu) c(J1$par, log(nu))),
                        lapply(c(.3, .8),    function(th) c(J0$par, log(th), log(.3)))), d, "J2")
  J2n <- polish_multi(list(J2$par[c(1:2, 5:10)],
                           c(J2$par[c(1:2, 5:8)], log(.5), log(.5))), d, "J2n")
  f <- list(J0 = J0, J1 = J1, J2 = J2, J2n = J2n)
  check_nesting(f, co)

  AIC <- sapply(f, aic)
  tab <- data.frame(cohort = co, model = names(f), k = sapply(f, function(m) length(m$par)),
                    logLik = sapply(f, function(m) -m$value), AIC = AIC, dAIC = AIC - min(AIC),
                    convergence = sapply(f, `[[`, "convergence"))
  print(tab, row.names = FALSE, digits = 6)
  best <- names(which.min(AIC))
  pt <- ci_tab(f[[best]], best)
  cat("best =", best, "\n"); print(pt, row.names = FALSE, digits = 3)

  # GOF null: simulate from the fitted joint model, refit W3, recompute the PIT KS distance.
  obs <- gof_D(d, w3$par)
  set.seed(303)
  survC <- sort(d$p$cens_yrs[d$D == 0])
  A <- d$p$cens_yrs; died <- d$D == 1
  A[died] <- vapply(d$p$cens_yrs[died], function(cc){   # admin censor for deaths: draw from survivors' censors
    cand <- survC[survC >= cc]
    if (length(cand)) cand[sample.int(length(cand), 1)] else max(survC) }, 1)
  Dn <- t(vapply(1:B, function(b){
    set.seed(1000 + b)
    ds <- sim_joint(unpack(f[[best]]$par, best), A, trunc4 = (co == "MrOS"))
    o  <- optim(w3$par, nll_w3, d = ds, method = "Nelder-Mead", control = list(maxit = 800))
    gof_D(ds, o$par)
  }, c(D = 0, umed = 0)))
  cat(sprintf("GOF: observed D = %.3f | null median %.3f, 95th pct %.3f\n",
              obs["D"], median(Dn[, 1]), quantile(Dn[, 1], .95)))

  res[[co]] <- list(fits = f, best = best, w3 = w3, gof = list(obs = obs, null = Dn))
  models[[co]] <- tab; params[[co]] <- cbind(cohort = co, model = best, pt)
}

saveRDS(res, "results/joint_fits.rds")
write.csv(bind_rows(models), "results/joint_models.csv", row.names = FALSE)
write.csv(bind_rows(params), "results/joint_params.csv", row.names = FALSE)
