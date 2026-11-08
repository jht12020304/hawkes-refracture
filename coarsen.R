# Sensitivity: does the MrOS self-exciting kernel survive when MrOS is observed the way SOF is?
# SOF only records the FIRST fracture per site after each visit anchor, so short same-site gaps
# (the information that identifies self-excitation) are lost. Mimic this on MrOS:
#   grid = SOF median anchor times; keep the first event per (person x grid cell x site).
# Two site resolutions bracket SOF: site_raw (4-digit code, lenient) and grp (8 groups, strict).
# Censor C is left unchanged. Requires results/joint_fits.rds (fit_joint.R) and SOF_DIR.
source("model.R")
suppressMessages(library(haven))

SOF <- Sys.getenv("SOF_DIR")
if (!dir.exists(SOF)) stop("Set SOF_DIR to the SOF SAS release folder.")
read_days <- function(k){
  f <- list.files(file.path(SOF, paste0("v", k, "demogr")), pattern = "sas7bdat$", ignore.case = TRUE,
                  full.names = TRUE, recursive = TRUE)[1]
  read_sas(f, col_select = all_of(c("ID", paste0("V", k, "DAYS")))) |> zap_labels()
}
v1 <- read_days(1)
anchor_med <- sapply(c(2:6, 8, 9), function(k){
  x <- inner_join(v1, read_days(k), by = "ID")
  median((x[[paste0("V", k, "DAYS")]] - x$V1DAYS)/365.25, na.rm = TRUE)
})
grid <- c(0, anchor_med, Inf)

coarsen <- function(rows, grid, sitecol){
  rows %>% mutate(cell = findInterval(t_yrs, grid)) %>%
    group_by(ID, cell, site = .data[[sitecol]]) %>% slice_min(t_yrs, n = 1, with_ties = FALSE) %>%
    ungroup() %>% arrange(ID, t_yrs) %>% group_by(ID) %>%
    mutate(eid = cumsum(t_yrs - lag(t_yrs, default = -Inf) > 1/365.25)) %>%   # same-day sites = 1 event
    group_by(ID, eid) %>% summarise(t = min(t_yrs), .groups = "drop")
}
# Toy check: HIP at 0.5 and 0.9 in (0,1] -> keep 0.5; WRIST 0.6 kept; HIP 1.2 in (1,2] kept.
toy <- data.frame(ID = "A", t_yrs = c(.5, .9, 1.2, .6), grp = c("HIP","HIP","HIP","WRIST"))
stopifnot(isTRUE(all.equal(sort(coarsen(toy, c(0, 1, 2, Inf), "grp")$t), c(.5, .6, 1.2))))

inp  <- load_inputs()
rows <- filter(inp$rows, cohort == "MrOS")
p    <- filter(inp$P, cohort == "MrOS")
f0   <- readRDS("results/joint_fits.rds")$MrOS$fits
n_orig <- nrow(merge_events(rows))

out <- list()
for (v in c("site_raw", "grp")){
  d <- build_d(p, coarsen(rows, grid, v))
  J0  <- polish(f0$J0$par, d, "J0")
  J2  <- polish_multi(list(f0$J2$par, c(J0$par, log(.3), log(.3))), d, "J2")
  J2n <- polish_multi(list(f0$J2n$par, J2$par[c(1:2, 5:10)],
                           c(J2$par[c(1:2, 5:8)], log(.5), log(.5))), d, "J2n")
  f <- list(J0 = J0, J2 = J2, J2n = J2n)
  check_nesting(f, v)
  pt <- ci_tab(J2, "J2")
  out[[v]] <- data.frame(resolution = v, events = d$nE, events_removed = n_orig - d$nE,
    logLik_J0 = -J0$value, logLik_J2 = -J2$value, logLik_J2n = -J2n$value,
    dAIC_J2n_vs_J2 = aic(J2n) - aic(J2),
    eta = pt$est[3], eta_lo = pt$lo[3], eta_hi = pt$hi[3], half_life = log(2)/pt$est[4])
}
tab <- bind_rows(out)
cat(sprintf("SOF grid (years): %s\n", paste(sprintf("%.1f", anchor_med), collapse = ", ")))
cat(sprintf("uncoarsened MrOS: dAIC(J2n - J2) = %.1f\n", aic(f0$J2n) - aic(f0$J2)))
print(tab, row.names = FALSE, digits = 4)
write.csv(tab, "results/coarsen.csv", row.names = FALSE)
