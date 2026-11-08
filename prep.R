# Build unified fracture-event and person tables for MrOS (men) and SOF (women).
# Input : raw SAS releases; set MROS_DIR (MrOSAUG24 folder) and SOF_DIR (folder holding v1demogr, v1endpt, ...)
# Output: data/events.rds (one row per bone/site), data/persons.rds (one row per participant)
# Time axis: years from baseline visit (MrOS V1, SOF V1).
suppressMessages({library(haven); library(dplyr); library(tidyr)})

MROS <- Sys.getenv("MROS_DIR"); SOF <- Sys.getenv("SOF_DIR")
if (!dir.exists(MROS) || !dir.exists(SOF)) stop("Set MROS_DIR and SOF_DIR to the raw SAS release folders.")
dir.create("data", showWarnings = FALSE)

read_ds <- function(base, folder, cols){
  f <- list.files(file.path(base, folder), pattern = "sas7bdat$", ignore.case = TRUE,
                  full.names = TRUE, recursive = TRUE)[1]
  nm <- names(read_sas(f, n_max = 0))
  read_sas(f, col_select = all_of(c("ID", intersect(cols, nm)))) |> zap_labels() |> zap_formats()
}

## ---- Site groups shared by both cohorts (SOF has no clinical spine field) ----
grp_mros <- function(lc){                       # FAFXLC 4-digit code -> first 2 digits
  p <- substr(sprintf("%04d", lc), 1, 2)
  case_when(p == "01" ~ "HIP", p == "02" ~ "WRIST",
            p %in% c("06","07") ~ "ARM",       p %in% c("12","14") ~ "SPINE",
            p %in% c("15","16") ~ "PELVIS",    p %in% c("10","11") ~ "RIBCHEST",
            p %in% c("17","18","19","20","21") ~ "LOWERLIMB",
            TRUE ~ "OTHER")                     # skull, face, cervical spine, hand, finger
}
grp_sof <- function(s) case_when(
  s == "HIP" ~ "HIP", s == "WRST" ~ "WRIST",
  s %in% c("HUM","ELB","CLAV") ~ "ARM",
  s == "PEL" ~ "PELVIS", s == "RIB" ~ "RIBCHEST",
  s %in% c("ULEG","LLEG","KNEE","ANK","FOOT","TOE","HEEL") ~ "LOWERLIMB",
  s %in% c("FACE","HAND","FING") ~ "OTHER", TRUE ~ "UNK")

## ================= MrOS: per-event times for the first 4 fractures =================
fa <- read_ds(MROS, "FAAUG24", c("FAFXN", paste0("FAFX", rep(c("FV","LC","TX"), each = 4), 1:4)))
ef <- read_ds(MROS, "EFAUG24", c("EFSTATUS", "FUCDTIME"))
m  <- left_join(fa, ef, by = "ID") %>% mutate(K = pmin(coalesce(FAFXN, 0), 4))

# FAFXFVk is filled with end of follow-up (not NA) when there is no k-th fracture -> real events are k <= FAFXN.
# FAFXTXk == 1: extreme trauma (adjudicated); excluded downstream, NA kept.
ev_m <- m %>% select(ID, K, matches("^FAFX(FV|LC|TX)[1-4]$")) %>%
  pivot_longer(-c(ID, K), names_to = c(".value", "k"), names_pattern = "FAFX(FV|LC|TX)([1-4])") %>%
  mutate(k = as.integer(k)) %>% filter(k <= K, !is.na(FV)) %>%
  mutate(t_yrs = FV/365.25) %>% arrange(ID, t_yrs) %>% group_by(ID) %>%
  mutate(event_id = cumsum(t_yrs - lag(t_yrs, default = -Inf) > 1/365.25)) %>% ungroup() %>%   # same-day bones = 1 event
  transmute(cohort = "MrOS", ID = as.character(ID), event_id, t_yrs,
            grp = grp_mros(LC), site_raw = as.character(LC), xtrauma = TX)

# Fracture surveillance censor: FAFXFV4 for people with events (also the 4th-event time), else FAFXFV1.
# Death = EFSTATUS == 1 (DADEAD is missing for some confirmed deaths).
pm <- m %>% transmute(cohort = "MrOS", ID = as.character(ID),
  cens_yrs = coalesce(if_else(as.character(ID) %in% ev_m$ID, FAFXFV4, NA_real_), FAFXFV1)/365.25,
  death = as.integer(EFSTATUS == 1), fu_yrs = FUCDTIME/365.25)

## ================= SOF: rebuild event times across visit anchors =================
# At anchor k, VkSITEF = days from Vk to the FIRST fracture at that site, VkSITEI = indicator.
# Event day (from enrolment) = VkDAYS + VkSITEF. The same event is seen from several anchors -> dedup within site.
# Limitation: a 2nd fracture at the same site within one anchor interval is invisible.
ANCHORS <- c(1:6, 8, 9)   # no v7 endpt
TOL     <- 30             # days; same-site records closer than this are one event
SITES   <- c("HIP","WRST","HUM","PEL","RIB","ELB","CLAV","FACE","ANK","HAND",
             "TOE","FOOT","HEEL","ULEG","LLEG","KNEE","FING")

cand <- list(); censk <- list()
for (k in ANCHORS){
  dg <- read_ds(SOF, paste0("v", k, "demogr"), paste0("V", k, "DAYS"))
  ep <- read_ds(SOF, paste0("v", k, "endpt"), paste0("V", k, c(SITES, "ANY"), rep(c("F","I"), each = length(SITES) + 1)))
  d  <- inner_join(dg, ep, by = "ID") %>% rename(DAYS = !!paste0("V", k, "DAYS")) %>% filter(!is.na(DAYS))
  for (s in c(SITES, "ANY")){
    Fc <- paste0("V", k, s, "F"); Ic <- paste0("V", k, s, "I")
    if (!all(c(Fc, Ic) %in% names(d))) next
    x <- d %>% select(ID, DAYS, F = all_of(Fc), I = all_of(Ic)) %>% filter(!is.na(F))
    cand[[paste(k, s)]] <- x %>% filter(I == 1) %>% transmute(ID, site = s, day = DAYS + F)
    if (s == "ANY") censk[[as.character(k)]] <- x %>% filter(I == 0) %>% transmute(ID, cday = DAYS + F)
  }
}
cand <- bind_rows(cand)

# Population: original cohort in v1endpt, excluding V1ANYI undetermined (they can never have events).
base <- inner_join(read_ds(SOF, "v1endpt", c("V1ANYI","V1DEATH","V1FOLALL")),
                   read_ds(SOF, "v1demogr", "V1DAYS"), by = "ID") %>% filter(!is.na(V1ANYI))
cand <- semi_join(cand, base, by = "ID")

dedup <- function(tbl, tol) tbl %>% arrange(ID, site, day) %>% group_by(ID, site) %>%
  mutate(eid = cumsum(day - lag(day, default = -Inf) > tol)) %>%
  group_by(ID, site, eid) %>% summarise(day = min(day), .groups = "drop")
ev_site <- dedup(filter(cand, site != "ANY"), TOL)
ev_any  <- dedup(filter(cand, site == "ANY"), TOL)

new_event_id <- function(x) x %>% arrange(ID, day) %>% group_by(ID) %>%
  mutate(event_id = cumsum(day - lag(day, default = -Inf) > 1)) %>% ungroup()   # <= 1 day apart = same event
ev <- new_event_id(ev_site)
# Some site fields disappear from v5 on: ANY-stream events not matched by a site event are added as UNK.
matched <- ev %>% group_by(ID, event_id) %>% summarise(day = min(day), .groups = "drop") %>%
  inner_join(ev_any %>% select(ID, aday = day), by = "ID", relationship = "many-to-many") %>%
  filter(abs(day - aday) <= TOL) %>% distinct(ID, day = aday)
ev <- new_event_id(bind_rows(select(ev, ID, site, day),
                             anti_join(ev_any, matched, by = c("ID","day")) %>% transmute(ID, site = "UNK", day)))

# Censor = latest event-free anchor follow-up (VkDAYS + VkANYF), last event, or V1 follow-up, whichever is later.
pers <- base %>%
  left_join(bind_rows(censk) %>% group_by(ID) %>% summarise(cens_anchor = max(cday)), by = "ID") %>%
  left_join(ev %>% group_by(ID) %>% summarise(last_ev = max(day)), by = "ID") %>%
  mutate(cens_day = pmax(coalesce(cens_anchor, -Inf), coalesce(last_ev, -Inf), V1DAYS + coalesce(V1FOLALL, 0)))

es <- ev %>% left_join(select(pers, ID, V1DAYS), by = "ID") %>%
  transmute(cohort = "SOF", ID = as.character(ID), event_id, t_yrs = (day - V1DAYS)/365.25,
            grp = grp_sof(site), site_raw = site, xtrauma = NA_integer_)   # SOF has no per-event trauma flag
ps <- pers %>% transmute(cohort = "SOF", ID = as.character(ID), cens_yrs = (cens_day - V1DAYS)/365.25,
                         death = as.integer(coalesce(V1DEATH, 0) == 1), fu_yrs = V1FOLALL/365.25)

## ================= Combine =================
events  <- bind_rows(ev_m, es)
persons <- bind_rows(pm, ps) %>%    # a fracture can be confirmed after last contact -> extend censor to it
  left_join(events %>% group_by(cohort, ID) %>% summarise(last_ev = max(t_yrs), .groups = "drop"),
            by = c("cohort","ID")) %>%
  mutate(cens_yrs = pmax(cens_yrs, last_ev, na.rm = TRUE),
         fu_yrs   = pmax(fu_yrs,   last_ev, na.rm = TRUE)) %>% select(-last_ev)

late <- events %>% inner_join(persons, by = c("cohort","ID")) %>% filter(t_yrs > cens_yrs + 1e-9)
stopifnot(nrow(late) == 0, !anyDuplicated(persons[c("cohort","ID")]))

saveRDS(events,  "data/events.rds")
saveRDS(persons, "data/persons.rds")
print(events %>% filter(coalesce(xtrauma, 0) != 1) %>% distinct(cohort, ID, event_id) %>% count(cohort, name = "events"))
print(count(persons, cohort, name = "persons"))
