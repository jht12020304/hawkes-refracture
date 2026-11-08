# Does a fracture trigger the next one?

Self-exciting (Hawkes) point-process analysis of recurrent fractures in two older-adult osteoporosis
cohorts: **MrOS** (Osteoporotic Fractures in Men Study, 5,994 men) and **SOF** (Study of
Osteoporotic Fractures, 9,182 women).

Recurrent fractures cluster within people. Two explanations give nearly the same event counts:

- **self-excitation**: each fracture raises the short-term risk of the next one;
- **heterogeneity**: some people are simply more fracture-prone (a gamma frailty).

The two differ only in the timing between events, so the analysis models the full event history.
Mortality is very high in these cohorts. Because death can depend on fracture proneness, fractures
and death are modelled jointly with a shared frailty.

## Model

For person *i* with frailty *Z_i* ~ Gamma(1/θ, θ):

```
fracture   λ_i(t) = Z_i    · ( μγ t^(γ-1)  + ηβ    Σ_j exp(-β  (t - t_j)) )
death      h_i(t) = Z_i^ν  · ( ρκ t^(κ-1)  + η_d β_d Σ_j exp(-β_d (t - t_j)) )
```

- η is the expected number of fractures triggered by each fracture; ln2/β is its half-life.
- η_d is the expected number of deaths triggered by each fracture.
- The likelihood is integrated over *Z* by Gauss-Legendre quadrature on the gamma quantile scale.

Four models are compared by AIC:

| Model | Description |
|---|---|
| J0 | independent (θ = 0) |
| J1 | shared frailty, ν = 1 (closed form) |
| J2 | shared frailty, ν free |
| J2n | J2 without the fracture kernel (η = 0): frailty only |

`ΔAIC(J2n − J2) > 0` means the self-exciting kernel is still needed once frailty is accounted for.

## Data

The data are **not included**. Both cohorts are distributed free of charge by the UCSF
coordinating center after registration and a data use agreement:

- MrOS: <https://mrosonline.ucsf.edu> (release AUG24; datasets `FAAUG24`, `EFAUG24`)
- SOF: <https://sofonline.ucsf.edu> (datasets `v{1-6,8,9}demogr`, `v{1-6,8,9}endpt`)

Point the scripts to the downloaded SAS folders:

```bash
export MROS_DIR=/path/to/MrOSAUG24
export SOF_DIR=/path/to/SOF/sas      # folder containing v1demogr/, v1endpt/, ...
```

## Run

Requires R (tested on 4.5) with `haven`, `dplyr` and `tidyr`.

```bash
Rscript prep.R        # raw SAS -> data/events.rds, data/persons.rds
Rscript fit_joint.R   # main analysis -> results/joint_*.csv, results/joint_fits.rds
Rscript coarsen.R     # sensitivity  -> results/coarsen.csv
```

| File | Purpose |
|---|---|
| `prep.R` | MrOS: per-event times of the first 4 fractures, trauma flag, same-day bones merged. SOF: event times rebuilt from visit anchors (`VkDAYS + VkSITEF`), deduplicated within site. Both: unified site groups and censor times. |
| `model.R` | Likelihood, optimiser with multiple starts, nesting check, joint simulator, goodness-of-fit, self-check. |
| `fit_joint.R` | Fits J0/J1/J2/J2n per cohort. Goodness of fit uses a parametric bootstrap from the fitted joint model. |
| `coarsen.R` | Coarsens MrOS to SOF-style observation (first event per site per visit interval) and refits J0/J2/J2n. |

`fit_joint.R` starts with a simulation self-check: it recovers known parameters, and J2 with ν = 1 must
reproduce J1. All fits also check that J2, which nests J0, J1 and J2n, has the highest log-likelihood.

## Known limitations

- MrOS releases times only for the first 4 fractures, and surveillance is censored at the 4th.
- SOF has no per-event trauma flag and no clinical spine fracture field.
- In SOF, a second fracture at the same site within one visit interval is not observed.
- Sex and cohort are completely confounded (MrOS is all men, SOF all women), so a sex difference
  cannot be separated from differences in study design or measurement.

## License

MIT
