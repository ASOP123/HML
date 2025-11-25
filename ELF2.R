# Install if needed:
# install.packages(c("brms", "dplyr", "purrr", "tidyr", "ggplot2"))

library(brms)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

set.seed(123)
#1. Prepare your data
# Make sure these are the right types
elf_data <- elf_data %>%
  mutate(
    State        = factor(State),
    HazardGroup  = factor(HazardGroup),
    InjuryType   = factor(InjuryType),
    AccidentYear = factor(AccidentYear)
  ) %>%
  filter(TrendedUltimate > 0)  # safety check
#2. Fit the Bayesian hierarchical severity model
priors <- c(
  set_prior("normal(10, 5)", class = "Intercept"),
  set_prior("exponential(1)", class = "sd"),
  set_prior("exponential(1)", class = "shape")
)

brm_elf <- brm(
  formula = TrendedUltimate ~ 1 +
    (1 | State) +
    (1 | HazardGroup) +
    (1 | InjuryType) +
    (1 | AccidentYear),
  data   = elf_data,
  family = Gamma(link = "log"),
  prior  = priors,
  chains = 2,
  cores  = 2,
  iter   = 2000,
  control = list(adapt_delta = 0.95)
)

summary(brm_elf)
#3. Function to compute ELF from simulated severities
limits <- c(25000, 50000, 100000, 250000, 500000, 1000000)

elf_curve_from_y <- function(y, limits) {
  tibble(
    Limit = limits,
    ELF   = sapply(limits, function(L) mean(pmin(y, L)) / mean(y))
  )
}

#4. Posterior predictive draws for each Hazard Group
haz_levels      <- levels(elf_data$HazardGroup)
baseline_state  <- levels(elf_data$State)[1]
baseline_injury <- levels(elf_data$InjuryType)[1]
baseline_ay     <- levels(elf_data$AccidentYear)[1]

# One row per hazard group
newdata_hg <- data.frame(
  State        = factor(baseline_state,  levels = levels(elf_data$State)),
  HazardGroup  = factor(haz_levels,      levels = levels(elf_data$HazardGroup)),
  InjuryType   = factor(baseline_injury, levels = levels(elf_data$InjuryType)),
  AccidentYear = factor(baseline_ay,     levels = levels(elf_data$AccidentYear))
)

# Posterior predictive draws *including random effects* (hierarchical model)
y_draws_hg <- posterior_predict(
  brm_elf,
  newdata  = newdata_hg,
  ndraws   = 2000         # <= total posterior draws
  # re_formula = NULL by default: include all group-level effects
)
# y_draws_hg is a matrix: ndraws x n_hazard_groups
dim(y_draws_hg)
#5. ELF curves for each hazard group
elf_hazard <- map_dfr(seq_along(haz_levels), function(j) {
  y <- as.numeric(y_draws_hg[, j])   # all draws for hazard group j
  
  elf_curve_from_y(y, limits) %>%
    mutate(
      HazardGroup = haz_levels[j],
      LineType    = "Hazard"
    )
})

#6. Population ELF curve as mixture of hazard groups
# Mix all hazard-group draws together (equal weight for each group)
y_pop <- as.numeric(y_draws_hg)  # flatten matrix → long vector of severities

elf_pop <- elf_curve_from_y(y_pop, limits) %>%
  mutate(
    HazardGroup = "Population",
    LineType    = "Population"
  )
#7. Combine + plot 7 hazard lines + 1 population line
elf_all <- bind_rows(elf_hazard, elf_pop)

# For labels on the right side
labels_df <- elf_all %>%
  group_by(HazardGroup) %>%
  filter(Limit == max(Limit))

ggplot(elf_all, aes(x = Limit, y = ELF, group = HazardGroup)) +
  geom_line(aes(color = LineType, size = LineType, alpha = LineType)) +
  geom_point(aes(color = LineType), size = 1) +
  
  # Label at the end of each line
  geom_text(
    data = labels_df,
    aes(label = HazardGroup, color = LineType),
    hjust = -0.1, vjust = 0.5,
    size  = 3.5,
    show.legend = FALSE
  ) +
  
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(xlim = c(min(elf_all$Limit), max(elf_all$Limit) * 1.05)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  
  scale_color_manual(
    values = c(Hazard = "red", Population = "turquoise3"),
    breaks = c("Hazard", "Population"),
    labels = c("Hazard groups", "Population")
  ) +
  scale_size_manual(values  = c(Hazard = 0.4, Population = 1.5)) +
  scale_alpha_manual(values = c(Hazard = 0.5, Population = 1.0)) +
  
  labs(
    title = "ELF Curves by HG",
    subtitle = paste(
      "Baseline:",
      baseline_state,
      "| InjuryType =", baseline_injury,
      "| AY =", baseline_ay
    ),
    x = "Per-Claim Limit",
    y = "ELF = E[min(X, L)] / E[X]",
    color = "Level of Effect"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(color = "darkblue")
  )
