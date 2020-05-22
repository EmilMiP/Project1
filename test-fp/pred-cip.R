library(ggplot2)
library(dplyr)

str(age_dist <- bigreadr::fread2("data-oleguer/RESULTS_age_distribution.txt"))

age_dist %>%
  filter(birth_year == 1950) %>%
  ggplot(aes(age_dx, pct, color = sex)) +
  geom_line() +
  facet_wrap(~ dx)

str(cip <- bigreadr::fread2("data-oleguer/RESULTS_cip.txt"))

cip %>%
  filter(birth_year == 1980) %>%
  ggplot(aes(age, cip, color = sex)) +
  geom_line() +
  facet_wrap(~ dx, scales = "free_y")


cip %>%
  filter(dx == "D33") %>%
  filter(birth_year %% 5 == 0) %>%
  ggplot(aes(age, cip, color = sex)) +
  geom_line() +
  facet_wrap(~ birth_year, scales = "free_y")

library(xgboost)
cip_adhd <- cip %>%
  filter(dx == "adhd") %>%
  select(-dx, -pct)
  
# y <- pmin(qnorm(cip_adhd$cip, lower.tail = FALSE), 5)
y <- cip_adhd$cip
X <- bigstatsr::covar_from_df(select(cip_adhd, -cip))
bst.params <- list(
  max_depth  = 10,
  base_score = 0,
  verbose    = 2,
  nthread    = 4,
  min_child_weight  = 10
)
xgb <- xgboost(X, y, nrounds = 200, params = bst.params)

xgboost::xgb.importance(model = xgb)


params <- expand.grid(age = 0:120, birth_year = 1930:2015, sex = c("F", "M"))
params$cip_pred <- pmax(predict(xgb, bigstatsr::covar_from_df(params)), 1e-4)
params <- params %>%
  group_by(birth_year, sex) %>%
  arrange(age) %>%
  mutate(cip_pred = cummax(cip_pred), 
         thr = qnorm(cip_pred, lower.tail = FALSE))
params[-4]  # would be provided to the user


bigstatsr::plot_grid(
  cip %>%
    filter(dx == "adhd") %>%
    filter(birth_year %% 5 == 0) %>%
    ggplot(aes(age, cip, color = sex)) +
    geom_line() +
    facet_wrap(~ birth_year, scales = "free_y", nrow = 3),
  
  params %>%
    filter(birth_year %% 5 == 0) %>%
    ggplot(aes(age, cip_pred, color = sex)) +
    geom_line() +
    facet_wrap(~ birth_year, scales = "free_y", nrow = 3),
  
  ncol = 1
)
