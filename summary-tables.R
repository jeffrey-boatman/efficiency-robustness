# ----------------------------------------
# creates simulation tables for manuscript
# ----------------------------------------
options(digits = 3)
library(xtable)
library(tidyverse)


ttl <- lapply(paste0("outfiles/", dir("outfiles")), 
  read.table, 
  header = TRUE,
  as.is = TRUE)
tt <- do.call(rbind, ttl)

tt <- as.tibble(tt)


# bias
bb <- tt %>%
  select(-ends_with("stderr"), -ends_with("covers")) %>%
  gather(ends_with("error"), key = "estimator", value = "error") %>%
  group_by(n, type, coef, known, mis.spec, estimator) %>%
  mutate(
    median = median(error, na.rm = TRUE),
    iqr = diff(quantile(error, c(0.25, 0.75), na.rm = TRUE)),
    lower = median - 9 / 4 * iqr,
    upper = median + 9 / 4 * iqr,
    outlier = as.numeric(error < lower | error > upper)) %>%
  ungroup() %>%
  separate(estimator, into = c("estimator", "junk"), 
    sep = ".error",
    convert = TRUE) %>%
  select(-junk)

# stderr 
ss <- tt %>%
  select(-ends_with("error"), -ends_with("covers")) %>%
  gather(ends_with("stderr"), key = "estimator", value = "stderr") %>%
  separate(estimator, into = c("estimator", "junk"), 
    sep = ".stderr",
    convert = TRUE) %>%
  select(-junk)

# coverage prob
cc <- tt %>%
  select(-ends_with("error"), -ends_with("stderr")) %>%
  gather(ends_with("covers"), key = "estimator", value = "cp") %>%
  separate(estimator, into = c("estimator", "junk"), 
    sep = ".covers",
    convert = TRUE) %>%
  select(-junk)

ll <- bb %>%
  left_join(ss) %>%
  left_join(cc)

grouped <- ll %>%
  group_by(n, type, auc, coef, known, mis.spec, estimator)

grouped.noout <- filter(grouped, outlier == 0)
##grouped.noout <- grouped

grouped.noout <- grouped.noout %>%
  group_by(n, type, auc, coef, known, mis.spec, estimator)

summ <- grouped.noout %>%
  summarize(
    bias = mean(error, na.rm = TRUE),
    mcsd = sd(error, na.rm = TRUE),
    mse = mean(error ^ 2, na.rm = TRUE),
    cp = mean(cp, na.rm = TRUE)
  ) 


bias <- summ %>%
  select(-mcsd, -mse, -cp) %>%
  spread(key = n, value = bias, sep = ".") %>%
  rename(
    bias.n.275 = n.275,
    bias.n.550 = n.550,
    bias.n.1100 = n.1100
  )

mcsd <- summ %>%
  select(-bias, -mse, -cp) %>%
  spread(key = n, value = mcsd, sep = ".") %>%
  rename(
    mcsd.n.275 = n.275,
    mcsd.n.550 = n.550,
    mcsd.n.1100 = n.1100
  )

mse <- summ %>%
  select(-mcsd, -bias, -cp) %>%
  spread(key = n, value = mse, sep = ".") %>%
  rename(
    mse.n.275 = n.275,
    mse.n.550 = n.550,
    mse.n.1100 = n.1100
  )

cp <- summ %>%
  select(-mcsd, -bias, -mse) %>%
  spread(key = n, value = cp, sep = ".") %>%
  rename(
    cp.n.275 = n.275,
    cp.n.550 = n.550,
    cp.n.1100 = n.1100
  )

summ.no <- grouped %>%
  summarize(
    num.outliers = sum(outlier, na.rm = TRUE)
  )

no <- summ.no %>%
  spread(key = n, value = num.outliers, sep = ".") %>%
  rename(
    num.outliers.n.275 = n.275,
    num.outliers.n.550 = n.550,
    num.outliers.n.1100 = n.1100
  )


est.order <- c("ipw", "reg", "aipw", "cure", 
  "cure.plus", "em.reg", "wreg", "acure")
est.names <- c("IPW", "REG", "A-IPW", "CURE", "CURE+", "EM-REG", 
  "W-REG", "A-CURE")

mm <- bias %>%
  left_join(mcsd) %>%
  left_join(mse) %>%
  left_join(cp) %>%
  left_join(no) 

cbind(1:ncol(mm), names(mm))

mm <- mm %>%
  select(c(1, 3:5, 2, 6, 7, 10, 13, 16, 19, 8, 11, 14, 17, 20, 9, 12, 15, 18, 21)) %>%
  mutate(
    ord = match(estimator, est.order),
    estimator = est.names[ord])

mm <- mm %>%
  ungroup()

digs <- c(0, 0, rep(c(rep(3, 4), 0), 3))
h.after <- c(3,8)
#,11,16)

# --------------- #
# --- Table 2 --- #
# --------------- #

t2 <- mm %>% 
  filter(coef == 0) %>%
  arrange(auc, ord) %>%
  select(-c(1:5), -ord)

print(xtable(t2, digits = digs), 
  file = "xtables/table2.tex",
  hline.after = h.after,
  only.contents = TRUE,
  include.colnames = FALSE,
  include.rownames = FALSE)

# --------------- #
# --- Table 3 --- #
# ---_----------- #

t3 <- mm %>% 
  filter(type == "Interaction",
         coef == 1) %>%
  arrange(mis.spec, ord) %>%
  select(-(1:5), -ord)

print(xtable(t3, digits = digs), 
      file = "xtables/table3.tex",
      hline.after = h.after,
      only.contents = TRUE,
      include.colnames = FALSE,
      include.rownames = FALSE)

# --------------- #
# --- Table 4 --- #
# --------------- #

t4 <- mm %>% 
  filter(type == "Quadratic",
         coef == 1) %>%
  arrange(mis.spec, ord) %>%
  select(-(1:5), -ord)

print(xtable(t4, digits = digs), 
      file = "xtables/table4.tex",
      hline.after = h.after,
      only.contents = TRUE,
      include.colnames = FALSE,
      include.rownames = FALSE)
