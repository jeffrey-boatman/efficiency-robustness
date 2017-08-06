# ----------------------------------------
# creates simulation tables for manuscript
# ----------------------------------------
options(digits = 3)
library(xtable)
library(dplyr)


# functions -------------------------------------------------
prop.na <- function(x) mean(is.na(x))


rm.outliers <- function(x) {
  mult <- 9 / 4
  iqr <- IQR(x, na.rm = TRUE)
  m <- median(x, na.rm = TRUE)
  lower <- m - (mult * iqr)
  upper <- m + (mult * iqr)
  x[(x < lower) | (upper < x)] <- NA
  x
}

rmdup <- function(x, dup) {
  x[dup] <- NA
  x 
}


sum.na <- function(x) sum(is.na(x))

last <- function(x) x[length(x)]

make.xtable <- function(df, add.consistency.col = FALSE,
  no.coef = FALSE, clean = FALSE) {
  ident.names <- c("n", "type", "coef")
  ident.pos <- match(ident.names, names(df))
  bias.pos <- grep("bias", names(df))
  bias.names <- names(df)[bias.pos]
  df <- df[, c(ident.pos, bias.pos)]
  new.ident.pos <- match(ident.names, names(df))

  iqr <- aggregate(df[bias.names], 
    df[ident.names],
    IQR,
    na.rm = TRUE)
  medians <- aggregate(df[bias.names], 
    df[ident.names],
    median,
    na.rm = TRUE)
  iqr.long <- reshape(iqr,
    idvar = ident.names,
    times = setdiff(names(iqr), ident.names),
    timevar = "Estimator",
    varying = list(setdiff(names(iqr), ident.names)),
    v.names = "iqr",
    direction = "long")
  row.names(iqr.long) <- NULL
  medians.long <- reshape(medians,
    idvar = ident.names,
    times = setdiff(names(medians), ident.names),
    timevar = "Estimator",
    varying = list(setdiff(names(medians), ident.names)),
    v.names = "median",
    direction = "long")
  row.names(medians.long) <- NULL
  df$dummy.seed <- seq_len(nrow(df))
  df.long <- reshape(df,
    idvar = c(ident.names, "dummy.seed"),
    times = setdiff(bias.names, ident.names),
    timevar = "Estimator",
    varying = list(setdiff(bias.names, ident.names)),
    v.names = "bias",
    direction = "long")
  row.names(df.long) <- NULL
  temp.df <- left_join(df.long, medians.long)
  temp.df <- left_join(temp.df, iqr.long)
  temp.df <- mutate(temp.df,
    lower = median - 9 / 4 * iqr,
    upper = median + 9 / 4 * iqr,
    outlier = bias < lower | upper < bias,
    bias = ifelse(outlier, NA, bias))
  temp.df$median <- temp.df$iqr <- temp.df$lower <- 
    temp.df$upper <- temp.df$outlier <- NULL
  temp.df <- reshape(temp.df,
    idvar = c("dummy.seed", ident.names),
    timevar = "Estimator",
    direction = "wide")
  temp.df$dummy.seed <- NULL
  df$dummy.seed <- NULL
  names(temp.df) <- names(df)

  num.outliers <- aggregate(temp.df[, -new.ident.pos], 
    temp.df[ident.names], 
    sum.na)
  num.outliers <- reshape(num.outliers,
    idvar = ident.names,
    times = setdiff(names(num.outliers), ident.names),
    timevar = "Estimator",
    varying = list(setdiff(names(num.outliers), ident.names)),
    v.names = "n.outliers",
    direction = "long")


  if(clean) {
    df <- temp.df
  }
  means = aggregate(df[, -new.ident.pos],
    df[ident.names],
    mean,
    na.rm = TRUE)
  means <- reshape(means, 
    idvar = ident.names, 
    times = setdiff(names(df), ident.names), 
    timevar = "Estimator",
    varying = list(setdiff(names(df), ident.names)),
    v.names = "bias",
    direction = "long")
  ses <- aggregate(df[, -new.ident.pos],
    df[ident.names],
    sd,
    na.rm = TRUE)
  ses <- reshape(ses,
    idvar = ident.names,
    times = setdiff(names(df), ident.names),
    timevar = "Estimator",
    varying = list(setdiff(names(df), ident.names)),
    v.names = "mcsd",
    direction = "long")
  tt <- left_join(means, ses)
  tt <- mutate(tt, mse = bias ^ 2 + mcsd ^ 2)
  tt <- left_join(tt, num.outliers)
  (ff <- reshape(tt,
    idvar = c("type", "coef", "Estimator"),
    v.names = c(c("bias", "mcsd", "mse", "n.outliers")),
    timevar = "n",
    direction = "wide"))
  ff <- ff[order(ff$type, ff$coef), ]
  rownames(ff) <- NULL

  ff$coef <- rmdup(ff$coef, duplicated(ff[c("type", "coef")]))
  ff$type <- rmdup(ff$type, duplicated(ff$type))
  rr <- c("IPW", "REG", "A-IPW", "CURE", "CURE+", 
    "EM-REG", "W-REG", "A-CURE")
  ff$Estimator <- rep(rr, nrow(ff) / length(rr))
  samp.sizes <- unique(df$n)
  nums.pos <- Map(grep, as.list(as.character(samp.sizes)), list(names(ff)))

  # add extra spaces to allow multiline commands in latex
  foo <- list()
  foo[seq_len(length(nums.pos)) * 2 - 1] <- nums.pos
  bar <- lapply(foo, function(x) ff[, x])
  bar[seq(2, length(bar) - 1, 2)] <- NA
  if(no.coef) {
    ff <- cbind(ff[, 3], do.call(cbind, bar))
  } else {
    ff <- cbind(ff[, 1:3], do.call(cbind, bar))
  }
  base.digits <- c(3, 3, 3, 0, 0)
  if(add.consistency.col) {
    cons <- c("Y", "N", "Y", "Y", "N", "N", "N", "Y")
    ff$consistent <- rep(cons, nrow(ff) / length(rr))
    dd <- c(0, 0, 0, 0, rep(base.digits, length(unique(df$n)) - 1), rev(rev(base.digits)[-1]), 0)
  } else {
    dd <- c(0, 0, 0, 0, rep(base.digits, length(unique(df$n)) - 1), rev(rev(base.digits)[-1]))
  } 
  if (no.coef) {
    dd <- dd[-c(1, 2)]
  }
  xff <- xtable(ff, digits = dd)
  print(xff,
    include.rownames = FALSE,
    include.colnames = FALSE,
    only.contents = TRUE)
}


ttl <- lapply(paste0("outfiles/out_", 1:5, ".txt"), 
  read.table, 
  header = TRUE,
  as.is = TRUE)
tt <- do.call(rbind, ttl)

# ------------------------------------ #
# --- results with known compliers --- #
# ------------------------------------ #
nocoef.out <- subset(tt, type == "Interaction" & coef == 0 & !mis.spec & known)
nomis.out <- subset(tt, coef != 0 & known & !mis.spec)
mis.out <- subset(tt, coef != 0 & known & mis.spec)


with(nocoef.out, ftable(n, type, coef))
with(nomis.out, ftable(type, n, coef))
with(mis.out, ftable(type, n, coef))

nrow(subset(nocoef.out, em.fail == 1)) # 0
nrow(subset(nomis.out, em.fail == 1))  # 0
nrow(subset(mis.out, em.fail == 1))    # 0

# tables with outliers removed ------------------------
make.xtable(nocoef.out, no.coef = TRUE, clean = TRUE)
make.xtable(nomis.out, clean = TRUE)
make.xtable(mis.out, clean = TRUE, add.consistency.col = FALSE)


# Tables for supplementary materials/appendix
make.xtable(nocoef.out, no.coef = TRUE, clean = FALSE)
make.xtable(nomis.out, clean = FALSE)
make.xtable(mis.out, clean = FALSE, add.consistency.col = FALSE)

# --------------------------------------- #
# --- results without known compliers --- #
# --------------------------------------- #
nocoef.out <- subset(tt, type == "Interaction" & coef == 0 & !mis.spec & !known)
nomis.out <- subset(tt, coef != 0 & !known & !mis.spec)
mis.out <- subset(tt, coef != 0 & !known & mis.spec)


with(nocoef.out, ftable(n, type, coef))
with(nomis.out, ftable(type, n, coef))
with(mis.out, ftable(type, n, coef))

nrow(subset(nocoef.out, em.fail == 1)) # 0
nrow(subset(nomis.out, em.fail == 1))  # 0
nrow(subset(mis.out, em.fail == 1))    # 0

# tables with outliers removed ------------------------
make.xtable(nocoef.out, no.coef = TRUE, clean = TRUE)
make.xtable(nomis.out, clean = TRUE)
make.xtable(mis.out, clean = TRUE, add.consistency.col = FALSE)

