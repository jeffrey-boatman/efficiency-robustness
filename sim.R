options(digits = 3)

# Simulation parameters -----------------------------------------------------------------
n.mc <- 1000                       # number of Monte Carlo Iterations
mc.cores <- 24                     # number of cores to use
sample.sizes <- c(250, 500, 1000)  # sample sies assuming no kown compliers
num.tries <- 4                     # number of sets of starting values for em algorithm
n.boot <- 50
# ---------------------------------------------------------------------------------------

library(parallel)

# designed to run on 5 servers
host <- system2("hostname", stdout = TRUE)
hosts <- paste0(c("carbon", "cesium", "chromium", 
  "potassium", "silicon"), 
  ".ccbr.umn.edu")
n.s <- length(hosts)
j <- match(host, hosts)

source("functions.R")

# Pr(C = 1)
pp <- 0.20 

# correlation parameters and mean vector for C', X1, X2
rho1 <- 0.2
rho2 <- 0.075
rho3 <- 0.20
m <- diag(3)
m[upper.tri(m)] <- c(rho1, rho2, rho3)
m[lower.tri(m)] <- c(rho1, rho2, rho3)
s <- sqrt(c(1, 1, 1))
V <- diag(s) %*% m %*% diag(s)    # covariance matrix
ts <- t(chol(V))
mu <- c(0, 2, 2)


# for debugging 
# seed <- 1
# n <- 1000
# coef <- 0
# type <- "Interaction"
# include.known <- FALSE
# mis.spec <- FALSE
# write <- FALSE
# n.boot <- 0
# sample.sizes <- 1000


sim <- function(seed, coef, type = c("Interaction", "Quadratic"), num.tries, 
  include.known = FALSE, mis.spec = TRUE, 
  auc = c("high", "low"),
  write) {
  # coef: coefficient for interaction or quadratic term
  # num.tries: number of sets of starting values for em algorithm. Must be at least 2
  # include.known: whether to include known compliers
  # mis.spec: whether model is misspecified
  # write: whether to write results to file. FALSE is used for debugging. If TRUE,
  #  will write results to outfile, which must be globally defined.

  type <- match.arg(type)
  auc <- match.arg(auc)

  if(missing(mis.spec)) stop("missing mis.spec arg")
  if(missing(write)) stop("missing write arg")

  y.coefs <- c(5, -2.25, 1, 1, coef)

  out <- NULL
  for(n in sample.sizes) {
    set.seed(seed)
    if(include.known) {
      n.known <- n * 0.1
    } else {
      n.known <- 0
    }
    known <- rep(c(TRUE, FALSE), times = c(n.known, n))
    n <- ceiling(n + n.known)
    dd <- t(replicate(n, mu + ts %*% rnorm(3), simplify = TRUE))
    cstar <- dd[, 1]
    x1 <- dd[, 2]
    x2 <- dd[, 3]
    comp <- cstar > qnorm(1 - pp)
    comp <- as.numeric(comp)
    comp[known] <- TRUE

    if(type == "Interaction") {
      x3 <- x1 * x2
      mean.x3 <- V[2, 3] + prod(mu[-1])
    } else if(type == "Quadratic") {
      x3 <- x1 ^ 2
      mean.x3 <- V[2, 2] + mu[2] ^ 2
    }
    (true.mu <- c(y.coefs %*% c(1, 1, mu[-1], mean.x3)))

    y.lp <- c(cbind(1, comp, x1, x2, x3) %*% y.coefs) 
    y <- y.lp + rnorm(n, sd = 1.25)


    # generate biomarker data
    gamma <- c(-9.3, -0.8, 0.7)    
    bio.lp <- c(cbind(1, comp, y) %*% gamma)
    biost <- rnorm(n)
    # b.sd <- 0.40 
    b.sd.vals <- list(
      high = 0.4400557,
      low = 0.6694395
    ) 
    bio <- bio.lp + b.sd.vals[[auc]] * biost

    estimators.l <- estimate(xx1 = x1, 
      xx2 = x2, 
      xx3 = x3,
      ccomp = comp, 
      bbio = bio, 
      yy = y,
      kknown = known,
      mmis.spec = mis.spec,
      coef = coef,
      num.tries = num.tries, 
      type = type,
      y.coefs = y.coefs,
      gamma = gamma)
    (estimators <- estimators.l$estimators)
    (em.fail <- as.numeric(estimators.l$em.fail))

    # --- bootstrap --- #
    # n.boot <- 10 # now defined globally
    if(coef < 2 | !mis.spec) {
      boot.dist <- array(dim = c(n.boot, length(estimators)))
      n.unknown <- n - n.known
      for(b in seq_len(n.boot)) {
        set.seed(72576 + seed * n.boot + b)
        known.samp <- sort(sample(seq_len(n.known), 
          size = n.known, replace = TRUE))
        unknown.samp <- sort(sample(n.known + seq_len(n.unknown), 
          size = n.unknown, 
          replace = TRUE))
        samp <- c(known.samp, unknown.samp)
        boot.estimators.l <- estimate(xx1 = x1[samp], 
          xx2 = x2[samp], 
          xx3 = x3[samp],
          ccomp = comp[samp], 
          bbio = bio[samp], 
          yy = y[samp],
          kknown = known[samp],
          mmis.spec = mis.spec,
          coef = coef,
          num.tries = num.tries, 
          type = type,
          y.coefs = y.coefs,
          gamma = gamma)
        (boot.estimators <- boot.estimators.l$estimators)
        boot.dist[b, ] <- boot.estimators
      }

      boot.sd <- apply(boot.dist, 2, sd, na.rm = TRUE)

      conf.level <- 0.95
      (z <- qnorm((1 + conf.level) / 2))
      (boot.conf.lims <- array(
        data = estimators + cbind(- z * boot.sd, z * boot.sd),
        dim = c(length(estimators), 2),
        dimnames = list(base.names, c("2.5%", "97.5%"))
      ))

      covers <- boot.conf.lims[, 1] < true.mu & 
        true.mu < boot.conf.lims[, 2]
      covers <- as.numeric(covers)
    } else{
      boot.sd <- covers <- rep(NA, length(estimators))
    }
    

    rr <- data.frame(t(c(n, seed, type, auc, coef, 
      estimators - true.mu,
      boot.sd,
      covers, 
      known = include.known,
      mis.spec = mis.spec,
      em.fail)))
    out <- rbind(out, rr)
  }

  names(out) <- col.names
  if(!write) return(out)
  write.table(out, 
    file = outfile, 
    append = TRUE, 
    quote = FALSE, 
    na = "NA",
    row.names = FALSE, 
    col.names = FALSE)
}
base.names <- c("ipw", "reg", "aipw", 
  "cure", "cure.plus",
  "em.reg", "wreg",
  "acure")
col.names <- c("n", "seed", "type", "auc", "coef", 
  paste0(base.names, ".error"),
  paste0(base.names, ".stderr"),
  paste0(base.names, ".covers"),
  "known",
  "mis.spec",
  "em.fail")

# for debugging
# debug(sim)
# n.boot <- 2
# sim(seed = 1,
  # coef = 2,
  # type = "Interaction",
  # num.tries = 4,
  # include.known = TRUE,
  # mis.spec = FALSE,
  # auc = "high",
  # write = FALSE)

sims <- ((j * n.mc / n.s) - (n.mc / n.s - 1)):(j * n.mc / n.s)

outfile <- paste0("outfiles/out_", j, ".txt")

write.table(t(col.names), 
  file = outfile, 
  quote = FALSE, 
  col.names = FALSE, 
  row.names = FALSE)

args.list <- list(
  # table 2 ----------------
  list(coef = 0, 
    type = "Interaction", 
    mis.spec = FALSE, 
    auc = "low"),
  list(coef = 0, 
    type = "Interaction",  
    mis.spec = FALSE, 
    auc = "high"),
  # table 3 ----------------
  list(coef = 1, 
    type = "Interaction",  
    mis.spec = FALSE, 
    auc = "high"),
  list(coef = 1, 
    type = "Interaction",  
    mis.spec = TRUE, 
    auc = "high"),
  # table 3 ----------------
  list(coef = 1, 
    type = "Quadratic",  
    mis.spec = FALSE, 
    auc = "high"),
  list(coef = 1, 
    type = "Quadratic",  
    mis.spec = TRUE, 
    auc = "high")
)

for(args in args.list) {
  cat(
    "Running with:\n",
      "mc.cores = ", mc.cores, "\n",
      "coef = ", args$coef, "\n",
      "type =", args$type, "\n",
      "num.tries =", num.tries, "\n",
      "include.known =", TRUE, "\n",
      "mis.spec = ", args$mis.spec, "\n",
      "auc =", args$auc, "\n",
      "write =", TRUE, "\n"
  )
  out <- mclapply(sims, 
    sim,
    mc.cores = mc.cores, 
    coef = args$coef, 
    type = args$type, 
    num.tries = num.tries, 
    include.known = TRUE, 
    mis.spec = args$mis.spec,
    auc = args$auc, 
    write = TRUE)
}
