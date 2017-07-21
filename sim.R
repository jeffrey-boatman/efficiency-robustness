options(digits = 3)

# Simulation parameters -----------------------------------------------------------------
n.mc <- 1000                       # number of Monte Carlo Iterations
mc.cores <- 24                     # number of cores to use
sample.sizes <- c(250, 500, 1000)  # sample sies assuming no kown compliers
num.tries <- 4                     # number of sets of starting values for em algorithm
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
  include.known = FALSE, mis.spec = TRUE, write) {
  # coef: coefficient for interaction or quadratic term
  # num.tries: number of sets of starting values for em algorithm. Must be at least 2
  # include.known: whether to include known compliers
  # mis.spec: whether model is misspecified
  # write: whether to write results to file. FALSE is used for debugging. If TRUE,
  #  will write results to outfile, which must be globally defined.

  type <- match.arg(type)

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
    b.sd <- 0.40 
    bio <- bio.lp + b.sd * biost

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
    rr <- data.frame(t(c(n, seed, type, coef, 
      estimators - true.mu, 
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
col.names <- c("n", "seed", "type", "coef", 
  paste0(base.names, ".bias"),
  "known",
  "mis.spec",
  "em.fail")

sims <- ((j * n.mc / n.s) - (n.mc / n.s - 1)):(j * n.mc / n.s)

outfile <- paste0("outfiles/out_", j, ".txt")

write.table(t(col.names), 
  file = outfile, 
  quote = FALSE, 
  col.names = FALSE, 
  row.names = FALSE)

# coef == 0
for(include.known in c(TRUE, FALSE)) {
  out <- mclapply(sims,
    sim,
    mc.cores = mc.cores,
    coef = 0,
    type = "Interaction",
    num.tries = num.tries,
    include.known = include.known,
    mis.spec = FALSE,
    write = TRUE)
}

# coef != 0
for(include.known in c(TRUE, FALSE)) {
  for(coef in c(1, 2)) {
    for(type in c("Interaction", "Quadratic")) {
      for(mis.spec in c(TRUE, FALSE)) {
        out <- mclapply(sims,
          sim,
          mc.cores = mc.cores,
          coef = coef,
          type = type,
          num.tries = num.tries,
          include.known = include.known,
          mis.spec = mis.spec,
          write = TRUE)
      }
    }
  }
}