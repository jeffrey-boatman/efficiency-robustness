mixfit <- function (bbio, yy, xx1, xx2, xx3, kknown, tol = 1e-06, alpha.start, gamma.start,
  lamda.start, reg, mis.spec, coef) 
{
  # Estimates the parameters of the mixture distribution.
  # Arguments:
  #  - bbio: vector with the biomarker
  #  - yy: vector with the outcome
  #  - xx1, xx2, xx3: vectors with the confounder
  #  - kknown: vector with indicator for whether compliance is known
  #  - tol: tolerance for changes in negative log likelihood.
  #  - gamma.start: starting coefficients for f(B | Y, C)
  #  - alpha.start: starting coefficients for Pr(C = 1 | X, Y)
  #  - lamda.start: starting coefficients for f(Y|X, X)
		
  bio <- bbio
  y <- yy
  x1 <- xx1
  x2 <- xx2
  x3 <- xx3
  known <- kknown

  alpha <- alpha.start
  gamma <- gamma.start
  lamda <- lamda.start

  alpha.trace <- alpha
  gamma.trace <- gamma
  lamda.trace <- lamda

  if(reg) {
    g.inv <- pnorm
    p.design <- cbind(1, x1, x2)
  } else {
    g.inv <- plogis
    if(coef == 0) {
      p.design <- cbind(1, x1, x2, y)
    } else {
      p.design <- cbind(1, x1, x2, y, x3)
    }
  }
  p.hat <- g.inv(c(p.design %*% alpha))

  # function to evaluate the density
  if(reg) {
    if(mis.spec | coef == 0) {
      x.design <- cbind(1, x1, x2)
    } else {
      x.design <- cbind(1, x1, x2, x3)
    }
    dens <- function(gamma, lamda){
      temp <- matrix(nrow = length(bio), ncol = 3)
      (gamma.c <- c(sum(gamma[1:2]), rev(rev(gamma)[-(1:2)])[-(1:2)]))
      (gamma.nc <- c(gamma[1], rev(rev(gamma)[-(1:2)])[-(1:2)]))
      (b.c.sd <- b.nc.sd <- tail(gamma, 1))

      (lamda.c <- c(sum(lamda[1:2]), rev(rev(lamda)[-(1:2)])[-(1:2)]))
      (lamda.nc <- c(lamda[1], rev(rev(lamda)[-(1:2)])[-(1:2)]))
      (y.c.sd <- y.nc.sd <- tail(lamda, 1))

      temp[, 1] <- dnorm(bio, c(cbind(1, y) %*% gamma.c), b.c.sd)
      temp[, 2] <- dnorm(bio, c(cbind(1, y) %*% gamma.c), b.c.sd) *
        dnorm(y, c(x.design %*% lamda.c), y.c.sd)
      temp[, 3] <- dnorm(bio, c(cbind(1, y) %*% gamma.nc), b.nc.sd) *
        dnorm(y, c(x.design %*% lamda.nc), y.nc.sd)
      temp
    }
  } else {
    dens <- function(gamma, lamda){
      temp <- matrix(nrow = length(bio), ncol = 3)
      gamma.c <- c(sum(gamma[1:2]), rev(rev(gamma)[-(1:2)])[-(1:2)])
      gamma.nc <- c(gamma[1], rev(rev(gamma)[-(1:2)])[-(1:2)])
      (b.c.sd <- b.nc.sd <- tail(gamma, 1))

      temp[, 1] <- temp[, 2] <- dnorm(bio, c(cbind(1, y) %*% gamma.c), b.c.sd)
      temp[, 3] <- dnorm(bio, c(cbind(1, y) %*% gamma.nc), b.nc.sd)
      temp
    }
  }

  # function to evaluate the incomplete negative log likelihood.
  incomp.nll <- function(gamma, lamda, p.hat){
    k <- cbind(known, 
      (1 - known) * p.hat, 
      (1 - known) * (1 - p.hat))
    -sum(log(apply(k * dens(gamma, lamda), 1, sum)))
  }

  k <- cbind(known, 
    (1 - known) * p.hat, 
    (1 - known) * (1 - p.hat))
  nll <- -sum(log(apply(k * dens(gamma, lamda), 1, sum)))

  # first E-step
  dens_ <- cbind(p.hat, 1- p.hat) * (dens(gamma, lamda)[, -1])
  w <- dens_ / apply(dens_, 1, sum)

  nll.diff <- 1
  iter <- 0		
  
  nll.total <- nll

  # weighted sum of squares function
  wss <- function(pars, yy, w1, w2, ...) {
    ll <- list(...)
    xx <- cbind(1, do.call(cbind, ll))
    cc <- w1 * (yy - c(xx %*% c(sum(pars[1:2]), pars[-(1:2)]))) ^ 2
    nc <- w2 * (yy - c(xx %*% c(pars[1], pars[-(1:2)]))) ^ 2
    sum(cc, nc)
  }

  while(nll.diff > tol) {
    gg <- numeric(length(gamma))
    aa <- numeric(length(alpha))
    if(reg) ll <- numeric(length(lamda))

    if(reg) {
      comp_model <- suppressWarnings(
        glm(w[, 1] ~ x1 + x2, weights = (1 - known),
        family = binomial(link = "probit")))
      if(mis.spec | coef == 0) {
        y.lm <- optim(rev(rev(lamda)[-(1:2)]), 
          fn = wss, 
          yy = y,
          w1 = (1 - known) * w[, 1],
          w2 = (1 - known) * w[, 2],
          x1 = x1,
          x2 = x2)
      } else {
        y.lm <- optim(rev(rev(lamda)[-(1:2)]), 
          fn = wss, 
          yy = y,
          w1 = (1 - known) * w[, 1],
          w2 = (1 - known) * w[, 2],
          x1 = x1,
          x2 = x2,
          x3 = x3)
      }
      ll[1:length(y.lm$par)] <- y.lm$par
      ll <- rev(ll)
      ll[1:2] <- weighted.sd(y, 
        mean.1 = c(x.design %*% c(sum(y.lm$par[1:2]), y.lm$par[-(1:2)])),
        mean.2 = c(x.design %*% c(y.lm$par[1], y.lm$par[-(1:2)])),
        weight.1 = (1 - known) * w[, 1],
        weight.2 = (1 - known) * w[, 2])
      ll <- rev(ll)
      lamda <- ll
    } else {
      if(coef == 0) {
        comp_model <- glm(w[, 1] ~ x1 + x2 + y, weights = (1 - known),
          family = binomial(link = logit),
          control = list(maxit = 100))
      } else {
        comp_model <- suppressWarnings(
          glm(w[, 1] ~ x1 + x2 + y + x3, weights = (1 - known),
          family = binomial(link = logit)))
      }
    }
    aa <- coef(comp_model)
  
    bio.lm <- optim(rev(rev(gamma)[-(1:2)]), 
      fn = wss, 
      yy = bio,
      w1 = known + (1 - known) * w[, 1],
      w2 = (1 - known) * w[, 2],
      y = y)

    gg[1:length(bio.lm$par)] <- bio.lm$par
    gg <- rev(gg)
    gg[1:2] <- weighted.sd(bio, 
      mean.1 = c(cbind(1, y) %*% c(sum(bio.lm$par[1:2]), bio.lm$par[-(1:2)])),
      mean.2 = c(cbind(1, y) %*% c(bio.lm$par[1], bio.lm$par[-(1:2)])),
      weight.1 = known + (1 - known) * w[, 1],
      weight.2 = (1 - known) * w[, 2])
    gg <- rev(gg)

    gamma <- gg
    alpha <- unname(aa) 

    alpha.trace <- rbind(alpha.trace, alpha)
    gamma.trace <- rbind(gamma.trace, gamma)
    if(reg) lamda.trace <- rbind(lamda.trace, lamda)

    p.hat <- g.inv(c(p.design %*% alpha))
  
    k <- cbind(known,
      (1 - known) * p.hat,
      (1 - known) * (1 - p.hat))
  
    # E-step
    dens_ <- cbind(p.hat, 1- p.hat) * (dens(gamma, lamda)[, -1])
    w <- dens_ / apply(dens_, 1, sum)

    nll.update <- incomp.nll(gamma, lamda, p.hat)
    nll.diff <- nll - nll.update

    nll <- nll.update
    iter <- iter + 1
    nll.total <- c(nll.total, nll)
  }

  fail <- any(diff(nll.total) >= 0) | nll.diff < 0


  dens_ <- cbind(p.hat, 1 - p.hat) * dens(gamma, lamda)[, -1]
  post.prob <- dens_[, 1]/apply(dens_, 1, sum)
  list(gamma.hat = gamma, 
    lamda.hat = lamda,
    alpha.hat = alpha,
    prob.compliance = weighted.mean(p.hat, 1 - known),
    post.prob = post.prob,
    num.its = length(nll.total),
    neg.log.lik = nll.total[length(nll.total)],
    fail = fail)
}

failwith <- function(default = NULL, f, quiet = TRUE) {
  # borrowed from dplyr package. 
  function(...) {
    out <- default
    try(out <- f(...), silent = quiet)
    out
  }
}

weighted.sd <- function(x, mean.1, mean.2, weight.1, weight.2) {
  # returns weighted sd
  sum.w <- sum(weight.1) + sum(weight.2)
  sqrt(sum(weight.1 * (x - mean.1) ^ 2 + weight.2 * (x - mean.2) ^ 2) / sum.w)
}


quiet <- function(f) {
  function(...) {
    suppressWarnings(f(...))
  }
}


estimate.mixture <- function(bbio, yy, xx1, xx2, xx3, kknown, ccomp, reg, mis.spec, coef, 
  num.tries = 10, type, y.coefs, gamma) {

  bio <- bbio
  y <- yy
  x1 <- xx1
  x2 <- xx2
  x3 <- xx3
  known <- kknown
  comp <- ccomp

  true.starting.vals <- get.true.starting.vals(reg = reg, 
    mis.spec = mis.spec, 
    type = type, 
    y.coefs = y.coefs, 
    gamma = gamma, 
    coef = coef)
  bio.lm <- lm(bio ~ comp + y)
  bio.coefs <- coef(bio.lm)
  bio.rmse <- sqrt(mean(residuals(bio.lm) ^ 2))
  gamma.start <- unname(c(coef(bio.lm), rep(bio.rmse, 2)))
  
  if(reg) {
    if(mis.spec) {
      y.form <- y ~ comp + x1 + x2
    } else {
      if(coef == 0) {
        y.form <- y ~ comp + x1 + x2
      } else {
        y.form <- y ~ comp + x1 + x2 + x3
      }
    }
    y.lm <- lm(y.form, weights = 1 - known)
    y.coefs <- coef(y.lm)

    if(mis.spec | coef == 0) {
      rr <- y - c(cbind(1, comp, x1, x2) %*% coef(y.lm))
    } else {
      rr <- y - c(cbind(1, comp, x1, x2, x3) %*% coef(y.lm))
    }
    y.rmse <- sqrt(weighted.mean(x = rr ^ 2, w = 1 - known))

    lamda.start <- unname(c(y.coefs, rep(y.rmse, 2)))
    alpha.glm <- glm(comp ~ x1 + x2, 
      family = binomial(link = "probit"),
      weights = 1 - known)
  }  else {
    lamda.start <- NA
    if(coef == 0) {
      alpha.glm <- glm(comp ~ x1 + x2 + y, 
        family = binomial(link = logit),
        weights = 1 - known)
      # design <- cbind(1, x1, x2, y)
    } else {
      alpha.glm <- glm(comp ~ x1 + x2 + y + x3, 
        family = binomial(link = logit),
        weights = 1 - known)
    }
  }
  alpha.start <- unname(coef(alpha.glm))



  it.num <- 1
  good.fits <- list()
  while(length(good.fits)  < num.tries & it.num <= 100) {
    if(it.num == 1) {
      a.start <- alpha.start
      g.start <- gamma.start
      l.start <- lamda.start
    } else if(it.num == 2) {
      a.start <- true.starting.vals$alpha.start
      g.start <- true.starting.vals$gamma.start
      l.start <- true.starting.vals$lamda.start
    } else {
      a.start <- add.noise(alpha.start)
      g.start <- add.noise(gamma.start)
      l.start <- add.noise(lamda.start)

      g.start[length(g.start)] <- g.start[length(g.start) - 1]
      if(reg) {
        l.start[length(l.start)] <- l.start[length(l.start) - 1]
      }
    }
    fit <- careful.mixfit(bbio = bio,
      yy = y,
      xx1 = x1,
      xx2 = x2,
      xx3 = x3,
      kknown = known,
      alpha.start = a.start,
      gamma.start = g.start,
      lamda.start = l.start,
      reg = reg,
      mis.spec = mis.spec,
      coef = coef)
    fail <- fit$fail
    if(!fail) good.fits[[length(good.fits) + 1]] <- fit
    it.num <- it.num + 1
  }

  if(length(good.fits) == 0) {
    return(fit)
  } else {
    nlls <- sapply(good.fits, '[[', 'neg.log.lik')
    best.fit <- good.fits[[which.min(nlls)]]
    return(best.fit)
  }
}

estimate <- function(xx1, xx2, xx3, ccomp, bbio, yy, kknown, mmis.spec, coef, num.tries, type, y.coefs, gamma) {
  x1 <- xx1
  x2 <- xx2
  x3 <- xx3
  comp <- ccomp
  bio <- bbio
  y <- yy
  known <- kknown
  mis.spec <- mmis.spec

  den.model <- glm(comp ~ x1 + x2, 
    family = binomial(link = probit),
    weights = 1 - known)

  pi.x <- pnorm(c(cbind(1, x1, x2) %*% coef(den.model)))

  (ipw <- weighted.mean(y, (1 - known) * comp / pi.x))

  if(mis.spec) {
    reg.fm <- y ~ x1 + x2
    x.design <- cbind(1, x1, x2)
  } else {
    if(coef == 0) {
      reg.fm <- y ~ x1 + x2
      x.design <- cbind(1, x1, x2)
    } else {
      reg.fm <- y ~ x1 + x2 + x3
      x.design <- cbind(1, x1, x2, x3)
    }
  }
  reg.model <- lm(reg.fm, weights = (1 - known) * comp)
  ey.x <- c(x.design %*% coef(reg.model))
  c.reg <- weighted.mean(ey.x, 1 - known)
  (aipw <- weighted.mean(
    x = (comp * y) / pi.x - (comp - pi.x) / pi.x * ey.x,
    w = 1 - known))


  fit <- Map(estimate.mixture, 
    bbio = list(bio), 
    yy = list(y), 
    xx1 = list(x1),
    xx2 = list(x2),
    xx3 = list(x3),
    kknown = list(known), 
    ccomp = list(comp),
    reg = as.list(c(FALSE, TRUE)),
    mis.spec = list(mis.spec),
    coef = as.list(rep(coef, 2)),
    num.tries = list(num.tries),
    type = list(type),
    y.coefs = list(y.coefs),
    gamma = list(gamma))

  fail <- sapply(fit, `[[`, "fail") 
  any.fail <- any(fail)
  # numerator for CURE
  cure.pi.bxy <- lapply(fit, `[[`, "post.prob")
  # denominator for CURE
  cure.den.model.func <- function(pi.bxy, weights) {
    glm(pi.bxy ~ x1 + x2, family = binomial(link = "probit"), 
    weights = weights)
  }
  cure.den.model <- lapply(cure.pi.bxy, quiet(cure.den.model.func), 
    weights = 1 - known)
  # lapply(cure.pi.bxy, cure.den.model.func)
  cure.den.coefs <- lapply(cure.den.model, coef)
  cure.pi.x <- lapply(cure.den.coefs, function(cc) c(pnorm(cbind(1, x1, x2) %*% cc)))
  cure.weights <- Map('/', cure.pi.bxy, cure.pi.x)
  cure.weights <- Map('*', cure.weights, list(1 - known))

  # run this to summarize the weights:
  # lapply(cure.pi.bxy, summary)
  # lapply(cure.weights, summary)
  

  na2zero <- function(x) {
    x[is.na(x)] <- 0
    x
  }
  cure.weights <- lapply(cure.weights, na2zero)

  (cure <- mapply(weighted.mean, 
    x = list(y),
    w = cure.weights))
  if(fail[1]) cure[1] <- NA
  if(fail[2]) cure[2] <- NA

  (lamda.comp <- rev(rev(fit[[2]]$lamda.hat)[-c(1, 2)]))
  (lamda.comp <- c(sum(lamda.comp[1:2]), lamda.comp[-c(1, 2)]))

  (reg.est <- weighted.mean(c(x.design %*% lamda.comp), 1 - known))
  if(fail[2]) reg.est <- NA   


  weighted.lm <- lm(y ~ x.design + 0, weights = (1 - known) * cure.pi.bxy[[1]])
  wreg.lp <- c(x.design %*% coef(weighted.lm))
  (wreg <- weighted.mean(wreg.lp, 1 - known))

  aug.fun <- function(pi.alpha, pi.beta, y, expectation.y, k.weights) {
    weighted.mean(
      x = pi.alpha * y / pi.beta - (pi.alpha - pi.beta) / pi.beta * expectation.y,
      w = k.weights)
  }

  (aug.est <- aug.fun(pi.alpha = cure.pi.bxy[[1]],
    pi.beta = cure.pi.x[[1]],
    y = y,
    expectation.y = wreg.lp,
    k.weights = 1 - known))
  if(fail[1]) aug.est <- NA

  estimators <- c(ipw, c.reg, aipw, cure, reg.est, wreg, aug.est)
  list(estimators = estimators, 
    em.fail = any.fail)
}

add.noise <- function(x, jitter.amount = 0.35) {
 x + jitter.amount * x * runif(length(x), -1 , 1)
}

careful.mixfit <- failwith(list(fail = TRUE), mixfit)


get.true.starting.vals <- function(reg, mis.spec, type,
  y.coefs, gamma, coef) {
  set.seed(101)

  nn <- 1e5
  ### large data set for starting vals
  dd <- t(replicate(nn, mu + ts %*% rnorm(3), simplify = TRUE))
  cstar <- dd[, 1]
  x1 <- dd[, 2]
  x2 <- dd[, 3]
  comp <- cstar > qnorm(1 - pp)
  comp <- as.numeric(comp)
  if(type == "Interaction") {
    x3 <- x1 * x2
    mean.x3 <- V[2, 3] + prod(mu[-1])
  } else if(type == "Quadratic") {
    x3 <- x1 ^ 2
  }
  y.lp <- c(cbind(1, comp, x1, x2, x3) %*% y.coefs) 
  y <- y.lp + rnorm(nn, sd = 1.25)
  bio.lp <- c(cbind(1, comp, y) %*% gamma)
  bio <- rnorm(nn, mean = bio.lp, sd = 0.40)
  known <- rep(FALSE, nn)


  bio.lm <- lm(bio ~ comp + y)
  bio.coefs <- coef(bio.lm)
  bio.rmse <- sqrt(mean(residuals(bio.lm) ^ 2))
  gamma.start <- unname(c(coef(bio.lm), rep(bio.rmse, 2)))
  
  if(reg) {
    if(mis.spec) {
      y.form <- y ~ comp + x1 + x2
    } else {
      if(coef == 0) {
        y.form <- y ~ comp + x1 + x2
      } else {
        y.form <- y ~ comp + x1 + x2 + x3
      }
    }
    y.lm <- lm(y.form, weights = 1 - known)
    y.coefs <- coef(y.lm)
    if(mis.spec | coef == 0) {
      rr <- y - c(cbind(1, comp, x1, x2) %*% coef(y.lm))
    } else {
      rr <- y - c(cbind(1, comp, x1, x2, x3) %*% coef(y.lm))
    }
    y.rmse <- sqrt(weighted.mean(x = rr ^ 2, w = 1 - known))
    lamda.start <- unname(c(y.coefs, rep(y.rmse, 2)))
    alpha.glm <- glm(comp ~ x1 + x2, 
      family = binomial(link = "probit"),
      weights = 1 - known)
  }  else {
    lamda.start <- NA
    if(coef == 0) {
      alpha.glm <- glm(comp ~ x1 + x2 + y, 
        family = binomial(link = logit),
        weights = 1 - known)
    } else {
      alpha.glm <- glm(comp ~ x1 + x2 + y + x3, 
        family = binomial(link = logit),
        weights = 1 - known)
    }
  }
  alpha.start <- unname(coef(alpha.glm))
  list(alpha.start = alpha.start,
    lamda.start = lamda.start,
    gamma.start = gamma.start)
}
