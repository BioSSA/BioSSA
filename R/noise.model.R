slwsapply <- function(v, w, FUN) {
  sapply(seq_len(length(v) - w + 1), function(i) FUN(v[i : (i+w-1)]))
}

noise.model.default <- function(x, trend,
                                offset = 0,
                                model = "estimate",
                                reg.type = c("winsor", "trim"),
                                reg.level = 0,
                                averaging.type = c("sliding-window", "quantile-break", "equal-break", "none"),
                                breaks = 2, window = 51,
                                FUN = median,
                                FUN.trend = FUN,
                                ...) {
  residuals <- x

  stopifnot(all(is.finite(residuals)))
  stopifnot(all(is.finite(trend)))
  residuals <- residuals[order(trend)]
  trend <- trend[order(trend)]

  residuals.original <- residuals
  trend.original <- trend

  if (is.character(offset)) {
    offset <- match.arg(offset, c("bottom-trend"))

    if (identical(offset, "bottom-trend")) {
      eps <- sqrt(.Machine$double.eps)
      offset <- -(min(trend) - eps)
    } else {
      stop("Unknown `offset' correction type")
    }
  }

  reg.level <- min(reg.level, 1 - reg.level)
  tb <- quantile(trend, probs = c(reg.level / 2, 1 - reg.level / 2))
  rb <- quantile(residuals, probs = c(reg.level / 2, 1 - reg.level / 2))

  reg.type <- match.arg(reg.type)
  if (identical(reg.type, "winsor")) {
    trend[trend < tb[1]] <- tb[1]
    trend[trend > tb[2]] <- tb[2]
    residuals[residuals < rb[1]] <- rb[1]
    residuals[residuals > rb[2]] <- rb[2]
  } else if (identical(reg.type, "trim")) {
    good <- trend >= tb[1] & trend <= tb[2] & residuals >= rb[1] & residuals <= rb[2]
    trend <- trend[good]
    residuals <- residuals[good]
  } else {
    stop("Unknown `reg.type'")
  }

  stopifnot(all(abs(residuals) > 0))
  stopifnot(all(abs(trend + offset) > 0))
  residuals <- log(abs(residuals))
  trend <- log(abs(trend + offset))

  averaging.type <- match.arg(averaging.type)
  if (identical(averaging.type, "none")) {
    trend.means <- trend
    residuals.means <- residuals
  } else if (identical(averaging.type, "equal-break") || identical(averaging.type, "quantile-break")) {
    p <- seq(0, 1, length.out = breaks + 1)
    p.mid <- p[-length(p)] + diff(p)
    p <- p[-c(1, length(p))]

    if (identical(averaging.type, "equal-break")) {
      borders <- log(min(exp(trend)) + p * diff(range(exp(trend))))

      borders.mid <- log(min(exp(trend)) + p.mid * diff(range(exp(trend))))
    } else if (identical(averaging.type, "quantile-break")) {
      borders <- quantile(trend, probs = p)
      borders.mid <- quantile(trend, probs = p.mid)
    }

    g <- findInterval(trend, borders) + 1

    trend.means <- as.vector(tapply(trend, g, FUN.trend))
    residuals.means <- as.vector(tapply(residuals, g, FUN))
  } else if (identical(averaging.type, "sliding-window")) {
    trend.means <- slwsapply(trend, window, FUN.trend)
    residuals.means <- slwsapply(residuals, window, FUN)
  } else {
    stop("Unknown `averaging.type'")
  }

  if (is.character(model))
    model <- match.arg(model, c("estimate", "additive", "multiplicative", "poisson"))
  if (identical(model, "estimate")) {
    R <- lm(residuals.means ~ trend.means)
    alpha <- as.numeric(coef(R)[2])
  } else {
    if (is.character(model)) {
      alpha <- switch(model,
                      additive = 0,
                      multiplicative = 1,
                      poisson = 0.5)
    } else {
      alpha <- model
    }
  }

  sigma <- sqrt(mean((exp(residuals.means) / exp(trend.means)^alpha)^2, na.rm = TRUE))
  residuals.means.fitted <- sigma * exp(trend.means * alpha)
  rresiduals <- residuals.original / (trend.original + offset)^alpha

  res <- list(alpha = alpha,
              sigma = sigma,
              sd = sqrt(mean(rresiduals^2)),
              residuals = residuals.original,
              trend = trend.original,
              residuals.means = exp(residuals.means),
              residuals.means.fitted = residuals.means.fitted,
              trend.means = exp(trend.means) - offset,
              offset = offset,
              averaging.type = averaging.type,
              call = match.call())

  class(res) <- "noise.model"

  res
}

noise.model <- function(x, ...)
  UseMethod("noise.model")

noise.model.BioSSA2d.reconstruction <- function(x, ...) {
  stopifnot(length(x) == 1)

  residuals <- residuals(x)$values
  trend <- x[[1]]$values

  finite <- is.finite(residuals) & is.finite(trend)
  residuals <- residuals[finite]
  trend <- trend[finite]


  noise.model(residuals, trend, ...)
}

noise.model.BioSSA2d <- function(x, groups, ...) {
  groups <- list(trend = if (missing(groups)) seq_len(min(nsigma(x$ssa), nu(x$ssa))) else unlist(groups))

  rec <- reconstruct(x, groups = groups)

  noise.model(rec, ...)
}

noise.model.BioSSA2d3d <- function(x, ...) {
  noise.model(x$bssa2d, ...)
}

noise.model.BioSSA3d <- noise.model.BioSSA2d
noise.model.BioSSA3d.reconstruction <- noise.model.BioSSA2d.reconstruction

print.noise.model <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Noise model:\n")
  cat("\tMultiplicity:", format(x$alpha, digits = digits), "\n")
  cat("\tsigma:", format(x$sigma, digits = digits), "\n")
  cat("\tNoise sd:", format(x$sd, digits = digits), "\n")

  invisible(x)
}

summary.noise.model <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  print(object, digits = digits, ...)
}

plot.noise.model <- function(x,
                             absolute = FALSE, #TODO Mb, remove it????
                             relative = TRUE,
                             draw.residuals = TRUE,
                             draw.means.fitted = FALSE,
                             print.alpha = TRUE,
                             ref = TRUE,
                             symmetric = !absolute,
                             ...,
                             dots.residuals = list(),
                             dots.means.fitted = list(),
                             digits = max(3, getOption("digits") - 3)) {
  if (!relative || x$alpha == 0) {
    ylab.default <- "residuals"
  } else {
    ylab.default <- sprintf("residuals / %s%s",
                            if (x$offset == 0)
                              "trend"
                            else
                              paste0("(trend + ", format(x$offset, digits = digits), ")"),
                            if (x$alpha == 1) "" else paste0("^", format(x$alpha, digits = digits)))
  }

  if (relative) {
    x$residuals <- x$residuals / (x$trend + x$offset) ^ x$alpha
    x$residuals.means <- x$residuals.means / (x$trend.means + x$offset) ^ x$alpha
    x$residuals.means.fitted <- x$residuals.means.fitted / (x$trend.means + x$offset) ^x$alpha
  }

  if (absolute) {
    x$residuals <- abs(x$residuals)
    ylab.default <- sprintf("|%s|", ylab.default)
  }

  xlim.default <- range(0, abs(x$trend))
  ylim.default <- range(0, x$residuals.means)
  if (draw.residuals)
    ylim.default <- range(ylim.default, x$residuals)
  if (symmetric)
    ylim.default <- range(ylim.default, -ylim.default)

  dots <- list(...)
  dots <- .defaults(dots,
                    xlab = "|trend|",
                    ylab = ylab.default)

  dots.means <- .defaults(dots,
                          type = "b",
                          par.settings = list(plot.symbol = list(pch = 15, col = "magenta"),
                                              plot.line = list(lwd = 1, col = "magenta")))

  dots.residuals <- .defaults(modifyList(dots, dots.residuals),
                              type = "p",
                              par.settings = list(plot.symbol = list(pch = 18, col = "blue")))

  dots.means.fitted <- .defaults(modifyList(dots, dots.means.fitted),
                                 type = "l",
                                 par.settings = list(plot.line = list(lwd = 0.5, col = "grey")))

  res <- do.call("xyplot", c(list(residuals.means ~ trend.means,
                                        data = x,
                                        panel = function(x, y, ...) {
                                          panel.abline(h = 0, ..., reference = TRUE)
                                          panel.xyplot(x, y, ...)
                                        }),
                                   dots.means))

  if (draw.residuals) {
    res <- do.call("xyplot", c(list(residuals ~ trend,
                                    data = x),
                               dots.residuals)) + res
  }

  if (draw.means.fitted) {
    res <- res + do.call("xyplot", c(list(residuals.means.fitted ~ trend.means,
                                          data = x),
                                     dots.means.fitted))
  }

  if (print.alpha) {
    res <- res + layer(panel.key(text = text,
                                 points = FALSE,
                                 lines = FALSE,
                                 border = TRUE),
                       data = list(text = sprintf("alpha = %s", format(x$alpha, digits = digits))))
  }

  res
}
