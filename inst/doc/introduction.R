## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4,
  dpi = 72
)
set.seed(1)

## ----library------------------------------------------------------------------
library(mhn)

## ----dmhn-basic---------------------------------------------------------------
dmhn(c(0.5, 1, 2), alpha = 2, beta = 1, gamma = 1)
dmhn(1, alpha = 2, beta = 1, gamma = 1, log = TRUE)

## ----dmhn-recycle-------------------------------------------------------------
dmhn(c(0.5, 1, 1.5), alpha = c(1, 2), beta = 1, gamma = c(0, 1, -1))

## ----dmhn-plot, fig.cap = "MHN densities for three parameter triples; the steelblue curve sits in the alpha < 1, gamma >> 0 regime where the density combines a boundary divergence at x to 0+ with an interior local maximum near x = 1.8 (Sun et al. 2023, Lemma 3c). The y-axis is clipped at 1; the divergent left tails of the tomato and steelblue curves continue upward beyond the plot."----
x <- c(seq(0.001, 0.3, length.out = 60), seq(0.3, 5, length.out = 140))
plot(x, dmhn(x, alpha = 2, beta = 1, gamma = 1),
     type = "l", lwd = 2, ylab = "density", ylim = c(0, 1),
     main = expression("MHN(" * alpha * ", " * beta * ", " * gamma * ")"))
lines(x, dmhn(x, alpha = 0.5, beta = 1, gamma = 1), lwd = 2, col = "tomato")
lines(x, dmhn(x, alpha = 0.3, beta = 1, gamma = 4), lwd = 2, col = "steelblue")
legend("topright", bty = "n",
       legend = c("(2, 1, 1)    log-concave, mode near 1",
                  "(0.5, 1, 1)  monotone decreasing",
                  "(0.3, 1, 4)  boundary spike + interior bump"),
       col = c("black", "tomato", "steelblue"), lwd = 2)

## ----pmhn-basic---------------------------------------------------------------
pmhn(c(0.5, 1, 1.5, 2), alpha = 2, beta = 1, gamma = 1)
pmhn(2, alpha = 2, beta = 1, gamma = 1, lower.tail = FALSE)
pmhn(2, alpha = 2, beta = 1, gamma = 1, log.p = TRUE)

## ----pmhn-check---------------------------------------------------------------
q <- 1.5
ref <- integrate(function(x) dmhn(x, 2, 1, 1), 0, q,
                 rel.tol = 1e-10)$value
all.equal(pmhn(q, 2, 1, 1), ref)

## ----pmhn-plot, echo = -c(1, 4), fig.cap = "Density and matching CDF for MHN(2, 1, 1)."----
op <- par(mfrow = c(1, 2))
plot(x, dmhn(x, 2, 1, 1), type = "l", lwd = 2, ylab = "f(x)", main = "Density")
plot(x, pmhn(x, 2, 1, 1), type = "l", lwd = 2, ylab = "F(x)", main = "CDF",
     ylim = c(0, 1))
par(op)

## ----qmhn-basic---------------------------------------------------------------
qmhn(c(0.1, 0.5, 0.9), alpha = 2, beta = 1, gamma = 1)

## ----qmhn-roundtrip-----------------------------------------------------------
p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
all.equal(pmhn(qmhn(p, 2, 1, 1), 2, 1, 1), p, tolerance = 1e-6)

## ----qmhn-flags---------------------------------------------------------------
qmhn(0.95, 2, 1, 1, lower.tail = FALSE)        # = qmhn(0.05)
qmhn(log(0.05), 2, 1, 1, log.p = TRUE)         # same value, log-input

## ----rmhn-basic---------------------------------------------------------------
rmhn(5, alpha = 2, beta = 1, gamma = 1)

# Vector parameters are recycled to length n.
rmhn(5, alpha = c(1, 2, 3, 4, 5))

## ----rmhn-overlay, fig.cap = "10,000 draws from MHN(2, 1, 1) with the true density overlaid."----
set.seed(42)
draws <- rmhn(10000, alpha = 2, beta = 1, gamma = 1)
hist(draws, breaks = 60, probability = TRUE, col = "grey90", border = "white",
     xlab = "x", main = "rmhn(10000, 2, 1, 1)")
lines(x, dmhn(x, 2, 1, 1), lwd = 2, col = "tomato")

## ----rmhn-method-equiv--------------------------------------------------------
set.seed(1); s_rtdr <- rmhn(5000, 2, 1, 1, method = "rtdr")
set.seed(1); s_sun  <- rmhn(5000, 2, 1, 1, method = "sun")
ks.test(s_rtdr, s_sun)$p.value

## ----moments-table------------------------------------------------------------
data.frame(
  Quantity = c("mean", "variance", "skewness", "excess kurtosis", "mode"),
  Function = c("mhn_mean", "mhn_var", "mhn_skewness",
               "mhn_kurtosis", "mhn_mode"),
  Value = c(
    mhn_mean(2, 1, 1),
    mhn_var(2, 1, 1),
    mhn_skewness(2, 1, 1),
    mhn_kurtosis(2, 1, 1),
    mhn_mode(2, 1, 1)
  )
)

## ----mode-na------------------------------------------------------------------
mhn_mode(0.5, 1, -1)

## ----special-sqrt-gamma-------------------------------------------------------
# gamma = 0: dmhn matches the change-of-variable sqrt-Gamma density.
xx <- c(0.5, 1, 1.5, 2)
mhn_d  <- dmhn(xx, alpha = 2, beta = 1, gamma = 0)
ref_d  <- dgamma(xx^2, shape = 1, rate = 1) * 2 * xx
all.equal(mhn_d, ref_d)

## ----special-overlay, fig.cap = "MHN(1, 1, 1) and the equivalent truncated normal density."----
mu_tn    <- 1 / (2 * 1)
sigma_tn <- 1 / sqrt(2 * 1)
tn_dens  <- function(x) {
  dnorm(x, mu_tn, sigma_tn) / pnorm(0, mu_tn, sigma_tn, lower.tail = FALSE)
}
plot(x, dmhn(x, 1, 1, 1), type = "l", lwd = 2, ylab = "density",
     main = "alpha = 1: MHN reduces to truncated normal")
points(x, tn_dens(x), pch = ".", col = "tomato", cex = 2)
legend("topright", bty = "n",
       legend = c("dmhn(x, 1, 1, 1)", "TN reference"),
       col = c("black", "tomato"), lwd = c(2, NA), pch = c(NA, 19))

