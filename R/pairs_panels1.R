#' Draw the SPLOM of several polychoric variables.
#' @description
#' Draw a scatter plot of matrices (SPLOM), with bivariate scatter plots below
#' the diagonal, histograms on the diagonal, and the polychoric correlation
#' coefficients with standard errors above the diagonal.
#' Correlation ellipses are drawn in the same graph. The red lines below the
#' diagonal are the LOESS smoothed lines, fitting a smooth curve between two variables.
#'
#' @param x A dataframe of matrix contains several polychoric variables.
#' @param MLE \code{TRUE} for MLE method, \code{FALSE} for IRLS method.
#' @param smooth \code{TRUE} for fitting a smooth curve between two variables.
#' @param density \code{TRUE} for adding density curves on histograms.
#' @param ellipses \code{TRUE} for adding correlation ellipses.
#' @param pch The points type.
#' @param jiggle \code{TRUE} for adding noise to x.
#' @param factor The noise intensity in jiggle.
#' @param hist.col The color of histograms.
#' @param show.points \code{TRUE} for showing the points.
#' @param rug \code{TRUE} for adding rugs to histograms.
#' @param breaks Similar to the \item{breaks} parameter in function \item{hist}.
#' @param cex.cor The size of coefficients.
#' @param smoother \code{TRUE} for adding gradient effect to points.
#'
#' @return SPLOM as described in the description.
#'
#' @export
#'
#' @examples
#' pairs_panels1(CECERS[,1:5])
#' pairs_panels1(CECERS[,1:5],MLE=T)

pairs_panels1 <- function (x, MLE=FALSE, smooth = TRUE, scale = FALSE, density = TRUE, ellipses = TRUE,
                          pch = 20, jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE,
                           rug = TRUE, breaks = "Sturges", cex.cor = 1, wt = NULL, smoother = FALSE,
                           stars = FALSE, ci = FALSE, alpha = 0.05, ...)
{
  "panel.hist.density" <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1], usr[2], 0, 1.5))
    tax <- table(x)
    if (length(tax) < 11) {
      breaks <- as.numeric(names(tax))
      y <- tax/max(tax)
      interbreak <- min(diff(breaks)) * (length(tax) - 1)/41
      rect(breaks - interbreak, 0, breaks + interbreak, y, col = hist.col)
    }
    else {
      h <- hist(x, breaks = breaks, plot = FALSE)
      breaks <- h$breaks
      nB <- length(breaks)
      y <- h$counts
      y <- y/max(y)
      rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
    }
    if (density) {
      tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd",
                               adjust = 1.2), silent = TRUE)
      if (class(tryd) != "try-error") {
        d$y <- d$y/max(d$y)
        lines(d)
      }
    }
    if (rug) {rug(x)}
  }
  "panel.cor" <- function(x, y, prefix = "", ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    if (is.null(wt)) {
      if(MLE){
        dat <- data.frame(x, y)
        dat <- na.omit(dat)
        dat = as.matrix(dat)
        out1 = polycor::polychor(dat[,1],dat[,2],ML=T,std.err = T)
        r = sprintf("%0.4f", out1$rho)
        std = sprintf("%0.4f", sqrt(out1$var[1,1]))
        rs = paste(r,'\n')
        rs = paste(rs,'(',sep = '')
        rs = paste(rs,std,sep = '')
        rs = paste(rs,')',sep = '')
      }else{
        dat <- cbind(x, y)
        dat <- dat[apply(is.na(dat), 1, sum) == 0, ]
        out1 = esti_polychoric(t(as.matrix(dat)))
        r = sprintf("%0.4f", out1$rho)
        std = sprintf("%0.4f", out1$std)
        rs = paste(r,'\n')
        rs = paste(rs,'(',sep = '')
        rs = paste(rs,std,sep = '')
        rs = paste(rs,')',sep = '')
      }
    }
    else {
      r <- cor.wt(data.frame(x, y), w = wt[, c(1:2)])$r[1,2]
    }
    txt = rs
    txt <- paste(prefix, txt, sep = "")
    if (stars) {
      pval <- r.test(sum(!is.na(x * y)), r)$p
      symp <- symnum(pval, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", " "), legend = FALSE)
      txt <- paste0(txt, symp)
    }
    cex <- cex.cor * 0.8/(max(strwidth("0.12***"), strwidth(txt)))
    if (scale) {
      cex1 <- cex * abs(r)
      if (cex1 < 0.25)
        cex1 <- 0.25
      text(0.5, 0.5, txt, cex = cex1)
    }
    else {
      text(0.5, 0.5, txt, cex = cex)
    }
  }
  "panel.smoother" <- function(x, y, pch = par("pch"), col.smooth = "red",
                               span = 2/3, iter = 3, ...) {
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    if(MLE){
      r = polycor::polychor(x,y,ML=T,std.err = F)
    }else{
      dat <- cbind(x, y)
      dat <- dat[apply(is.na(dat), 1, sum) == 0, ]
      x <- dat[,1]
      y <- dat[,2]
      r <- esti_polychoric(t(as.matrix(dat)))$rho
    }
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (smoother) {
      smoothScatter(x, y, add = TRUE, nrpoints = 0)
    }
    else {
      if (show.points)
        points(x, y, pch = pch, ...)
    }
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      if (smooth & ci) {
        lml <- loess(y ~ x, degree = 1, family = "symmetric")
        tempx <- data.frame(x = seq(min(x, na.rm = TRUE),
                                    max(x, na.rm = TRUE), length.out = 47))
        pred <- predict(lml, newdata = tempx, se = TRUE)
        if (ci) {
          upperci <- pred$fit + confid * pred$se.fit
          lowerci <- pred$fit - confid * pred$se.fit
          polygon(c(tempx$x, rev(tempx$x)), c(lowerci, rev(upperci)),
                  col = adjustcolor("light grey", alpha.f = 0.8), border = NA)
        }
        lines(tempx$x, pred$fit, col = col.smooth, ...)
      }
      else {
        if (smooth)
          lines(stats::lowess(x[ok], y[ok], f = span,
                              iter = iter), col = col.smooth)
      }
    }
    if (ellipses)
      draw.ellipse(xm, ym, xs, ys, r, col.smooth = col.smooth,
                   ...)
  }
  "draw.ellipse" <- function(x = 0, y = 0, xs = 1, ys = 1,
                             r = 0, col.smooth, add = TRUE, segments = 51, ...) {
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    if (!is.na(r)) {
      if (abs(r) > 0)
        theta <- sign(r)/sqrt(2)
      else theta = 1/sqrt(2)
      shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(
        c(theta, theta, -theta, theta), ncol = 2, byrow = TRUE)
      ellipse <- unit.circle %*% shape
      ellipse[, 1] <- ellipse[, 1] * xs + x
      ellipse[, 2] <- ellipse[, 2] * ys + y
      if (show.points)
        points(x, y, pch = 19, col = col.smooth, cex = 1.5)
      lines(ellipse, ...)
    }
  }
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (missing(cex.cor))
    cex.cor <- 1
  for (i in 1:ncol(x)) {
    if (is.character(x[[i]])) {
      x[[i]] <- as.numeric(as.factor(x[[i]]))
      colnames(x)[i] <- paste(colnames(x)[i], "*", sep = "")
    }
  }
  n.obs <- nrow(x)
  confid <- qt(1 - alpha/2, n.obs - 2)
  pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor,
        lower.panel = panel.smoother, pch = pch)
}
