
############### Documentation ###############

#' Toydata
#' 
#' This dataset allows to reproduce
#' the examples shown in the vignette.
#' 
#' @docType data
#' @keywords datasets
#' @name toydata
#' @usage data(toydata)
#' @return All entries are numeric.
#' @format A list of numeric vectors and numeric matrices.
NULL

# \itemize{
#   \item first example
#   \item second example
# }

#' Internal functions
#' 
#' @description
#' This page lists and describes all internal functions of
#' the R package \code{\link{globalSeq}}.
#' 
#' \strong{Preparation}
#' \cr \code{\link{intern.estim}}
#' estimates the parameters of the negative binomial distribution
#' by maximum likelihood.
#' \cr \code{\link{intern.permu}}
#' permutes values across samples,
#' either across all samples or across samples within subgroups.
#' \cr \code{\link{intern.score}}
#' computes the score test statistic.
#' 
#' \strong{Testing}
#' \cr \code{\link{intern.crude}}
#' calculates p-values by permutation.
#' \cr \code{\link{intern.focus}}
#' calculates p-values by permutation,
#' focusing on a region of interest.
#' \cr \code{\link{intern.conva}}
#' calculates p-values by permutation,
#' using the method of control variates.
#' 
#' \strong{Decomposition}
#' \cr \code{\link{intern.cov}}
#' decomposes the test statistic to show the influence of covariates.
#' \cr \code{\link{intern.sam}}
#' decomposes the test statistic to show the influence of samples.
#' \cr \code{\link{intern.plot}}
#' plots the contributions of covariates or samples.
#' 
#' \strong{Communication}
#' \cr \code{\link{intern.chromo}}
#' runs through all genes on a chromosome.
#' \cr \code{\link{intern.select}}
#' identifies local covariates.
#' \cr \code{\link{intern.matrix}}
#' transforms data to a numeric matrix.
#' 
#' @name internal
#' @keywords documentation
#' @seealso
#' The user functions of the R package \code{\link{globalSeq}} are
#' \code{\link{cursus}}, \code{\link{omnibus}} and \code{\link{proprius}}.
NULL

# at export internal <- function(){?internal} # old

# globalSeq <- function(){ print(packageDescription('globalSeq'))
# ?globalSeq }


############### Preparation ###############

#' Internal function
#' 
#' This functions estimates the parameters of the negative binomial
#' distribution by maximum likelihood. It is called by the functions
#' \code{\link{omnibus}} and \code{\link{proprius}}.
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams omnibus
#' 
#' @param y
#' random variable: numeric vector of length \code{n}
#' 
#' @details
#' We assume the negative binomial distribution \code{y_i ~ NB(mu,phi)},
#' where the samples are indexed by \code{i} (\code{i=1,...,n}).
#' Our parametrisation leads to \code{E[y]= mu} 
#' and \code{Var[y]= mu + phi mu^2}.
#' With the offset \code{a} the model becomes \code{y_i ~ NB(a_i*mu,phi)},
#' where the \code{a_i} are known.
#' 
#' @return
#' The function returns a list of numeric vectors.
#' 
#' @references
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' @seealso
#' This is an \code{\link{internal}} function. The user functions
#' are \code{\link{cursus}}, \code{\link{omnibus}},
#' and \code{\link{proprius}}.
#' 
#' @examples
#' set.seed(1)
#' y <- rnbinom(n=1000,mu=10,size=1/0.2)
#' intern.estim(y)
#' 
intern.estim <- function(y, offset = NULL) {
    if (is.null(offset)) {
        mu <- rep(mean(y), length(y))
    } else {
        mu <- offset * sum(y)/sum(offset)
    }
    loglik <- function(phi) sum(lgamma(y + 1/phi) - lgamma(1/phi) - lgamma(y + 
        1) - 1/phi * log(1 + mu * phi) + y * log(mu) - y * log(1/phi + 
        mu))
    phi <- suppressWarnings(stats::optimize(loglik, interval = c(0, 1000), tol = 10^{
        -10
    }, maximum = TRUE)$maximum)
    list(mu = mu, phi = phi)
}

#' Internal function
#' 
#' The number of permutations of \code{n} elements is \code{n!}.
#' This function randomly rearranges the elements \code{it} times,
#' and then deletes all duplicates.
#' Thus it finds always less than \code{it} and \code{n!} permutations.
#' If a confounding variable is provided,
#' the function uses stratified permutation.
#' This function is called by the functions \code{\link{omnibus}}
#' and \code{\link{proprius}}.
#' 
#' @export
#' @keywords misc
#' 
#' @param n
#' Number of samples.
#' @param it
#' Number of repetitions.
#' @param group
#' Either \code{NULL} or a factor
#' of length \code{n}.
#' 
#' @return 
#' The function returns a matrix.
#' 
#' @references
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' @seealso
#' This is an \code{\link{internal}} function. The user functions
#' are \code{\link{cursus}}, \code{\link{omnibus}},
#' and \code{\link{proprius}}.
#' 
#' @examples
#' group <- as.factor(c('A','A','B','B','B'))
#' set.seed(1)
#' intern.permu(n=5,it=1000,group=group)
#' 
intern.permu <- function(n, it, group) {
    if (is.null(group)) {
        temp <- matrix(NA, nrow = n, ncol = it + 1)
        temp[, 1] <- 1:n
        temp[, -1] <- replicate(it, sample(1:n))
        unique(temp, MARGIN = 2)
    } else {
        levels <- unique(group)
        temp <- matrix(NA, nrow = n, ncol = it + 1)
        temp[, 1] <- 1:n
        for (i in 1:length(levels)) {
            which <- group == levels[i]
            temp[which, -1] <- replicate(it, sample((1:n)[which]))
        }
        unique(temp, MARGIN = 2)
    }
}

#' Internal function
#' 
#' This function calculates the test statistic.
#' It is called by the function \code{\link{omnibus}}.
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams omnibus
#' 
#' @param y
#' response variable: numeric vector of length \code{n}
#' @param R
#' numeric matrix of dimensions \code{n*n} (see example)
#' 
#' @return
#' The function returns a real number.
#' 
#' @references
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' @seealso
#' This is an \code{\link{internal}} function. The user functions
#' are \code{\link{cursus}}, \code{\link{omnibus}},
#' and \code{\link{proprius}}.
#' 
#' @examples
#' # simulate high-dimensional data
#' n <- 30
#' p <- 100
#' set.seed(1)
#' y <- rnbinom(n,mu=10,size=1/0.25)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' 
#' # calculate test statistic
#' R <- X %*% t(X) / ncol(X)
#' mu <- mean(y)
#' phi <- (var(y)-mu)/mu^2
#' intern.score(y,R,mu,phi)
#' 
intern.score <- function(y, R, mu, phi) {
    0.5 * matrix((y - mu)/(1 + phi * mu), nrow = 1) %*% R %*% matrix((y - 
        mu)/(1 + phi * mu), ncol = 1) - 0.5 * matrix((mu + y * phi * mu)/(1 + 
        phi * mu)^2, nrow = 1) %*% matrix(diag(R), ncol = 1)
}

############### Testing ###############

#' Internal function
#' 
#' Using the parameter estimates \code{mu} and \code{phi}
#' and the permutation matrix \code{perm}, these functions
#' tests for global association between \code{y} and \code{X}.
#' The function \code{\link{intern.crude}} calculates
#' p-values by permutation (without repetitions).
#' The functions \code{\link{intern.focus}} and
#' \code{\link{intern.conva}} use different tricks
#' to increase precision and decrease computational expense.
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams omnibus
#' 
#' @param y
#' response variable: numeric vector of length \code{n}
#' @param X
#' covariate set: numeric matrix with \code{n} rows (samples)
#' and \code{p} columns (covariates)
#' @param mu
#' mean parameters: numeric vector of length \code{n}
#' @param phi
#' dispersion parameter: non-negative real number
#' @param perm
#' permutations: matrix with \code{n} rows (see example)
#' @param focus
#' number between 0 and 1
#' 
#' @details
#' 
#' The function \code{\link{intern.focus}}
#'  uses permutations in chunks.
#' If the remaining permutations do not allow
#' to reach a specified significance level,
#' it stops and rounds the p-value to one.
#' 
#' The function \code{\link{intern.conva}}
#' uses the method of control variates
#' from Senchaudhuri et al. (1995).
#' Roughly speaking, if the test statistics
#' from Rauschenberger et al. (2015)
#' and Goeman et al. (2004) are highly correlated,
#' it returns the asymptotic p-value from Goeman et al. (2004).
#' 
#' @return
#' Each function returns a dataframe,
#' with the p-value in the first row,
#' and the test statistic in the second row.
#' 
#' @references
#' 
#' P Senchaudhuri, CR Mehta, and NR Patel.
#' Estimating exact p values by the method of control variates
#' or Monte Carlo rescue.
#' Journal of the American Statistical Association
#' 1995;90:640-648
#' \href{http://www.jstor.org/stable/2291077}{abstract}
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' JJ Goeman, SA van de Geer, F de Kort, and HC van Houwelingen.
#' A global test for groups of genes:
#' testing association with a clinical outcome.
#' Bioinformatics 2004;20:93-99.
#' \href{http://bioinformatics.oxfordjournals.org/content/20/1/93}{full text}
#' 
#' @seealso
#' 
#' These are an \code{\link{internal}} functions. The user functions 
#' of the R package \code{\link{globalSeq}} are \code{\link{cursus}},
#' \code{\link{omnibus}}, and \code{\link{proprius}}.
#' 
#' @examples
#' 
#' # simulate high-dimensional data
#' n <- 30
#' p <- 100
#' set.seed(1)
#' y <- rnbinom(n,mu=10,size=1/0.25)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' 
#' # prepare arguments
#' mu <- rep(mean(y),n)
#' phi <- (var(y)-mu)/mu^2
#' perm <- intern.permu(n=n,it=100,group=NULL)
#' 
#' # perform tests
#' intern.crude(y,X,mu,phi,perm)
#' intern.focus(y,X,mu,phi,perm,0.05)
#' intern.conva(y,X,mu,phi,perm,NULL)
#' # globaltest::gt(y~X) gives pstar
#' 
intern.crude <- function(y, X, mu, phi, perm) {
    if (ncol(X) == 0) {
        pvalue <- NA
        teststat <- NA
    } else if (sum(y) == 0) {
        pvalue <- 1
        teststat <- NA
    } else {
        R <- X %*% t(X)/ncol(X)
        nb_sim <- apply(perm, 2, function(perm) intern.score(y = y[perm], 
            R = R, mu = mu[perm], phi = phi))
        pvalue <- sum(nb_sim >= nb_sim[1])/ncol(perm)
        teststat <- nb_sim[1]
    }
    data.frame(pvalue = pvalue, teststat = teststat)
}
#' @export
#' @keywords misc
#' @rdname intern.crude
intern.focus <- function(y, X, mu, phi, perm, focus) {
    # focus <- 0.01 # SHOULD BECOME AN ARGUMENT !!!
    if (ncol(X) == 0) {
        pvalue <- NA
        teststat <- NA
    } else if (sum(y) == 0) {
        pvalue <- 1  # was: pvalue <- focus
        teststat <- NA
    } else {
        R <- X %*% t(X)/ncol(X)
        it <- ncol(perm)
        target <- focus * it
        i <- -Inf
        z <- 0
        sim <- rep(NA, it)
        pos <- c(2^(0:floor(log(it, base = 2))), it + 1)
        for (j in 1:(length(pos) - 1)) {
            if (z < target & i <= it) {
                for (i in pos[j]:(pos[j + 1] - 1)) {
                  sim[i] <- intern.score(y = y[perm[, i]], R = R, mu = mu[perm[, 
                    i]], phi = phi)
                }
                z <- z + sum(sim[pos[j]:(pos[j + 1] - 1)] >= sim[1])
            } else {
                z <- i  # new
                break
            }
        }
        # pvalue <- pmin(focus,z/i) # old
        pvalue <- z/i  # new
        teststat <- sim[1]
    }
    data.frame(pvalue = pvalue, teststat = teststat)
}
#' @export
#' @keywords misc
#' @rdname intern.crude
intern.conva <- function(y, X, mu, phi, perm, offset) {
    if (ncol(X) == 0) {
        out <- data.frame(pvalue = NA, teststat = NA, rausch = NA, goeman = NA, 
            cor = NA, pstar = NA)
    } else if (length(unique(y))<=1) { # was sum(y)==0
        out <- data.frame(pvalue = 1, teststat = 0, rausch = 1, goeman = 1, 
            cor = 1, pstar = NA)
    } else {
        R <- X %*% t(X)/ncol(X)
        it <- ncol(perm)
        ### rausch ###
        nb_sim <- apply(perm, 2, function(perm) intern.score(y = y[perm], 
            R = R, mu = mu[perm], phi = phi))
        nb <- nb_sim >= nb_sim[1]
        p_nb <- sum(nb)/it
        ### goeman ###
        gt_sim <- apply(perm, 2, function(perm) t(y[perm] - mean(y)) %*% 
            R %*% (y[perm] - mean(y))/var(y))
        gt <- gt_sim >= gt_sim[1]
        p_gt <- sum(gt)/it
        ### pstar ### pstar <- globaltest::p.value(globaltest::gt(y,~X)) # old
        if(!is.null(offset)){y <- offset*y}
        n <- length(y)
        mu2 <- var(y)
        H <- 1/n * matrix(1, nrow = n) %*% matrix(1, ncol = n)
        I <- diag(n)
        RT <- t(I - H) %*% R %*% (I - H)
        exp <- sum(diag(RT))
        var <- 2/(n + 1) * ((n - 1) * sum(diag(RT %*% RT)) - sum(diag(RT))^2)
        c <- var/(2 * exp)
        v <- 2 * exp^2/var
        pstar <- 1 - stats::pchisq(gt_sim[1]/c, df = v)
        ### monte carlo rescue ###
        D <- nb - gt + pstar
        p <- sum(D)/it
        if (p < 0 | p > 1) {
            if (p < 0) {
                p <- 0
            }
            if (p > 1) {
                p <- 1
            }
        }
        cor <- cor(nb, gt)
        out <- data.frame(pvalue = p, teststat = nb_sim[1], rausch = p_nb, 
            goeman = p_gt, cor = round(cor, 2), pstar = pstar)
        #SE_pvalue = sqrt(sum((D - p)^2)/(it - 1))
        #SE_rausch = sqrt(p_nb * (1 - p_nb)/it)
        #SE_goeman = sqrt(p_gt * (1 - p_gt)/it)
    }
    out
}


############### Decomposition ###############

#' Internal function
#' 
#' These functions calculate the contribution of covariate
#' or samples to the test statistic.
#' They are called by the function \code{\link{proprius}}.
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams intern.crude
#' 
#' @return 
#' Both functions return a numeric vector.
#' 
#' @references
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' JJ Goeman, SA van de Geer, F de Kort, and HC van Houwelingen.
#' A global test for groups of genes:
#' testing association with a clinical outcome.
#' Bioinformatics 2004;20:93-99.
#' 
#' @seealso
#' 
#' This is an \code{\link{internal}} function. The user functions 
#' of the R package \code{\link{globalSeq}} are \code{\link{cursus}},
#' \code{\link{omnibus}}, and \code{\link{proprius}}.
#' 
#' @examples
#' 
#' # simulate high-dimensional data
#' n <- 30
#' p <- 100
#' set.seed(1)
#' y <- rnbinom(n,mu=10,size=1/0.25)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' 
#' # prepare arguments
#' mu <- rep(mean(y),n)
#' phi <- (var(y)-mean(y))/mean(y)^2 
#' 
#' # decompose test statistic
#' intern.sam(y,X,mu,phi)
#' intern.cov(y,X,mu,phi)
#' 
intern.sam <- function(y, X, mu, phi) {
    n <- nrow(X)  # number of samples
    R <- X %*% t(X)/ncol(X)
    u <- rep(NA, n)
    for (i in 1:n) {
        u[i] <- sum(0.5 * (y[i] - mu[i])/(1 + phi * mu[i]) * R[i, ] * (y - 
            mu)/(1 + phi * mu)) - 0.5 * R[i, i] * (mu[i] + y[i] * phi * 
            mu[i])/(1 + phi * mu[i])^2
    }
    names(u) <- rownames(X)
    u
}
#' @export
#' @keywords misc
#' @rdname intern.sam
intern.cov <- function(y, X, mu, phi) {
    p <- ncol(X)  # number of covariates
    u <- rep(NA, p)
    for (i in 1:p) {
        R <- 1/p * matrix(X[, i], ncol = 1) %*% matrix(X[, i], nrow = 1)
        u[i] <- 0.5 * matrix((y - mu)/(1 + phi * mu), nrow = 1) %*% R %*% 
            matrix((y - mu)/(1 + phi * mu), ncol = 1) - 0.5 * matrix((mu + 
            y * phi * mu)/(1 + phi * mu)^2, nrow = 1) %*% matrix(diag(R), 
            ncol = 1)
    }
    names(u) <- colnames(X)
    u
}

#' Internal function
#' 
#' This function plots the individual contributions
#' to the test statistic.
#' It is called by the function \code{\link{proprius}}.
#' 
#' @export
#' @keywords misc
#' 
#' @param u
#' influence:
#' numeric vector of length \code{n}
#' @param upper
#' critical values:
#' numeric vector of length \code{n}
#' @param xlab
#' label of horizontal axis:
#' character string
#' 
#' @return
#' The function plots the arguments.
#' 
#' @references
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' @seealso
#' This is an \code{\link{internal}} function. The user functions
#' are \code{\link{cursus}}, \code{\link{omnibus}},
#' and \code{\link{proprius}}.
#' 
#' @examples
#' 
#' # simulate influences
#' set.seed(1)
#' u <- rchisq(n=100,df=2)
#' 
#' # influence plot
#' upper <- rep(qchisq(p=0.95,df=2),times=100)
#' intern.plot(u,upper)
#' 
intern.plot <- function(u, upper = NULL, xlab = "indices") {
    lwd <- log(1000/length(u))
    lwd <- max(lwd, 0.1)
    lwd <- min(lwd, 5)
    n <- length(u)
    # min <- floor(min(u,upper)) max <- ceiling(max(u,upper))
    min <- min(u, upper)
    max <- max(u, upper)
    if (is.null(upper)) {
        col <- ifelse(u >= 0, "black", "grey")
    } else {
        col <- ifelse(u > upper, "black", "black")
    }
    graphics::par(mar = c(5, 4, 1, 1))
    graphics::plot.new()
    graphics::plot.window(xlim = c(1, n), ylim = c(min - abs(0.1 * min), max + abs(0.1 * 
        max)))
    graphics::box()
    graphics::abline(a = 0, b = 0, lty = 2)
    for (i in 1:n) {
        graphics::segments(x0 = i, y0 = 0, x1 = i, y1 = u[i], col = col[i], lwd = lwd)
    }
    if (is.null(upper)) {
        cond <- TRUE
    } else {
        cond <- ifelse(u > upper, TRUE, FALSE)
        h <- c(0, 0, rep(1:n, each = 2)) + 0.5
        v <- c(min - abs(max), rep(upper, each = 2), min - abs(max))
        graphics::polygon(x = h, y = v, density = 40, col = "grey")
    }
    if (n <= 25) {
        graphics::axis(side = 1, at = (1:n)[!cond], labels = names(u)[!cond], las = 2, 
            cex.axis = 0.8, col.axis = "grey")
        graphics::axis(side = 1, at = (1:n)[cond], labels = names(u)[cond], las = 2, 
            cex.axis = 0.8, col.axis = "black")
        graphics::axis(side = 2)
        graphics::title(ylab = "contribution")
    } else {
        graphics::axis(side = 1)
        graphics::axis(side = 2)
        graphics::title(ylab = "contribution", xlab = xlab)
    }
}


############### Communication ###############

#' Internal function
#' 
#' Communicates between \code{\link{cursus}}
#' and \code{\link{omnibus}}.
#'
#' @export 
#' @keywords misc
#' 
#' @inheritParams intern.chromo
#' @inheritParams cursus
#' @inheritParams omnibus
#' 
#' @param i index
#' @param Ystart location (or start location)
#' @param Yend location (or end location)
#' 
#' @return
#' The function returns a dataframe,
#' with the p-value in the first column,
#' and the test statistic in the second column.
#' 
#' @references
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes.
#' Testing for association between RNA-Seq
#' and high-dimensional data. Manuscript in preparation.
#' \link[=globalSeq]{abstract}.
#' 
#' @seealso
#' 
#' This is an \code{\link{internal}} function. The user functions
#' are \code{\link{cursus}}, \code{\link{omnibus}},
#' and \code{\link{proprius}}.
#' 
#' @examples
#' # simulate high-dimensional data
#' n <- 30
#' q <- 10
#' p <- 100
#' set.seed(1)
#' Y <- matrix(rnbinom(q*n,mu=10,
#'     size=1/0.25),nrow=q,ncol=n)
#' X <- matrix(rnorm(p*n),nrow=p,ncol=n)
#' Yloc <- seq(0,1,length.out=q)
#' Xloc <- seq(0,1,length.out=p)
#' window <- 1
#' 
#' # hypothesis testing
#' cursus(Y,Yloc,X,Xloc,window)
#' 
#' @usage 
#' intern.select(i, Y, Ystart, Yend, X, Xloc,
#'              window, offset, group,
#'              perm, phi, kind)
#' 
intern.select <- function(i, Y, Ystart, Yend, X, Xloc, window, offset, 
    group, perm, phi, kind) {
    Y <- intern.matrix(Y)
    if (class(X) == "matrix" | class(X) == "data.frame") {
        sel <- Ystart[i] - window <= Xloc & Xloc <= Yend[i] + window
        y <- Y[i, ]
        Xsel <- t(X[sel, , drop = FALSE])
        if (nrow(Xsel) == 0) {
            out <- NA
        } else {
            out <- omnibus(y = y, X = Xsel, offset = offset, group = group, 
                perm = perm, phi = phi, kind = kind)
        }
    } else {
        Xsel <- list()
        for (j in 1:length(X)) {
            sel <- Ystart[i] - window[[j]] <= Xloc[[j]] & Xloc[[j]] <= 
                Yend[i] + window[[j]]
            y <- Y[i, ]
            Xsel[[j]] <- t(X[[j]][sel, , drop = FALSE])
        }
        if (any(sapply(X, nrow) == 0)) {
            NA
        } else {
            out <- omnibus(y = y, X = Xsel, offset = offset, group = group, 
                perm = perm, phi = phi, kind = kind)
            # instead of out was temp, out <- c(temp$joint,temp$single) new
            # names(out) <- c('joint',paste('single',
            # 1:length(temp$single),sep='')) # new
        }
    }
    out
}

#' Internal function
#' 
#' Chromosome-wide analysis
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams cursus
#' @inheritParams omnibus
#' 
#' @param Y
#' RNA-Seq data\strong{:}
#' numeric matrix with \code{q} rows (genes)
#' and \code{n} columns (samples);
#' or a SummarizedExperiment object
#' @param Ystart
#' start location of genes\strong{:}
#' numeric vector of length \code{q}
#' @param Yend
#' end location of genes\strong{:}
#' NULL or numeric vector of length \code{q}
#' @param X
#' genomic profile\strong{:}
#' numeric matrix with \code{p} rows (covariates)
#' and \code{n} columns (samples)
#' @param Xloc
#' location covariates\strong{:}
#' numeric vector of length \code{p}
#' @param window
#' maximum distance\strong{:}
#' non-negative real number
#' 
#' @return
#' 
#' The function returns a dataframe,
#' with the p-value in the first column,
#' and the test statistic in the second column.
#' 
#' @examples
#' # simulate high-dimensional data
#' n <- 30
#' q <- 10
#' p <- 100
#' set.seed(1)
#' Y <- matrix(rnbinom(q*n,mu=10,
#'     size=1/0.25),nrow=q,ncol=n)
#' X <- matrix(rnorm(p*n),nrow=p,ncol=n)
#' Yloc <- seq(0,1,length.out=q)
#' Xloc <- seq(0,1,length.out=p)
#' window <- 1
#' 
#' # hypothesis testing
#' cursus(Y,Yloc,X,Xloc,window)
#' 
#' @usage
#' intern.chromo(Y, Ystart, Yend, X, Xloc,
#'              window, offset, group, perm,
#'              nodes, disp, kind) 
#' 
intern.chromo <- function(Y, Ystart, Yend, X, Xloc, window, offset, group, 
    perm, nodes, disp, kind) {
    Y <- intern.matrix(Y)
    if (disp == FALSE) {
        phi <- 0
    } else {
        phi <- NULL
    }
    if (nodes == 1) {
        out <- sapply(1:nrow(Y), function(i) intern.select(i = i, Y = Y, 
            Ystart = Ystart, Yend = Yend, X = X, Xloc = Xloc, window = window, 
            offset = offset, group = group,
            perm = perm, phi = phi, kind = kind))
    } else {
        cluster <- parallel::makeCluster(nodes)
        parallel::clusterExport(cluster, "intern.select")
        parallel::clusterExport(cluster, c("Ystart", "Y", "Yend", "X", 
            "Xloc", "window", "offset", "group", "perm", "phi", "kind"), 
            envir = environment())
        out <- parallel::parSapply(cluster, 1:nrow(Y),
            function(i) intern.select(i = i, Y = Y, Ystart = Ystart,
            Yend = Yend, X = X, Xloc = Xloc, window = window, 
            offset = offset, group = group, perm = perm,
            phi = phi, kind = kind))
        parallel::stopCluster(cluster)
        rm(cluster)
    }
    out
}

#' Internal function
#' 
#' Convert RNA-Seq data to a numeric matrix
#' 
#' @export
#' @keywords misc
#' 
#' @inheritParams intern.chromo
#' 
#' @return
#' 
#' The function returns a matrix.
#' 
#' @examples
#' # simulate RNA-Seq data
#' Y <- matrix(rnbinom(30,mu=10,size=1/0.2),nrow=10,ncol=3)
#' rownames(Y) <- paste("gene",1:nrow(Y),sep="")
#' colnames(Y) <- paste("cell",1:ncol(Y),sep="")
#' 
#' # create data structure
#' # Z <- SummarizedExperiment::SummarizedExperiment(
#' #      S4Vectors::SimpleList(counts=Y))
#' 
#' # conversion to matrix
#' # all.equal(Y,intern.matrix(Z))
#' 
#' @usage
#' intern.matrix(Y) 
#' 
intern.matrix <- function(Y){
    if(class(Y)!="matrix"){
        if(class(Y) %in% c("RangedSummarizedExperiment","SummarizedExperiment0","SummarizedExperiment")){
            if(!is.element("SummarizedExperiment",utils::installed.packages()[,1])){
                stop("Please transform Y to a matrix, or type:
                     source(\"http://bioconductor.org/biocLite.R\")
                     BiocInstaller::biocLite(\"SummarizedExperiment\")")
            } else {
               Y <- SummarizedExperiment::assays(Y)$counts 
            }
        #} else if(class(Y)=="DGEList"){
        #    Y <- Y$counts
        #} else if(class(Y)=="data.frame"){
        #    Y <- as.matrix(Y)
        }
    }
    Y
}



# CHECKING INPUTS. DO NOT DELETE. #' Internal function #' #' This
# function tests whether the inputs satisfy some #' formal
# requirements. It is called by the functions #'
# \code{\link{omnibus}} and \code{\link{proprius}}.  #' @param y #'
# Response variable.  #' @param X #' Explanatory variables.  #' @param
# mu #' Mean parameters.  #' @param phi #' Dispersion parameter.  #'
# @param group #' Confounding variable.  #' @param offset #' Offset
# values.  #' @param it #' Number of iterations.  #' @param type #'
# Type of decomposition.  #' @keywords misc #' @export #' @examples
# #' intern.check(y=c(1,0)) #' A Rauschenberger, MA Jonker, MA van de
# Wiel, and RX Menezes.  #' Testing for association between RNA-Seq #'
# and high-dimensional data. Manuscript in preparation.  #'
# \link[=globalSeq]{abstract}..  #' \link[=globalSeq]{abstract}.
# intern.check <- function(y,X=NULL,mu=NULL,phi=NULL,group=NULL,
# offset=NULL,it=NULL,type=NULL){ ####################### ### default
# values #### ####################### if(is.null(X)){X <-
# matrix(1,nrow=length(y),ncol=100)} if(is.null(mu)){mu <-
# rep(1,length(y))} if(is.null(phi)){phi <- 1} if(is.null(group)){group
# <- rep(1,max(2,length(y)))} if(is.null(it)){it <- 1000}
# if(is.null(type)){type <- 'covariates'} if(is.null(offset)){offset <-
# rep(1,max(2,length(y)))} ######################## ### nested
# functions ### ######################## intern.nan <- function(names){
# n <- 0 for(i in 1:length(names)){ x <- eval(parse(text=names[i]))
# if(any(is.na(x))){ cat(paste('Object \'',names[i],'\' must not
# contain missing values.\n',sep='')) n <- n + 1 } } n } intern.num <-
# function(names){ n <- 0 for(i in 1:length(names)){ x <-
# eval(parse(text=names[i])) if(!is.numeric(x)){ cat(paste('Object
# \'',names[i],'\' must be numeric.\n',sep='')) n <- n + 1 } } n }
# intern.vec <- function(names){ n <- 0 for(i in 1:length(names)){ x <-
# eval(parse(text=names[i])) if(!is.vector(x)){ if(all(is.matrix(x) |
# is.data.frame(x)) & (ncol(x)==1 | nrow(x)==1)==FALSE){
# cat(paste('Object \'',names[i],'\' must be a vector.\n',sep='')) n
# <- n + 1 } } } n } intern.mat <- function(names){ n <- 0 for(i in
# 1:length(names)){ x <- eval(parse(text=names[i])) if(!is.matrix(x) &
# !is.data.frame(x)){ cat(paste('Object \'',names[i],'\' must be a
# matrix.\n',sep='')) n <- n + 1 } } n } intern.len <-
# function(ref,names){ n <- 0 for(i in 1:length(names)){ x <-
# eval(parse(text=names[i])) if(length(x)!=ref){ cat(paste('Object
# \'',names[i],'\' must have length ',ref,'.\n',sep='')) n <- n + 1
# } } n } intern.ext <- function(y,X,mu,phi,group,offset,it,type){ n <-
# 0 # response variable if(length(y) < 2){ cat('At least two samples
# are necessary:\n') cat('Vector \'y\' must be longer than one.\n')
# n <- n + 1 } if(any(y<0) | any(y==Inf)){ cat('Values of \'y\' are
# outside the negative binomial support [0,Inf).\n') n <- n + 1 } #
# covariate matrix if(class(X)=='matrix'){ if(length(y)!=nrow(X)){
# if(length(y)==ncol(X)){ cat('Maybe you just have to transpose the
# matrix \'X\',\n') cat('such that samples are in rows and
# covariates in columns.\n') } cat('The length of vector \'y\' and
# the number of rows of \'X\' must be equal.\n') n <- n + 1 } } #
# mean if(any(mu<0) | any(mu==Inf)){ cat('Values of \'mu\' are
# outside the parameter space [0,Inf).\n') n <- n + 1 } # dispersion
# if(any(phi<=0) | any(phi==Inf)){ cat('Value of \'phi\' is outside
# the parameter space (0,Inf).\n') n <- n + 1 } # confounding m <-
# sapply(unique(group),function(x) sum(x==group)) if(min(m)<2){
# cat('Object \'group\' must attribute at least two samples to each
# subgroup.\n') n <- n + 1 } # offset if(any(offset<=0) |
# any(offset==Inf)){ cat('All elements in \'group\' must be inside
# (0,Inf).\n') } # iterations if(round(it)!=it){ cat('Object \'it\'
# must be an integer.\n') n <- n + 1 } if(it < 2){ cat('Object
# \'it\' must be larger than one.\n') n <- n + 1 } # decomposition
# if(type!='samples' & type!='covariates'){ cat('Object \'type\' must
# either be \'samples\' or \'covariates\'.\n') n <- n + 1 } n }
# ################ ### checking ### ################ n <- 0 # number of
# errors n <- n + intern.nan(c('y','X','mu','phi',
# 'group','offset','it','type')) n <- n +
# intern.num(c('y','X','mu','phi','offset','it')) n <- n +
# intern.vec(c('y','mu','phi','group', 'offset','it','type')) n <- n +
# intern.mat(c('X')) n <- n + intern.len(1,c('phi','it','type')) n <- n
# + intern.len(length(y),c('mu','group','offset')) n <- n +
# intern.ext(y=y,X=X,mu=mu,phi=phi,group=group,
# offset=offset,it=it,type=type) if(n==0){ cat('Inputs accepted!\n') }
# else { if(n==1){ stop('Please correct the error mentioned
# above.',call.=FALSE) } else { stop(paste('Please correct
# the',n,'errors mentioned above.'),call.=FALSE) } } } 
