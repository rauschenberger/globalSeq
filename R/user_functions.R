############### Genome-wide analysis ###############

#' Genome-wide analysis
#' 
#' This function tests for associations between gene expression
#' or exon abundance (\code{Y})
#' and genetic or epigenetic alterations (\code{X}).
#' Using the locations of genes (\code{Yloc}),
#' and the locations of genetic
#' or epigenetic alterations (\code{Xloc}),
#' the expression of each gene is tested for associations with 
#' alterations on the same chromosome that are closer to the gene
#' than a given distance (\code{window}).
#' 
#' @export
#' @keywords methods
#' 
#' @inheritParams omnibus
#' @param Y
#' \strong{RNA-Seq data}\strong{:}
#' numeric matrix with \code{q} rows (genes)
#' and \code{n} columns (samples);
#' or a SummarizedExperiment object
#' @param Yloc
#' \strong{location RNA-Seq}\strong{:}
#' numeric vector of length \code{q}
#' (point location)\strong{;}
#' numeric matrix with \code{q} rows
#' and two columns (start and end locations)
#' @param X
#' \strong{genomic profile}\strong{:}
#' numeric matrix with \code{p} rows (covariates)
#' and \code{n} columns (samples)
#' @param Xloc
#' \strong{location covariates}\strong{:}
#' numeric vector of length \code{p}
#' @param Ychr
#' chromosome RNA-Seq\strong{:}
#' factor of length \code{q}
#' @param Xchr
#' chromosome covariates\strong{:}
#' factor of length \code{p}
#' @param window
#' \strong{maximum distance}\strong{:}
#' non-negative real number
#' @param nodes
#' number of cluster nodes for parallel computation
#' @param phi
#' dispersion parameters\strong{:} vector of length \code{q}
#' 
#' @details
#' Note that \code{Yloc}, \code{Xloc} and \code{window} must
#' be given in the same unit, usually in base pairs.
#' If \code{Yloc} indicates interval \strong{locations},
#' and \code{window} is zero,
#' then only covariates between the start and end location
#' of the gene are of interest.
#' Typically \code{window} is larger than one million base pairs.
#' 
#' If \code{Y} and \code{X} include data from a single chromosome,
#' \code{Ychr} and \code{Xchr} are redundant.
#' If \code{Y} or \code{X} include data
#' from \strong{multiple chromosomes},
#' \code{Ychr} and \code{Xchr} should be specified
#' in order to prevent confusion between chromosomes. 
#' 
#' For the simultaneous analysis of
#' \strong{multiple genomic profiles}
#' \code{X} should be a list of numeric matrices with
#' \code{n} columns (samples),
#' \code{Xloc} a list of numeric vectors,
#' and \code{window} a list of non-negative real numbers.
#' If provided, \code{Xchr} should be alist of of numeric vectors.
#' 
#' The \code{offset} is meant to account for
#' different \strong{libary sizes}.
#' By default the \code{offset} is calculated based on \code{Y}.
#' Different library sizes can be ignored by 
#' setting the \code{offset} to \code{rep(1,n)}.
#' 
#' The user can provide the \strong{confounding} variable \code{group}.
#' Note that each level of \code{group} must appear at least twice
#' in order to allow stratified permutations.
#' 
#' Efficient alternatives to classical \strong{permutation} (\code{kind=1})
#' are the method of control variates (\code{kind=0})
#' and permutation in chunks (0 < \code{kind} < 1)
#' \link[=intern.crude]{details}.
#' 
#' @return
#'
#' The function returns a dataframe,
#' with the p-values in the first row
#' and the test statistics in the second row.
#' 
#' @references
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016).
#' "Testing for association between RNA-Seq and high-dimensional data",
#' \emph{BMC Bioinformatics}. 17:118.
#' \href{http://dx.doi.org/10.1186/s12859-016-0961-5}{html}
#' \href{http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf}{pdf}
#' (open access)
#' 
#' RX Menezes, M Boetzer, M Sieswerda, GJB van Ommen, and JM Boer (2009).
#' "Integrated analysis of DNA copy number
#' and gene expression microarray data using gene sets",
#' \emph{BMC Bioinformatics}. 10:203.
#' \href{http://dx.doi.org/10.1186/1471-2105-10-203}{html}
#' \href{http://www.biomedcentral.com/content/pdf/1471-2105-10-203.pdf}{pdf}
#' (open access)
#' 
#' @seealso The function \code{\link{omnibus}} tests for associations
#' between an overdispersed response variable
#' and a high-dimensional covariate set.
#' The function \code{\link{proprius}} calculates the contributions
#' of individual samples or covariates to the test statistic.
#' All other function of the R package
#' \code{\link{globalSeq}} are \code{\link{internal}}.
#' 
#' @examples
#' # simulate high-dimensional data
#' n <- 30; q <- 10; p <- 100
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
#' cursus(Y, Yloc, X, Xloc, window,
#'         Ychr = NULL, Xchr = NULL,
#'         offset = NULL, group = NULL,
#'         perm = 1000, nodes = 2,
#'         phi = NULL, kind = 0.01)
#'           
cursus <- function(Y, Yloc, X, Xloc, window, Ychr = NULL, Xchr = NULL, 
    offset = NULL, group = NULL, perm = 1000,
    nodes = 2, phi = NULL, kind = 0.01) { # instead of phi = NULL was disp = TRUE
    Y <- globalSeq::intern.matrix(Y)
    if (is.vector(Yloc)) {
        Ystart <- Yend <- Yloc
    } else {
        Ystart <- Yloc[, 1]
        Yend <- Yloc[, 2]
    }
    if (length(perm) == 1) {
        perm <- globalSeq::intern.permu(n = ncol(Y), it = perm - 1, group = group, kind=kind) # new: kind=kind
    }
    if (is.null(offset)) {
        ls <- colSums(Y)
        gm <- exp(1/length(ls) * sum(log(ls)))
        offset <- ls/gm
    }
    # #  alternative with edgeR
    # offset <- edgeR::calcNormFactors(Y)
    # utils::capture.output({phi <- edgeR::estimateDisp(edgeR::DGEList(counts=Y,norm.factors=offset))})
    # phi <- phi$tagwise.dispersion
    if (is.null(Ychr) | is.null(Xchr)) {
        ########## Analysing a single chromosome ##########
        cat("Analysing a single chromosome.\n")
        out <- globalSeq::intern.chromo(Y = Y, Ystart = Ystart, Yend = Yend, X = X, 
            Xloc = Xloc, window = window, offset = offset, group = group, 
            perm = perm, nodes = nodes, phi = phi, kind = kind) # was disp = disp instead of phi = phi
    } else if (!is.null(Ychr) & !is.null(Xchr)) {
        ########## Analysing multiple chromosomes ##########
        out <- list()
        chr <- unique(Ychr)
        cat("Analysing multiple chromosomes:\n")
        for (i in 1:length(chr)) {
            cat(paste("chromosome", chr[i],"\n"))
            ###### START NEW ########
            if (class(X) == "list") {
                Xpass <- lapply(1:length(X), function(j) X[[j]][Xchr[[j]] == 
                  chr[i], , drop = FALSE])
            } else {
                Xpass <- X[Xchr == chr[i], , drop = FALSE]
            }
            if (class(Xloc) == "list") {
                Xlocpass <- lapply(1:length(Xloc),
                function(j) Xloc[[j]][Xchr[[j]] == chr[i]])
            } else {
                Xlocpass <- Xloc[Xchr == chr[i]]
            }
            out[[i]] <- globalSeq::intern.chromo(Y = Y[Ychr == chr[i], , drop = FALSE], 
                Ystart = Ystart[Ychr == chr[i]], Yend = Yend[Ychr == chr[i]], 
                X = Xpass, Xloc = Xlocpass, window = window, perm = perm, 
                offset = offset, group = group, nodes = nodes, phi = phi, # was disp = disp instead of phi = phi
                kind = kind)
            #### END NEW ######## The following line was active !  out[[i]] <-
            #### globalSeq::intern.chromo(Y=Y[Ychr==chr[i],,drop=FALSE],
            #### Ystart=Ystart[Ychr==chr[i]],Yend=Yend[Ychr==chr[i]],
            #### X=X[Xchr==chr[i],,drop=FALSE],Xloc=Xloc[Xchr==chr[i]],
            #### window=window,perm=perm,offset=offset,group=group,
            #### nodes=nodes,disp=disp,kind=kind)
        }
        # if(class(out[[1]])=="numeric"){out <- unlist(out)}else{
        out <- do.call(cbind, out)
        # }
    } else {
        ### Ambiguous arguments ###
        stop("Please provide either\n both or none of \"Ychr\" and \"Xchr\".")
    }
    # everything from here is new (before only: out)
    col <- matrix(NA, nrow = ncol(out), ncol = nrow(out))
    for (i in 1:nrow(out)) {
        col[, i] <- unlist(out[i, ])
    }
    colnames(col) <- rownames(out)
    rownames(col) <- rownames(Y)
    as.data.frame(col)
}






############### OMNIBUS ###############

#' Omnibus test
#' 
#' Test of association between a count response and
#' one or more covariate sets.
#' This test may be conceptualised as
#' a test of overall significance in regression analysis,
#' where the response variable is overdispersed, and where
#' the number of explanatory variables (\code{p})
#' exceeds the sample size (\code{n}).
#' The negative binomial distribution accounts for overdispersion
#' and a random effect model accounts for high dimensionality
#' (\code{p}>>\code{n}).
#' 
#' @export
#' @keywords methods
#' 
#' @param y
#' \strong{response variable}\strong{:}
#' numeric vector of length \code{n}
#' @param X
#' \strong{one covariate set}\strong{:}
#' numeric matrix with \code{n} rows (samples)
#' and \code{p} columns (covariates);
#' \cr \strong{multiple covariate sets}\strong{:}
#' list of numeric matrices with \code{n} rows (samples)
#' @param offset
#' numeric vector of length \code{n}
#' @param group
#' confounding variable\strong{:}
#' factor of length \code{n}
#' @param perm
#' number of iterations\strong{:}
#' positive integer
#' @param mu
#' mean parameters\strong{:}
#' numeric vector of length \code{1} or \code{n}
#' @param phi
#' dispersion parameter\strong{:}
#' non-negative real number
#' @param kind
#' computation \strong{:}
#' number between 0 and 1
#' 
#' @details
#' 
#' The user can provide a common \code{mu} for all samples
#' or sample-specific \code{mu}, and a common \code{phi}.
#' Setting \code{phi} equal to zero is equivalent
#' to using the Poisson model.
#' If \code{mu} is missing, then \code{mu} is estimated from \code{y}.
#' If \code{phi} is missing, then \code{mu} and \code{phi}
#' are estimated from \code{y}.
#' The \code{offset} is only taken into account
#' for estimating \code{mu} or \code{phi}.
#' By default the offset is \code{rep(1,n)}.
#' 
#' The user can provide the \strong{confounding} variable \code{group}.
#' Note that each level of \code{group} must appear at least twice
#' in order to allow stratified permutations.
#' 
#' Efficient alternatives to classical \strong{permutation} (\code{kind=1})
#' are the method of control variates (\code{kind=0})
#' and permutation in chunks (0 < \code{kind} < 1)
#' \link[=intern.crude]{details}.
#' 
#' @return
#' 
#' The function returns a dataframe,
#' with the p-value in the first column,
#' and the test statistic in the second column.
#' 
#' @references
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016).
#' "Testing for association between RNA-Seq and high-dimensional data",
#' \emph{BMC Bioinformatics}. 17:118.
#' \href{http://dx.doi.org/10.1186/s12859-016-0961-5}{html}
#' \href{http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf}{pdf}
#' (open access)
#' 
#' RX Menezes, L Mohammadi, JJ Goeman, and JM Boer (2016).
#' "Analysing multiple types of molecular profiles simultaneously:
#' connecting the needles in the haystack",
#' \emph{BMC Bioinformatics}. 17:77.
#' \href{http://dx.doi.org/10.1186/s12859-016-0926-8}{html}
#' \href{http://www.biomedcentral.com/content/pdf/s12859-016-0926-8.pdf}{pdf}
#' (open access)
#' 
#' S le Cessie, and HC van Houwelingen (1995).
#' "Testing the fit of a regression model 
#' via score tests in random effects models",
#' \emph{Biometrics}. 51:600-614.
#' \href{http://dx.doi.org/10.2307/2532948}{html}
#' \href{http://www.jstor.org/stable/pdf/2532948.pdf?acceptTC=true}{pdf}
#' (restricted access)
#' 
#' @seealso
#' 
#' The function \code{\link{proprius}} calculates
#' the contributions of individual samples or covariates
#' to the test statistic.
#' The function \code{\link{cursus}} tests for association
#' between RNA-Seq and local genetic or epigenetic alternations
#' across the whole genome.
#' All other functions of the R package \code{\link{globalSeq}}
#' are \code{\link{internal}}.
#' 
#' @examples
#' 
#' # simulate high-dimensional data
#' n <- 30; p <- 100
#' y <- rnbinom(n,mu=10,size=1/0.25)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#'
#' # hypothesis testing
#' omnibus(y,X)
#' 
#' @usage
#' omnibus(y, X, offset = NULL, group = NULL,
#'         mu = NULL, phi = NULL,
#'         perm = 1000, kind = 1)
#'         
omnibus <- function(y, X, offset = NULL, group = NULL, mu = NULL, phi = NULL, 
    perm = 1000, kind = 1) {
    ########## initialisation ##########
    n <- length(y)  # number of samples
    if (is.null(mu) | is.null(phi)) {
        est <- globalSeq::intern.estim(y = y, offset = offset)
        mu <- est$mu
        if (is.null(phi)) {
            phi <- est$phi
        }
    } else {
        if (length(mu) == 1) {
            mu <- rep(mu, n)
        }
    }
    if (length(perm) == 1) {
        it <- perm
        perm <- globalSeq::intern.permu(n = n, it = perm - 1, group = group, kind = kind) # new: kind=kind
        it <- ncol(perm)
    } else {
        it <- ncol(perm)
    }
    ########## testing one covariate set ##########
    if (class(X) == "matrix" | class(X) == "data.frame") {
        if (kind == 1) {
            globalSeq::intern.crude(y = y, X = X, mu = mu, phi = phi, perm = perm)
        } else if (kind == 0) {
            globalSeq::intern.conva(y = y, X = X, mu = mu, phi = phi, perm = perm,
                offset = offset)
        } else if (kind > 0 & kind < 1) {
            globalSeq::intern.focus(y = y, X = X, mu = mu, phi = phi, perm = perm, 
                focus = kind)
        } else {
            stop("Argument \"kind\" must be within 0 and 1.")
        }
    } else {
        ########## testing multiple sets ##########
        sets <- X
        k <- length(sets)  # number of sets
        # # # # # # # # # # # separate testing # # # # # # # # # # #
        single <- covs <- rep(NA, times = k)  # single p-values                        # new: covs <- 
        sim <- matrix(NA, nrow = k, ncol = it)  # simulated test statistics
        for (i in 1:k) {
            X <- sets[[i]]
            R <- X %*% t(X)/ncol(X)
            temp <- apply(perm, 2, function(perm) globalSeq::intern.score(y = y[perm], 
                R = R, mu = mu[perm], phi = phi))
            single[i] <- sum(temp >= temp[1])/it  # single p-values
            covs[i] <- ncol(X)                                                         # new
            sim[i, ] <- temp
        }
        # # # # # # # # # joint testing # # # # # # # # #
        sim_mu <- rowMeans(sim)
        sim_sd <- apply(sim, 1, stats::sd)
        com <- rep(NA, times = it)  # simulated test statistics
        for (i in 1:it) {
            score <- sim[, i, drop = FALSE]
            com[i] <- sum((score - sim_mu)/sim_sd)  # standardisation
        }
        joint <- sum(com >= com[1])/it  # joint p-value
        data.frame(joint = joint, teststat = com[1], single = matrix(single, 
            nrow = 1),covs=matrix(covs,nrow=1))                                                              # new: p=p
    }
}



############### PROPRIUS ###############

#' Decomposition
#' 
#' Even though the function \code{\link{omnibus}} tests
#' a single hypothesis on a whole covariate set,
#' this function allows to calculate
#' the individual contributions of \code{n} samples or
#' \code{p} covariates to the test statistic.
#' 
#' @export
#' @keywords methods
#' 
#' @inheritParams omnibus
#' @param X \strong{covariate set}\strong{:}
#' numeric matrix with \code{n} rows (samples) 
#' and \code{p} columns (covariates)
#' @param type
#' character '\strong{covariates}' or '\strong{samples}'
#' @param plot
#' plot of results\strong{:} logical
#' @param alpha
#' significance level\strong{:} real number between 0 and 1
#' 
#' @details
#' 
#' The user can provide a common \code{mu} for all samples
#' or sample-specific \code{mu}, and a common \code{phi}.
#' Setting \code{phi} equal to zero
#' is equivalent to using the Poisson model.
#' If \code{mu} is missing, then \code{mu} is estimated from \code{y}.
#' If \code{phi} is missing, then \code{mu} and \code{phi}
#' are estimated from \code{y}.
#' The \code{offset} is only taken into account
#' for estimating \code{mu} or \code{phi}.
#' 
#' The user can provide the confounding variable \code{group}.
#' Note that each level of \code{group} must appear at least twice
#' in order to allow stratified permutations.
#' 
#' @return
#' 
#' If \code{alpha=NULL}, then the function returns a numeric vector,
#' and else a list of numeric vectors.
#' 
#' @references
#' 
#' A Rauschenberger, MA Jonker, MA van de Wiel, and RX Menezes (2016).
#' "Testing for association between RNA-Seq and high-dimensional data",
#' \emph{BMC Bioinformatics}. 17:118.
#' \href{http://dx.doi.org/10.1186/s12859-016-0961-5}{html}
#' \href{http://www.biomedcentral.com/content/pdf/s12859-016-0961-5.pdf}{pdf}
#' (open access)
#' 
#' JJ Goeman, SA van de Geer, F de Kort, and HC van Houwelingen (2004).
#' "A global test for groups of genes:
#' testing association with a clinical outcome",
#' \emph{Bioinformatics}. 20:93-99.
#' \href{http://dx.doi.org/10.1093/bioinformatics/btg382}{html}
#' \href{http://bioinformatics.oxfordjournals.org/content/20/1/93.full.pdf}{pdf}
#' (open access)
#' 
#' @seealso
#' 
#' The function \code{\link{omnibus}} tests for associations
#' between an overdispersed response variable and a high-dimensional
#' covariate set.
#' The function \code{\link{cursus}} tests for association
#' between RNA-Seq and local genetic or epigenetic alternations
#' across the whole genome.
#' All other functions of the R package \code{\link{globalSeq}}
#' are \code{\link{internal}}.
#' 
#' @examples
#' 
#' # simulate high-dimensional data
#' n <- 30; p <- 100
#' y <- rnbinom(n,mu=10,size=1/0.25)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#'
#' # decomposition
#' proprius(y,X,type="samples")
#' proprius(y,X,type="covariates")
#' 
#' @usage
#' 
#' proprius(y, X, type, offset = NULL, group = NULL,
#'         mu = NULL, phi = NULL,
#'         alpha = NULL, perm = 1000, plot = TRUE)
#' 
proprius <- function(y, X, type, offset = NULL, group = NULL, mu = NULL, 
    phi = NULL, alpha = NULL, perm = 1000, plot = TRUE) {
    if (is.null(mu) | is.null(phi)) {
        est <- globalSeq::intern.estim(y = y, offset = offset)
        mu <- est$mu
        if (is.null(phi)) {
            phi <- est$phi
        }
    } else {
        if (length(mu) == 1) {
            mu <- rep(mu, length(y))
        }
    }
    ### decomposition ###
    if (is.null(alpha)) {
        if (type == "samples") {
            u <- globalSeq::intern.sam(y, X, mu, phi)
        }
        if (type == "covariates") {
            u <- globalSeq::intern.cov(y, X, mu, phi)
        }
        upper <- NULL
    } else {
        if (length(perm) == 1) 
            {
                perm <- globalSeq::intern.permu(n = length(y),
                it = perm - 1, group = group, kind = Inf) # new: kind= Inf
            }  # new
        # perm <- globalSeq::intern.permu(length(y),1000-1,group) # old
        if (type == "covariates") {
            sim <- apply(perm, 2, function(perm) globalSeq::intern.cov(y = y[perm], 
                X = X, mu = mu[perm], phi = phi))
        }
        if (type == "samples") {
            sim <- apply(perm, 2, function(perm) globalSeq::intern.sam(y = y[perm], 
                X = X, mu = mu[perm], phi = phi))
        }
        u <- sim[, 1]
        upper <- apply(sim, 1, function(x) stats::quantile(x, p = 1 - alpha))
    }
    ### plot ###
    if (plot == "TRUE") {
        globalSeq::intern.plot(u = u, upper = upper, xlab = paste("indices of", type))
    }
    ### output ###
    if (is.null(alpha)) {
        u
    } else {
        list(u = u, upper = upper)
    }
} 
