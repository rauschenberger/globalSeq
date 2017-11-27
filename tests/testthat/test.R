# Avoid parallel processing by setting the argument nodes equal to one.

testthat::test_that("estimates are inside parameter space", {
    estimates <- globalSeq::intern.estim(rnbinom(n = 1000, mu = 0.1, size = 1/0.01))
    testthat::expect_more_than(min(estimates$mu), -10^(-10))
    testthat::expect_more_than(estimates$phi, -10^(-10))
})

testthat::test_that("permutations are without repetitions", {
    perms <- globalSeq::intern.permu(n = 3, it = 1000, group = NULL,kind=rnorm(1))
    testthat::expect_less_than(ncol(perms), factorial(3) + 1)
    testthat::expect_equal(perms[, 1], 1:3)
    testthat::expect_false(any(duplicated(t(perms))))
})

testthat::test_that("permutations are only within groups \n and not between groups", 
    {
        group <- c(1, 2, 2, 1, 2)
        perm <- globalSeq::intern.permu(n = length(group), it = 1000, group = group,kind=rnorm(1))
        obs <- apply(perm, 1, function(x) sort(unique(x)))
        testthat::expect_identical(length(unique(obs)), length(unique(group)))
        testthat::expect_true(identical(obs[[1]], obs[[4]]))
        testthat::expect_true(identical(obs[[2]], obs[[3]]))
        testthat::expect_true(identical(obs[[2]], obs[[5]]))
        testthat::expect_false(identical(obs[[1]], obs[[2]]))
    })

testthat::test_that("score is calculated correctly", {
    # simulate data
    y <- rnbinom(30, mu = 10, size = 1/0.25)
    X <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
    # calculate scores
    R <- X %*% t(X)/ncol(X)
    score.1 <- 0.5 * y %*% R %*% y
    score.2 <- globalSeq::intern.score(y, R, mu = 0, phi = rnorm(1))
    # compare
    testthat::expect_equal(score.1, score.2)
})

testthat::test_that("test statistic equals sum over individual contributions \n (using internal functions)", 
    {
        # simulate data
        y <- rnbinom(30, mu = 10, size = 1/0.25)
        X <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
        # calculate test statistic
        R <- X %*% t(X)/ncol(X)
        mu <- rep(mean(y), length(y))
        phi <- (var(y) - mean(y))/mean(y)^2
        tstat <- globalSeq::intern.score(y, R, mu, phi)
        # calculate individual contributions
        cont.sam <- globalSeq::intern.sam(y, X, mu, phi)
        cont.cov <- globalSeq::intern.cov(y, X, mu, phi)
        # compare
        testthat::expect_equal(sum(cont.sam), as.numeric(tstat))
        testthat::expect_equal(sum(cont.cov), as.numeric(tstat))
    })

testthat::test_that("test statistic equals sum over individual contributions \n (using interface functions)", 
    {
        # simulate data
        Y <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
        X <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
        # calculate test statistic and individual contributions
        tstat.1 <- globalSeq::cursus(Y = Y, Yloc = 1, X = X, Xloc = rep(1, 
            nrow(X)), window = 0, offset = rep(1, ncol(Y)), nodes = 1)$teststat
        tstat.2 <- globalSeq::omnibus(y = Y, X = t(X))$teststat
        cont.sam <- globalSeq::proprius(y = Y, X = t(X), type = "samples")
        cont.cov <- globalSeq::proprius(y = Y, X = t(X), type = "covariates")
        # compare
        testthat::expect_equal(tstat.1, sum(cont.sam))
        testthat::expect_equal(tstat.1, sum(cont.cov))
        testthat::expect_equal(tstat.1, tstat.1)
    })

testthat::test_that("trivial p-values are given correctly", {
    # initialise
    y <- rnbinom(30, mu = 10, size = 1/0.25)
    offset <- rpois(30, lambda=100)
    X <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
    mu <- rep(mean(y), 30)
    phi <- (var(y) - mu)/mu^2
    set.seed(1)
    perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=runif(1))
    # If there is no differential expression, the p-value should be equal
    # to one.
    y0 <- rep(0, length(y))
    crude <- globalSeq::intern.crude(y0, X, mu, phi, perm)
    focus <- globalSeq::intern.focus(y0, X, mu, phi, perm, focus = 0.01)
    conva <- globalSeq::intern.conva(y0, X, mu, phi, perm, offset = offset)
    testthat::expect_identical(1, crude$pvalue)
    testthat::expect_identical(1, focus$pvalue)
    testthat::expect_identical(1, conva$pvalue)
    # If there are no covariates to test, no p-value should be given.
    X0 <- matrix(nrow = nrow(X), ncol = 0)
    crude <- globalSeq::intern.crude(y, X0, mu, phi, perm)
    focus <- globalSeq::intern.focus(y, X0, mu, phi, perm, focus = 0.01)
    conva <- globalSeq::intern.conva(y, X0, mu, phi, perm, offset = offset)
    testthat::expect_identical(NA, crude$pvalue)
    testthat::expect_identical(NA, focus$pvalue)
    testthat::expect_identical(NA, conva$pvalue)
})

testthat::test_that("it does not matter whether a common mean is supplied \n as a scalar or as a vector", 
    {
        # simulate data
        y <- rnbinom(30, mu = 10, size = 1/0.25)
        X <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
        # mean scalar versus mean vector
        mu_scalar <- mean(y)
        mu_vector <- rep(mean(y), length(y))
        # provide phi such that mu is not estimated inside function
        phi <- (var(y) - mu_scalar)/mu_scalar^2
        # test statistic
        ts1 <- globalSeq::omnibus(y = y, mu = mu_scalar, phi = phi, X = t(X))$teststat
        ts2 <- globalSeq::omnibus(y = y, mu = mu_vector, phi = phi, X = t(X))$teststat
        # individual contributions
        ts3 <- sum(globalSeq::proprius(y = y, mu = mu_scalar, phi = phi, X = t(X), 
            type = "samples"))
        ts4 <- sum(globalSeq::proprius(y = y, mu = mu_vector, phi = phi, X = t(X), 
            type = "samples"))
        ts5 <- sum(globalSeq::proprius(y = y, mu = mu_scalar, phi = phi, X = t(X), 
            type = "covariates"))
        ts6 <- sum(globalSeq::proprius(y = y, mu = mu_vector, phi = phi, X = t(X), 
            type = "covariates"))
        # compare
        testthat::expect_true(all.equal(ts1, ts2, ts3, ts4, ts5, ts6))
    })

testthat::test_that("different permutation tests do not contradict each other \n (using internal functions)", 
    {
        # simulate data
        y <- rnbinom(30, mu = 10, size = 1/0.25)
        X <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
        mu <- rep(mean(y), 30)
        phi <- (var(y) - mu)/mu^2
        # different permutation tests
        perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=rnorm(1))
        crude <- globalSeq::intern.crude(y, X, mu, phi, perm)
        focus <- globalSeq::intern.focus(y, X, mu, phi, perm, focus = 0.01)
        conva <- globalSeq::intern.conva(y, X, mu, phi, perm, NULL)
        # various checks
        testthat::expect_less_than(crude$pvalue, focus$pvalue + 10^(-10))
        testthat::expect_less_than(min(crude$pvalue, focus$pvalue, conva$pvale), 
            1 + 10^(-10))
        testthat::expect_more_than(min(crude$pvalue, focus$pvalue, conva$pvale), 
            0 - 10^(-10))
    })

testthat::test_that("different permutation tests do not contradict each other \n (using interface functions)", 
    {
        # simulate data
        y <- rnbinom(30, mu = 10, size = 1/0.25)
        X <- matrix(rnorm(30 * 100), nrow = 30, ncol = 100)
        mu <- rep(mean(y), 30)
        phi <- (var(y) - mu)/mu^2
        # different permutation tests
        perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=rnorm(1))
        crude <- globalSeq::omnibus(y, X, mu = mu, phi = phi, perm = perm, 
            kind = 1)
        focus <- globalSeq::omnibus(y, X, mu = mu, phi = phi, perm = perm, 
            kind = 0.01)
        conva <- globalSeq::omnibus(y, X, mu = mu, phi = phi, perm = perm, 
            kind = 0)
        # various tests
        testthat::expect_less_than(crude$pvalue, focus$pvalue + 10^(-10))
        testthat::expect_less_than(min(crude$pvalue, focus$pvalue, conva$pvale), 
            1 + 10^(-10))
        testthat::expect_more_than(min(crude$pvalue, focus$pvalue, conva$pvale), 
            0 - 10^(-10))
        testthat::expect_error(globalSeq::omnibus(y, X, perm = perm, kind = -abs(rnorm(1))))
        testthat::expect_error(globalSeq::omnibus(y, X, perm = perm, kind = 1 + 
            abs(rnorm(1))))
    })

testthat::test_that("significance levels are intepreted correctly", {
    # simulate data
    y <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    X <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    # different significance levels
    perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=rnorm(1))
    alpha0.00 <- globalSeq::proprius(y, t(X), type = "covariates", perm = perm, 
        alpha = 0)
    alpha0.01 <- globalSeq::proprius(y, t(X), type = "covariates", perm = perm, 
        alpha = 0.01)
    alpha0.05 <- globalSeq::proprius(y, t(X), type = "covariates", perm = perm, 
        alpha = 0.05)
    alpha1.00 <- globalSeq::proprius(y, t(X), type = "covariates", perm = perm, 
        alpha = 1)
    # contributions must be the same
    testthat::expect_equal(alpha0.00$u, alpha0.01$u)
    testthat::expect_equal(alpha0.00$u, alpha0.05$u)
    testthat::expect_equal(alpha0.00$u, alpha1.00$u)
    # observations must be within extreme boundary
    testthat::expect_true(all(alpha0.00$u <= alpha0.00$upper))
    testthat::expect_true(all(alpha0.00$u >= alpha1.00$upper))
    # ordering
    testthat::expect_true(all(alpha0.00$upper >= alpha0.01$upper))
    testthat::expect_true(all(alpha0.01$upper >= alpha0.05$upper))
    testthat::expect_true(all(alpha0.05$upper >= alpha1.00$upper))
})


testthat::test_that("permutating inside or outside omnibus is equivalent",{
    # simulate data
    y <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    X <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    # permuting outside
    set.seed(1)
    perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=rnorm(1))
    res.out <- globalSeq::omnibus(y,t(X),perm=perm)
    # permuting inside
    set.seed(1)
    res.in <- globalSeq::omnibus(y,t(X),perm=101)
    # compare
    testthat::expect_identical(res.out,res.in)
})

testthat::test_that("permutating inside or outside proprius is equivalent",{
    # simulate data
    y <- matrix(rnbinom(24, mu = 10, size = 1/0.25), nrow = 1)
    X <- matrix(rnorm(100 * 24), nrow = 100, ncol = 24)
    # permuting outside
    set.seed(1)
    perm <- globalSeq::intern.permu(n = 24, it = 100, group = NULL,kind=rnorm(1))
    res.out <- globalSeq::proprius(y,t(X),type="samples",perm=perm,alpha=0.01)$upper
    # permuting inside
    set.seed(1)
    res.in <- globalSeq::proprius(y,t(X),type="samples",perm=101,alpha=0.01)$upper
    # compare
    testthat::expect_identical(res.out,res.in)
})

testthat::test_that("multiple covariate sets can be analysed", 
    {
        # simulate data
        y <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
        X1 <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
        X2 <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
        # different tests
        perm <- globalSeq::intern.permu(n = 30, it = 100, group = NULL,kind=rnorm(1))
        single.1 <- globalSeq::omnibus(y, t(X1), perm = perm)$pvalue
        single.2 <- globalSeq::omnibus(y, t(X2), perm = perm)$pvalue
        joint1 <- globalSeq::omnibus(y, list(t(X1), t(X2)), perm = perm)
        joint2 <- globalSeq::cursus(Y = as.matrix(y, nrow = 1), Yloc = 1, 
            X = list(X1, X2), Xloc = list(1, 1), window = list(0, 0), perm = perm, 
            nodes = 1, offset = rep(1, length(y)))
        # first covariate set
        testthat::expect_identical(single.1, joint1$single.1)
        testthat::expect_identical(single.1, joint2$single.1)
        # second covariate set
        testthat::expect_identical(single.2, joint1$single.2)
        testthat::expect_identical(single.2, joint2$single.2)
    })

testthat::test_that("multiple chromosomes can be analysed", {
    # simulate data
    Y1 <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    Y2 <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    X1 <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    X2 <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    # analysing chromosomes separately or simultaneously
    chr.1 <- globalSeq::cursus(Y = Y1, Yloc = 1, X = X1, Xloc = rep(1, 
        nrow(X1)), window = 0, offset = rep(1, ncol(Y1)), nodes = 1)$teststat
    chr.2 <- globalSeq::cursus(Y = Y2, Yloc = matrix(c(1, 1), nrow = 1), 
        X = X2, Xloc = rep(1, nrow(X2)), window = 0, offset = rep(1, ncol(Y1)), 
        nodes = 1)$teststat
    joint <- globalSeq::cursus(Y = rbind(Y1, Y2), Yloc = c(1, 1), Ychr = c(1, 
        2), X = rbind(X1, X2), Xloc = rep(1, nrow(X1) + nrow(X2)), Xchr = rep(c(1, 
        2), times = c(nrow(X1), nrow(X2))), window = 0, offset = rep(1, 
        ncol(Y1)), nodes = 1)$teststat
    # compare results
    testthat::expect_identical(chr.1, joint[1])
    testthat::expect_identical(chr.2, joint[2])
})

testthat::test_that("multiple chromosomes and multiple sets can be analysed",{
    # simulate data
    Y1 <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    Y2 <- matrix(rnbinom(30, mu = 10, size = 1/0.25), nrow = 1)
    X1A <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    X1B <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    X2A <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    X2B <- matrix(rnorm(100 * 30), nrow = 100, ncol = 30)
    # analysing chromosomes separately or simultaneously
    set.seed(1)
    chr.1 <- globalSeq::cursus(Y = Y1, Yloc = 1, X = list(X1A,X1B),
                Xloc = list(rep(1, nrow(X1A)),rep(1,nrow(X1B))), window = list(0,0),
                offset = rep(1, ncol(Y1)), nodes = 1)$joint
    set.seed(1)
    chr.2 <- globalSeq::cursus(Y = Y2, Yloc = matrix(c(1, 1), nrow = 1), 
                X = list(X2A,X2B), Xloc = list(rep(1, nrow(X2A)),rep(1,nrow(X2B))),
                window = list(0,0),
                offset = rep(1, ncol(Y1)), nodes = 1)$joint
    set.seed(1)
    joint <- globalSeq::cursus(Y = rbind(Y1, Y2), Yloc = c(1, 1),
                Ychr = c(1,2), X = list(rbind(X1A, X2A),rbind(X1B,X2B)),
                Xloc = list(rep(1, nrow(X1A) + nrow(X2A)),rep(1, nrow(X1B) + nrow(X2B))),
                Xchr = list(rep(c(1,2),times = c(nrow(X1A), nrow(X2A))),rep(c(1,2),times = c(nrow(X1A), nrow(X2A)))),
                window = list(0,0),
                offset = rep(1,ncol(Y1)), nodes = 1)
    # compare results
    testthat::expect_identical(chr.1, joint[1,1])
    testthat::expect_identical(chr.2, joint[2,1])
})

testthat::test_that("conversion to matrix works",{
    # simulate RNA-Seq data
    Y <- matrix(rnbinom(30,mu=10,size=1/0.2),nrow=10,ncol=3)
    rownames(Y) <- paste("gene",1:nrow(Y),sep="")
    colnames(Y) <- paste("cell",1:ncol(Y),sep="")
    # create data structure
    Z <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(counts=Y))
    # conversion
    testthat::expect_identical(Y,globalSeq::intern.matrix(Z))
})

testthat::test_that("calculating offset inside or outside cursus is equivalent",{
    # simulate high-dimensional data
    n <- 30; q <- 10; p <- 100
    Y <- matrix(rnbinom(q*n,mu=10,size=1/0.25),nrow=q,ncol=n)
    X <- matrix(rnorm(p*n),nrow=p,ncol=n)
    Yloc <- seq(0,1,length.out=q)
    Xloc <- seq(0,1,length.out=p)
    window <- 1
    # calculate inside function
    set.seed(1)
    inside <- globalSeq::cursus(Y,Yloc,X,Xloc,window,nodes=1)
    # calculate outside function
    set.seed(1)
    ls <- colSums(Y)
    gm <- exp(1/length(ls) * sum(log(ls)))
    offset <- ls/gm
    outside <- globalSeq::cursus(Y,Yloc,X,Xloc,window,offset=offset,nodes=1)
    testthat::expect_identical(inside,outside)
})
