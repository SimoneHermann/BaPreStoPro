

#' S4 class for the jump diffusion process
#' @slot theta
#' @slot phi ...

setClass(Class = "jumpDiffusion", representation = representation(theta = "numeric", phi = "numeric", gamma2 = "numeric", xi = "numeric",
                                                                  b.fun = "function", s.fun = "function", h.fun = "function",
                                                                  Lambda = "function", priorRatio = "list", start = "list"))

setClass(Class = "Merton", representation = representation(thetaT = "numeric", phi = "numeric", gamma2 = "numeric", xi = "numeric",
                                                           Lambda = "function", prior = "list", start = "list", priorRatio = "function"))

setClass(Class = "Diffusion", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                              b.fun = "function", sT.fun = "function", prior = "list", start = "list"))

setClass(Class = "mixedDiffusion", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                   y0.fun = "function", b.fun = "function", sT.fun = "function", prior = "list", start = "list"))

setClass(Class = "hiddenDiffusion", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                                    sigma2 = "numeric",
                                                                    b.fun = "function", sT.fun = "function", y0.fun = "function",
                                                                    prior = "list", start = "list"))

setClass(Class = "hiddenmixedDiffusion", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                    sigma2 = "numeric",
                                                                    b.fun = "function", sT.fun = "function", y0.fun = "function",
                                                                    prior = "list", start = "list"))

setClass(Class = "jumpRegression", representation = representation(theta = "numeric", gamma2 = "numeric", xi = "numeric", fun = "function",
                                                                   Lambda = "function", sT.fun = "function", prior = "list", start = "list"))

setClass(Class = "NHPP", representation = representation(xi = "numeric", Lambda = "function", priorRatio = "function", start = "numeric"))

setClass(Class = "Regression", representation = representation(phi = "numeric", gamma2 = "numeric",
                                                              fun = "function", sT.fun = "function", prior = "list", start = "list"))

setClass(Class = "mixedRegression", representation = representation(phi = "matrix", mu = "numeric", Omega = "numeric", gamma2 = "numeric",
                                                                   fun = "function", sT.fun = "function", prior = "list", start = "list"))

####################
# estimation classes
####################

setClass(Class = "est.jumpDiffusion", representation = representation(theta = "numeric", phi = "numeric", gamma2 = "numeric", xi = "matrix",
                                                                      model = "list", N.est = "matrix", t = "numeric",
                                                                      Y = "numeric", N = "numeric",
                                                                      burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Merton", representation = representation(thetaT = "numeric", phi = "numeric", gamma2 = "numeric", xi = "matrix",
                                                               model = "list", N.est = "matrix", t = "numeric",
                                                               Y = "numeric", N = "numeric",
                                                               burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Diffusion", representation = representation(phi = "matrix", gamma2 = "numeric",
                                                                  model = "list", t = "numeric", Y = "numeric",
                                                                  burnIn = "numeric", thinning = "numeric"))
# Y, t as list ?!?
setClass(Class = "est.mixedDiffusion", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                       model = "list", t = "numeric", Y = "matrix", t.list = "list", Y.list = "list",
                                                                       burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.hiddenDiffusion", representation = representation(phi = "matrix", gamma2 = "numeric", sigma2 = "numeric", Y.est = "matrix",
                                                                        model = "list", t = "numeric", Z = "numeric",
                                                                        burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.hiddenmixedDiffusion", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                         sigma2 = "numeric", Y.est = "list",
                                                                         model = "list", t = "numeric", Z = "matrix", t.list = "list", Z.list = "list",
                                                                         burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.jumpRegression", representation = representation(theta = "matrix", gamma2 = "numeric", xi = "matrix", N.est = "matrix",
                                                                       model = "list", t = "numeric", Y = "numeric", N = "numeric",
                                                                       burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.NHPP", representation = representation(xi = "matrix", model = "list", N = "numeric", t = "numeric", jumpTimes = "numeric",
                                                             burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.Regression", representation = representation(phi = "matrix", gamma2 = "numeric",
                                                                   model = "list", t = "numeric", Y = "numeric",
                                                                   burnIn = "numeric", thinning = "numeric"))

setClass(Class = "est.mixedRegression", representation = representation(phi = "list", mu = "matrix", Omega = "matrix", gamma2 = "numeric",
                                                                        model = "list", t = "numeric", Y = "matrix", t.list = "list", Y.list = "list",
                                                                        burnIn = "numeric", thinning = "numeric"))






