aldex_bis <- function (clr, conditions, mc.samples = 128, test = "t", effect = TRUE, 
          include.sample.summary = FALSE, verbose = FALSE, denom = "all") 
{
  print("aldex.clr: generating Monte-Carlo instances and clr values")
  # x <- aldex.clr(reads = reads, conds = conditions, mc.samples = mc.samples, 
  #                denom = denom, verbose = verbose, useMC = FALSE)
  x <- clr
  if (test == "iterative") {
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test = FALSE)
    print("aldex.ttest: seeding a second t-test")
    nonDE.i <- which(rownames(reads) %in% rownames(x.tt[x.tt$wi.eBH > 
                                                          0.05 | x.tt$we.eBH > 0.05, ]))
    if (length(nonDE.i) == 0) 
      stop("no non-DE references found")
    x.tt <- aldex(reads, conditions, mc.samples = mc.samples, 
                  test = "t", effect = effect, include.sample.summary = include.sample.summary, 
                  verbose = verbose, denom = nonDE.i)
  }
  else if (test == "t") {
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test = FALSE)
  }
  else if (test == "glm") {
    print("aldex.glm: doing Kruskal Wallace and glm test")
    x.tt <- aldex.glm(x, conditions)
  }
  else {
    stop("argument 'test' not recognized")
  }
  if (effect == TRUE && test == "t") {
    print("aldex.effect: calculating effect sizes")
    x.effect <- aldex.effect(x, conditions, include.sample.summary = include.sample.summary, 
                             verbose = verbose)
    z <- data.frame(x.effect, x.tt)
  }
  else {
    z <- data.frame(x.tt)
  }
  return(z)
}
