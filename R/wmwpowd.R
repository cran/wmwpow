#' @title Power Calculation with wmwpowd
#' @name wmwpowd
#' @description \emph{wmwpowd} has two purposes:
#'
#' 1. Determine the power for a one-sided or two-sided Wilcoxon-Mann-Whitney test
#' with an exact p-value given two user specified distributions.
#'
#' 2. Determine p, the P(X<Y), where X represents random draws from one continuous
#' probability distribution and Y represents random draws from another distribution;
#' p is useful for quantifying the effect size that the Wilcoxon-Mann-Whitney test is
#' assessing.
#'
#' Both 1. and 2. are calculated empirically using simulated data and output automatically.
#'
#' @importFrom stats integrate nlm pnorm qnorm rnorm wilcox.test
#' 
#' @usage wmwpowd(n, m, distn, distm, sides, alpha = 0.05, nsims = 10000)
#' @param n Sample size for the first distribution (numeric)
#' @param m Sample size for the second distribution (numeric)
#' @param alpha Type I error rate or significance level (numeric)
#' @param distn Base R’s name for the first distribution and any required parameters 
#' ("norm", "beta", "cauchy", "f", "gamma", "lnorm", "unif", "weibull","exp", "chisq", "t", "doublex")
#' @param distm Base R’s name for the second distribution and any required parameters 
#' ("norm", "beta", "cauchy", "f", "gamma", "lnorm", "unif", "weibull","exp", "chisq", "t", "doublex")
#' @param sides Options are “two.sided”, “less”, or “greater”. “less” means the
#' alternative hypothesis is that distn is less than distm (string)
#' @param nsims Number of simulated datasets for determining power; 10,000 is the default (numeric)
#' @note Example of distn, distm: “norm(1,2)” or “exp(1)”
#'
#' In addition to all continuous distributions supported in Base R, \emph{wmwpowd} also supports the 
#' double exponential distribution from the smoothmest package
#'
#' The output WMWOdds is p expressed as odds (p/(1-p))
#'
#' Use $ notation to select specific output parameters
#'
#' The function has been optimized to run through simulations quickly; long wait times are unlikely 
#' for n and m of 50 or fewer

#' @examples
#' # 1. We want to determine the statistical power to compare body length measured on two groups of
#' # rabbits. Each group (X and Y) has 7 rabbits. We assume that body length will be normally 
#' # distributed and have a constant standard deviation of 2 cm among groups. We assume that Group X 
#' # will have a mean of 35 cm and Group Y will have a mean of 32 cm; the desired type I error is 0.05.
#'
#' \donttest{wmwpowd(n = 7, m = 7, distn = "norm(35,2)", distm = "norm(32,2)", sides = "two.sided", 
#'         alpha = 0.05, nsims=10000)}
#'         
#' \dontshow{wmwpowd(n = 5, m = 5, distn = "norm(35,2)", distm = "norm(32,2)", sides = "two.sided", 
#'         alpha = 0.05, nsims=100)}
#'
#' # 2. We are interested in determining the statistical power (with type I error = 0.05) for a
#' # comparison of the use of ornamentation among fiddle players living in two regions of the United 
#' # States: X county, Texas and Y county, North Carolina. A random sample of 18 fiddlers will be 
#' # collected within each state. The fiddlers will practice and perform a standardized version of the 
#' # Tennessee Waltz. The proportion of melody notes that are ornamented (including vibrato) will be 
#' # calculated. We assume that the proportion will follow a beta distribution with a mean of 0.40 and 
#' # a shape well described by alpha = 8 and beta = 12 among the Texas fiddlers. We assume that 
#' # the distribution will be shifted to a lower mean of 0.25 and have the shape alpha = 2, 
#' # beta = 6 for the North Carolina fiddlers.
#'
#' \donttest{wmwpowd(n=18, m=18, distn = "beta(8,12)", distm = "beta(2,6)", sides = "two.sided", 
#'         alpha = 0.05, nsims=10000)}
#'         
#' \dontshow{wmwpowd(n=5, m=5, distn = "beta(8,12)", distm = "beta(2,6)", sides = "two.sided", 
#'         alpha = 0.05, nsims=100)}
#' @export
#'
#'
#'

###########################################################################################################################
#Name: wmwpowd.R
#Programmer: Camden Bay
#Purpose: A flexible function to perform a power analysis for an exact WMW test through simulation AND output p''
#Notes: p'' is determined empirically. smoothmest must be installed.
#Date Completed: 11/22/2017
###########################################################################################################################

library(smoothmest)

wmwpowd <- function(n,m,distn,distm,sides="two.sided",alpha = 0.05,nsims=10000)
{
  dist1<-distn
  dist2<-distm
  n1=n
  n2=m
  if(is.numeric(n1) == F | is.numeric(n2) == F)
  {
    stop("n1 and n2 must be numeric")
  }

  if(is.character(dist1) == F | is.character(dist2) == F |
     !(sub("[^A-z]+", "", dist1) %in% c("norm", "beta", "cauchy", "f", "gamma", "lnorm", "unif", "weibull","exp", "chisq", "t", "doublex")) |
     !(sub("[^A-z]+", "", dist2) %in% c("norm", "beta", "cauchy", "f", "gamma", "lnorm", "unif", "weibull","exp", "chisq", "t", "doublex")))
  {
    stop("distn and distm must be characters in the form of distribution(parmater1) or distribution(paramter1, parameter2).
         See documentation for details.")
  }

  if(is.numeric(alpha) == F)
  {
    stop("alpha must be numeric")
  }

  if(is.numeric(nsims) == F)
  {
    stop("nsims must be numeric")
  }

  if (!(sides %in% c("less","greater","two.sided"))){
    stop("sides must be character value of less, greater, or two.sided")
  }

  if(sides == "two.sided"){test_sides <- "Two-sided"}
  if(sides %in% c("less", "greater")){test_sides <- "One-sided"}

  #Power Simulation
  dist1_func_char <- paste("r",sub("\\(.*", "", dist1),"(",n1,",",sub(".*\\(", "", dist1),sep="")
  dist2_func_char <- paste("r",sub("\\(.*", "", dist2),"(",n2,",",sub(".*\\(", "", dist2),sep="")

  power_sim_func <- function()
  {
    wilcox.test(eval(parse(text = dist1_func_char)),eval(parse(text = dist2_func_char)),paired=F,correct=F,alternative=sides,
                exact=T)$p.value
  }
  pval_vect <- replicate(nsims,power_sim_func())
  empirical_power <- round(sum(pval_vect<alpha)/length(pval_vect),3)

  #p
  dist1_func_ppp <- paste("r",sub("\\(.*", "", dist1),"(",10000000,",",sub(".*\\(", "", dist1),sep="")
  dist2_func_ppp <- paste("r",sub("\\(.*", "", dist2),"(",10000000,",",sub(".*\\(", "", dist2),sep="")

  p <- round(sum(eval(parse(text = dist1_func_ppp)) < eval(parse(text = dist2_func_ppp)))/10000000,3)
  wmw_odds <- round(p/(1-p),3)

  #Output
  cat("Supplied distribution 1: ", dist1, "; n = ", n1, "\n",
      "Supplied distribution 2: ", dist2, "; m = ", n2, "\n\n",
      "p: ", p, "\n",
      "WMW odds: ", wmw_odds, "\n",
      test_sides, " exact WMW test (alpha = ", alpha, ")\n\n",
      "Empirical power: ", empirical_power,sep = "")

  output_list <- list(empirical_power = empirical_power, alpha = alpha, test_sides = test_sides, p = p,
                      wmw_odds = wmw_odds, distn = dist1, distm = dist2, n = n1, m= n2)
  }



