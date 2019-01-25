#' @title Exact Monte Carlo Power Calculation by Inputting P (wmwpowp)
#' @name wmwpowp
#' @import smoothmest
#' @description \emph{wmwpowp} has two purposes:
#'
#' 1. Calculate the power for a
#' one-sided or two-sided Wilcoxon-Mann-Whitney test with an exact Monte Carlo p-value given
#' one user specified distribution and p (defined as P(X<Y)).
#'
#' 2. Calculate the parameters of the second
#' distribution. It is assumed that the second population is from the same type of
#' continuous probability distribution as the first population.
#'
#' Power is calculated empirically using simulated data and the parameters are calculated using derived
#' mathematical formulas for P(X<Y).
#'
#' @importFrom stats integrate nlm pnorm qnorm rnorm wilcox.test
#' 
#' @usage wmwpowp(n, m, distn, k = 1, p = NA, wmwodds = NA, sides, alpha = 0.05, nsims = 10000)
#' @param n Sample size for the first distribution (numeric)
#' @param m Sample size for the second distribution (numeric)
#' @param p The effect size, i.e., the probability that the first random variable is less than the 
#' second random variable (P(X<Y)) (numeric)
#' @param alpha Type I error rate or significance level (numeric)
#' @param distn Base R’s name for the first distribution (known as X in the above notation) and any
#' required parameters. Supported distributions are normal, exponential, and double exponential 
#' ("norm","exp", "doublex"). User may enter distribution without parameters, and default parameters will
#' be set (i.e., "norm" defaults to "norm(0,1)"), or user may specify both distribution and parameters 
#' (i.e., "norm(0,1)"). 
#' @param sides Options are “two.sided”, “less”, or “greater”. “less” means the alternative
#' hypothesis is that distn is less than distm (string)
#' @param k Standard deviation (SD) scalar for use with the normal or double 
#' exponential distribution options. The SD for distm is computed as k multiplied by
#'  the SD for distn. Equivalently, k is the ratio of the SDs of the second and first 
#'  distribution (k = SDm/SDn). Default is k=1 (equal SDs) (numeric)
#' @param wmwodds The effect size expressed as odds = p/(1-p). Either p or wmwodds must be
#' input (numeric)
#' @param nsims Number of simulated datasets for calculating power; 10,000 is the default.
#' For exact power to the hundredths place (e.g., 0.90 or 90\%) around 100,000 simulated
#' datasets is recommended (numeric)
#'
#' @references 
#' Mollan K.R., Trumble I.M., Reifeis S.A., Ferrer O., Bay C.P., Baldoni P.L.,
#' Hudgens M.G. Exact Power of the Rank-Sum Test for a Continuous Variable, 
#' arXiv:1901.04597 [stat.ME], Jan. 2019.
#' 
#' @examples
#' # We want to calculate the statistical power to compare the distance between mutations on a DNA 
#' # strand in two groups of people. Each group (X and Y) has 10 individuals. We assume that the 
#' # distance between mutations in the first group is exponentially distributed with rate 3. We assume
#' # that the probability that the distance in the first group is less than the distance in the second 
#' # group (i.e., P(X<Y)) is 0.8. The desired type I error is 0.05.
#'
#' wmwpowp(n = 10, m = 10, distn = "exp(3)", p = 0.8, sides = "two.sided", alpha = 0.05)
#'
#' @export

###########################################################################################################################
#Name: wmwpowp.R
#Programmer: Ilana Trumble
#Purpose: Write a flexible function to perform a power analysis for an exact Monte Carlo WMW test,
# given p'' and one distribution, through simulation. Also return shifted distribution
#Notes: p''=P(X<Y) given by the user; Works with continuous pdfs: norm, exp, double exponential
#Date Completed: 22NOV2017
###########################################################################################################################
library(smoothmest)

wmwpowp<- function(n,m,distn,k=1,p=NA,wmwodds=NA,sides="two.sided",alpha=0.05,nsims=10000)
{
  dist1<-distn
  if (dist1=="norm"){
    dist1<-"norm(0,1)"
  }
  if (dist1=="exp"){
    dist1<-"exp(1)"
  }
  if (dist1=="doublex"){
    dist1<-"doublex(0,1)"
  }
  n1=n
  n2=m
  #warnings
  if(is.numeric(n1) == F | is.numeric(n2) == F)
  {
    stop("n1 and n2 must be numeric")
  }

  if(is.numeric(k) == F)
  {
    stop("k must be numeric")
  }

  if(is.character(dist1) == F |
     !(sub("[^A-z]+", "", dist1) %in% c("norm","exp", "doublex")))
  {
    stop("distn must be characters in the form of distribution(parmater1) or distribution(paramter1, parameter2).
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

  if(!(sides %in% c("less","greater","two.sided")))
  {
    stop("sides must be character value of less, greater, or two.sided")
  }

  if(sub("[^A-z]+", "", dist1) == "exp" & k != 1)
  {
    stop("the parameter k is not applicable to the exponential distribution")
  }

  #wmwodds and p warnings
  if(is.na(p) == F & is.na(wmwodds) == F)
  {
    stop("please enter either p or wmwodds, not both")
  }

  if(is.na(p) == F)
  {
    if (is.numeric(p) == F | p<0 | p>1)
    {
      stop("p must be a probability between 0 and 1")
    }
    wmwodds <- round(p/(1-p),3)
  }

  else if(is.na(wmwodds) == F)
  {
    if (is.numeric(wmwodds) == F | wmwodds < 0)
    {
      stop("wmwodds must be a positive number")
    }
    p <- wmwodds/(1+wmwodds)
  }

  # calculate dist2
  create_dist2 <- function(input_dist_char,p)
  {
    dist_name <- sub("[^A-z]+", "", input_dist_char)


    if(dist_name =="norm")
    {
      # calculate mu2
      mu1<- as.numeric(gsub(".*\\(|,.*","",input_dist_char))
      sd1 <- as.numeric(gsub(".*,|\\)","",input_dist_char))
      sd2<-k*sd1

      mu2<- (sqrt(sd1^2+sd2^2)*qnorm(p))+mu1

      return(paste(dist_name,"(",mu2,",",sd2,")",sep=""))
    }

    else if(dist_name == "exp")
    {
      mu <- as.numeric(gsub(".*\\(|\\).*","",input_dist_char))
      lambda<- (mu/p)*(1-p)
      return(paste(dist_name,"(",lambda,")",sep=""))
    }

    else if(dist_name=="doublex")
    {
      mu1<- as.numeric(gsub(".*\\(|,.*","",input_dist_char))
      sigma1 <- as.numeric(gsub(".*,|\\)","",input_dist_char))
      sigma2<- k*sigma1

      integral<-function(mu1,mu2,sigma1,sigma2){
        a <- integrate(function(y)
          (1/(2*sigma2))*exp(-abs(y-mu2)/sigma2) * (1/2)*exp((y-mu1)/sigma1),
          -Inf,mu1
        )$value

        b <- integrate(function(y)
          (1/(2*sigma2))*exp(-abs(y-mu2)/sigma2) * (1-(1/2)*exp(-(y-mu1)/sigma1)),
          mu1, Inf
        )$value
        return(a+b)
      }

      rootfunc<-function(mu1,mu2,sigma1,sigma2,pr){
        s<-integral(mu1,mu2,sigma1,sigma2)-pr
        m<-s^2 #square so can use minimization
        return(m)
      }

      muy<-function(mu1,sigma1,sigma2,pr){
        start=mu1
        # non linear minimization
        mu2<-nlm(rootfunc,p=start,
                 mu1=mu1,sigma1=sigma1,sigma2=sigma2,pr=pr)$estimate
        return(mu2)
      }

      mu2<-muy(mu1=mu1,sigma1=sigma1,sigma2=sigma2,pr=p)

      return(paste(dist_name,"(",mu2,",",sigma2,")",sep=""))
    }
  }

  dist2<-create_dist2(dist1,p)

  #Power Simulation
  power_sim_func<-function(){
    dist1_func_char <- paste("r",sub("\\(.*", "", dist1),"(",n1,",",sub(".*\\(", "", dist1),sep="")
    dist2_func_char <- paste("r",sub("\\(.*", "", dist2),"(",n2,",",sub(".*\\(", "", dist2),sep="")
    return(wilcox.test(eval(parse(text = dist1_func_char)),eval(parse(text = dist2_func_char)),
                       paired=F,correct=F,exact=T,alternative=sides)$p.value)
  }
  #vectorize
  pval_vect<-replicate(nsims,power_sim_func())
  empirical_power <- round(sum(pval_vect<alpha)/length(pval_vect),3)

  # shorten decimals in dist2
  dist2print <- function(input_dist_char)
  {
    dist_name <- sub("[^A-z]+", "", input_dist_char)

    if(dist_name =="norm")
    {
      muorig<- as.numeric(gsub(".*\\(|,.*","",input_dist_char))
      munew<-round(muorig,3)
      var<- as.numeric(gsub(".*,|\\)","",input_dist_char))
      return(paste(dist_name,"(",munew,",",var,")",sep=""))
    }

    if(dist_name =="doublex")
    {
      muorig<- as.numeric(gsub(".*\\(|,.*","",input_dist_char))
      munew<-round(muorig,3)
      sigma<- as.numeric(gsub(".*,|\\)","",input_dist_char))
      return(paste(dist_name,"(",munew,",",sigma,")",sep=""))
    }

    if(dist_name == "exp")
    {
      muorig <- as.numeric(gsub(".*\\(|\\).*","",input_dist_char))
      munew<-round(muorig,3)
      return(paste(dist_name,"(",munew,")",sep=""))
    }
  }

  #Output
  if (sides=="two.sided"){test_sides<-"Two-sided"}
  if (sides %in% c("less","greater")){test_sides<-"One-sided"}

  # Round p and wmwodds as last step for printing output
  
  cat("Supplied distribution: ", dist1, "; n = ", n1, "\n",
      "Shifted distribution: ", dist2print(dist2), "; m = ", n2, "\n\n",
      "p: ", round(p,3), "\n",
      "WMW odds: ", round(wmwodds,3), "\n",
      "Number of simulated datasets: ", nsims, "\n",
      test_sides, " exact WMW test (alpha = ", alpha, ")\n\n",
      "Empirical power: ", empirical_power,sep = "")

  output_list <- list(empirical_power = empirical_power, alpha = alpha, test_sides = test_sides, p = p,
                      wmw_odds = wmwodds, distn = dist1, distm = dist2print(dist2), n = n1, m = n2)
  }


#wmwpowp(n1=10,n2=9,dist1="exp(2)",k=1,p=.67,wmwodds=NA,sides="two.sided",alpha=0.05,nsims=10000)
