#' @title Power Calculation Using the Shieh et. al. Approach
#' @name shiehpow
#' @import lamW
#' @description The purpose of \emph{shiehpow} is to perform a power analysis for a one
#' or two-sided Wilcoxon-Mann-Whitney test using the method developed by Shieh and
#' colleagues.
#'
#' @importFrom stats integrate nlm pnorm qnorm rnorm wilcox.test
#' 
#' @param n Sample size of first sample (numeric)
#' @param m Sample size of second sample (numeric)
#' @param p Effect size, P(X<Y) (numeric)
#' @param alpha Type I error rate (numeric)
#' @param dist The distribution type for the two groups (“exp”, “dexp”, or “norm”) (string)
#' @param sides Options are “two.sided” and “one.sided” (string)
#' @note When calculating power for dist=”norm”, \emph{shiehpow} uses 100,000 draws from a Z ~ N(0,1)
#' distribution for the internal calculation of p2 and p3 from Shieh et al. (2006); thus 
#' \emph{shiehpow} normal distribution power results may vary in the thousandths place from one run 
#' to the next.
#'
#' @references 
#' Shieh, G., Jan, S. L., Randles, R. H. (2006). On power and sample size
#' determinations for the Wilcoxon–Mann–Whitney test. Journal of Nonparametric
#' Statistics, 18(1), 33-43.
#' 
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
#' shiehpow(n = 10, m = 10, p = 0.80, alpha = 0.05, dist = "exp", sides = "two.sided")
#'
#' @export



#################################################################################
# Author: Orlando Ferrer
# Edited by: Camden Bay, Katie Mollan
# Date: 29JUNE2017, last edited: 20JAN2018
# Notes:
# Ha: theta is not equal to 0, theta is abbreviated as 'tht'
# N is the total sample size, n and m are the group sizes
# calculated parameters: tht, muo, sigo, mu, sig
# p = P(X_1 < Y_1), alpha is set at a default of 0.05,
# "dist" will take inputs: "exp","norm", or "dexp". And the "sides" input will
# take either "one.sided" or "two.sided"
#################################################################################

library(lamW)
shiehpow <- function(n, m, p, alpha=.05, dist, sides="two.sided")
{
  if(!(dist %in% c("exp", "norm", "dexp")))
  {
    stop('dist must be equal to "exp", "norm", or "dexp."')
  }

  if(sides == "two.sided"){test_sides <- "Two-sided"}
  if(sides == "one.sided"){test_sides <- "One-sided"}


  # Below we define internal parameters needed for calculating power
  N=m+n
  np<-length(p)
  nn<-length(n)
  tht = vector(length=length(p))
  muo=vector(length=length(n))
  mu=matrix(0,nrow=length(p),ncol=length(n))
  sigo=vector(length=length(n))
  p2=vector(length=length(p))
  p3=vector(length=length(p))
  sig=matrix(0,nrow=length(p),ncol=length(n))
  p1=vector(length=length(p))
  power <- matrix(0,nrow=np,ncol=nn)

  # Within this double for-loop we calculate power for 3 different cases: Exponential
  # distributions, Standard Normal distributrions, and Double Exponentials distributions
  for(i in 1:np) {
    for(j in 1:nn) {

      # Exponential Case
      if (dist == "exp"){

        # Here we fix for when theta is negative, which happens when p < 0.5
        if(p[i] >= .5) {
          tht[i]= -log(2*(1-p[i]))
          p1[i] <- p[i]
        }
        else if (p[i] < .5) {
          p1[i] <- 1 - p[i] #When p < 0.5, we evaluate theta[1-p] instead of theta[p]
          tht[i]= -log(2*(1-p1[i]))

        }
        p2[i]= 1-(2/3)*exp(-tht[i]) #P2 value for exp dist
        p3[i]= 1-exp(-tht[i]) + (1/3)*exp(-2*tht[i]) #P3 value for exp dist
      }

      # Standard normal dist case
      if(dist=="norm"){
        if(p[i] >= .5) {
          tht[i]= sqrt(2)*qnorm(p[i])
          p1[i] <- p[i]
        }
        else if (p[i] < .5) {
          p1[i] <- 1 - p[i]
          tht[i]= sqrt(2)*qnorm(p1[i])

        }
        p2[i]= mean((pnorm(rnorm(100000)+tht[i]))^2) #P2 value for standard norm
        p3[i]= p2[i] #P3 value for standard norm
      }

      # Double Exp dist case
      if(dist=="dexp"){
        if(p[i] >= .5) {

          # Here we use a "Lambert W" function to be able to write theta in terms of p
          # Need to install "lamW" package
          tht[i]= (-lamW::lambertWm1(4*(p[i]-1)/(exp(1)^2))) -2
          p1[i] <- p[i]
        }
        else if (p[i] < .5) {
          p1[i] <- 1 - p[i]
        }
        tht[i]= (-lamW::lambertWm1(4*(p1[i]-1)/(exp(1)^2)))-2


        p2[i]= 1-(7/12 + tht[i]/2)*exp(-tht[i]) -(1/12)*exp(-2*tht[i])
        p3[i]= p2[i]
      }


      N[j]= m[j]+n[j]
      muo[j]= m[j]*n[j]/2
      mu[i,j]= (m[j])*(n[j])*p1[i]

      sigo[j]= sqrt(m[j]*n[j]*(m[j]+n[j]+1)/12)


      sig[i,j]=sqrt(m[j]*n[j]*p1[i]*(1-p1[i]) + m[j]*n[j]*(n[j]-1)*(p2[i]-(p1[i])^2) + m[j]*n[j]*(m[j]-1)*(p3[i]-(p1[i])^2))

      za1<-qnorm(alpha,lower.tail=F) #one-sided
      za2<-qnorm(alpha/2,lower.tail=F) #two-sided

      power[i,j]<-pnorm((mu[i,j]-muo[j]-za2*sigo[j])/sig[i,j]) + pnorm((-mu[i,j]+muo[j]-za2*sigo[j])/sig[i,j])
      # if one-sided, we use: power[i,j]<-pnorm((mu[i,j]-muo[j]-za2*sigo[j])/sig[i,j])
      # and we would use za1 instead of za2


      # one-sided case
      if (sides=="one.sided") {
        power[i,j]<-pnorm((mu[i,j]-muo[j]-za1*sigo[j])/sig[i,j])
      }
    }

  }
  wmw_odds <- round(p/(1-p),3)

  cat("Distribution: ", dist, "\n",
      "Sample sizes: ", n, " and ", m, "\n",
      "p: ", p, "\n",
      "WMW odds: ", wmw_odds, "\n",
      "sides: ", test_sides, "\n",
      "alpha: ", alpha, "\n\n",
      "Shieh Power: ", round(power,3),sep = "")

  output_list <- list(distribution = dist,n = n, m = m, p = p, wmw_odds = wmw_odds,
                      alpha = alpha, test_sides = test_sides, power = round(power,3))
}

#shiehpow(n = 9, m = 9, p = 0.76, alpha=.05, dist = "exp", sides="two.sided")



