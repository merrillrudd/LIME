#' Use Newtons method to find instantaneous F
#' \code{getFt} Approximates F using Newtons method based on catch, life history, and selectivity
#'
#' @author S.J.D. Martell
#' @param ct catch
#' @param m natural mortality
#' @param sa selectivity at age
#' @param wa weight at age
#' @param na abundance at age
#' 
#' @return instantaneous fishing mortality rate F
#' @export
getFt <- function(ct,m,sa,wa,na)
{	#use newtons root finding method to get ft
	ft 	<- ct/(sum(na*exp(-m/2)*wa*sa))	
	for(iter in 1:9)
	{
		T1	<- wa*na
		T2	<- exp(-m-ft*sa)
		T3	<- (1-T2)
		T4	<- m+ft*sa
		c1	<- sum(ft*sa*T1*T3/T4)
		c2	<- sum(sa*T1*T3/T4 - ft*sa^2*T1*T3/T4^2 + ft*sa^2*T1*T2/T4)
		ft	<- ft - (c1-ct)/c2	#newton step.
	}
	
	return (ft)
}

