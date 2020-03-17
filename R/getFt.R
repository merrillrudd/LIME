#' use newtons method to find instantaneous f
#' \code{getft} approximates f using newtons method based on catch, life history, and selectivity
#'
#' @author s.j.d. martell
#' @param ct catch
#' @param m natural mortality
#' @param sa selectivity at age
#' @param wa weight at age
#' @param na abundance at age
#' 
#' @return instantaneous fishing mortality rate f
#' @export
getFt <- function(ct,m,sa,wa,na)
{	#use newtons root finding method to get ft
	ft 	<- ct/(sum(na*exp(-m/2)*wa*sa))	
	for(iter in 1:10)
	{
		t1	<- wa*na
		t2	<- exp(-m-ft*sa)
		t3	<- (1-t2)
		t4	<- m+ft*sa
		c1	<- sum(ft*sa*t1*t3/t4)
		c2	<- sum(sa*t1*t3/t4 - ft*sa^2*t1*t3/t4^2 + ft*sa^2*t1*t2/t4)
		ft	<- ft - (c1-ct)/c2	#newton step.
	}
	
	return (ft)
}

