calc_msy <- function(ages, S_a, M, R0, W_a, Mat_a, F_vec){

	Nage <- calc_equil_abund(ages=ages, S_a=S_a, M=M, F=F, R0=R0)
	YPR <- sum(Nage * (1 - exp(-M - F*S_a)) * (F * S_a) / (M + F*S_a))
	
	return(YPR)
}