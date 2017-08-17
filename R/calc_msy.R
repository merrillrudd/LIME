calc_msy <- function(F, ages, S_a, M, R0, W_a, Mat_a){

	Nage <- calc_equil_abund(ages=ages, S_a=S_a, M=M, F=F, R0=R0)
	YPR <- sum(Nage * W_a * (1 - exp(-M - F*S_a)) * (F * S_a) / (M + F*S_a))

	return(YPR)
}