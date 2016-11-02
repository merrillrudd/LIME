#' Example data from Costa Rican spotted rose snapper.
#'
#' A dataset containing the length composition, biological information, and abundance index for Costa Rican spotted rose snapper.
#'
#' @format A list with 39 items:
#' \describe{
#'   \item{vbk}{number; von Bertanffy k Brody growth coefficient}
#'   \item{linf}{number; asymptotic length in von Bertalanffy curve}
#'   \item{t0}{number, default=-0.01; theoretical age at length=0}
#'   \item{binwidth}{integer, default=1; width of length bin}
#'   \item{CVlen}{number; default=0.1; coefficient of variation for growth curve}
#'   \item{SigmaC}{number; default=0.2; standard deviation for observation error in catch data}
#'   \item{SigmaI}{number; default=0.2; standard deviation for observation error in abundance index}
#'   \item{SigmaR}{number; default=0.6; standard deviation for process error in recruitment}
#'   \item{SigmaF}{number; default=0.2; standard deviation for process error in fishing mortality}
#'   \item{R0}{number; default=1; equilibrium recruitment - default of 1 represents relative recruitment deviations without measure of scale}
#'   \item{lwa}{number; scaling parameter in length-weight curve}
#'   \item{lwb}{number; allometric parameter in length-weight curve}
#'   \item{S50}{number; age or length at 50% selectivity}
#'   \item{h}{steepness parameter; default=1}
#'   \item{qcoef}{catchability coefficient relevant when abundance index is included; default=1e-5}
#'   \item{M}{natural mortality rate; default = 1.5*vbk when not specified directly}
#'   \item{F1}{starting value for fishing mortality rate estimation, default=0.2}
#'   \item{AgeMax}{maximum age in population, used as plus group. default is 1% of remaining population based on natural mortality M}
#'   \item{mids}{midpoints of length bins.}
#'   \item{highs}{upper end of length bins. May be different than other argument 'lbins'. 'lbins' is based on the observed length frequency data, this argument "highs" is based on 1.5*linf. will default to smaller bins in create_input function.}
#'   \item{lows}{lower end of length bins.}
#'   \item{S_a}{vector of selectivity at age calculated from inputs}
#'   \item{L_a}{vector of length at age calculated from inputs}
#'   \item{W_a}{vector of weight at age calculated from inputs}
#'   \item{M50}{number; age or length at 50% maturity}
#'   \item{Mat_a}{vector of maturity at age calculated from inputs}
#'   \item{Fequil}{equilibrium fishing mortality (used in simulation; not relevant here)}
#'   \item{Frate}{rate of expected increase for fishing mortality (used in simulation; not relevant here)}
#'   \item{Fmax}{maximum fishing mortality (used in simulation; not relevant here)}
#'   \item{I_t}{vector of abundance index over time, with names indicating the year in the time series that they come from}
#'   \item{LF}{matrix with the proportion of the catch in each length bin (columsn) in each year (rows) with the upper limit of each bin naming the columns and the year of the length composition naming the rows}
#' 	 \item{LFprop}{the same as LF in this case, but indicating that it is possible to specify the length composition LF in counts and the LFprop in proportions}
#'   \item{years}{a vector of the total years to model; e.g. 1997-2015}
#'   \item{years_i}{a vector of the total years to model, as an index; e.g. 1-19}
#'   \item{lbins}{a vector with the upper limit of each length bin based on the observed length composition (may be different than other argument 'highs' based on 1.5*linf). The final length bins are output in the create_inputs function.}
#'   \item{ML_t}{a vector with the mean length in each year where length information is available, with the vector names the years as an index}
#'   \item{Nyears}{an integer of the total number of years to model}
#'   \item{Nyears_comp}{an integer of the total number of years where length data is available}
#'   \item{obs_per_year}{vector of the effective sample size of length measurements in each year}
#'   ...
#' }
#' @source {Conservation International}
"snapper"