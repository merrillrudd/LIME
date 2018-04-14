# 2017 Blue/Deacon Rockfish Assessment for California
# Point Conception to the CA-OR border
# E.J. Dick, NMFS/SWFSC, edward.dick@noaa.gov
# Stock Synthesis V3.30.03.07

0 # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
#
2 # recr_dist_method for parameters:  1=like 3.24; 2=main effects for GP, Settle timing, Area; 3=each Settle entity; 4=none when N_GP*Nsettle*pop==1
1 # Recruitment: 1=global; 2=by area (future option)
1 #  number of recruitment settlement assignments
0 # year_x_area_x_settlement_event interaction requested (only for recr_dist_method=1)
#GP month area age
  1     7    1   0 # settlement month set to July (close approximation, 7cm lowest pop'n size bin matches this)
#
2 #_Nblock_Patterns
2 1 #_blocks_per_pattern
# begin and end years of blocks
# BLOCK 1: 1971 and 2000 bag limit changes
1899 1971
1972 1999
# BLOCK 2: 2000 bag limit change only
1899 1999
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
0 0 0 0 0 # autogen
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
0    #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
1    # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K; 4=not implemented
2    #_Growth_Age_for_L1
30   #_Growth_Age_for_L2 (999 to use as Linf)
0.13 #_exponential decay for growth above maxage (fixed at 0.2 in 3.24; value should approx initial Z; -999 replicates 3.24)
0    #_placeholder for future growth feature
0    #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
1    #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1    #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2    #_First_Mature_Age
2    #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0    #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
2    #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
###########
# FEMALES #
###########
#     LO        HI      INIT     PRIOR  PR_SD  PR_type  PHASE env_var dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# natural mortality
   0.001       0.4     0.119   -2.0272  0.438        3      2       0        0         0         0      0     0         0 # NatM_p_1_Fem_GP_1
# length-at-age
      10        30      17.4        21   1000        0      2       0        0         0         0      0     0         0 # L_at_Amin_Fem_GP_1
      25        45      37.4        38   1000        0      2       0        0         0         0      0     0         0 # L_at_Amax_Fem_GP_1
    0.01       0.3     0.117      0.12   1000        0      2       0        0         0         0      0     0         0 # VonBert_K_Fem_GP_1
    0.01       0.5       0.1       0.1   1000        0      2       0        0         0         0      0     0         0 # CV_young_Fem_GP_1
    0.01       0.5       0.1       0.1   1000        0      2       0        0         0         0      0     0         0 # CV_old_Fem_GP_1
# weight-length
3.4e-005  3.4e-005  3.4e-005  3.4e-005   1000        0     -2       0        0         0         0      0     0         0 # Wtlen_1_Fem
       1         3      2.87      2.87   1000        0     -2       0        0         0         0      0     0         0 # Wtlen_2_Fem
# maturity ogive
      22        32        26        26   1000        0     -2       0        0         0         0      0     0         0 # Mat50%_Fem
    -0.7      -0.5      -0.6      -0.6   1000        0     -2       0        0         0         0      0     0         0 # Mat_slope_Fem
# fecundity F=aL^b; 'a' has been scaled so spawning output has units of millions of eggs
       0         1 1.143e-08 1.143e-08   1000        0     -2       0        0         0         0      0     0         0 # 'a' in F=a*L^b
       2         5     4.816     4.816   1000        0     -2       0        0         0         0      0     0         0 # 'b' in F=a*L^b
#########
# MALES # -- M and vB parameters are exponential offsets; offset = ln(M/F); M = F*exp(offset)
#########
#     LO        HI      INIT     PRIOR  PR_SD  PR_type  PHASE env_var dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# natural mortality
      -3         3       0.3         0   1000        0      2       0        0         0         0      0     0         0 # NatM_p_1_Mal_GP_1
# length-at-age
      -3         3         0         0   1000        0      2       0        0         0         0      0     0         0 # L_at_Amin_Mal_GP_1
      -3         3      -0.2         0   1000        0      2       0        0         0         0      0     0         0 # L_at_Amax_Mal_GP_1
      -3         3      0.22         0   1000        0      2       0        0         0         0      0     0         0 # VonBert_K_Mal_GP_1
      -3         3         0         0   1000        0      2       0        0         0         0      0     0         0 # CV_young_Mal_GP_1
      -3         3       0.5         0   1000        0      2       0        0         0         0      0     0         0 # CV_old_Mal_GP_1
# weight-length
       0         1  2.9e-005  2.9e-005   1000        0     -2       0        0         0         0      0     0         0 # Wtlen_1_Mal
       1         3      2.89      2.89   1000        0     -2       0        0         0         0      0     0         0 # Wtlen_2_Mal
#
# recruitment apportionment
       0         0         0         0      0        0     -1       0        0         0         0      0     0         0 # RecrDist_GP_1
       0         0         0         0      0        0     -1       0        0         0         0      0     0         0 # RecrDist_Area_1
       0         0         0         0      0        0     -1       0        0         0         0      0     0         0 # RecrDist_Bseas_1
# cohort growth deviation
       1         1         1         1   1000        0     -1       0        0         0         0      0     0         0 # CohortGrowDev
# fraction female 
   0.001     0.999       0.5       0.5   1000        0     -1       0        0         0         0      0     0         0 # FracFemale_GP_1
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#
#
#_Spawner-Recruitment
3  #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepard_3Parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#   LO    HI     INIT  PRIOR  PR_SD  PR_type  PHASE  env-var  use_dev  dev_mnyr  dev_mxyr  dev_PH  Block  Blk_Fxn #  parm_name
     5    12        9     10     10        0      1        0        0         0         0       0      0        0 # SR_LN(R0)
 0.201 0.999     0.65  0.718  0.158        2      2        0        0         0         0       0      0        0 # SR_std_B-H
   0.1     1      0.5    0.5   1000        0     -1        0        0         0         0       0      0        0 # SR_sigmaR
    -5     5        0      0   1000        0     -1        0        0         0         0       0      0        0 # SR_regime
     0   0.5        0      0   1000        0     -1        0        0         0         0       0      0        0 # SR_autocorr
#
1    #do_recdev:  0=none; 1=devvector; 2=simple deviations
1960 # first year of main recr_devs; early devs can preceed this era
2015 # last year of main recr_devs; forecast devs start in following year
3    #_recdev phase
1    # (0/1) to read 13 advanced options
 1950      #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 4         #_recdev_early_phase
 0         #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1         #_lambda for Fcast_recr_like occurring before endyr+1
 1966.0   #_last_early_yr_nobias_adj_in_MPD 
 1969.8   #_first_yr_fullbias_adj_in_MPD 
 2013.7   #_last_yr_fullbias_adj_in_MPD 
 2017.8   #_first_recent_yr_nobias_adj_in_MPD 
 0.462   #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models)   
 0         #_period of cycles in recruitment (N parms read below)
-5         #_min rec_dev
 5         #_max rec_dev
 0         #_read_recdevs
#_end of advanced SR options
#
#Fishing Mortality info 
0.5 # F ballpark
-2002 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
4 # Number of tuning iterations in hybrid F method
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
#_Q_setup
# 1. Fleet
# 2. Link type: 1 = simple Q (proportional); 2 = mirror simple Q; 3 = Q with power, 2 parameters
# 3. Extra input for link (i.e. mirror fleet); >0 = mirror the Q from another (lower numbered survey designated by abs(value))
# 4. Do extra SD; 0 = skip (typical); 1 = estimate constant to be added to the input stddev of survey
# 5. Bias adjustment; 0 = no bias adjustment applied, 1 = apply bias adjustment
# 6. Q float; 0 = no float (parameter is estimated); 1 = float (analytical solution, parameter line still required)
#_   fleet  link link_info  extra_se  biasadj  float  #  fleetname
         1     1         0         1        0      1  #  DOCKSIDE MRFSS CPFV INDEX
         3     1         0         1        0      1  #  DOCKSIDE CRFS PR INDEX
         9     1         0         1        0      1  #  ONBOARD_CPFV_88_98_N
        10     1         0         1        0      1  #  JUV_SWFSC_N
        11     1         0         1        0      1  #  ONBOARD_CPFV_01_16_N
        14     1         0         0        0      0  #  SSB_Survey_2017
-9999 0 0 0 0 0
#
#_Q_parms(if_any); Qunits_are_ln(q)
# Q PARAMETER LINES ARE REQUIRED, EVEN WHEN USING ANALYTICAL SOLUTION
# extra_se parameter lines only needed when activated (extra_se = 1) above
# LO   HI      INIT  PRIOR  PR_SD  PR_type  PHASE  env-var  use_dev  dev_mnyr  dev_mxyr  dev_PH  Block  Blk_Fxn  #  parm_name
 -15   15        -9      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_DOCKSIDE_MRFSS_CPFV_N
  0  0.75      0.22   0.05      1        0      2        0        0         0         0       0      0        0  #  extra_se_DOCKSIDE_MRFSS_CPFV_N
 -15   15        -9      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_DOCKSIDE_CRFS_PR_N
  0  0.75      0.55   0.05      1        0      2        0        0         0         0       0      0        0  #  extra_se_DOCKSIDE_CRFS_PR_N
 -15   15        -9      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_ONBOARD_CPFV_88_98_N
  0  0.75      0.17   0.05      1        0      2        0        0         0         0       0      0        0  #  extra_se_ONBOARD_CPFV_88_98_N
 -15   15        -9      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_JUV_SWFSC_N
  0  0.75      0.19   0.05      1        0      2        0        0         0         0       0      0        0  #  extra_se_JUV_SWFSC_N
 -15   15        -9      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_ONBOARD_CPFV_01_16_N
  0  0.75      0.68   0.05      1        0      2        0        0         0         0       0      0        0  #  extra_se_ONBOARD_CPFV_01_16_N
# fix q=1 for decision table analysis (survey lambda=0 in base)
 -15   15         0      0      1        0     -1        0        0         0         0       0      0        0  #  LnQ_base_SSB_Survey_2017
#
#_size_selex_types
#discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead
#_Pattern Discard Male Special
       24       0    0       0 # 1 REC_CPFV_N
       24       0    0       0 # 2 REC_PRIV_N_MRFSS
        5       0    0       2 # 3 REC_PRIV_N_CRFS -- mirror MRFSS
       24       0    0       0 # 4 REC_DISC_N
       24       0    0       0 # 5 COM_HKL_N
       24       0    0       0 # 6 COM_NET_N
        5       0    0       5 # 7 COM_OTH_N -- mirror HKL
       24       0    0       0 # 8 COM_DISC_N
       24       0    0       0 # 9 ONBOARD_CPFV_88_98_N
        0       0    0       0 #10 JUV_SWFSC_N
        5       0    0       1 #11 ONBOARD_CPFV_01_16_N -- mirror CPFV
       24       0    0       0 #12 SCHMIDT_N
       24       0    0       0 #13 ABRAMS_N
        0       0    0       0 #14 SSB_Survey_2017
#
#_age_selex_types
#_Pattern Discard Male Special
       10       0    0       0 # 1 REC_CPFV_N
       10       0    0       0 # 2 REC_PRIV_N_MRFSS
       10       0    0       0 # 3 REC_PRIV_N_CRFS
       10       0    0       0 # 4 REC_DISC_N
       10       0    0       0 # 5 COM_HKL_N
       10       0    0       0 # 6 COM_NET_N
       10       0    0       0 # 7 COM_OTH_N
       10       0    0       0 # 8 COM_DISC_N
       10       0    0       0 # 9 ONBOARD_CPFV_CenCA_N
       11       0    0       0 #10 JUV_SWFSC_N -- selects only age-0 fish (YOY index)
       10       0    0       0 #11 ONBOARD_CPFV_01_16_N
       10       0    0       0 #12 SCHMIDT_N
       10       0    0       0 #13 ABRAMS_N
       10       0    0       0 #14 SSB_Survey_2017
#
# LO  HI      INIT  PRIOR  PR_SD  PR_type  PHASE  env-var  use_dev  dev_mnyr  dev_mxyr  dev_PH  Block  Blk_Fxn  #  parm_name
  20  50        34     35    0.5        0      4        0        0         0         0       0      1        2  #  SizeSel_PEAK_REC_CPFV_N (1)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_REC_CPFV_N (1)
   1  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_REC_CPFV_N (1)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_REC_CPFV_N (1)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_REC_CPFV_N (1)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_COM_REC_CPFV_N (1)

  20  50        35     35    0.5        0      4        0        0         0         0       0      2        2  #  SizeSel_PEAK_REC_PRIV_N_MRFSS (2)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_REC_PRIV_N_MRFSS (2)
   1  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_REC_PRIV_N_MRFSS (2)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_REC_PRIV_N_MRFSS (2)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_REC_PRIV_N_MRFSS (2)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_COM_REC_PRIV_N_MRFSS (2)

  -2   0        -1      1    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P1_REC_PRIV_N_CRFS (3)
  -2   0        -1     31    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P2_REC_PRIV_N_CRFS (3)

  14  30        22     21    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_REC_DISC_N (4)
 -12   0       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_REC_DISC_N (4)
   1  10       3.5    3.5    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_REC_DISC_N (4)
   1  10       3.5    3.5    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_REC_DISC_N (4)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_REC_DISC_N (4)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_REC_DISC_N (4)

  20  50        39     39    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_COM_HKL_N (5)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_COM_HKL_N (5)
   1  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_COM_HKL_N (5)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_COM_HKL_N (5)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_COM_HKL_N (5)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_COM_HKL_N (5)

  20  50        42     41    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_COM_NET_N (6)
 -12   0       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_COM_NET_N (6)
   1  10         3      3    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_COM_NET_N (6)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_COM_NET_N (6)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_COM_NET_N (6)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_COM_NET_N (6)

  -2   0        -1      1    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P1_COM_OTH_N (7)
  -2   0        -1     31    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P2_COM_OTH_N (7)

  14  50        27     25    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_COM_DISC_N (8)
 -12   0       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_COM_DISC_N (8)
0.01  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_COM_DISC_N (8)
0.01  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_COM_DISC_N (8)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_DISC_N (8)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_COM_DISC_N (8)

  20  50        30     29    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_ONBOARD_CPFV_CenCA_N (9)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_ONBOARD_CPFV_CenCA_N (9)
   1  10         4      4    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_ONBOARD_CPFV_CenCA_N (9)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_ONBOARD_CPFV_CenCA_N (9)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_ONBOARD_CPFV_CenCA_N (9)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_ONBOARD_CPFV_CenCA_N (9)

  -2   0        -1      1    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P1_ONBOARD_CPFV_01_16_N (11)
  -2   0        -1     31    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_P2_ONBOARD_CPFV_01_16_N (11)

  20  50        33     31    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_SCHMIDT_N (12)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_SCHMIDT_N (12)
   1  10         5      5    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_SCHMIDT_N (12)
   1  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_SCHMIDT_N (12)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_SCHMIDT_N (12)
 -11  -9        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_SCHMIDT_N (12)

  20  50        34     31    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_PEAK_ABRAMS_N (13)
 -12   0        -9     -9    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_TOP_ABRAMS_N (13)
   1  10         5      5    0.5        0      4        0        0         0         0       0      0        0  #  SizeSel_ASC-WIDTH_ABRAMS_N (13)
0.01  10         8      8    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_DSC-WIDTH_ABRAMS_N (13)
 -11  -9       -10    -10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_INIT_ABRAMS_N (13)
 -10  10        10     10    0.5        0     -4        0        0         0         0       0      0        0  #  SizeSel_FINAL_ABRAMS_N (13)

### AGE SELEX PARAMETERS ###
# NOTE: LENGTH-BASED SELEX EXCLUDE AGE-0 FISH BY DEFAULT; ADD AGE-BASED PATTERN 11 TO INCLUDE YOY IN JUV. INDICES
# LO  HI      INIT  PRIOR  PR_SD  PR_type  PHASE  env-var  use_dev  dev_mnyr  dev_mxyr  dev_PH  Block  Blk_Fxn  #  parm_name
   0   0         0      0     10        0     -4        0        0         0         0       0      0        0  #  AgeSel_P1_JUV_SWFSC_N(10)
   0   0         0      0     10        0     -4        0        0         0         0       0      0        0  #  AgeSel_P2_JUV_SWFSC_N(10)

0  #_ 0/1 to request experimental 2D_AR selectivity smoother options.
0  # TG_custom:  0=no read; 1=read if tags exist
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
###
### WEIGHTS FOR LENGTH COMPS BASED ON FRANCIS METHOD
### CAAL WEIGHTS BASED ON HARMONIC MEAN
###
#_Factor  Fleet  Value
 4        1      0.13
 4        2      0.07
 4        3      0.035 # CRFS sample sizes are MUCH larger than MRFSS
 4        4      0.199
 4        5      0.326
 4        6      1
 4        7      1
 4        8      0.35
 4        9      0.62
 4       10      1
 4       11      1
 4       12      0.27
 4       13      0.81
 5        1      0.24
 5       12      0.32
 5       13      0.21
 -9999    1      0 # terminator
#
1 #_maxlambdaphase
1 #_sd_offset
#
#### LAMBDAS (LIKELIHOOD MULTIPLIERS)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age;
#                   6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 10=recrdev;
#                   11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp;
#                   16=Tag-negbin; 17=F_ballpark
#like_comp fleet  phase  value  sizefreq_method
# base model excludes SSB_Survey_2017
    1 14  1  0  1
-9999  1  1  1  1  #  terminator
#
0 # (0/1) read specs for more stddev reporting 
# 0 1 -1 5 1 5 1 -1 5 # placeholder for selex type, len/age, year, N selex bins, Growth pattern, N growth ages, NatAge_area(-1 for all), NatAge_yr, N Natages
# placeholder for vector of selex bins to be reported
# placeholder for vector of growth ages to be reported
# placeholder for vector of NatAges ages to be reported
999
