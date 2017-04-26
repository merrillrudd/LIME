// Estimates more similar parameters to LB-SPR
// with more similar inputs
//
// Requires:
//   1) estimates of VB growth parameters, including CV on Linf
//   2) maturity at size
//   3) length composition of catch
// 
//  Estimates:
//   1) selectivity at length
//   2) annual F as fixed effect
//   3) annual R as random effect 

# include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ========== Inputs ============================

    // Indices
    DATA_INTEGER(n_t); // total number of years
    DATA_INTEGER(n_lb); // number of length bins
    DATA_INTEGER(n_c); // number of years of catch data available
    DATA_INTEGER(n_i); // number of years of index data available
    DATA_INTEGER(n_lc); // number of years of length composition data
    DATA_INTEGER(n_ml); // number of years of mean length data
    DATA_VECTOR(T_yrs); // vector of years of all data types
    DATA_VECTOR(C_yrs); // vector of years of catch data available
    DATA_VECTOR(I_yrs); // vector of years of index data available
    DATA_VECTOR(LC_yrs); // vector of years of length comp data
    DATA_VECTOR(ML_yrs); // vector of years of mean length data available
    DATA_VECTOR(obs_per_yr); // number of independent observation times annually, likely between 1 and C_t

    // Data in likelihood
    DATA_VECTOR(I_t); // CPUE for each year
    DATA_VECTOR(C_t); // catch each year
    DATA_INTEGER(C_opt); // if C_opt=0, no catch data, if C_opt=1, numbers, if C_opt=2, biomass
    DATA_VECTOR(ML_t); // mean length each year
    DATA_MATRIX(LF); // length composition
    DATA_INTEGER(LFdist); // distribution for length composition; 0=multinomial (and requires obs_per_year), 1=dirichlet-multinomial
    DATA_INTEGER(theta_type); //0= annual theta, 1=single theta
    DATA_VECTOR(lbhighs); // upper length bins
    DATA_VECTOR(lbmids);
    DATA_INTEGER(binwidth);

    // Known values
    DATA_INTEGER(n_a);
    DATA_VECTOR(ages);
    DATA_VECTOR(L_a);
    DATA_VECTOR(W_a);
    DATA_SCALAR(M);
    DATA_SCALAR(h);
    DATA_VECTOR(Mat_a); // maturity
    // DATA_SCALAR(ML50);

    // penalties
    DATA_INTEGER(Fpen);
    DATA_INTEGER(SigRpen);
    DATA_VECTOR(SigRprior);

    // option for fixed time series for selectivity
    DATA_VECTOR(S_l_input);
    // DATA_VECTOR(M_l_input);

    // option for shorter time-step than years
    DATA_IVECTOR(S_yrs); // matching each time step with a year
    DATA_INTEGER(n_s); // number of time steps within a year
    DATA_INTEGER(n_y); // number of years

  // ======== Parameters =================================
    // Fixed, estimated parameters
    PARAMETER(log_sigma_F); // SD of F
    PARAMETER_VECTOR(log_F_t_input);  // F(t)
    PARAMETER(log_q_I); // catachability associated with index
    PARAMETER(beta);
    PARAMETER(log_sigma_R); // magnitude of temporal variation
    PARAMETER(logS50); // Length at 50% selectivity
    PARAMETER(log_sigma_C); // log sigma catch
    PARAMETER(log_sigma_I); // log sigma index
    PARAMETER(log_CV_L); // log sigma length comp
    PARAMETER_VECTOR(log_theta); // dirichlet-multinomial parameter


    // Random effects
    PARAMETER_VECTOR(Nu_input); // temporal variation in recruitment


  // ============ Global values ============================

  using namespace density;
  Type jnll=0;
  vector<Type> jnll_comp(7);
  jnll_comp.setZero();

  // ======= Transform parameters =========================
  Type q_I = exp(log_q_I);
  Type sigma_F = exp(log_sigma_F);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_C = exp(log_sigma_C);
  Type sigma_I = exp(log_sigma_I);
  Type CV_L = exp(log_CV_L);
  Type S50 = exp(logS50);
  int amax;
  amax = n_a/n_s;

  vector<Type> theta(n_lc);
  int l;
  if(theta_type==0){
    for(int l=0;l<n_lc;l++){
      theta(l) = exp(log_theta(l));
    }
  }
  if(theta_type==1){
    for(int l=0;l<n_lc;l++){
      theta(l) = exp(log_theta(0));
    }
  }


  // Transform vectors
  vector<Type> F_t(n_t); //number of total time steps
  vector<Type> F_y(n_y); //number of years
  Type F_equil;
  F_equil = exp(log_F_t_input(0));
  int t;
  int y;
  int tmp; 
  for(int y=0;y<n_y;y++){
    F_y(y) = exp(log_F_t_input(y));
  }
  for(int t=0;t<n_t;t++){
    tmp = S_yrs(t) - 1;
    // match fishing mortality in the total number of time steps to the annual parameter F
    F_t(t) = F_y(tmp)/n_s;
  }


  // ========= Convert inputs  =============================

  /////probability of being in a length bin given INPUT age
  matrix<Type> plba(n_a,n_lb);
  Type sum_sublast = 0;
  for(int a=0;a<n_a;a++){
    for(int l=0;l<n_lb;l++){
      if(l==0){
        plba(a,l) = pnorm(lbhighs(l), L_a(a), L_a(a)*CV_L);
        sum_sublast += plba(a,l);
      }
      if(l>=1){
        if(l<(n_lb-1)){
          plba(a,l) = pnorm(lbhighs(l), L_a(a), L_a(a)*CV_L) - pnorm(lbhighs(l-1), L_a(a), L_a(a)*CV_L);
          sum_sublast += plba(a,l);
        }
        if(l==(n_lb-1)) plba(a,l) = Type(1.0) - sum_sublast;
      }
    }
    sum_sublast = 0;
  }

 
  //selectivity at length
  vector<Type> S_l(n_lb);
  S_l.setZero();
  for(int l=0;l<n_lb;l++){
    if(S_l_input(0)<0) S_l(l) = 1 / (1 + exp(S50 - lbmids(l)));
    if(S_l_input(0)>=0) S_l(l) = S_l_input(l);
  }

  // input age
  vector<Type> S_a(n_a);
  S_a.setZero();
  vector<Type> sub_plba(n_lb);
  sub_plba.setZero();
  for(int a=0;a<n_a;a++){
    sub_plba = plba.row(a);
    for(int l=0;l<n_lb;l++){
      S_a(a) += sub_plba(l)*S_l(l);
    }
    sub_plba.setZero();
  }           

  // ============ Probability of random effects =============
  jnll_comp(0) = 0;
  for(int y=0;y<n_y;y++){
    jnll_comp(0) -= dnorm(Nu_input(y), Type(0.0), sigma_R, true);
  }

  // ============ equilibrium spawning biomass ===============
  Type SB0 = 0;
  Type TB0 = 0;
  for(int a=1;a<n_a;a++){
    SB0 += exp(beta) * exp(-M*Type(a)) * W_a(a) * Mat_a(a);
    TB0 += exp(beta) * exp(-M*Type(a)) * W_a(a);
  }
  
  // ============ Project dynamics ===========================

  vector<Type> R_t(n_t); // recruitment for each time step
  matrix<Type> N_ta(n_t,n_a);
  matrix<Type> SB_ta(n_t,n_a);
  matrix<Type> TB_ta(n_t,n_a);
  matrix<Type> Cn_ta(n_t,n_a);
  vector<Type> N_t(n_t);
  vector<Type> SB_t(n_t);
  vector<Type> TB_t(n_t);
  vector<Type> C_t_hat(n_t);
  vector<Type> Cw_t_hat(n_t);
  R_t.setZero();
  N_ta.setZero();
  SB_ta.setZero();
  TB_ta.setZero();
  Cn_ta.setZero();
  N_t.setZero();
  SB_t.setZero();
  TB_t.setZero();
  C_t_hat.setZero();
  Cw_t_hat.setZero();

  // initialize
  R_t(0) = exp(beta) * exp(Nu_input(0) - pow(sigma_R,2)/Type(2));

  for(int a=0;a<n_a;a++){
    // Population abundance
    if(a==0) N_ta(0,a) = R_t(0);
    if(a>=1 & a<(n_a-1)) N_ta(0,a) = N_ta(0,a-1) * exp(-M - F_t(0) * S_a(a-1));
    if(a==(n_a-1)) N_ta(0,a) = (N_ta(0,a-1) * exp(-M - F_t(0) * S_a(a))) / (1 - exp(-M - F_t(0) * S_a(a)));

    // Spawning biomass
    SB_ta(0,a) = N_ta(0,a) * Mat_a(a) * W_a(a);

    // Total biomass
    TB_ta(0,a) = N_ta(0,a) * W_a(a);

    // Catch
    Cn_ta(0,a) = N_ta(0,a) * (Type(1.0) - exp(-M - F_t(0) * S_a(a))) * (F_t(0) * S_a(a))/(M + F_t(0) * S_a(a));

    //Annual values
    if(a>0) N_t(0) += N_ta(0,a);
    SB_t(0) += SB_ta(0,a);
    TB_t(0) += TB_ta(0,a);
    C_t_hat(0) += Cn_ta(0,a);
    Cw_t_hat(0) += Cn_ta(0,a)*W_a(a);
  }

  // Project forward in time
  for(int t=1;t<n_t;t++){
    // Recruitment
    R_t(t) = ((4 * h * exp(beta) * SB_t(t-1)) / (SB0 * (1-h) + SB_t(t-1) * (5*h-1))) * exp(Nu_input(S_yrs(t)-1) - pow(sigma_R,2)/Type(2));
    
    // Age-structured dynamics
    for(int a=0;a<n_a;a++){
      
      // Population abundance
      if(t>=1 & a==0) N_ta(t,a) = R_t(t);
      if(t>=1 & a>=1 & a<(n_a-1)) N_ta(t,a) = N_ta(t-1,a-1) * exp(-M - F_t(t-1) * S_a(a-1));
      if(t>=1 & a==(n_a-1)) N_ta(t,a) = (N_ta(t-1,a-1) * exp(-M - F_t(t-1) * S_a(a-1))) + (N_ta(t-1,a) * exp(-M - F_t(t-1) * S_a(a)));

      // Spawning biomass
      SB_ta(t,a) = N_ta(t,a)*Mat_a(a)*W_a(a);


      // Total biomass
      TB_ta(t,a) = N_ta(t,a)*W_a(a);
      
      // Catch
      Cn_ta(t,a) = N_ta(t,a) * (Type(1.0) - exp(-M - F_t(t) * S_a(a))) * (F_t(t)*S_a(a)) / (M + F_t(t) * S_a(a));

      //Annual values
      if(a>0) N_t(t) += N_ta(t,a);
      SB_t(t) += SB_ta(t,a);
      TB_t(t) += TB_ta(t,a);
      C_t_hat(t) += Cn_ta(t,a);
      Cw_t_hat(t) += Cn_ta(t,a)*W_a(a);
    }
  }

  // ========== Length composition ================================

  // probability of being harvested at an age
  matrix<Type> page(n_t,n_a);
  for(int t=0;t<n_t;t++){
    for(int a=0;a<n_a;a++){
      page(t,a) = (N_ta(t,a)*S_a(a))/(N_t(t));
    }
  }

  // probability of sampling given length bin
  matrix<Type> plb_init(n_t,n_lb);
  plb_init = page*plba;
  vector<Type> plb_sums(n_t);
  for(int t=0;t<n_t;t++){
    plb_sums(t) = 0;
    for(int l=0;l<n_lb;l++){
      if(plb_init(t,l)==0) plb_init(t,l) = 1e-20;
      plb_sums(t) += plb_init(t,l);
    }
  }
  matrix<Type> plb(n_t,n_lb);
  for(int t=0;t<n_t;t++){
    for(int l=0;l<n_lb;l++){
      plb(t,l) = plb_init(t,l)/plb_sums(t);
    }
  }

  // expected mean length each year
  vector<Type> L_t_hat(n_t);
  vector<Type> Vul_pop(n_t);
  Type temp = 0;
  Type temp2 = 0;
  for(int t=0;t<n_t;t++){
    for(int a=0;a<n_a;a++){
      temp2 += N_ta(t,a)*S_a(a);
    }
    Vul_pop(t) = temp2;
    temp2 = 0;
    for(int l=0;l<n_lb;l++){
       temp += Vul_pop(t)*plb(t,l)*lbmids(l);
    }
    L_t_hat(t) = temp/Vul_pop(t);
    temp = 0;
  }

  // // ========= spawning potential ratio ==============================
  matrix<Type> Na0(n_t,n_a);
  matrix<Type> Naf(n_t,n_a);

  // SPR
  vector<Type> SB0_t(n_t);
  vector<Type> SBf_t(n_t);
  vector<Type> SPR_t(n_t);
  SB0_t.setZero();
  SBf_t.setZero();
  SPR_t.setZero();

  for(int t=0;t<n_t;t++){
    for(int a=0;a<n_a;a++){
      if(a==0){
        Na0(t,a) = 1;
        Naf(t,a) = 1;
      }
      if(a>0 & a<(n_a-1)){
        Na0(t,a) = Na0(t,a-1)*exp(-M);
        Naf(t,a) = Naf(t,a-1)*exp(-M-S_a(a-1)*F_t(t));
      }
      if(a==(n_a-1)){
        Na0(t,a) = (Na0(t,a-1)*exp(-M))/(1-exp(-M));
        Naf(t,a) = (Naf(t,a-1)*exp(-M-S_a(a)*F_t(t)))/(1-exp(-M-S_a(a)*F_t(t)));
      }

      if(a>0){
        SB0_t(t) += Na0(t,a)*Mat_a(a)*W_a(a);
        SBf_t(t) += Naf(t,a)*Mat_a(a)*W_a(a);
      }
    }
    SPR_t(t) = SBf_t(t)/SB0_t(t);
  }

  // ========= Build likelihood ==============================

  // Likelihood contribution from observations
  vector<Type> log_pL_t(n_t);
  log_pL_t.setZero();
  vector<Type> neff(n_lc);
  neff.setZero();
    vector<Type> LFprob(n_lb);
    vector<Type> LFraw(n_lb);
    vector<Type> prob(n_lb);
    vector<Type> sum1(n_lc);
    vector<Type> sum2(n_lc);
    sum1.setZero();
    sum2.setZero();
  int lc;
  if(n_lc>0){
    for(int t=0;t<n_t;t++){
      log_pL_t(t) = 0;
      for(int lc=0;lc<n_lc;lc++){
        if(LC_yrs(lc)==T_yrs(t)){
          LFraw = LF.row(lc);
          prob = plb.row(t);
          if(LFdist==0){
            LFprob = obs_per_yr(t)*(LFraw/LFraw.sum()); 
            log_pL_t(t) += dmultinom(LFprob, prob, true); // check log
          }
          if(LFdist==1){
            LFprob = (LFraw/LFraw.sum());
            for(int l=0;l<n_lb;l++){
              sum1(lc) += lgamma(LFraw.sum()*LFprob(l)+1);
              sum2(lc) += lgamma(LFraw.sum()*LFprob(l) + theta(lc)*LFraw.sum()*prob(l)) - lgamma(theta(lc)*LFraw.sum()*prob(l));
            }
            log_pL_t(t) += lgamma(LFraw.sum()+1) - sum1(lc) + lgamma(theta(lc)*LFraw.sum()) - lgamma(LFraw.sum()+theta(lc)*LFraw.sum()) + sum2(lc);
          }
        }
      }
    }
  }

  if(LFdist==1){
    for(int lc=0;lc<n_lc;lc++){
      LFraw = LF.row(lc);
      neff(lc) = (1+theta(lc)*LFraw.sum())/(1+theta(lc));
    }
  }

  vector<Type> I_t_hat(n_t);
  for(int t=0;t<n_t;t++){
      I_t_hat(t) = q_I*TB_t(t);
  }


  vector<Type> log_pI_t(n_t);
  log_pI_t.setZero();
  if(n_i>0){
    for(int t=0;t<n_t;t++){
      log_pI_t(t) = 0;
      for(int i=0;i<n_i;i++){
        if(I_yrs(i)==T_yrs(t)){
          // probability of index at that sample
          log_pI_t(t) += dlognorm( I_t(i), log(I_t_hat(t)), sigma_I, true);
          }
        }
      }
    }

  vector<Type> log_pC_t(n_t);
  log_pC_t.setZero();
  if(n_c>0){
    for(int t=0;t<n_t;t++){  
      log_pC_t(t) = 0;
      for(int c=0;c<n_c;c++){
        if(C_yrs(c)==T_yrs(t)){
          // probability of index at that sample
          if(C_opt==1) log_pC_t(t) += dlognorm( C_t(c), log(C_t_hat(t)), sigma_C, true);
          if(C_opt==2) log_pC_t(t) += dlognorm( C_t(c), log(Cw_t_hat(t)), sigma_C, true);
          }
        }
      }
    }



  vector<Type> log_pML_t(n_t);
  log_pML_t.setZero();
  int ml;
  if(n_ml>0){
    for(int t=0;t<n_t;t++){
      log_pML_t(t) = 0;
      for(int ml=0;ml<n_ml;ml++){
         if(ML_yrs(ml)==T_yrs(t)){
            log_pML_t(t) += dnorm( ML_t(ml), L_t_hat(t), L_t_hat(t)*CV_L, true);
         }
      }
    }
  }

  jnll_comp(1) = 0;
  if(n_lc>0) jnll_comp(1) = Type(-1)*sum( log_pL_t );
  jnll_comp(2) = 0;
  if(n_ml>0) jnll_comp(2) = Type(-1)*sum( log_pML_t ); 

  jnll_comp(3) = 0;
  if(n_i>0) jnll_comp(3) = Type(-1)*sum( log_pI_t );
  jnll_comp(4) = 0;
  if(n_c>0) jnll_comp(4) = Type(-1)*sum( log_pC_t );

  // Likelihood contributions from site-aggregated observations
    vector<Type> N_t_hat(n_t);
    vector<Type> SB_t_hat(n_t);
    vector<Type> TB_t_hat(n_t);
    vector<Type> R_t_hat(n_t);
    vector<Type> D_t(n_t);
    for(int t=0;t<n_t;t++){
       N_t_hat(t) = N_t(t);
       R_t_hat(t) = R_t(t);
       SB_t_hat(t) = SB_t(t);
       TB_t_hat(t) = TB_t(t);
       D_t(t) = SB_t_hat(t)/SB_t_hat(0);
    }

    vector<Type> lN_t(n_t);
    vector<Type> lSB_t(n_t);
    vector<Type> lTB_t(n_t);
    vector<Type> lR_t(n_t);
    vector<Type> lF_t(n_t);
    vector<Type> lC_t(n_t);
    vector<Type> lI_t(n_t);
    vector<Type> lD_t(n_t);
    for(t=0;t<n_t;t++){
      lN_t(t) = log(N_t_hat(t));
      lSB_t(t) = log(SB_t_hat(t));
      lTB_t(t) = log(TB_t_hat(t));
      lR_t(t) = log(R_t_hat(t));
      lF_t(t) = log(F_t(t));
      lC_t(t) = log(C_t_hat(t));
      lI_t(t) = log(I_t_hat(t));
      lD_t(t) = log(D_t(t));
    }

    vector<Type> lF_y(n_y);
    for(t=0;t<n_y;t++){
      lF_y(t) = log(F_y(t));
    }



  // F
    jnll_comp(5) = 0;
    if(Fpen==1){
        for(int y=1;y<n_y;y++) jnll_comp(5) -= dnorm(F_y(y), F_y(y-1), sigma_F, true);
    }

    // // SigmaR
    Type sigrp;
    sigrp = 0;
    jnll_comp(6) = 0;
    sigrp = dlognorm(sigma_R, log(SigRprior(0)), SigRprior(1), true);
    if(SigRpen==1) jnll_comp(6) = Type(-1)*sigrp;

    jnll = sum(jnll_comp);

  // ============ Reporting section ======================================

  ADREPORT( lC_t );
  ADREPORT( lN_t );
  ADREPORT( lR_t );
  ADREPORT( lI_t );
  ADREPORT( L_t_hat );
  ADREPORT( lSB_t );
  ADREPORT( lF_t );
  ADREPORT( lF_y );
  ADREPORT( lD_t );
  ADREPORT( SPR_t );
  ADREPORT( S50 );
  ADREPORT( S_l );
  ADREPORT( S_A );

  // Parameters
  REPORT( F_equil );
  REPORT( q_I );
  REPORT( sigma_F );
  REPORT( beta );
  REPORT( sigma_R );
  REPORT( log_sigma_R );
  REPORT( S50 );
  REPORT( S_a );
  REPORT( S_A );
  REPORT( S_l );
  REPORT( S_l_input );
  REPORT( sigma_C );
  REPORT( sigma_I );
  REPORT( CV_L );
  REPORT( SPR_t );
  REPORT( SB0_t );
  REPORT( SBf_t );

  // Random effects
  REPORT( Nu_input );

   // State variables
  REPORT( R_t_hat );
  REPORT( F_t );
  REPORT( F_y );
  REPORT( N_t_hat );

  // Predicted quantities
  REPORT( L_t_hat );
  REPORT( I_t_hat );
  // Derived quantities
  REPORT( C_t_hat );
  REPORT( Cw_t_hat );
  REPORT( SB_t_hat );
  REPORT( TB_t_hat );
  REPORT( SB0 );
  REPORT(D_t);
  REPORT(N_ta);
  REPORT(Cn_ta);
  REPORT(plba);
  REPORT(page);
  REPORT(plb);
  REPORT(W_a);
  REPORT(L_a);
  REPORT(Mat_a);
  REPORT(M);
    // Likelihoods
  REPORT(log_pC_t);
  REPORT(log_pI_t);
  REPORT(log_pL_t);
  REPORT(neff);
  REPORT(theta);
  REPORT(log_pML_t);
  REPORT(sigrp);
  REPORT(jnll_comp);
  REPORT(jnll); 
  return(jnll);
}
