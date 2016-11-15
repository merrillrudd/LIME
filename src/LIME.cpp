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
    DATA_INTEGER(start_f); // year to start estimating fishing mortality
    DATA_VECTOR(T_yrs); // vector of years of all data types
    DATA_VECTOR(C_yrs); // vector of years of catch data available
    DATA_VECTOR(I_yrs); // vector of years of index data available
    DATA_VECTOR(LC_yrs); // vector of years of length comp data
    DATA_VECTOR(ML_yrs); // vector of years of mean length data available
    DATA_VECTOR(obs_per_yr); // number of independent observation times annually, likely between 1 and C_t
    DATA_INTEGER(RecType); //0==global median recruitment, 1==Beverton-Holt stock recruit curve, 2==single value

    // Data in likelihood
    DATA_VECTOR(I_t); // CPUE for each year
    DATA_VECTOR(C_t); // catch each year
    DATA_VECTOR(ML_t); // mean length each year
    DATA_MATRIX(LF); // length composition
    DATA_INTEGER(rel_c); // relative catch, 0=no, 1=yes
    DATA_INTEGER(rel_i); // relative index, 0=no, 1=yes


    // Known values
    DATA_SCALAR(linf);
    DATA_SCALAR(vbk);
    DATA_SCALAR(t0);
    DATA_SCALAR(M);
    DATA_SCALAR(h);
    DATA_INTEGER(AgeMax);
    DATA_VECTOR(lbhighs); // upper length bins
    DATA_VECTOR(lbmids);
    DATA_VECTOR(Mat_a); // maturity
    DATA_SCALAR(lwa);
    DATA_SCALAR(lwb);

    // penalties
    DATA_INTEGER(Fpen);
    DATA_INTEGER(SigRpen);
    DATA_VECTOR(SigRprior);


  // ======== Parameters =================================
    // Fixed, estimated parameters
    PARAMETER(log_F_sd); // SD of F
    PARAMETER_VECTOR(log_F_t_input);  // F(t)
    PARAMETER(log_q_I); // catachability associated with index
    PARAMETER(beta);
    PARAMETER(log_sigma_R); // magnitude of temporal variation
    PARAMETER(logS50); // Age at 50% selectivity
    PARAMETER(log_sigma_C); // log sigma catch
    PARAMETER(log_sigma_I); // log sigma index
    PARAMETER(log_CV_L); // log sigma length comp


    // Random effects
    PARAMETER_VECTOR(Nu_input); // temporal variation in recruitment


  // ============ Global values ============================

  using namespace density;
  int a,t,lc,c,i,ml;
  Type jnll=0;
  vector<Type> jnll_comp(7);
  jnll_comp.setZero();

  // ======= Transform parameters =========================
  Type q_I = exp(log_q_I);
  Type F_sd = exp(log_F_sd);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_C = exp(log_sigma_C);
  Type sigma_I = exp(log_sigma_I);
  Type CV_L = exp(log_CV_L);
  Type S50 = exp(logS50);

  // Transform vectors
  vector<Type> F_t(n_t);
  Type F_equil = exp(log_F_t_input(start_f));
  for(int t=0;t<start_f;t++){
    if(RecType!=2) F_t(t) = exp(log_F_t_input(start_f));
    if(RecType==2) F_t(t) = F_equil;
  }
  for(int t=start_f;t<n_t;t++){
    if(RecType!=2) F_t(t) = exp(log_F_t_input(t));
    if(RecType==2) F_t(t) = F_equil;    
  }


  // ========= Convert inputs  =============================
  vector<Type> L_a(AgeMax+1);
  vector<Type> W_a(AgeMax+1);
  for(int a=0;a<=AgeMax;a++){
    L_a(a) = linf*(1-exp(-vbk*(a - t0)));
    W_a(a) = lwa*pow(L_a(a), lwb);
  }

  vector<Type> S_a(AgeMax+1);
  S_a(0) = 1e-20;
  for(int a=1;a<=AgeMax;a++){
    S_a(a) = 1 / (1 + exp(S50 - a));
  }
  

  // ============ Probability of random effects =============
  jnll_comp(0) = 0;
  for(int t=0;t<n_t;t++){
    jnll_comp(0) -= dnorm(Nu_input(t), Type(0.0), sigma_R, true);
  }

  // ============ equilibrium spawning biomass ===============
  Type SB0 = 0;
  for(int a=1;a<=AgeMax;a++){
    SB0 += exp(beta) * exp(-M*Type(a)) * W_a(a) * Mat_a(a);
  }
  
  // ============ Project dynamics ===========================

  vector<Type> R_t(n_t);
  matrix<Type> N_ta(n_t,AgeMax+1);
  matrix<Type> SB_ta(n_t,AgeMax+1);
  matrix<Type> Cn_ta(n_t,AgeMax+1);
  vector<Type> N_t(n_t);
  vector<Type> SB_t(n_t);
  vector<Type> C_t_hat(n_t);
  vector<Type> Cw_t_hat(n_t);

  // initialize
  R_t(0) = exp(beta) * exp(Nu_input(0) - pow(sigma_R,2)/Type(2));

  N_t(0) = 0;
  SB_t(0) = 0;
  C_t_hat(0) = 0;
  Cw_t_hat(0) = 0;
  for(int a=0;a<=AgeMax;a++){
    // Population abundance
    if(a==0) N_ta(0,a) = R_t(0);
    if(a>=1 & a<AgeMax) N_ta(0,a) = N_ta(0,a-1)*exp(-M-F_t(0)*S_a(a-1));
    if(a==AgeMax) N_ta(0,a) = (N_ta(0,a-1)*exp(-M-F_t(0)*S_a(a-1)))/(1-exp(-M-F_t(0)*S_a(a)));

    // Spawning biomass
    SB_ta(0,a) = N_ta(0,a)*Mat_a(a)*W_a(a);

    // Catch
    Cn_ta(0,a) = N_ta(0,a) * (Type(1.0)-exp(-M-F_t(0)*S_a(a))) * (F_t(0)*S_a(a))/(M+F_t(0)*S_a(a));

    //Annual values
    if(a>0) N_t(0) += N_ta(0,a);
    if(a>0) SB_t(0) += SB_ta(0,a);
    C_t_hat(0) += Cn_ta(0,a);
    Cw_t_hat(0) += Cn_ta(0,a)*W_a(a);
  }

  // Project forward in time
  for(int t=1;t<n_t;t++){
    // Recruitment
      // single global mean recruitment
    // if(RecType==2) R_t(t) = exp(beta) * exp(Nu_input(t) - RecDev_biasadj(t)*pow(sigma_R,2)/Type(2));
      // Beverton-Holt
    R_t(t) = ((4*h*exp(beta)*SB_t(t-1)) / (SB0*(1-h)+SB_t(t-1)*(5*h-1))) * exp(Nu_input(t) - pow(sigma_R,2)/Type(2));
    
    // Age-structured dynamics
    N_t(t) = 0;
    SB_t(t) = 0;
    C_t_hat(t) = 0;
    Cw_t_hat(t) = 0;
    for(int a=0;a<=AgeMax;a++){
      
      // Population abundance
      if(t>=1 & a==0) N_ta(t,a) = R_t(t);
      if(t>=1 & a>=1 & a<AgeMax) N_ta(t,a) = N_ta(t-1,a-1)*exp(-M-F_t(t-1)*S_a(a-1));
      if(t>=1 & a==AgeMax) N_ta(t,a) = (N_ta(t-1,a-1)*exp(-M-F_t(t-1)*S_a(a-1))) + (N_ta(t-1,a)*exp(-M-F_t(t-1)*S_a(a)));

      // Spawning biomass
      SB_ta(t,a) = N_ta(t,a)*Mat_a(a)*W_a(a);
      
      // Catch
      Cn_ta(t,a) = N_ta(t,a) * (Type(1.0)-exp(-M-F_t(t)*S_a(a))) * (F_t(t)*S_a(a))/(M+F_t(t)*S_a(a));

      //Annual values
      if(a>0) N_t(t) += N_ta(t,a);
      if(a>0) SB_t(t) += SB_ta(t,a);
      C_t_hat(t) += Cn_ta(t,a);
      Cw_t_hat(t) += Cn_ta(t,a)*W_a(a);
    }
  }

  // ========== Length composition ================================

  /////probability of being in a length bin given age
  matrix<Type> plba(AgeMax+1,n_lb);
  Type sum_sublast = 0;
  for(int a=0;a<=AgeMax;a++){
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

  // probability of being harvested at an age
  matrix<Type> page(n_t,AgeMax+1);
  for(int t=0;t<n_t;t++){
    for(int a=0;a<=AgeMax;a++){
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
    for(int a=0;a<=AgeMax;a++){
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

  // ========= spawning potential ratio ==============================
  matrix<Type> Na0(n_t,AgeMax+1);
  matrix<Type> Naf(n_t,AgeMax+1);
  vector<Type> SB0_t(n_t);
  vector<Type> SBf_t(n_t);
  SB0_t.setZero();
  SBf_t.setZero();
  vector<Type> SPR_t(n_t);
  for(int t=0;t<n_t;t++){
    for(int a=0;a<=AgeMax;a++){
      if(a==0){
        Na0(t,a) = 1;
        Naf(t,a) = 1;
      }
      if(a>0 & a<AgeMax){
        Na0(t,a) = Na0(t,a-1)*exp(-M);
        Naf(t,a) = Naf(t,a-1)*exp(-M-S_a(a-1)*F_t(t));
      }
      if(a==AgeMax){
        Na0(t,a) = (Na0(t,a-1)*exp(-M))/(1-exp(-M));
        Naf(t,a) = (Naf(t,a-1)*exp(-M-S_a(a-1)*F_t(t)))/(1-exp(-M-S_a(a-1)*F_t(t)));
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
  if(n_lc>0){
    vector<Type> dat(n_lb);
    vector<Type> dat1(n_lb);
    vector<Type> prob(n_lb);
    for(int t=0;t<n_t;t++){
      log_pL_t(t) = 0;
      for(int lc=0;lc<n_lc;lc++){
        if(LC_yrs(lc)==T_yrs(t)){
          dat1 = LF.row(lc);
          dat = obs_per_yr(t)*(dat1/dat1.sum());
          prob = plb.row(t);
          log_pL_t(t) += dmultinom(dat, prob, true);
        }
      }
    }
  }

  vector<Type> I_t_hat(n_t);
  for(int t=0;t<n_t;t++){
      I_t_hat(t) = q_I*SB_t(t);
  }

  // relative
  // vector<Type> I_t_hat_rel(n_t);
  // vector<Type> C_t_hat_rel(n_t);
  // Type I_last = I_t_hat(n_t-1);
  // Type C_last = C_t_hat(n_t-1); 
  // for(int t=0;t<n_t;t++){
  //     I_t_hat_rel(t) = I_t_hat(t)/I_last;
  //     C_t_hat_rel(t) = C_t_hat(t)/C_last;
  // }



  vector<Type> log_pI_t(n_t);
  log_pI_t.setZero();
  if(n_i>0){
    for(int t=0;t<n_t;t++){
      log_pI_t(t) = 0;
      for(int i=0;i<n_i;i++){
        if(I_yrs(i)==T_yrs(t)){
          // probability of index at that sample
          if(rel_i==0) log_pI_t(t) += dnorm( I_t(i), I_t_hat(t), sigma_I, true);
          // if(rel_i==1) log_pI_t(t) += dnorm( I_t(i), I_t_hat_rel(t), I_t_hat_rel(t)*CV_c, true);
          if(rel_i==2){
            log_pI_t(t) += dnorm( I_t(i), I_t_hat(t), sigma_I, true);
            // log_pI_t(t) += dnorm( I_t(i), I_t_hat_rel(t), I_t_hat_rel(t)*CV_c, true);
          }
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
          if(rel_c==0) log_pC_t(t) += dnorm( C_t(c), C_t_hat(t), sigma_C, true);
          // if(rel_c==1) log_pC_t(t) += dnorm( C_t(c), C_t_hat_rel(t), CV_c, true);
          if(rel_c==2){
            log_pC_t(t) += dnorm( C_t(c), C_t_hat(t), sigma_C, true);
            // log_pC_t(t) += dnorm( C_t(c), C_t_hat_rel(t), CV_c, true);
          }
        }
      }
    }
  }


  vector<Type> log_pML_t(n_t);
  log_pML_t.setZero();
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
    vector<Type> R_t_hat(n_t);
    vector<Type> Depl(n_t);
    for(int t=0;t<n_t;t++){
       N_t_hat(t) = N_t(t);
       R_t_hat(t) = R_t(t);
       SB_t_hat(t) = SB_t(t);
       Depl(t) = SB_t_hat(t)/SB0;
    }

    vector<Type> lN_t(n_t);
    vector<Type> lSB_t(n_t);
    vector<Type> lR_t(n_t);
    vector<Type> lF_t(n_t);
    vector<Type> lC_t(n_t);
    vector<Type> lI_t(n_t);
    vector<Type> lD_t(n_t);
    for(t=0;t<n_t;t++){
      lN_t(t) = log(N_t_hat(t));
      lSB_t(t) = log(SB_t_hat(t));
      lR_t(t) = log(R_t_hat(t));
      lF_t(t) = log(F_t(t));
      lC_t(t) = log(C_t_hat(t));
      lI_t(t) = log(I_t_hat(t));
      lD_t(t) = log(Depl(t));
    }


  // F
    jnll_comp(5) = 0;
    if(Fpen==1){
      jnll_comp(5) -= dnorm(F_t(start_f), F_equil, F_sd, true);
      for(int t=(start_f+1);t<n_t;t++) jnll_comp(5) -= dnorm(F_t(t), F_t(t-1), F_sd, true);
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
  ADREPORT( lD_t );
  ADREPORT( SPR_t );

  // Parameters
  REPORT( F_equil );
  REPORT( q_I );
  REPORT( F_sd );
  REPORT( beta );
  REPORT( sigma_R );
  REPORT( log_sigma_R );
  REPORT( S50 );
  REPORT( S_a );
  REPORT( sigma_C );
  REPORT( sigma_I );
  REPORT( CV_L );
  REPORT( SPR_t );

  // Random effects
  REPORT( Nu_input );

   // State variables
  REPORT( R_t_hat );
  REPORT( F_t );
  REPORT( N_t_hat );

  // Predicted quantities
  REPORT( L_t_hat );
  REPORT( I_t_hat );
  // REPORT( I_t_hat_rel );
  // Derived quantities
  REPORT( C_t_hat );
  REPORT( Cw_t_hat );
  // REPORT( C_t_hat_rel );
  REPORT( SB_t_hat );
  REPORT( SB0 );
  REPORT(Depl);
  REPORT(N_ta);
  REPORT(Cn_ta);
  REPORT(plba);
  REPORT(page);
  REPORT(plb);
  REPORT(W_a);
  REPORT(L_a);
  REPORT(Mat_a);
  REPORT(AgeMax);
  REPORT(M);
    // Likelihoods
  REPORT(log_pC_t);
  REPORT(log_pI_t);
  REPORT(log_pL_t);
  REPORT(log_pML_t);
  REPORT(sigrp);
  REPORT(jnll_comp);
  REPORT(jnll); 
  return(jnll);
}