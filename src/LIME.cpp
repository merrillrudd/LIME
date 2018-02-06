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

#define TMB_LIB_INIT R_init_LIME
#include <TMB.hpp>

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
    DATA_INTEGER(n_fl); // number of fleets
    DATA_INTEGER(n_a); // number of ages

    // Data in likelihood
    DATA_ARRAY(LF_tlf); // length composition
    DATA_MATRIX(n_lc_ft); // number of independent observation times annually, likely between 1 and C_t
    DATA_MATRIX(I_ft); // CPUE for each year
    DATA_MATRIX(C_ft); // catch each year
    DATA_INTEGER(C_opt); // if C_opt=0, no catch data, if C_opt=1, numbers, if C_opt=2, biomass
    DATA_MATRIX(ML_ft); // mean length each year

    // Known values
    DATA_VECTOR(ages); // ages
    DATA_VECTOR(L_a); // length-at-age
    DATA_VECTOR(W_a); // weight-at-age
    DATA_SCALAR(M); // natural mortality
    DATA_SCALAR(h); // steepness
    DATA_VECTOR(Mat_a); // maturity-at-age
    DATA_VECTOR(lbhighs); // upper length bins
    DATA_VECTOR(lbmids); // mid length bins

    // penalties
    DATA_INTEGER(Fpen); // penalty on annual fishing mortality
    DATA_INTEGER(SigRpen); // penalty on sigmaR
    DATA_VECTOR(SigRprior); //mean and standard deviation for sigmaR prior

    // option for fixed time series for selectivity
    DATA_IVECTOR(selex_type_f); //0 =fixed, 1=estimated

    // option for likelihood distribution for length comps
    DATA_INTEGER(LFdist); // 0=multinomial, 1=dirichlet-multinomial

    // option for shorter time-step than years
    DATA_IVECTOR(S_yrs); // matching each time step with a year
    DATA_INTEGER(n_s); // number of time steps within a year
    DATA_INTEGER(n_y); // number of years

  // ======== Parameters =================================
    // Fixed, estimated parameters
    PARAMETER_MATRIX(log_F_ft);  // fishing mortality by fleet
    PARAMETER_VECTOR(log_q_f); // catachability by fleet associated with index
    PARAMETER(beta); // equilibrium recruitment
    PARAMETER(log_sigma_R); // recruitment variation
    PARAMETER_VECTOR(log_S50_f); // Length at 50% selectivity
    PARAMETER_VECTOR(log_Sdelta_f); // log(L95 - L50)
    PARAMETER(log_sigma_F) // fishing mortality standard deviation
    PARAMETER(log_sigma_C); // observation error - catch
    PARAMETER(log_sigma_I); // observation error - index
    PARAMETER(log_CV_L); // coefficient of variation in the age-length curve
    PARAMETER_VECTOR(log_theta); // dirichlet-multinomial parameter


    // Random effects
    PARAMETER_VECTOR(Nu_input); // temporal variation in recruitment


  // ============ Global values ============================

  using namespace density;
  Type jnll=0;
  vector<Type> jnll_comp(7);
  jnll_comp.setZero();

  // ======= Transform parameters =========================

  // variation terms
  Type sigma_F = exp(log_sigma_F);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_C = exp(log_sigma_C);
  Type sigma_I = exp(log_sigma_I);
  Type CV_L = exp(log_CV_L);

  // selectivity
  vector<Type> S50_f(n_fl);
  vector<Type> Sdelta_f(n_fl);
  vector<Type> S95_f(n_fl);
  for(int f=0;f<n_fl;f++){
    S50_f(f) = exp(log_S50_f(f));
    Sdelta_f(f) = exp(log_Sdelta_f(f));
    S95_f(f) = S50_f(f) + Sdelta_f(f);
  }

  // dirichlet-multinomial parameter
  vector<Type> theta(n_fl);
  for(int f=0;f<n_fl;f++){
    theta(f) = exp(log_theta(f));
  }

  // catchability coefficient
  vector<Type> q_f(n_fl);
  for(int f=0;f<n_fl;f++){
    q_f(f) = exp(log_q_f(f));
  } 

  // ========= Probability of being length at age  =============================

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

  // ========= By fleet  =============================

  // Transform vectors
  matrix<Type> F_ft(n_fl,n_t); //number of total time steps
  matrix<Type> F_fy(n_fl,n_y); //number of years
  for(int f=0;f<n_fl;f++){
    for(int y=0;y<n_y;y++){
      F_fy(f,y) = exp(log_F_ft(f,y));
    }
  }

  int tmp; 
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      tmp = S_yrs(t) - 1;
      // match fishing mortality in the total number of time steps to the annual parameter F
      F_ft(f,t) = F_fy(f,tmp)/n_s;
    }
  }

  //total F
  vector<Type> F_t(n_t);
  F_t.setZero();
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      F_t(t) += F_ft(f,t);
    }
  }
  vector<Type> F_y(n_y);
  F_y.setZero();
  for(int f=0;f<n_fl;f++){
    for(int y=0;y<n_y;y++){
      F_y(y) += F_fy(f,y);
    }
  }

  //selectivity at length
  matrix<Type> S_fl(n_fl,n_lb);
  S_fl.setZero();
  for(int f=0;f<n_fl;f++){
    for(int l=0;l<n_lb;l++){
      if(selex_type_f(f)==1) S_fl(f,l) = 1 / (1 + exp(-log(Type(19))*(lbmids(l) - S50_f(f))/(S95_f(f) - S50_f(f))));
    }   
  }


  // input age
  matrix<Type> S_fa(n_fl,n_a);
  S_fa.setZero();
  vector<Type> sub_plba(n_lb);
  sub_plba.setZero();
  for(int f=0;f<n_fl;f++){
    for(int a=0;a<n_a;a++){
      sub_plba = plba.row(a);
      for(int l=0;l<n_lb;l++){
        S_fa(f,a) += sub_plba(l)*S_fl(f,l);
      }
      sub_plba.setZero();
    }     
  }
          
  // ============ Probability of random effects =============
  jnll_comp(0) = 0;
  for(int y=0;y<n_y;y++){
    jnll_comp(0) -= dnorm(Nu_input(y), Type(0.0), sigma_R, true);
  }

  // ============ equilibrium spawning biomass ===============
  Type SB0 = 0;
  for(int a=1;a<n_a;a++){
    SB0 += exp(beta) * exp(-M*Type(a)) * W_a(a) * Mat_a(a);
  }
  
  // // ============ joint F rate including selectivity ===========================

  // calculate fishing mortality at age over time by fleet
  // includes fishing pressure by fleet and selectivity
  // may need to include fleet percentage?
  array<Type> F_atf(n_a,n_t,n_fl);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        F_atf(a,t,f) = F_ft(f,t) * S_fa(f,a);
      }
    }
  }


  // =========================== Population dynamics ===============================  

  // ============ initialize =============  
  vector<Type> R_t(n_t); //recruitment
  R_t.setZero();
  R_t(0) = exp(beta) * exp(Nu_input(0) - pow(sigma_R,2)/Type(2));

  // combine total fishing mortality at age over time across fleets
  matrix<Type> F_ta(n_t,n_a);
  F_ta.setZero();
  for(int f=0;f<n_fl;f++){
    for(int a=0;a<n_a;a++){
      F_ta(0,a) += F_atf(a,0,f);     
    }
  }

  //over time by age
  matrix<Type> N_ta(n_t,n_a); // abundance
  matrix<Type> SB_ta(n_t,n_a); // spawning biomass
  matrix<Type> TB_ta(n_t,n_a); // total biomass
  matrix<Type> Cn_ta(n_t,n_a); //catch in numbers
  matrix<Type> Cw_ta(n_t,n_t); //catch in biomass
  N_ta.setZero();
  SB_ta.setZero();
  TB_ta.setZero();
  Cn_ta.setZero();
  Cw_ta.setZero();

  //over time
  vector<Type> N_t(n_t); // abundance
  vector<Type> SB_t(n_t); //spawning biomass
  vector<Type> TB_t(n_t); //total biomass
  vector<Type> Cn_t_hat(n_t); //catch in numbers
  vector<Type> Cw_t_hat(n_t); //catch in biomass
  N_t.setZero();
  SB_t.setZero();
  TB_t.setZero();
  Cn_t_hat.setZero();
  Cw_t_hat.setZero();

  for(int a=0;a<n_a;a++){
    // Population abundance
    if(a==0) N_ta(0,a) = R_t(0);
    if((a>=1) & (a<(n_a-1))) N_ta(0,a) = N_ta(0,a-1) * exp(-M - F_ta(0,a));
    if(a==(n_a-1)) N_ta(0,a) = (N_ta(0,a-1) * exp(-M - F_ta(0,a))) / (1 - exp(-M - F_ta(0,a)));

    // Spawning biomass
    SB_ta(0,a) = N_ta(0,a) * Mat_a(a) * W_a(a);

    // Total biomass
    TB_ta(0,a) = N_ta(0,a) * W_a(a);

    // Catch
    Cn_ta(0,a) = N_ta(0,a) * (Type(1.0) - exp(-M - F_ta(0,a))) * (F_ta(0,a))/(M + F_ta(0,a));
    Cw_ta(0,a) = Cn_ta(0,a) * W_a(a);

    //Annual values
    if(a>0) N_t(0) += N_ta(0,a);
    SB_t(0) += SB_ta(0,a);
    TB_t(0) += TB_ta(0,a);
    Cn_t_hat(0) += Cn_ta(0,a);
    Cw_t_hat(0) += Cw_ta(0,a);
  }

  //catch by fleet
  array<Type> Cn_taf(n_t,n_a,n_fl);
  array<Type> Cw_taf(n_t,n_a,n_fl);
  matrix<Type> Cn_ft(n_fl,n_t);
  matrix<Type> Cw_ft(n_fl,n_t);
  Cn_taf.setZero();
  Cw_taf.setZero();
  Cn_ft.setZero();
  Cw_ft.setZero();

  for(int f=0;f<n_fl;f++){
    for(int a=0;a<n_a;a++){
      Cn_taf(0,a,f) = N_ta(0,a) * (Type(1.0) - exp(-M - F_atf(a,0,f))) * (F_atf(a,0,f))/(M + F_atf(a,0,f));
      Cw_taf(0,a,f) = W_a(a) * N_ta(0,a) * (Type(1.0) - exp(-M - F_atf(a,0,f))) * (F_atf(a,0,f))/(M + F_atf(a,0,f));

      Cn_ft(f,0) += Cn_taf(0,a,f);
      Cw_ft(f,0) += Cw_taf(0,a,f);
    }
  }


  // ============ project forward in time =============  
  for(int t=1;t<n_t;t++){
    // Recruitment
    R_t(t) = ((4 * h * exp(beta) * SB_t(t-1)) / (SB0 * (1-h) + SB_t(t-1) * (5*h-1))) * exp(Nu_input(S_yrs(t)-1) - pow(sigma_R,2)/Type(2));
    
    // Age-structured dynamics
    for(int a=0;a<n_a;a++){
      
      // Population abundance
      if((t>=1) & (a==0)) N_ta(t,a) = R_t(t);
      if((t>=1) & (a>=1) & (a<(n_a-1))) N_ta(t,a) = N_ta(t-1,a-1) * exp(-M - F_ta(t-1,a-1)); 
      if((t>=1) & (a==(n_a-1))) N_ta(t,a) = (N_ta(t-1,a-1) * exp(-M - F_ta(t-1,a-1))) + (N_ta(t-1,a) * exp(-M - F_ta(t-1,a)));

      // Spawning biomass
      SB_ta(t,a) = N_ta(t,a)*Mat_a(a)*W_a(a);

      // Total biomass
      TB_ta(t,a) = N_ta(t,a)*W_a(a);
      
      // Catch
      Cn_ta(t,a) = N_ta(t,a) * (Type(1.0) - exp(-M - F_ta(t,a))) * (F_ta(t,a) / (M + F_ta(t,a)));
      Cw_ta(t,a) = Cn_ta(t,a) * W_a(a);

      //Annual values
      if(a>0) N_t(t) += N_ta(t,a);
      SB_t(t) += SB_ta(t,a);
      TB_t(t) += TB_ta(t,a);
      Cn_t_hat(t) += Cn_ta(t,a);
      Cw_t_hat(t) += Cw_ta(t,a);
    }
  }

  //catch by fleet
  for(int f=0;f<n_fl;f++){
    for(int t=1;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        Cn_taf(t,a,f) = N_ta(t,a) * (Type(1.0) - exp(-M - F_atf(a,t,f))) * (F_atf(a,t,f))/(M + F_atf(a,t,f));
        Cw_taf(t,a,f) = Cn_taf(t,a,f) * W_a(a);  

        Cn_ft(f,t) += Cn_taf(t,a,f);
        Cw_ft(f,t) += Cw_taf(t,a,f);
      }      
    }
  }

  // // ========== Length composition ================================

  // probability of being harvested at an age
  array<Type> page(n_t,n_a,n_fl);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        page(t,a,f) = (N_ta(t,a)*S_fa(f,a))/(N_t(t));
      }
    }    
  }

  // probability of sampling given length bin
  matrix<Type> plb_init(n_t,n_lb);
  matrix<Type> page_temp(n_t,n_a);
  matrix<Type> plb_sums(n_fl,n_t);
  array<Type> plb(n_t,n_lb,n_fl);

  for(int f=0;f<n_fl;f++){
    //find page for fleet f saved as page_temp
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
          page_temp(t,a) = page(t,a,f);
      }
    }
    //for each fleet f:
    plb_init = page_temp*plba;
    for(int t=0;t<n_t;t++){
      plb_sums(f,t) = 0;
      for(int l=0;l<n_lb;l++){
        if(plb_init(t,l)==0) plb_init(t,l) = 1e-20;
        plb_sums(f,t) += plb_init(t,l);
      }
    }
    for(int t=0;t<n_t;t++){
      for(int l=0;l<n_lb;l++){
        plb(t,l,f) = plb_init(t,l)/plb_sums(f,t);
      }
    }    
  }

  // expected mean length each year
  matrix<Type> ML_ft_hat(n_fl,n_t);
  matrix<Type> Vul_pop(n_fl,n_t);
  Type temp = 0;
  Type temp2 = 0;
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      for(int a=0;a<n_a;a++){
        temp2 += N_ta(t,a)*S_fa(f,a);
      }
      Vul_pop(f,t) = temp2;
      temp2 = 0;
      for(int l=0;l<n_lb;l++){
         temp += Vul_pop(f,t)*plb(t,l,f)*lbmids(l);
      }
      ML_ft_hat(f,t) = temp/Vul_pop(f,t);
      temp = 0;
    }    
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
      if((a>0) & (a<(n_a-1))){
        Na0(t,a) = Na0(t,a-1)*exp(-M);
        Naf(t,a) = Naf(t,a-1)*exp(-M-F_ta(t,a-1));
      }
      if(a==(n_a-1)){
        Na0(t,a) = (Na0(t,a-1)*exp(-M))/(1-exp(-M));
        Naf(t,a) = (Naf(t,a-1)*exp(-M-F_ta(t,a)))/(1-exp(-M-F_ta(t,a)));
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
  matrix<Type> log_pL_ft(n_fl,n_t);
  log_pL_ft.setZero();
  vector<Type> neff(n_fl);
  neff.setZero();
  Type checklc;
    vector<Type> LFprob(n_lb);
    vector<Type> LFraw(n_lb);
    vector<Type> prob(n_lb);
    vector<Type> sum1(n_t);
    vector<Type> sum2(n_t);
    sum1.setZero();
    sum2.setZero();

  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
      checklc = 0;
      for(int l=0;l<n_lb;l++){
        prob(l) = plb(t,l,f);
        LFraw(l) = LF_tlf(t,l,f);
        checklc += LFraw(l);
      }
      if(checklc > 0){
        if(LFdist==0){
          LFprob = n_lc_ft(f,t)*(LFraw/LFraw.sum());
          log_pL_ft(f,t) += dmultinom(LFprob, prob, true);
        }
        if(LFdist==1){
          LFprob = (LFraw/LFraw.sum());
          for(int l=0;l<n_lb;l++){
            sum1(t) += lgamma(LFraw.sum()*LFprob(l)+1);
            sum2(t) += lgamma(LFraw.sum()*LFprob(l) + theta(f)*LFraw.sum()*prob(l)) - lgamma(theta(f)*LFraw.sum()*prob(l));
          }
          log_pL_ft(f,t) += lgamma(LFraw.sum()+1) - sum1(t) + lgamma(theta(f)*LFraw.sum()) - lgamma(LFraw.sum() + theta(f)*LFraw.sum()) + sum2(t);
        }
      }
    }
  }

  if(LFdist==1){
    for(int f=0;f<n_fl;f++){
      for(int t=0;t<n_t;t++){
        for(int l=0;l<n_lb;l++){
          LFraw(l) = LF_tlf(t,l,f);
        }
      }
      neff(f) = (1 + theta(f)*LFraw.sum())/(1+theta(f));
    }
  }

  matrix<Type> I_ft_hat(n_fl,n_t);
  for(int f=0;f<n_fl;f++){
    for(int t=0;t<n_t;t++){
        I_ft_hat(f,t) = q_f(f)*TB_t(t);
    }
  }

  // matrix<Type> log_pI_ft(n_fl,n_t);
  // log_pI_ft.setZero();
  // Type n_i;
  // n_i = I_ft.size()
  // if(n_i > 1){
  //   for(int f=0;f<n_fl;f++){
  //     for(int t=0;t<n_t;t++){
  //         // probability of index at that sample
  //         log_pI_ft(t) += dlognorm( I_ft(i), log(I_ft_hat(t)), sigma_I, true);
  //     }
  //   }
  // }

  // vector<Type> log_pC_t(n_t);
  // log_pC_t.setZero();
  // if(n_c>0){
  //   for(int t=0;t<n_t;t++){  
  //     log_pC_t(t) = 0;
  //     for(int c=0;c<n_c;c++){
  //       if(C_yrs(c)==T_yrs(t)){
  //         // probability of index at that sample
  //         if(C_opt==1) log_pC_t(t) += dlognorm( C_t(c), log(C_t_hat(t)), sigma_C, true);
  //         if(C_opt==2) log_pC_t(t) += dlognorm( C_t(c), log(Cw_t_hat(t)), sigma_C, true);
  //         }
  //       }
  //     }
  //   }



  // vector<Type> log_pML_t(n_t);
  // log_pML_t.setZero();
  // // int ml;
  // if(n_ml>0){
  //   for(int t=0;t<n_t;t++){
  //     log_pML_t(t) = 0;
  //     for(int ml=0;ml<n_ml;ml++){
  //        if(ML_yrs(ml)==T_yrs(t)){
  //           log_pML_t(t) += dnorm( ML_t(ml), L_t_hat(t), L_t_hat(t)*CV_L, true);
  //        }
  //     }
  //   }
  // }

  jnll_comp(1) = 0;
  jnll_comp(1) = Type(-1)*sum( log_pL_ft );
  // jnll_comp(2) = 0;
  // if(n_ml>0) jnll_comp(2) = Type(-1)*sum( log_pML_t ); 

  // jnll_comp(3) = 0;
  // if(n_i>0) jnll_comp(3) = Type(-1)*sum( log_pI_t );
  // jnll_comp(4) = 0;
  // if(n_c>0) jnll_comp(4) = Type(-1)*sum( log_pC_t );

  // Likelihood contributions from site-aggregated observations
    vector<Type> D_t(n_t);
    for(int t=0;t<n_t;t++){
       D_t(t) = SB_t(t)/SB0;
    }

    vector<Type> lN_t(n_t);
    vector<Type> lSB_t(n_t);
    vector<Type> lTB_t(n_t);
    vector<Type> lR_t(n_t);
    vector<Type> lF_t(n_t);
    // vector<Type> lC_t(n_t);
    // vector<Type> lI_t(n_t);
    vector<Type> lD_t(n_t);
    for(int t=0;t<n_t;t++){
      lN_t(t) = log(N_t(t));
      lSB_t(t) = log(SB_t(t));
      lTB_t(t) = log(TB_t(t));
      lR_t(t) = log(R_t(t));
      lF_t(t) = log(F_t(t));
      // if(C_opt==0) lC_t(t) = log(Cw_t_hat(t));
      // if(C_opt==1) lC_t(t) = log(Cn_t_hat(t));
      // if(C_opt==2) lC_t(t) = log(Cw_t_hat(t));
      // lI_t(t) = log(I_ft_hat(t));
      lD_t(t) = log(D_t(t));
    }

    vector<Type> lF_y(n_y);
    for(int t=0;t<n_y;t++){
      lF_y(t) = log(F_y(t));
    }

    // F
    jnll_comp(5) = 0;
    if(Fpen==1){
        for(int y=1;y<n_y;y++) jnll_comp(5) -= dnorm(F_y(y), F_y(y-1), sigma_F, true);
    }

    // SigmaR
    Type sigrp;
    sigrp = 0;
    jnll_comp(6) = 0;
    sigrp = dlognorm(sigma_R, log(SigRprior(0)), SigRprior(1), true);
    if(SigRpen==1) jnll_comp(6) = Type(-1)*sigrp;

    jnll = sum(jnll_comp);

  // // ============ Reporting section ======================================

  // ADREPORT( lC_t );
  ADREPORT( lN_t );
  ADREPORT( lR_t );
  // ADREPORT( lI_t );
  ADREPORT( ML_ft_hat );
  ADREPORT( lSB_t );
  ADREPORT( lF_t );
  ADREPORT( lF_y );
  ADREPORT( lD_t );
  ADREPORT( SPR_t );
  ADREPORT( S50_f );
  ADREPORT( S_fl );
  ADREPORT( F_ft );
  ADREPORT( F_fy );

  // Parameters
  REPORT( q_f );
  REPORT( sigma_F );
  REPORT( beta );
  REPORT( sigma_R );
  REPORT( log_sigma_R );
  REPORT( S50_f );
  REPORT( S95_f );
  REPORT( S_fa );
  REPORT( S_fl );
  REPORT( sigma_C );
  REPORT( sigma_I );
  REPORT( CV_L );
  REPORT( SPR_t );
  REPORT( SB0_t );
  REPORT( SBf_t );

  // Random effects
  REPORT( Nu_input );

   // State variables
  REPORT( R_t );
  REPORT( F_t );
  REPORT( F_y );
  REPORT( N_t );
  REPORT( F_ft );
  REPORT( F_fy );

  // Predicted quantities
  REPORT( ML_ft_hat );
  REPORT( I_ft_hat );
  // Derived quantities
  REPORT( Cn_t_hat );
  REPORT( Cw_t_hat );
  REPORT( SB_t );
  REPORT( TB_t );
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
  // REPORT(log_pC_t);
  // REPORT(log_pI_t);
  REPORT(log_pL_ft);
  REPORT(neff);
  REPORT(theta);
  // REPORT(log_pML_ft);
  REPORT(sigrp);
  REPORT(jnll_comp);
  REPORT(jnll); 
  return(jnll);
}
