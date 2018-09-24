// this is a simplified version of Stan's ADMB model which was used to
// match the TMB version

DATA_SECTION
  init_int ntows // init indicates that values are in .dat file, it declares variable and gets the data as well
  init_int nlayers
  init_matrix SA(1,ntows,1,nlayers) //matrix with ES-60 sa data
  vector sum_SA1(1,ntows); // this just declares variable
  vector sum_SA2(1,ntows);
  
  !!cout <<sum_SA1<<endl;//exit(1);
  
  init_matrix Y(1,ntows,1,19)  // matrix of data by haul
  vector BSA(1,ntows) // BT sa
  vector BD(1,ntows) // bottom depth
  vector ST(1,ntows) // surface temperature
  vector BT(1,ntows) // bottom temperature
  vector CC(1,ntows) // cross current
  vector PC(1,ntows) // parallel current
  vector PP(1,ntows) // grain size, predicted phi
  vector BL(1,ntows) // near_bottom light
  vector FL(1,ntows) // mean fish length
  vector TC(1,ntows) //total current
  vector SL(1,ntows) //local slope
  ivector h(1,nlayers) // ???
  
  !! ad_comm::change_datafile_name("htrial.dat"); // ???
  init_int h1 // ???
  init_int h2 // ???

 LOCAL_CALCS
  h.fill_seqadd(1,nlayers);   // h_trial = 6;
  for (int i=1;i<=ntows;i++) sum_SA1(i) = sum(SA(i)(3,h1)); //???
  for (int i=1;i<=ntows;i++) sum_SA2(i) = sum(SA(i)(3,h2)); //???
  //for (int i=1;i<=ntows;i++) SL(i) = log(SA(i)(5)+.001)-log(SA(i)(3)+.001); //???
  BSA = column(Y,18);  // Extracted bottom sA
  BD = column(Y,11);  // Extracted bottom depth data
  ST = column(Y,12);  // Extracted surface temperatue data
  BT = column(Y,13);  // Extracted bottom temperature data
  CC = column(Y,14);  // Extracted 
  PC = column(Y,15);  // Extracted 
  PP = column(Y,16);  // Extracted 
  BL = column(Y,17);  // Extracted 
  FL = column(Y,19);
  TC = pow(elem_prod(square(CC),square(PC)),0.5);
 END_CALCS


INITIALIZATION_SECTION
  //log_q -1.0 // ???
  log_q 2.5
  log_a 8.5
  log_c 3.1
  d2 -2.1
  b_BD 0
  cb_BD 0
  logSigma -.3 

PARAMETER_SECTION
  //init_number log_q(1)
  init_number log_q
  number q
  init_bounded_number log_a(1,10)
  number a
  init_number d2
  init_number b_BD
  init_bounded_number b_ST(-10,10,-1)
  init_bounded_number b_BT(-10,10,-1)
  init_bounded_number b_CC(-10,10,-1)
  init_bounded_number b_PC(-10,10,-1)
  init_bounded_number b_PP(-10,10,-1)
  init_bounded_number b_BL(-10,10,-1)
  init_bounded_number b_FL(-10,10,-1)
  init_bounded_number log_c(-1,25,1)
  number c
  init_number logSigma
  sdreport_number sigmasq
  init_bounded_number b_TC(-10,10,-1)
  init_number cb_BD
  init_bounded_number cb_ST(-10,10,-1)
  init_bounded_number cb_BT(-10,10,-1)
  init_bounded_number cb_CC(-10,10,-1)
  init_bounded_number cb_PC(-10,10,-1)
  init_bounded_number cb_PP(-10,10,-1)
  init_bounded_number cb_BL(-10,10,-1)
  init_bounded_number cb_FL(-10,10,-1)
  init_bounded_number cb_TC(-10,10,-1)

  sdreport_number q2
  sdreport_number a2
  sdreport_number c2

  vector BSA_hat(1,ntows)  // ???
  objective_function_value nll

PROCEDURE_SECTION
   q = mfexp(log_q);
   c = mfexp(log_c);
   a = mfexp(log_a);

   q2 = q;
   a2 = a;
   c2 = c;
   
   dvar_vector eta(1,ntows);
   dvar_vector eta_c(1,ntows);
  // dvar_vector eta_q(1,ntows);
   eta = d2+b_BD*BD + b_ST*ST + b_BT*BT + b_CC*CC + b_PC*PC + 
         b_PP*PP + b_BL*BL + b_FL*FL + b_TC*TC;
   eta_c = cb_BD*BD + cb_ST*ST + cb_BT*BT + cb_CC*CC + cb_PC*PC + 
         cb_PP*PP + cb_BL*BL + cb_FL*FL + cb_TC*TC;
   //eta_q = d1 + qb_BD*BD + qb_ST*ST + qb_BT*BT + qb_CC*CC + qb_PC*PC + 
   //      qb_PP*PP + qb_BL*BL + qb_FL*FL;
   // BSA_hat =
   // q*sum_SA1 + elem_prod( mfexp(eta),sum_SA2);//+mfexp(eta_c)*c);//+1/a;
   BSA_hat =1/(1/
   (q*sum_SA1 + elem_prod( mfexp(eta),sum_SA2) +  mfexp(eta_c)*c)+1/a);
  // cout << "c=" << c << endl;
     // cout<<"eta: "<<endl<<eta<<endl;
     // cout<<"eta_c: "<<endl<<eta_c<<endl;	  
   //nll = regression(log(BSA+.001),log(BSA_hat+.001));

   sigmasq=exp(2*logSigma);
   nll=0.5*(ntows*log(2*M_PI*sigmasq)+norm2(log(BSA)-log(BSA_hat))/sigmasq);   


REPORT_SECTION
 cout<< "nll=" << nll << endl;
  // //REPORT(sigmasq);
  // REPORT(BSA);
  // REPORT(BSA_hat);

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	//#define REPORT(object) report << #object "\n" << object << endl;
       	#define REPORT(object) report << object << endl;
