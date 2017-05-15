*! version 3.0.0 10May2017
program define rdperm, rclass
version 13

syntax varlist [if] [in]                                                        ///
               [, Cutoff(real 0)                                                ///
				  Qobservations(integer 0)                                      ///
				  CVM                                                           ///
				  Permutations(integer 999)] 
 
//DESCRIPTION
//VARLIST				  
//First variable is the running variable. 
//Second variable onwards are covariates.
//OPTIONS
//Cutoff value has default 0 and must be real.
//Q is the effective number observations with default of 0, which is rot.
//CVM specifies that cvm test statistic used for the joint test.
//Permutations is the number of permutations performed with default of 999.
//OUPUT
//P.Val is the p value

marksample touse 
gettoken running covariates :varlist                                            //parsing the running variable and covariates
tempname pval N N_l N_r h_l h_r ts q_new                                        //introducing tempname to save output of mata function

local cvm_count : word count `cvm' 
local var_count : word count `covariates' 
if `cvm_count'==1 | `var_count' == 1 {
  scalar `ts' = 0                                                               //If max specified and no. of covariates above 1, then max statistic used
}
else {
  scalar `ts' = 1                                                               //Otherwise CvM statistic is default
}

mata: rdperm_work("`covariates'", "`running'", "`touse'", `cutoff',             ///
                   `qobservations', `permutations', "`ts'", "`pval'",           ///
				   "`N'", "`N_l'", "`N_r'", "`h_l'", "`h_r'", "`q_new'")        //running the mata function using the inputs from data

return scalar pval = `pval'                                                     //storing pvalue
return scalar h_r = `h_r'                                                       //storing effective neighbourhood to the right
return scalar h_l = `h_l'                                                       //storing effective neighbourhood to the left
return scalar N_r = `N_r'                                                       //storing no. of obs. to the right
return scalar N_l = `N_l'                                                       //storing no. of obs to the left
return scalar N = `N'                                                           //storing number of observations
return scalar P = `permutations'                                                //storing number of permutations
return scalar q = `q_new'                                                       //storing effective observations used
return scalar c = `cutoff'                                                      //storing threshold value

return local cmd "rdperm"                                                       //storing command name
if `cvm_count'==1 | `var_count' == 1 {
  return local teststat "CvM"                                                   //storing test statistic used
}
else {
  return local teststat "max"                                                   //storing test statistic used
}
return local covariates `covariates'                                            //storing covariate names
return local runninvar `running'                                                //storing running variable name
				   
end


mata:

void rdperm_work(string scalar covariates_s, string scalar running_s,           ///
                 string scalar touse_s, real scalar cutoff_s,                   ///
				 real scalar q_s, real scalar P_s, string scalar ts_s,          /// 
				 string scalar pval_s, string scalar N_s, string scalar N_l_s,  /// 
				 string scalar N_r_s, string scalar h_l_s, string scalar h_r_s,  ///
				 string scalar q_new_s) {


   rseed(12345)                                                                 //Setting seed for the permutations

   ts_s              = st_numscalar(ts_s) 	                                    //  
   w                 = st_data(., tokens(covariates_s), touse_s)                //This is the vector of covariates
   z                 = st_data(., running_s, touse_s)                           //This is vector if running variables
   z_bar             = cutoff_s                                                 //This is threshold value
   q                 = ceil(q_s)                                                //Number of observations used from each side of threshold
   z                 = z :- z_bar                                               //Normalized running variable
   n_w               = rows(w)                                                  //Sample size of w
   n_z               = rows(z)                                                  //Sample size of z
   K                 = cols(w)                                                  //
   N                 = n_w                                                      //
   P                 = P_s                                                      //Number of permutations

   if (n_w != n_z) {                                                            //Error message for size of w and z not equal
     printf("size of variables do not match")
	 exit(error(3498))
   }
   
   //---------------------------------------------------------------------------
   //----Computing rule of thumb q value when specified-------------------------
   //---------------------------------------------------------------------------

   if (q == 0) {
     sz              = sqrt(variance(z))                                        //
	 mz              = mean(z)                                                  //
	 h               = 1.84*sz*(N^(-1/5))                                       //
	 f               = sum((abs(z :/ h) :<= 1) :* (1 :- abs(z :/ h)))/(N*h)     //
	 q               = J(K,1,1)                                                 //
	 for (k=1; k <= K; k++) {
	   w_temp        = w[.,k]                                                   //    
	   sw            = sqrt(variance(w_temp))                                   //
	   mw            = mean(w_temp)                                             //
	   p             = sum((z :- mz) :* (w_temp :- mw))/(N*sw*sz)               //
	   q[k,.]        = sqrt((1-p^2))*f*sz*(N^(0.9)/log(N))                    //
	 }
	 q               = ceil(colmin(q))                                          //
	 if (q < 10) {
       q             = 10                                                       //
	 }
	 else if (q > N^(0.9)/log(N)) {
	   q             = ceil(N^(0.9)/log(N))                                     //
     }
	 else {
	   q             = q                                                        //
	 }
   }

   //---------------------------------------------------------------------------
   
   z_left            = select(z, z[.,1] :< 0)                                   //Z obs to left of threshold
   w_left            = select(w, z[., 1] :< 0)                                  //W Obs to left of threshold
   n_left            = rows(z_left)                                             //No. of obs to left of threshold
   z_right           = select(z, z[.,1] :>= 0)                                  //Z obs to right of threshold
   w_right           = select(w, z[., 1] :>= 0)                                 //W Obs to right of threshold
   n_right           = rows(z_right)                                            //No. of obs to right of threshold
   
   if (n_left < q) {                                                            //Error message for size of w and z not equal
     printf("insufficient sample size for q")
	 exit(error(3498))
   }
   if (n_right < q) {                                                           //Error message for size of w and z not equal
     printf("insufficient sample size for q")
	 exit(error(3498))
   }
   
   w_left            = w_left[order(z_left,1), .]                               //Induced order of W obs to left of threshold
   z_left            = sort(z_left,1)                                           //Sorting Z obs to left of threshold
   w_right           = w_right[order(z_right,1), .]                             //Induced order of W obs to right of threshold
   z_right           = sort(z_right,1)                                          //Sorting Z obs to right of threshold

   left_max          = z_left[(n_left+1-q), 1]                                  //Effective neighbourhood from left
   right_max         = z_right[q,1]                                             //Effective neighbourhood from right

   w_left            = w_left[((n_left-q+1)..n_left), .]                        //q closest W from the left of threshold
   w_right           = w_right[(1..q), .]                                       //q closest W from the right of threshold
   w_q               = w_left \ w_right                                         //W obs used in permutations
   
   //---------------------------------------------------------------------------
   //----Computing Test Statistic-----------------------------------------------
   //---------------------------------------------------------------------------

   if (ts_s == 0) {
     A_left          = J(q,2*q,1)                                               //
	 A_right         = J(q,2*q,1)                                               //
     for (i=1; i <= K; i++) {
	   w_3           = J(q,1,w_q[.,i]')                                         //
	   w_left_1      = J(1,2*q,w_left[.,i])                                     //
	   A_temp        = w_left_1 :<= w_3                                         //
	   A_left        = A_left :* A_temp                                         //
	   w_right_1     = J(1,2*q,w_right[.,i])                                    //
	   A_temp        = w_right_1 :<= w_3                                        //
	   A_right       = A_right :* A_temp                                        //
	 }
	 m_minus         = mean(A_left)                                             //
	 m_plus          = mean(A_right)                                            //
	 TS              = mean(((m_minus - m_plus):^2)')                           //
   }
   else {
     C               = 100 - K                                                  //
	 c               = J(K,100,0) //
	 c[.,1..C]       = runiform(K,C)                                            //
	 c               = c :/ J(K,1,colsum(c))                                    //
	 for (i=1; i <= K; i++) {
	   temp          = J(K,1,0)                                                 //
	   temp[i,.]     = 1                                                        //
	   c[.,C+i]      = temp                                                     //
	 }
	 C               = 100                                                      //
	 w_temp          = w_q*c                                                    //
	 w_3             = J(q,1,rowshape(w_temp',1))                               //
	 w_left_1        = colshape(J(2*q,1,w_left*c)',q)'                          //
	 w_right_1       = colshape(J(2*q,1,w_right*c)',q)'                         //
	 A_left          = w_left_1 :<= w_3                                         //
	 A_right         = w_right_1 :<= w_3                                        //
	 m_minus         = mean(A_left)                                             //
	 m_plus          = mean(A_right)                                            //
	 TS              = colshape((m_minus - m_plus):^2,2*q)'                     //
	 TS              = max(mean(TS))                                            //
   }

   //---------------------------------------------------------------------------
   //----Computing Permutation Test---------------------------------------------
   //---------------------------------------------------------------------------

   TS_perm           = J(P,1,0)                                                 //
   if (ts_s == 0) {
     for (j=1; j <= P; j++) {
       w_perm        = jumble(w_q)                                              //Random permutation of W obs
	   w_perm_left   = w_perm[(1..q), .]                                        //
	   w_perm_right  = w_perm[((q+1)..2*q), .]                                  //   
	   A_left        = J(q,2*q,1)                                               //
	   A_right       = J(q,2*q,1)                                               //
       for (i=1; i <= K; i++) {
	     w_3         = J(q,1,w_q[.,i]')                                         //
	     w_left_1    = J(1,2*q,w_perm_left[.,i])                                //
	     A_temp      = w_left_1 :<= w_3                                         //
	     A_left      = A_left :* A_temp                                         //
	     w_right_1   = J(1,2*q,w_perm_right[.,i])                               //
	     A_temp      = w_right_1 :<= w_3                                        //
	     A_right     = A_right :* A_temp                                        //
	   }
	   m_minus       = mean(A_left)                                             //
	   m_plus        = mean(A_right)                                            //
	   TS_perm[j,1]  = mean(((m_minus - m_plus):^2)')                           //
     }
   }
   else {
     for (j=1; j <= P; j++) {
	   w_perm        = jumble(w_q)                                              //Random permutation of W obs
	   w_perm_left   = w_perm[(1..q), .]                                        //
	   w_perm_right  = w_perm[((q+1)..2*q), .]                                  //
	   w_left_1      = colshape(J(2*q,1,w_perm_left*c)',q)'                     //
	   w_right_1     = colshape(J(2*q,1,w_perm_right*c)',q)'                    //
	   A_left        = w_left_1 :<= w_3                                         //
	   A_right       = w_right_1 :<= w_3                                        //
	   m_minus       = mean(A_left)                                             //
	   m_plus        = mean(A_right)                                            //
	   TS_perm_temp  = colshape((m_minus - m_plus):^2,2*q)'                     //
	   TS_perm[j,1]  = max(mean(TS_perm_temp))                                  //
	 }
   }

   p_val             = mean((TS_perm :>= TS))                                   //Computed p value

   //---------------------------------------------------------------------------
   //----Results to view in Stata-----------------------------------------------
   //---------------------------------------------------------------------------
   
   st_numscalar(pval_s, p_val)                                                  //Storing pvalue to be viewed in stata
   st_numscalar(N_s, n_z)                                                       //Storing sample size N
   st_numscalar(N_l_s, n_left)                                                  //Storing sample size to the left
   st_numscalar(N_r_s, n_right)                                                 //Storing sample size to the right
   st_numscalar(h_l_s, left_max)                                                //Storing effective neighbourhood from left
   st_numscalar(h_r_s, right_max)                                               //Storing effective neighbourhood from right
   st_numscalar(q_new_s, q)                                                     //Storing effective number of observations
   
   //Producing output table below
   printf("RD Distribution Test using permutations.\n")
   printf("Cutoff c =%10.00g  {c |}   Left of c   Right of c     Number of obs   = %10.0g\n", z_bar, n_w)
   printf("{hline 22}{c +}{hline 26}    Fixed q         = %10.0g\n", q)
   printf("        Number of obs {c |}  %10.0g   %10.0g     Number of perms = %10.0g\n", n_left, n_right, P)
   if (ts_s == 0) {
     printf("   Eff. number of obs {c |}  %10.0g   %10.0g     Test statistic  =        CvM\n", q, q)
   }
   else {
     printf("   Eff. number of obs {c |}  %10.0g   %10.0g     Test statistic  =        max\n", q, q)
   }
   printf("   Eff. neighbourhood {c |}  %10.0g   %10.0g \n", left_max, right_max)
   printf("Running variable : {txt}%5s\n", running_s)
   printf("{bf:Covariates}       :{txt}%5s\n", covariates_s)
   printf("{hline 22}{c TT}{hline 26}\n")
   printf("{txt}{space 22}{c |}         {bf:P.value}\n")
   printf("{hline 22}{c +}{hline 26}\n")
   printf("               {bf:Result} {c |} {res}     %10.0g\n", p_val)
   printf("{hline 22}{c BT}{hline 26}\n")
   
   }

end
