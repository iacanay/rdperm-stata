*! version 1.0.0 08October2015
program define rdperm, rclass
version 13

syntax varlist [if] [in]                                  ///
               [, Cutoff(real 0)                          ///
				  Qobservations(integer 50)               ///
				  Permutations(integer 999)] 
 
//DESCRIPTION
//VARLIST				  
//First variable is the running variable. 
//Second variable onwards are covariates.
//OPTIONS
//Cutoff value has default 0 and must be real.
//Qobservations is the effective number observations with default of 50.
//Permutations is the number of permutations performed with default of 999.
//OUPUT
//P.Val is the p value

marksample touse 
gettoken running covariates :varlist                                       //parsing the running variable and covariates
tempname pval N N_l N_r h_l h_r                                            //introducing tempname to save output of mata function

mata: rdperm_work("`covariates'", "`running'", "`touse'",                  ///
                   `cutoff', `qobservations', `permutations',              ///
				   "`pval'", "`N'", "`N_l'", "`N_r'", "`h_l'", "`h_r'")    //running the mata function using the inputs from data

return scalar pval = `pval'                   //storing pvalue
return scalar h_r = `h_r'                     //storing effective neighbourhood to the right
return scalar h_l = `h_l'                     //storing effective neighbourhood to the left
return scalar N_r = `N_r'                     //storing no. of obs. to the right
return scalar N_l = `N_l'                     //storing no. of obs to the left
return scalar N = `N'                         //storing number of observations
return scalar P = `permutations'              //storing number of permutations
return scalar q = `qobservations'             //storing effective observations used
return scalar c = `cutoff'                    //storing threshold value

return local cmd "rdperm"                     //storing command name
return local teststat "CvM"                   //storing test statistic used
return local covariates `covariates'          //storing covariate names
return local runninvar `running'              //storing running variable name
				   
end


mata:

void rdperm_work(string scalar covariates_s, string scalar running_s,               ///
                 string scalar touse_s, real scalar cutoff_s,                       ///
				 real scalar q_s, real scalar P_s, string scalar pval_s,            ///
				 string scalar N_s, string scalar N_l_s, string scalar N_r_s,       ///
				 string scalar h_l_s, string scalar h_r_s) {


   rseed(12345)                                            //Setting seed for the permutations
				  
   w                 = st_data(., tokens(covariates_s),    ///
                                              touse_s)     //This is the vector of covariates
   z                 = st_data(., running_s, touse_s)      //This is vector if running variables
   c                 = cutoff_s                            //This is threshold value
   q                 = ceil(q_s)                           //Number of observations used from each side of threshold
   z                 = z :- c                              //Normalized running variable
   n_w               = rows(w)                             //Sample size of w
   n_z               = rows(z)                             //Sample size of z
   P                 = P_s                                 //Number of permutations

   if (n_w != n_z) {                                       //Error message for size of w and z not equal
     printf("size of variables do not match")
	 exit(error(3498))
   }

   z_left            = select(z, z[.,1] :< 0)             //Z obs to left of threshold
   w_left            = select(w, z[., 1] :< 0)            //W Obs to left of threshold
   n_left            = rows(z_left)                       //No. of obs to left of threshold
   z_right           = select(z, z[.,1] :>= 0)            //Z obs to right of threshold
   w_right           = select(w, z[., 1] :>= 0)           //W Obs to right of threshold
   n_right           = rows(z_right)                      //No. of obs to right of threshold
   
   if (n_left < q) {                                       //Error message for size of w and z not equal
     printf("insufficient sample size for q")
	 exit(error(3498))
   }
   if (n_right < q) {                                       //Error message for size of w and z not equal
     printf("insufficient sample size for q")
	 exit(error(3498))
   }
   
   w_left            = w_left[order(z_left,1), .]         //Induced order of W obs to left of threshold
   z_left            = sort(z_left,1)                     //Sorting Z obs to left of threshold
   w_right           = w_right[order(z_right,1), .]       //Induced order of W obs to right of threshold
   z_right           = sort(z_right,1)                    //Sorting Z obs to right of threshold

   left_max          = z_left[(n_left+1-q), 1]            //Effective neighbourhood from left
   right_max         = z_right[q,1]                       //Effective neighbourhood from right

   w_left            = w_left[((n_left-q+1)..n_left), .]  //q closest W from the left of threshold
   w_right           = w_right[(1..q), .]                 //q closest W from the right of threshold
   w_q               = w_left \ w_right                   //W obs used in permutations
   
   m_minus           = J(2*q,1,0)                         //
   m_plus            = J(2*q,1,0)                         //
   mperm_minus       = J(2*q,1,0)                         //
   mperm_plus        = J(2*q,1,0)                         //
   CvM_perm          = J(P,1,0)                           //

   for (i=1; i <= 2*q; i++) {
     w_3               = J(q, 1, w_q[i, .])               //
     A                 = w_left :<= w_3                   //
     m_minus[i,.]      = mean(rowmin(A))                  //
     A                 = w_right :<= w_3                  //
     m_plus[i,.]       = mean(rowmin(A))                  //
   }
   CvM_stat          = mean((m_minus - m_plus):^2)        //CvM TS in sample
   
   for (j=1; j <= P; j++) {
     w_perm            = jumble(w_q)                           //Random permutation of W obs
   for (i=1; i <= 2*q; i++) {
       w_3               = J(q, 1, w_q[i, .])                  //
       A                 = w_perm[(1..q), .] :<= w_3           //
       mperm_minus[i,.]  = mean(rowmin(A))                     //
       A                 = w_perm[((q+1)..2*q), .] :<= w_3     //
       mperm_plus[i,.]   = mean(rowmin(A))                     //
   }
     CvM_perm[j,1]      = mean((mperm_minus - mperm_plus):^2) //CvM TS in each permutation
   }
   p_val             = mean((CvM_perm :>= CvM_stat))       //Computed p value
   
   st_numscalar(pval_s, p_val)                             //Storing pvalue to be viewed in stata
   st_numscalar(N_s, n_z)                                  //Storing sample size N
   st_numscalar(N_l_s, n_left)                             //Storing sample size to the left
   st_numscalar(N_r_s, n_right)                            //Storing sample size to the right
   st_numscalar(h_l_s, left_max)                           //Storing effective neighbourhood from left
   st_numscalar(h_r_s, right_max)                          //Storing effective neighbourhood from right
   
   //Producing output table below
   printf("RD Distribution Test using permutations.\n")
   printf("Cutoff c =%10.00g  {c |}   Left of c   Right of c     Number of obs   = %10.0g\n", c, n_w)
   printf("{hline 22}{c +}{hline 26}    Fixed q         = %10.0g\n", q)
   printf("        Number of obs {c |}  %10.0g   %10.0g     Number of perms = %10.0g\n", n_left, n_right, P)
   printf("   Eff. number of obs {c |}  %10.0g   %10.0g     Test statistic  =        CvM\n", q, q)
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
