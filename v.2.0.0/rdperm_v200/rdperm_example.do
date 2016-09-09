clear
set more off
capture log close
cd "..."                                     //Change to working directory
log using "rdperm_example.smcl", replace 

********************************************************************************
//Canay & Kamat(2016) RDD Distribution Test using permutations on Lee(2008)
********************************************************************************

use table_two_final.dta, clear    //Loading data
capture program drop rdperm       //Installing rdperm program

//Test for `demshareprev'
rdperm difdemshare demshareprev if use==1, c(0)

//Test for `demwinprev'
rdperm difdemshare demwinprev if use==1, c(0)

//Test for `demofficeexp'
rdperm difdemshare demofficeexp if use==1, c(0)

//Test for `othofficeexp'
rdperm difdemshare othofficeexp if use==1, c(0)

//Test for `demelectexp'
rdperm difdemshare demelectexp if use==1, c(0)

//Test for `othelectexp'
rdperm difdemshare othelectexp if use==1, c(0)

//Joint test using max test statistic
rdperm difdemshare demshareprev demwinprev demofficeexp othofficeexp  ///
       demelectexp othelectexp if use==1, c(0) 

//Joint test using CvM test statistic
rdperm difdemshare demshareprev demwinprev demofficeexp othofficeexp  ///
       demelectexp othelectexp if use==1, c(0) cvm

log off 
