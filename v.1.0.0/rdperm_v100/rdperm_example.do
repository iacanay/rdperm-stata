clear
set more off
capture log close
cd "..."                                     //Change to working directory
log using "rdperm_example.smcl", replace 

********************************************************************************
//Canay & Kamat(2015) RDD Distribution Test using permutations on Lee(2008)
********************************************************************************

use table_two_final.dta, clear    //Loading data
capture program drop rdperm       //Installing rdperm program

//Testing exogeneity of `demshareprev':Col3 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  demshareprev if use==1, c(0) q(`x') p(999)
}

//Testing exogeneity of `demwinprev':Col4 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  demwinprev if use==1, c(0) q(`x') p(999)
}

//Testing exogeneity of `demofficeexp':Col5 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  demofficeexp if use==1, c(0) q(`x') p(999)
}

//Testing exogeneity of `othofficeexp':Col6 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  othofficeexp if use==1, c(0) q(`x') p(999)
}

//Testing exogeneity of `demelectexp':Col7 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  demelectexp if use==1, c(0) q(`x') p(999)
}

//Testing exogeneity of `othelectexp':Col8 of Table 1 Lee(2008)
foreach x in 25 50 100 { 
rdperm  difdemshare  othelectexp if use==1, c(0) q(`x') p(999)
}

//Testing for joint exogeneity
foreach x in 25 50 100 { 
rdperm  difdemshare demshareprev demwinprev demofficeexp othofficeexp  ///
        demelectexp othelectexp if use==1, c(0) q(`x') p(999)
}

log off 
