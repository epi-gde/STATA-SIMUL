// Simulation for article on bias of COVID-19 VE related to self-testing
// December 2023
// EpiConcept
// Esther Kissling
// Part of VEBIS project
// This script is for a simulation for bias mechanism 1 (please see associated protocol)
// Please see also the plan of analysis for "Simulation study for article on self-testing and COVID-19 VE"
* Essentially we are creating a dataset of a population of ARI that have specific properties of vaccination being a COVID-19 case, self-testing, seeing a GP and a true VE
* Then we select those seeing the GP and calculate the "estimated VE"
* This script creates the simulated datasets along with the VE estimated from the simulated data. 
* We would like to have 3 VE estimates: the VE in the total population of the simulated data, which should be essentially equal to the "true VE",
* the VE among those consulting the GP, unadjusted for self-testing, and the VE among those consulting the GP adjusted for self-testing.
* This script also creates the "states" dataset (very large) that collects the parameters needed to replicate exactly each analysis.


// Version History
// v2
//     - New version August 2024



***************************************************************************
* Preliminaries ***********************************************************
***************************************************************************

// We drop all existing programmes and set the seed. 
prog drop _all
set rngstream 1		// This ensures that we can reconstruct the dataset 
set seed 576819506
capture postclose simcheck1
capture postclose rngstates1


*********************************************************
* Starting a program that will be replicated many times *
* We distinguish between a program that generates data 	*
* and a program that does the analysis					*
*********************************************************

//Here we define the simulation program, called "datagen"
clear


program define datagen
 
 syntax [, ve(real 0.2)	st(real 0.1) s_rr(real 1) possee(real 0.5) negsee(real 0.5)]

 // Here are our 5 varying parameters. We cycle the simulation through each value of parameter.
	* ve		//  Parameter for true vaccine effectiveness
	* st		// 	Parameter for self-testing in the unvaccinated population (here: among unvaccinated, as we assume an association between self-testing and vaccination)
	* s_rr		// 	Parameter for the association between vaccination and self-testing (it is an RR)
	* possee 	// 	Association (RR) between a positive self-test result and seeing a GP
	* negsee 	// 	Association (RR) between a negative self-test result and seeing a GP



// Here we have as an included file the main part of our data generation
include "simulation self test mechanism datagen analysis_include"	// This randomly selects patients to allocate the various parameters to.
													// Open the included do-file and have a look.
end
* End of data generation programme


********************************************
* Here is the programe to run the analysis *
********************************************

program define analysis_data


	syntax, [rep(int 0) post(string) ve1(string) st1(string)  s_rr1(string) possee1(string) negsee1(string)]

// We collect some N to check that nothing is being too strange (e.g. negative people)
// The following three sections are in the general population (before restricting to those seeing the GP)
* Overall pop
 tab case vacc,  matcell(A)
 scalar b4_a = A[2,1]
 scalar b4_b = A[2,2]
 scalar b4_c = A[1,1]
 scalar b4_d = A[1,2]
 
*Overall pop selftesting
 tab case vacc if selftest==1,  matcell(A)
 scalar b4st_a = A[2,1]
 scalar b4st_b = A[2,2]
 scalar b4st_c = A[1,1]
 scalar b4st_d = A[1,2]

* Overall pop not selftesting
  tab case vacc if selftest==0,  matcell(A)
 scalar b4notst_a = A[2,1]
 scalar b4notst_b = A[2,2]
 scalar b4notst_c = A[1,1]
 scalar b4notst_d = A[1,2]
 
 * Here we are calculating the estimated VE in the total simulated population, which should be more or less equal to the "true VE"
capture noisily logistic case vacc 
 // Here we prepare to output the information - the N replications for each combination of ve, st, s_rr, possee, negsee, the logistic values, the checking values
 // We also collect information if the regression does not run: _rc>0
 	if _rc==0 & !mi("`post'") post `post'  (`rep') ("Total") /// 
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")  /// 
		(_b[vacc]) (_se[vacc]) (e(N)) (b4_a) (b4_b) (b4_c) (b4_d)  (b4st_a) (b4st_b) (b4st_c) (b4st_d)  (b4notst_a) (b4notst_b) (b4notst_c) (b4notst_d) 
	if _rc>0 & !mi("`post'") post `post' ((`rep') ("Total") ///
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   /// 
		(.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) (.) (.) (.)
 
// Then we do the analysis											
* We need to keep those seeing the GP, this emulates the TND study
keep if seeGP==1

* Here we are calculating the estimated VE, the VE in the population seeing the GP, which could be influenced by self-testing
capture noisily logistic case vacc 

* And we collect some N among those seeing the GP to check that nothing is being too strange (e.g. negative people)
 tab case vacc,  matcell(A)
 scalar a11 = A[2,1]
 scalar b11 = A[2,2]
 scalar c11 = A[1,1]
 scalar d11 = A[1,2]
 
*Overall pop selftesting
 tab case vacc if selftest==1,  matcell(A)
 scalar st_a = A[2,1]
 scalar st_b = A[2,2]
 scalar st_c = A[1,1]
 scalar st_d = A[1,2]
 di A[1,1]

* Overall pop not selftesting
 tab case vacc if selftest==0,  matcell(A)
 scalar notst_a = A[2,1]
 scalar notst_b = A[2,2]
 scalar notst_c = A[1,1]
 scalar notst_d = A[1,2]

 // Here we output the information after restricting to those consulting the GP - the N replications for each combination of ve, st, s_rr, possee, negsee, the logistic values, the checking values
 // We also collect information if the regression does not run: _rc>0
 	if _rc==0 & !mi("`post'") post `post'  (`rep') ("Noadj") /// 
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")  /// 
		(_b[vacc]) (_se[vacc]) (e(N)) (a11) (b11) (c11) (d11) (st_a) (st_b) (st_c) (st_d)  (notst_a) (notst_b) (notst_c) (notst_d)
	if _rc>0 & !mi("`post'") post `post' ((`rep') ("Noadj") ///
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   /// 
		(.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) (.) (.) (.)

// We then carry out a VE analysis adjusted for selftest use: 		
 capture noisily logistic case vacc selftest



 // Here we prepare to output the information - the N replications for each combination of ve, st, s_rr, possee, negsee, the logistic values, the checking values
 // We also collect information if the regression does not run: _rc>0
  	if _rc==0 & !mi("`post'") post `post'  (`rep') ("Adj")  ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   (_b[vacc]) (_se[vacc]) (e(N)) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) 
	if _rc>0 & !mi("`post'") post `post' (`rep') ("Adj")  ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) 



end
* End of analysis programme


**************************************
* Performing the simulation **********
**************************************

timer clear 1

// Perform simulation (1000 repetitions) - TO BE CHANGED TO 10000
local reps 10
local repsplus1 = `reps'+1

* Here is our output file of the simulation - we are giving the variable names here for what we collect.
* As a reminder: "Method" refers to the VE in the VE in the population that consults the GP, unadjusted and adjusted for self-testing
postfile simcheck1  int(rep) str8(method) str8(ve1 st1 s_rr1 possee1 negsee1) float(b se) int(N) float(case1vacc0 case1vacc1 case0vacc0 case0vacc1 st_case1vacc0 st_case1vacc1 st_case0vacc0 st_case0vacc1 nost_case1vacc0 nost_case1vacc1 nost_case0vacc0 nost_case0vacc1)   ///
	using simcheck1_postfile, replace
	
* With this output file we collect information that we need to reconstruct  the dataset for an ith repetition		
postfile rngstates1 str8(ve1 st1 s_rr1 possee1 negsee1) int(rep) str2000(rngstate1 rngstate2 rngstate3) ///
	using rngstates1_postfile, replace
	
	local timer1 = 0
	
// For each value we want for each parameter (we have a very large number of scenarios!)
	foreach j of numlist 0.2 0.4 0.6 {
			foreach k of numlist 0.1 0.2 0.3 {
				foreach l of numlist 1 1.5 2 2.5 {
				foreach m of numlist 0.5 0.7 1 1.5 2 {
					foreach n of numlist  0.5 0.7 1 1.5 2 {

	    timer on 1

			* For each of the scenarios above we would like 1000 repetitions (to be changed to 10000 repetitions)								
						forvalues i=1/`repsplus1' {
						    
		if `i'==1 _dots 0 , title("Simulation running (`reps' repetitions)  `timer1' ")
		_dots `i' 0
		// Here we collect info that can help us reconstruct the dataset (we need to partition to 3, as they are quite large)
		local rngstate1 = substr(c(rngstate),1,2000)
		local rngstate2 = substr(c(rngstate),2001,2000)
		local rngstate3 = substr(c(rngstate),4001,.)
		
		* Here is our output file to capture the states, so we can replicate the simulation
		post rngstates1   ("`j'") ("`k'") ("`l'") ("`m'") ("`n'")  (`i')  ("`rngstate1'") ("`rngstate2'") ("`rngstate3'")
		
		if `i'>`reps' continue, break
		
		// This is where we run the data generation and analysis programmes, to obtain the results:
		 quietly datagen, ve(`j') st(`k') s_rr(`l') possee(`m') negsee(`n')
		 quietly analysis_data, rep(`i') post(simcheck1) ve1("`j'") st1("`k'")  s_rr1("`l'") possee1("`m'") negsee1("`n'")
				

	                      }  // repetitions 

		timer off 1				  
		quietly timer list 
		local timer1 = r(t1)  
	    timer clear 1

					}   // n negsee
				}    //  m possee
			}    //  l  RR
		}   //  k  st
	}    //   j  ve 


// And we close the files that capture the outputs:
postclose simcheck1
postclose rngstates1

