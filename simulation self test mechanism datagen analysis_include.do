// Simulation for article on bias of COVID-19 VE related to self-testing
// December 2023
// EpiConcept
// Esther Kissling
// Part of VEBIS project
// This script is for a simulation for bias mechanism 1 (please see associated protocol)
// Please see also the plan of analysis for "Simulation study for article on self-testing and COVID-19 VE"
// Gilles Desve Review in september 2024


// Version History
// v3.2
//     - Original working version
//     - Restructured 12 September 2024
//     - use of binomial for binary variables instead of uniform/sort mechanism
//     - final reviewed by Gilles Desvé => v3.2


// Question, how do I loop around these different values?



// this part is called by the loop to generate one subset of data  
// Parameters are passed as variable to the command
program define datagen
 
 syntax [, ve(real 0.2)	st(real 0.1) s_rr(real 1) possee(real 0.5) negsee(real 0.5)]

 // Here are our 5 varying parameters. We cycle the simulation through each value of parameter.
	* ve		//  Parameter for true vaccine effectiveness
	* st		// 	Parameter for self-testing in the unvaccinated population (here: among unvaccinated, as we assume an association between self-testing and vaccination)
	* s_rr		// 	Parameter for the association between vaccination and self-testing (it is an RR)
	* possee 	// 	Association (RR) between a positive self-test result and seeing a GP
	* negsee 	// 	Association (RR) between a negative self-test result and seeing a GP

* OR from the varying VE
local OR = (1-`ve')
	

***************************************************************************
* Clearing the Stata console and setting the number of cases and controls *
***************************************************************************

drop _all	// drops all observations and variables
set obs $TOT // Setting the total number of observations


// Cases	 - the ratio of cases to controls was 1:5
gen case = 0 in 1/$CONTR
replace case = 1 if case==.


*********************************************
* Setting the parameters for the simulation *
*********************************************

***** Vaccine coverage
// vaccine coverage - 45% among controls - this is fixed.

gen vacc = rbinomial(1,$VC) if case == 0

// vaccine coverage among cases - depends on true VE, 20, 40, 60 
* VC among cases
* There is probably an easier simplification, but the number of vacc cases is OR*(N vacc controls/N unvacc controls)* N cases/(1 + (OR*(N vacc controls/N unvacc controls)))
count if case==0 & vacc==1
local N_VACCCONTROL = r(N)
count if case==0 & vacc==0
local N_UNVACCCONTROL = r(N)

local N_VACCCASE = (`OR'*(`N_VACCCONTROL'/`N_UNVACCCONTROL')* $CASE/(1 + (`OR'*(`N_VACCCONTROL'/`N_UNVACCCONTROL'))))
local VC_CASE = `N_VACCCASE'/$CASE
di "`VC_CASE'"

replace vacc = rbinomial(1,`VC_CASE') if case == 1

***** Self-testing
// We have information on the proportion self-testing overall (among vaccinated and unvaccinated), and assume 10%, 20% or 30% in an unvaccinated population
// We assume vaccination is positively associated with self-testing, 1, 1.5, 2, 2.5
* First we create the proportion self-testing among unvaccinated:
gen selftest=  rbinomial(1,`st') if vacc == 0


// Association (odds ratio) between vaccination and self-testing
// 1, 1.5, 2, 2.5
local PROP_VACCSELFTEST = `st'*`s_rr'
// Default selftest to 0 for vaccinated 
replace selftest = rbinomial(1,`PROP_VACCSELFTEST') if vacc == 1


***** Result of the self-test
// Correlation between self-test result and PCR SARS-CoV-2 result
// Sensitivity: 60% (but we could run sensitivity analyses) 
// Specificity: 99%
// Default selftestresult to 0 in that group 
gen selftestresult = rbinomial(1,$SENS) if case==1  &  selftest==1	


// Among controls
// Default to 0 in that group 
replace selftestresult = rbinomial(1,$INVSPEC_SELFTEST) if case==0 &  selftest==1 


****** Consulting a GP
// We assume the probability of non-self-testing people seeing the GP is 10%
gen seeGP =rbinomial(1,$SEEGP) if selftest==0


// Now we change "seeing the GP" among those doing a self-test and the test is positive
// Default to 0 in that group 
replace seeGP = rbinomial(1,`possee'*$SEEGP) if selftestresult==1 &  selftest==1 

// Now we change "seeing the GP" among those doing a self-test and the test is negative
// Default to 0 in that group 
replace seeGP =rbinomial(1,`negsee'*$SEEGP) if selftestresult==0 &  selftest==1 	// The rest are 0


end
* End of data generation programme



********************************************************************************
* Here is the programe to run the analysis *
********************************************
// This part is called for each subset of data (each repetition) 
// and append to the results file 3 records.
// one record for each type of analysis :  Total, NoAdj and Adj  analysis
program define analysis_data

	syntax, [rep(int 0) post(string) ve1(string) st1(string)  s_rr1(string) possee1(string) negsee1(string) loop(int 0)]

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
 
 
 // --------------------------  Total -----------------------------
 * Here we are calculating the estimated VE in the total simulated population, which should be more or less equal to the "true VE"
capture noisily logistic case vacc 
 // Here we prepare to output the information - the N replications for each combination of ve, st, s_rr, possee, negsee, the logistic values, the checking values
 // We also collect information if the regression does not run: _rc>0
 	if _rc==0 & !mi("`post'") post `post'  (`rep') (`loop') ("Total") /// 
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")  /// 
		(_b[vacc]) (_se[vacc]) (e(N)) (b4_a) (b4_b) (b4_c) (b4_d)  (b4st_a) (b4st_b) (b4st_c) (b4st_d)  (b4notst_a) (b4notst_b) (b4notst_c) (b4notst_d) 
	if _rc>0 & !mi("`post'") post `post' ((`rep') (`loop') ("Total") ///
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   /// 
		(.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) (.) (.) (.)
 
 
 
// --------------------------  NoAdj -----------------------------
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
 	if _rc==0 & !mi("`post'") post `post'  (`rep')  (`loop') ("Noadj") /// 
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")  /// 
		(_b[vacc]) (_se[vacc]) (e(N)) (a11) (b11) (c11) (d11) (st_a) (st_b) (st_c) (st_d)  (notst_a) (notst_b) (notst_c) (notst_d)
	if _rc>0 & !mi("`post'") post `post' ((`rep') (`loop') ("Noadj") ///
	 ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   /// 
		(.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) (.) (.) (.)

		
// --------------------------  Adj -----------------------------
		
// We then carry out a VE analysis adjusted for selftest use: 		
 capture noisily logistic case vacc selftest



 // Here we prepare to output the information - the N replications for each combination of ve, st, s_rr, possee, negsee, the logistic values, the checking values
 // We also collect information if the regression does not run: _rc>0
  	if _rc==0 & !mi("`post'") post `post'  (`rep') (`loop') ("Adj")  ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   (_b[vacc]) (_se[vacc]) (e(N)) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) 
	if _rc>0 & !mi("`post'") post `post' (`rep') (`loop') ("Adj")  ("`ve1'") ("`st1'") ("`s_rr1'") ("`possee1'") ("`negsee1'")   (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)  (.) 



end
* End of analysis programme






