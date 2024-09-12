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


													
// Perform simulation (1000 repetitions) - TO BE CHANGED TO 10000
// Number of repetitions for that stream
local reps 10


***************************************************************************
* Setting fixed global parameters ******************************************
***************************************************************************


* N total (this is fixed)
global TOT = 60000
* N controls (this is fixed) - should it be?
global CONTR = 50000
* N case 
global CASE = $TOT - $CONTR


* Vaccine coverage  (this is fixed)
global VC = 0.45
* Sensitivity of self-test results - 60% (note that we could do sensitivity analyses varying this between 50 and 70%, for example, if we like )
global SENS = 0.6
* Specificity of self-test results 99% (fixed)  (seem to not be in use ? )
global SPEC_SELFTEST = 0.99     
global INVSPEC_SELFTEST = (1-0.99)
* Proportion of ARI who see GP (this is fixed)
global SEEGP = 0.1


// Here we have as an included file the main part of our data generation
// This randomly selects patients to allocate the various parameters to.
// And we also define the data analysis program which will be used in the simulation
include "simulation self test mechanism datagen analysis_include"	
													// Open the included do-file and have a look.

// The dosubset programm is designed to allow redo of a subset with fixed parameters 
// Before calling this part you can retrieve the rngstate which was used into the rngstate file 
													
program define dosubset
		 syntax [, reps(int 10) ve(real 0.2)  st(real 0.1) s_rr(real 1) possee(real 0.5) negsee(real 0.5) timer(real 0)]

        local repsplus1 = `reps'+1
		* For each of the scenarios above we would like 1000 repetitions (to be changed to 10000 repetitions)								
		forvalues i=1/`repsplus1' {
						    
			if `i'==1 _dots 0 , title("Simulation running (`reps' repetitions) `ve' `st' `s_rr' `poseee' `negsee' (`timer') ")
			_dots `i' 0
			
			
			if `i'>`reps' continue, break
			
			// This is where we run the data generation and analysis programmes, to obtain the results:
			 quietly datagen, ve(`ve') st(`st') s_rr(`s_rr') possee(`possee') negsee(`negsee')
			 quietly analysis_data, rep(`i') post(simcheck1) ve1("`ve'") st1("`st'")  s_rr1("`s_rr'") possee1("`possee'") negsee1("`negsee'")
		
	     }
end


**************************************
* Performing the simulation **********
**************************************

timer clear 1


local timer1 = 0

* Here is our output file of the simulation - we are giving the variable names here for what we collect.
* As a reminder: "Method" refers to the VE in the VE in the population that consults the GP, unadjusted and adjusted for self-testing
postfile simcheck1  int(rep) str8(method) str8(ve1 st1 s_rr1 possee1 negsee1) float(b se) int(N) float(case1vacc0 case1vacc1 case0vacc0 case0vacc1 st_case1vacc0 st_case1vacc1 st_case0vacc0 st_case0vacc1 nost_case1vacc0 nost_case1vacc1 nost_case0vacc0 nost_case0vacc1)   ///
	using simcheck1_postfile, replace
	
* With this output file we collect information that we need to reconstruct  the dataset for an ith repetition		
postfile rngstates1 str8(ve1 st1 s_rr1 possee1 negsee1) int(rep) str2000(rngstate1 rngstate2 rngstate3) ///
	using rngstates1_postfile, replace
	

// For each value we want for each parameter (we have a very large number of scenarios!)
	foreach j of numlist 0.2 0.4 0.6 {             // VE
			foreach k of numlist 0.1 0.2 0.3 {         // st
				foreach l of numlist 1 1.5 2 2.5 {        // RR  
				foreach m of numlist 0.5 0.7 1 1.5 2 {        // possee  
					foreach n of numlist  0.5 0.7 1 1.5 2 {        //  negsee

					
		timer on 1	
		// Here we collect info that can help us reconstruct the dataset (we need to partition to 3, as they are quite large)
		local rngstate1 = substr(c(rngstate),1,2000)
		local rngstate2 = substr(c(rngstate),2001,2000)
		local rngstate3 = substr(c(rngstate),4001,.)
		
		* Here is our output file to capture the states, so we can replicate the simulation for starting from first repetition  
		post rngstates1   ("`j'") ("`k'") ("`l'") ("`m'") ("`n'")  (1)  ("`rngstate1'") ("`rngstate2'") ("`rngstate3'")				
		
		dosubset, reps(`reps')  ve(`j') st(`k') s_rr(`l') poseee(`m') negsee(`n') timer(`timer1')		
	  
		timer off 1 
		quietly timer list 1
		local timer1 = r(t1)
		timer clear 1

					}   // n negsee
				}    //  m poseee
			}    //  l  RR
		}   //  k  st
	}    //   j  VE


// And we close the files that capture the outputs:
postclose simcheck1
postclose rngstates1

