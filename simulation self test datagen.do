// Simulation for article on bias of COVID-19 VE related to self-testing
// December 2023
// EpiConcept
// Esther Kissling
// Gilles Desve Review in September 2024
// This script creates a simulation to understand bias introduced by self-testing in primary care test-negative design COVID-19 vaccine effectiveness studies
// Part of VEBIS project
* We are creating a dataset of a population of ARI that have specific properties of COVID-19 vaccination, being a COVID-19 case, self-testing, seeing a GP and a true VE
* Then we select those seeing the GP and calculate the "estimated VE"
* This script creates the simulated datasets along with the log odds ratio estimated from the simulated data. 
* We obtain 3 VE estimates: the VE in the total population of the simulated data (just as a validation exercise), which should be essentially equal to the "true VE",
* the VE among those consulting the GP, unadjusted for self-testing, and the VE among those consulting the GP adjusted for self-testing.
* This script also creates the "states" dataset that collects the parameters needed to replicate exactly each analysis.


// Version History
// v3.3
//   - New version 2 August 2024
//   - Restructured 12 September 2024 by Gilles Desvé  => v3.2 
//   - use of binomial for binary variables instead of uniform
//   - Final version updated 2024-11-4


***************************************************************************
* Preliminaries ***********************************************************
***************************************************************************

// We drop all existing programmes and set the seed. 
prog drop _all
set seed 576819506
set rngstream 1 


capture postclose simcheck1
capture postclose rngstates1

clear


***************************************************************************
* Setting simulation parameters  ******************************************
***************************************************************************
// Global flag to identify the result file for further analysis 
// if empty name fixed to simcheck1_postfile
global FILEFLAG 
												
// Perform simulation (5000 repetitions) 
// Number of repetitions for that stream
local reps 5000




***************************************************************************
* Setting fixed global parameters ******************************************
***************************************************************************

* N total (this is fixed)
global TOT = 60000
* N controls (this is fixed) 
global CONTR = 50000
* N case 
global CASE = $TOT - $CONTR


* Vaccine coverage  (this is fixed)
global VC = 0.45
* Sensitivity of self-test results - 60% 
global SENS = 0.6
* Specificity of self-test results 99% (fixed)
global SPEC_SELFTEST = 0.99     
global INVSPEC_SELFTEST = (1-0.99)
* Proportion of ARI who see GP (this is fixed)
global SEEGP = 0.1

// List of values for each parameter to be tested in loops
local VE_LIST   0.2 0.4 0.6   
local ST_LIST  0.1 0.2 0.3
local RR_LIST   1 1.5 2 2.5
local POSSEE_LIST  0.5 0.7 1 1.5 2
local NEGSEE_LIST  0.5 0.7 1 1.5 2


*********************************************************
* Include datagen and analysis routines                 *
* Datagen generates data for a set of parameters        *
* and analysis_data does the analysis					*
*********************************************************

// Here we have as an included file the main part of our data generation
// This randomly selects patients to allocate the various parameters to.
// And we also define the data analysis program which will be used in the simulation
include "simulation self test datagen_include"	
													// Open the included do-file and have a look.


***************************************************************************
* Defines the core routine to be executed repetively with ONE set of parameters 
* This routine is call repetitively by the main loop for EACH set of parameters
***************************************************************************
// The dosubset programm is designed to allow redo of a subset with fixed parameters 
// Before calling this part independently, you can retrieve the rngstate 
// which was used originaly saved into the rngstate file 
													
program define dosubset
	 syntax [, reps(int 1) ve(real 0.2)  st(real 0.1) s_rr(real 1) possee(real 0.5) negsee(real 0.5) timer(real 0) loop(int 0) ]

        local repsplus1 = `reps'+1
		* For each scenarios tested we would like xxxx repetitions (xxxx defined in local "reps", see above)
		* All parameters values for the scenario are transmited to that procedure 
		* in addition to the fixed parameters defined at the beginning of this script
		forvalues i=1/`repsplus1' {
						    
			if `i'==1 _dots 0 , title("Simulation running (`reps' repetitions) `ve' `st' `s_rr' `possee' `negsee' (`timer')/`loop' ")
			_dots `i' 0
			
			
			if `i'>`reps' continue, break
			
			// This is where we run the data generation and analysis programmes, to obtain the results:
			 quietly datagen , ve(`ve') st(`st') s_rr(`s_rr') possee(`possee') negsee(`negsee')
			 quietly analysis_data, rep(`i') post(simcheck1) ve1("`ve'") st1("`st'") ///
			            s_rr1("`s_rr'") possee1("`possee'") negsee1("`negsee'") loop(`loop') 
		      
	     }
		 
end


**************************************
* Performing the simulation 
* Here we start the job! 
**************************************

timer clear 1   // Timer of each loop (allows to estimate global time)
timer clear 2   // Store global time of the process
timer on 2


* Here is our output file of the simulation - we are giving the variable names here for what we collect.
* As a reminder: "Method" refers to the VE in the VE in the population that consults the GP, unadjusted and adjusted for self-testing
postfile simcheck1  int(rep) int(loop) str8(method) str8(ve1 st1 s_rr1 possee1 negsee1) float(b se) int(N) float(case1vacc0 case1vacc1 case0vacc0 case0vacc1 st_case1vacc0 st_case1vacc1 st_case0vacc0 st_case0vacc1 nost_case1vacc0 nost_case1vacc1 nost_case0vacc0 nost_case0vacc1)   ///
	using simcheck1_postfile$FILEFLAG, replace

	
* With this output file we collect information that we need to reconstruct the dataset for an ith repetition		
postfile rngstates1 str8(ve1 st1 s_rr1 possee1 negsee1) int(rep) str2000(rngstate1 rngstate2 rngstate3 )  ///
	using rngstates1_postfile$FILEFLAG, replace

	
**************************************
* Executing  the main loop  **********
**************************************

local timer1 = 0   // used to diplay loop time
local counter = 1  // used to diplay loop number

		
// For each value we want for each parameter (we have a very large number of scenarios!)
	foreach j of numlist `VE_LIST' {             // VE
			foreach k of numlist `ST_LIST' {         // st
				foreach l of numlist `RR_LIST' {        // RR  
				foreach m of numlist `POSSEE_LIST' {        // possee  
					foreach n of numlist  `NEGSEE_LIST' {        //  negsee

					
		timer on 1	

		// Here we collect info that can help us reconstruct the dataset (we need to partition to 3, as they are quite large)
		local rngstate1 = substr(c(rngstate),1,2000)
		local rngstate2 = substr(c(rngstate),2001,2000)
		local rngstate3 = substr(c(rngstate),4001,.)
		* Here is our output file to capture the states, so we can replicate the simulation for starting from first repetition  
		post rngstates1   ("`j'") ("`k'") ("`l'") ("`m'") ("`n'")  (1) ("`rngstate1'") ("`rngstate2'") ("`rngstate3'")				
		
		
		dosubset, reps(`reps')  ve(`j') st(`k') s_rr(`l') possee(`m') negsee(`n') timer(`timer1') loop(`counter')	
				
		
		// The counter gives feedback on where we are in the iterative process
	    local ++counter
		timer off 1 
		// The timer gives an indiction of time needed for each scenario (first display is 0)
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


timer off 2
timer list 
