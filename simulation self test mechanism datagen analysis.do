// Simulation for article on bias of COVID-19 VE related to self-testing
// December 2023
// EpiConcept
// Esther Kissling
// Gilles Desve Review in september 2024
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
// v3.1
//   - New version 2 August 2024
//   - Restructured 12 September 2024  
//   - use of binomial for binary variables instead of uniform



***************************************************************************
* Preliminaries ***********************************************************
***************************************************************************

// We drop all existing programmes and set the seed. 
prog drop _all
set seed 576819506
capture postclose simcheck1
capture postclose rngstates1

clear
timer clear

***************************************************************************
* Setting simulation parameters  ******************************************
***************************************************************************
// Global flag to identify the result file for further analysis 
// if empty name fixed to simcheck1_postfile
global FILEFLAG _GDRUN2
												
// Perform simulation (1000 repetitions) - TO BE CHANGED TO 10000
// Number of repetitions for that stream
// For Individual manual run uncomment next line
local reps 1000

// Automatic scripted run, un comment this (see run simul.do)
// local reps = $ireps

// Comment this to execute run simul.do, uncomment for manual run 
global rmethod = "datagenbin" 

// Set the stream to distribute the simulation. One stream for each computer
global stream_number = 11

***************************************************************************
* Setting fixed global parameters ******************************************
***************************************************************************

// This ensures that we can reconstruct the dataset 
set rngstream $stream_number	



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

local VE_LIST   0.2 0.4 0.6
local ST_LIST   0.1 0.2 0.3
local RR_LIST   1 1.5 2 2.5
local POSSEE_LIST  0.5 0.7 1 1.5 2
local NEGSEE_LIST  0.5 0.7 1 1.5 2


*********************************************************
* include datagen and analys routines                   *
* Datagen generates data for a set of parameters        *
* and analysis_data does the analysis					*
*********************************************************

// Here we have as an included file the main part of our data generation
// This randomly selects patients to allocate the various parameters to.
// And we also define the data analysis program which will be used in the simulation
include "simulation self test mechanism datagen analysis_include"	
													// Open the included do-file and have a look.


***************************************************************************
* define the core routine to be executed repetively with one set of parameters 
***************************************************************************
// The dosubset programm is designed to allow redo of a subset with fixed parameters 
// Before calling this part independently, you can retrieve the rngstate 
// which was used originaly saved into the rngstate file 
													
program define dosubset
	 syntax [, reps(int 1) ve(real 0.2)  st(real 0.1) s_rr(real 1) possee(real 0.5) negsee(real 0.5) timer(real 0) loop(int 0) genonly ]

        local repsplus1 = `reps'+1
		* For each scenarios tested we would like xxxx repetitions (xxxx defined in local "reps", see above)
		* All parameters values for the scenario are transmited to that procedure 
		* in addition to the fixed parameters defined at the beginning of this script
		forvalues i=1/`repsplus1' {
						    
			if `i'==1 _dots 0 , title("Simulation running (`reps' repetitions) `ve' `st' `s_rr' `possee' `negsee' (`timer')/`loop' ")
			_dots `i' 0
			
			
			if `i'>`reps' continue, break
			
			// This is where we run the data generation and analysis programmes, to obtain the results:
			 quietly $rmethod , ve(`ve') st(`st') s_rr(`s_rr') possee(`possee') negsee(`negsee')
			 if "`genonly'" == "" { 
			     quietly analysis_data, rep(`i') post(simcheck1) ve1("`ve'") st1("`st'")  s_rr1("`s_rr'") possee1("`possee'") negsee1("`negsee'")                         loop(`loop') 
		     } 
	     }
		 
end


**************************************
* Performing the simulation 
* Here we start the job ! 
**************************************

timer clear 1
timer clear 2
timer on 2

local timer1 = 0
local counter = 1


* Here is our output file of the simulation - we are giving the variable names here for what we collect.
* As a reminder: "Method" refers to the VE in the VE in the population that consults the GP, unadjusted and adjusted for self-testing
postfile simcheck1  int(rep) int(loop stream) str8(method) str8(ve1 st1 s_rr1 possee1 negsee1) float(b se) int(N) float(case1vacc0 case1vacc1 case0vacc0 case0vacc1 st_case1vacc0 st_case1vacc1 st_case0vacc0 st_case0vacc1 nost_case1vacc0 nost_case1vacc1 nost_case0vacc0 nost_case0vacc1)   ///
	using simcheck1_postfile$FILEFLAG, replace

	
* With this output file we collect information that we need to reconstruct  the dataset for an ith repetition		
postfile rngstates1 str8(ve1 st1 s_rr1 possee1 negsee1) int(rep) int(streamnb) str2000(rngstate1 rngstate2 rngstate3 )  ///
	using rngstates1_postfile$FILEFLAG, replace

/*	
// Simulation testing
// Uncomment this block to run only one step of the simulation
dosubset, reps(1) genonly	
postclose simcheck1
postclose rngstates1	
// End simulation testing
*/
	
**************************************
* Executing  the main loop  **********
**************************************
	
	
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
		post rngstates1   ("`j'") ("`k'") ("`l'") ("`m'") ("`n'")  (1) ($stream_number) ("`rngstate1'") ("`rngstate2'") ("`rngstate3'")				
		
		dosubset, reps(`reps')  ve(`j') st(`k') s_rr(`l') possee(`m') negsee(`n') timer(`timer1') loop(`counter')	
				
		// Counter give feedback on where we are in the iterative process
	    local ++counter
	    
		timer off 1 
		// Timer give an indiction of time needed for each scenario (first display is 0)
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


use simcheck1_postfile$FILEFLAG, clear
keep if method == "Noadj"
bysort loop : egen meanb = mean(b) 
by loop : egen SD = sd(b)
by loop : egen iter = max(rep) 
by loop : gen SE = SD/sqrt(iter)
gen OR = exp(meanb)
gen VE = (1 - OR) *100 
keep if rep == 1

keep meanb OR VE SD  SE  stream iter loop ve1 st1 s_rr1 possee1 negsee1
gen rmethod = "$rmethod"

if fileexists("simresult.dta") append using "simresult"
save simresult, replace 


timer off 2
timer list 
