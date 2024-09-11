// Simulation for article on bias of COVID-19 VE related to self-testing
// December 2023
// EpiConcept
// Esther Kissling
// Part of VEBIS project
// This script is for a simulation for bias mechanism 1 (please see associated protocol)
// Please see also the plan of analysis for "Simulation study for article on self-testing and COVID-19 VE"
// GD Review 


// Version History
// v1
//     - Original working version

// Question, how do I loop around these different values?


***************************************************************************
* Setting fixed local parameters ******************************************
***************************************************************************


* N total (this is fixed)
local TOT = 60000
* N controls (this is fixed) - should it be?
local CONTR = 50000
* Vaccine coverage  (this is fixed)
local VC = 0.45
* OR from the varying VE
local OR = (1-`ve')
* Sensitivity of self-test results - 60% (note that we could do sensitivity analyses varying this between 50 and 70%, for example, if we like )
local SENS = 0.6
* Specificity of self-test results 99% (fixed)
local SPEC_SELFTEST = 0.99
local INVSPEC_SELFTEST = (1-0.99)
* Proportion of ARI who see GP (this is fixed)
local SEEGP = 0.1



***************************************************************************
* Clearing the Stata console and setting the number of cases and controls *
***************************************************************************

drop _all	// drops all observations and variables
set obs `TOT'	// Setting the total number of observations


// Cases	 - the ratio of cases to controls was 1:5
gen case = 0 in 1/`CONTR'
replace case = 1 if case==.


*********************************************
* Setting the parameters for the simulation *
*********************************************

***** Vaccine coverage
// vaccine coverage - 45% among controls - this is fixed.
gen random=runiform() if case==0 // Creates a random variable only for controls
sort random	
count if case==0
local coverage = int(r(N)*`VC')	// Here we obtain the number of controls that should be vaccinated (45% of controls, 2250)
di `coverage'
gen vacc=1 in 1/`coverage' 	// We create the vaccination variable and add 1 in the randomly sorted controls from 1 to 2250 (45% of controls)
recode vacc .=0 if case==0 	// The rest are 0, unvaccinated.

// vaccine coverage among cases - depends on true VE, 20, 40, 60 
* VC among cases
* There is probably an easier simplification, but the number of vacc cases is OR*(N vacc controls/N unvacc controls)* N cases/(1 + (OR*(N vacc controls/N unvacc controls)))
count if case==0 & vacc==1
local N_VACCCONTROL = r(N)
count if case==0 & vacc==0
local N_UNVACCCONTROL = r(N)
count if case==1
local N_CASE = r(N)
local N_VACCCASE = (`OR'*(`N_VACCCONTROL'/`N_UNVACCCONTROL')* `N_CASE'/(1 + (`OR'*(`N_VACCCONTROL'/`N_UNVACCCONTROL'))))
local VC_CASE = `N_VACCCASE'/`N_CASE'
di "`VC_CASE'"
drop random
gen random=runiform() if case==1 // Creates a random variable only for cases
sort random	
count if case==1
local coverage = int(r(N)*`VC_CASE')
di `coverage'
replace vacc=1 in 1/`coverage' if case==1
recode vacc .=0 if case==1 	// The rest are 0

***** Self-testing
// We have information on the proportion self-testing overall (among vaccinated and unvaccinated), and assume 10%, 20% or 30% in an unvaccinated population
// We assume vaccination is positively associated with self-testing, 1, 1.5, 2, 2.5
* First we create the proportion self-testing among unvaccinated:
drop random
gen random=runiform() if vacc==0 // Creates a random variable only for unvaccinated
sort random	
count if vacc==0
local prop = int(r(N)*`st')
di `prop'
gen selftest=1 in 1/`prop' if vacc==0 
recode selftest .=0 if vacc==0  	// The rest are 0


// Association (odds ratio) between vaccination and self-testing
// 1, 1.5, 2, 2.5
drop random
local PROP_VACCSELFTEST = `st'*`s_rr'
gen random=runiform() if vacc==1 // Creates a random variable only for vaccinated
sort random	
count if vacc==1
local prop = int(r(N)*`PROP_VACCSELFTEST')
di `prop'
replace selftest=1 in 1/`prop' if vacc==1 
recode selftest .=0 if vacc==1  	// The rest are 0


***** Result of the self-test
// Correlation between self-test result and PCR SARS-CoV-2 result
// Sensitivity: 60% (but we could run sensitivity analyses) 
// Specificity: 99%
drop random
gen random=runiform() if case==1 & selftest==1 
sort random	
count if case==1 & selftest==1
local prop = int(r(N)*`SENS')
di `prop'
gen selftestresult=1 in 1/`prop' if case==1 &  selftest==1
recode selftestresult .=0 if case==1  &  selftest==1	// The rest are 0

// Among controls
drop random
gen random=runiform() if case==0 &  selftest==1 
sort random	
count if case==0 &  selftest==1
local prop = int(r(N)*`INVSPEC_SELFTEST')
di `prop'
replace selftestresult=1 in 1/`prop' if case==0   &  selftest==1
recode selftestresult .=0 if case==0 &  selftest==1 	// The rest are 0


****** Consulting a GP
// We assume the probability of non-self-testing people seeing the GP is 10%

drop random
gen random=runiform() if selftest==0
sort random	
count if selftest==0
local prop = int(r(N)*`SEEGP')
di `prop'
gen seeGP=1 in 1/`prop' if selftest==0
recode seeGP .=0  if selftest==0


// Now we change "seeing the GP" among those doing a self-test and the test is positive
drop random
gen random=runiform() if selftestresult==1 &  selftest==1 
sort random	
count if selftestresult==1 &  selftest==1 
local prop = int(r(N)*`possee'*`SEEGP')	// Among those normally seeing a GP there is a different association
di `prop'
replace seeGP=1 in 1/`prop' if selftestresult==1 &  selftest==1 
recode seeGP .=0 if selftestresult==1 &  selftest==1 	// The rest are 0

// Now we change "seeing the GP" among those doing a self-test and the test is negative
drop random
gen random=runiform() if selftestresult==0 &  selftest==1 
sort random	
count if selftestresult==0 &  selftest==1
local prop = int(r(N)*`negsee'*`SEEGP')	// Among those normally seeing a GP there is a different association
di `prop'
replace seeGP=1 in 1/`prop' if selftestresult==0 &  selftest==1 
recode seeGP .=0 if selftestresult==0 &  selftest==1 	// The rest are 0






