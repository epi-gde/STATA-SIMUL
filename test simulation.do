

prog drop _all
timer clear 

program define testit


drop _all
set obs 100000

timer on 1

gen case = 0 in 1/90000
replace case = 1 if case==.

gen random=runiform() if case==1 // Creates a random variable only for cases
sort random	
count if case==1
local coverage = int(r(N)*0.4)
di `coverage'
gen vacc=1 in 1/`coverage' 
recode vacc .=0 if case==1 	// The rest are 0
tab case vacc

timer off 1


drop _all
set obs 100000

timer on 2

gen case = 1
replace case = 0 in 1/90000

gen random=runiform() if case==1 // Creates a random variable only for cases
sort random	
count if case==1
local coverage = int(r(N)*0.4)
di `coverage'
gen vacc = 0 if case ==1
replace vacc=1 in 1/`coverage' 
tab case vacc

timer off 2


drop _all
set obs 100000

timer on 3

gen case = 0 in 1/90000
replace case = 1 if case==.

gen vacc=rbinomial(1,0.4) if case==1 // Creates a random variable only for cases
tab case vacc

timer off 3



end


testit
timer list 
