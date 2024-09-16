

prog drop _all
timer clear 

program define testit

local pop = 60000
local control = `pop'/10
drop _all

set obs `pop'


timer on 1

gen case = 0 in 1/`control'
replace case = 1 if case==.

gen random=runiform() if case==1 // Creates a random variable only for cases
sort random	
count if case==1
local coverage = int(r(N)*0.4)
di `coverage'
gen vacc=1 in 1/`coverage' 
recode vacc .=0 if case==1 	// The rest are 0
tab case vacc
di "recode"

timer off 1


drop _all
set obs `pop'

timer on 2

gen case = 1
replace case = 0 in 1/`control'

gen random=runiform() if case==1 // Creates a random variable only for cases
sort random	
count if case==1
local coverage = int(r(N)*0.4)
di `coverage'
gen vacc = 0 if case ==1
replace vacc=1 in 1/`coverage' 
tab case vacc
di "replace"
timer off 2


drop _all
set obs `pop'

timer on 3

gen case = 0 in 1/`control'
replace case = 1 if case==.

gen vacc=rbinomial(1,0.4) if case==1 // Creates a random variable only for cases
di "rbinomial"
tab case vacc

timer off 3


drop _all
set obs `pop'

timer on 4

gen case = 1 
replace case = 0 in 1/`control'

gen random=runiformint(1,`pop') if case==1 // Creates a random variable only for cases
sort random 
count if case==1
local coverage = int(r(N)*0.4)
di `coverage'
gen vacc = 0 if case ==1
replace vacc=1 in 1/`coverage' 
di "runiformint"
tab case vacc

timer off 4


drop _all
set obs `pop'

timer on 5

gen case = 1
replace case = 0 in 1/`control'

gen vacc = runiform() <= 0.4000 if case ==1
di "runiform<"
tab case vacc

timer off 5



end


testit
timer list 



