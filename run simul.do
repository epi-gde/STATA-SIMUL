global rmethod = "datagenuniform"

foreach iloop of numlist 10 100 500 1000 2000 3000 4000 5000 {
	
	global ireps = `iloop'
	do "simulation self test mechanism datagen analysis.do"
	
	
	
}


global rmethod = "datagenbin"

foreach iloop of numlist 10 100 500 1000 2000 3000 4000 5000 {
	
	global ireps = `iloop'
	do "simulation self test mechanism datagen analysis.do"
	
	
	
}
