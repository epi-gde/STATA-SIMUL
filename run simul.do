/*
global rmethod = "datagenuniform"

foreach iloop of numlist 3000 {
	
	global ireps = `iloop'
	global FILEFLAG = "_uni"
	do "simulation self test mechanism datagen analysis.do"
	
	
	
}
*/

global rmethod = "datagenbin"

fovalues iloop = 1/10  {
	
	global ireps = 100
	global FILEFLAG = "_bin_`iloop'"
	do "simulation self test mechanism datagen analysis.do"
	
	
	
}
*/
