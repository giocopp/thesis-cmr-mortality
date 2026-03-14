clear all
set more off

* Set maxiter - maxvar
set maxiter 16000
set maxvar 32767


ssc install mypkg, replace
mypkg 

*Data cleaning
ssc install egenmore, replace
ssc install carryforward, replace
ssc install charlist, replace
ssc install fsum, replace
ssc install fre, replace
ssc install ftools, replace
ssc install moremata, replace

*Graphs
ssc install palettes, replace
ssc install colrspace, replace
ssc install coefplot, replace
ssc install shp2dta, replace
ssc install mif2dta, replace 
ssc install spmap, replace
ssc install heatplot, replace
ssc install hmap, replace

*Export output
ssc install outreg2, replace
ssc install estout, replace


*Estimations
ssc install ranktest, replace
ssc install xtile2, replace  
ssc install egenmore, replace 
ssc install reghdfe, replace
ssc install ppmlhdfe, replace
ssc install ivreg2, replace
ssc install ivreg2hdfe, replace
ssc install xtqreg, replace
ssc install psacalc, replace


* Uncategorized
ssc install gtools, replace 
ssc install parmest, replace 
ssc install geodist, replace 
ssc install geoinpoly, replace


