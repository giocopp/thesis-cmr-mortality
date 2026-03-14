***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
* Replication material for:
* Political and environmental risks influence migration and human smuggling across the Mediterranean Sea
* Camarena, Claudy, Wang, Wright
*
* Note: loading the data using Stata will show variable labels
* Corresponding input file: monthly_trends.dta
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig 3. Impact of Italian Intervention
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
	
	tsset edate	

// PANEL C
	
	tw (tsline arrival_total if year==2017, lc(blue)) || (tsline arrival_total_g if year==2017, lc(green)) || (tsline arrival_total_s if year==2017, lc(black)), ///
	legend(size(vsmall) pos(1) ring(0) col(2) order(1 "Italy" 2 "Greece" 3 "Spain")) ///
	scheme(lean1) yscale(range(0,25000)) ylabel(0(5000)25000) /// xline(21015) /// 
	xlab(20825(30)21160 ,nogrid angle(45) format(%tdm-cy)) ylab( ,nogrid) ///
	ytitle("Total Arrivals, by Month") xtitle("Date") ysize(10) xsize(12)
		
	gr export "~/Results/TS_LibyaIntervention_Displacement.pdf", replace
	
// PANEL D
	
	tw (tsline arrival_total if year==2017, lc(blue) yaxis(1)) || (tsline total_deathmiss if year==2017, lc(red) yaxis(2)), ///
	legend(size(vsmall) pos(1) ring(0) col(2) order(1 "Arrivals" 2 "Deaths")) ///
	scheme(lean1) yscale(range(0,25000)) ylabel(0(5000)25000) ///
	xlab(20825(30)21160 ,nogrid angle(45) format(%tdm-cy)) ylab( ,nogrid) ///
	ytitle("Total Arrivals, by Month", axis(1)) ytitle("Death/Missing Totals, by Month", axis(2)) xtitle("Date") 
		
	gr export "~/Results/TS_LibyaIntervention_Deaths.pdf", replace
	
