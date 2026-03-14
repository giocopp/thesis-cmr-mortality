***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
* Replication material for:
* Political and environmental risks influence migration and human smuggling across the Mediterranean Sea
* Camarena, Claudy, Wang, Wright
*
* Note: loading the data using Stata will show variable labels
* Corresponding input file: event_study.dta
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig 3. Impact of Italian Intervention
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

// PANEL A
	
	tw  (tsline arrival_total if year==2017) (tsline arrival_total_pred if year==2017, lpattern(dash) lcolor(gs7)), ///
	xline(28, lcolor(blue)) xtitle("Week") ytitle(Weekly Arrivals) ///
	legend(size(vsmall) pos(1) ring(0) col(1) order(1 "Actual" 2 "Predicted")) ///
	scheme(lean1) 	text(2500 28 "Intervention {&rarr}", placement(left)) text(1850 0 "Pre-treat", placement(right)) text(1750 0 "R{superscript:2} = .35", placement(right)) text(1850 52 "Post-treat", placement(left)) text(1750 52 "R{superscript:2} = .02", placement(left))
	
	gr export "~/MiM_Replication/OUTPUT/TS_LibyaIntervention_SeasonalPred_Levels.pdf", replace

// PANEL B	
	
	matrix intresults = J(4,4,.)
	
	quietly ivreghdfe diff post if year==2017, absorb(absorb) cluster(week) bw(4)

	matrix intresults[1,1] = 1
	matrix intresults[1,2] = _b[post]
	matrix intresults[1,3] = _b[post]-1.96*_se[post]
	matrix intresults[1,4] = _b[post]+1.96*_se[post]
	
	quietly ivreghdfe diff post L.diff if year==2017, absorb(absorb) cluster(week) bw(4)
	
	matrix intresults[2,1] = 2
	matrix intresults[2,2] = _b[post]
	matrix intresults[2,3] = _b[post]-1.96*_se[post]
	matrix intresults[2,4] = _b[post]+1.96*_se[post]
	
	quietly ivreghdfe diff post_early if year==2017, absorb(absorb) cluster(week) bw(4)
	
	matrix intresults[3,1] = 3
	matrix intresults[3,2] = _b[post_early]
	matrix intresults[3,3] = _b[post_early]-1.96*_se[post_early]
	matrix intresults[3,4] = _b[post_early]+1.96*_se[post_early]
	
	quietly ivreghdfe diff post_early  L.diff if year==2017, absorb(absorb) cluster(week) bw(4)

	matrix intresults[4,1] = 4
	matrix intresults[4,2] = _b[post_early]
	matrix intresults[4,3] = _b[post_early]-1.96*_se[post_early]
	matrix intresults[4,4] = _b[post_early]+1.96*_se[post_early]
	
	svmat intresults, names(ir_)

gen y1var_base = ir_3
gen y2var_base = ir_4

gen x1var_base = ir_1 // - .1
gen x2var_base = ir_1 // - .1

gen coefest_base = ir_2

	tw (pccapsym y1var_base x1var_base y2var_base x2var_base, msymbol(none) lcolor(gray) mcolor(gray)) (scatter coefest_base x1var_base, msymbol(o) mcolor(gray)), legend(off) xlabel(1 "Treatment" 2 "+ Lag" 3 "Early Treatment" 4 "+ Lag", valuelabels) ///
	xscale(range(.75 4.25)) yline(0, lcolor(black)) yscale(range(-800 50)) ///
	text( -175 1.1 ///
        "`=ustrunescape("\u23AB")'" /* RCB UPPER HOOK   */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AC")'" /* RCB MIDDLE PIECE */ ///    
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AD")'" /* RCB LOWER HOOK   */ , size(5)  color("`RCBcol'") orient(horizontal)) ///
	text(-175 1.15 ".6 SD", placement(right))	

	gr export "~/MiM_Replication/OUTPUT/Libya_CoefPlots.pdf", replace

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig SI-1. Impact of Italian Intervention, accounting for post-intervention Battle of Sabratha
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
	
drop ir_1 ir_2 ir_3 ir_4 y1var_base y2var_base x1var_base x2var_base coefest_base
	
	matrix intresults = J(4,4,.)
	
	quietly ivreghdfe diff post if year==2017&week<=36, absorb(absorb) cluster(week) bw(4)

	matrix intresults[1,1] = 1
	matrix intresults[1,2] = _b[post]
	matrix intresults[1,3] = _b[post]-1.96*_se[post]
	matrix intresults[1,4] = _b[post]+1.96*_se[post]
	
	quietly ivreghdfe diff post L.diff if year==2017&week<=36, absorb(absorb) cluster(week) bw(4)
	
	matrix intresults[2,1] = 2
	matrix intresults[2,2] = _b[post]
	matrix intresults[2,3] = _b[post]-1.96*_se[post]
	matrix intresults[2,4] = _b[post]+1.96*_se[post]
	
	quietly ivreghdfe diff post_early if year==2017&week<=36, absorb(absorb) cluster(week) bw(4)
	
	matrix intresults[3,1] = 3
	matrix intresults[3,2] = _b[post_early]
	matrix intresults[3,3] = _b[post_early]-1.96*_se[post_early]
	matrix intresults[3,4] = _b[post_early]+1.96*_se[post_early]
	
	quietly ivreghdfe diff post_early  L.diff if year==2017&week<=36, absorb(absorb) cluster(week) bw(4)

	matrix intresults[4,1] = 4
	matrix intresults[4,2] = _b[post_early]
	matrix intresults[4,3] = _b[post_early]-1.96*_se[post_early]
	matrix intresults[4,4] = _b[post_early]+1.96*_se[post_early]
	
	svmat intresults, names(ir_)

gen y1var_base = ir_3
gen y2var_base = ir_4

gen x1var_base = ir_1 // - .1
gen x2var_base = ir_1 // - .1

gen coefest_base = ir_2

	tw (pccapsym y1var_base x1var_base y2var_base x2var_base, msymbol(none) lcolor(gray) mcolor(gray)) (scatter coefest_base x1var_base, msymbol(o) mcolor(gray)), legend(off) xlabel(1 "Treatment" 2 "+ Lag" 3 "Early Treatment" 4 "+ Lag", valuelabels) ///
	xscale(range(.75 4.25)) yline(0, lcolor(black)) yscale(range(-800 50)) ///
	text( -210 1.1 ///
        "`=ustrunescape("\u23AB")'" /* RCB UPPER HOOK   */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AC")'" /* RCB MIDDLE PIECE */ ///    
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AA")'" /* RCB EXTENSION    */ ///
        "`=ustrunescape("\u23AD")'" /* RCB LOWER HOOK   */ , size(5)  color("`RCBcol'") orient(horizontal)) ///
	text(-210 1.15 ".7 SD", placement(right))	

	gr export "~/MiM_Replication/OUTPUT/Libya_CoefPlots_Sabratha.pdf", replace
