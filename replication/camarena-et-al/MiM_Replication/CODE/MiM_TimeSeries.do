***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
* Replication material for:
* Political and environmental risks influence migration and human smuggling across the Mediterranean Sea
* Camarena, Claudy, Wang, Wright
*
* Note: loading the data using Stata will show variable labels
* Corresponding input file: time_series.dta
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//


***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table 1. Riots, sea conditions, and migrant flows to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly reg ln_arrival_total ln_riots_prevweek if sample==1, vce(robust)
	eststo: quietly reg ln_arrival_total ln_wave_height_prevweek if sample==1, vce(robust)
	eststo: quietly reg ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, vce(robust)
	eststo: quietly reghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month) vce(robust)
	eststo: quietly reghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) vce(robust)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main.tex", style(tex) numbers nomtitles order(ln_riots_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$")) varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-1. Summary statistics for time series analysis
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

sutex arrival_total riots_prevweek wave_height_prevweek ln_arrival_total ln_arrivalmissing_total ln_arrivalmissing1L_total ln_arrivalmissing2L_total ln_arrivalmissing3L_total asinh_arrival_total asinh_arrivalmissing_total ln_arrival_3dayMA ln_arrivalmissing_3dayMA deathrate ln_riots_prevweek asinh_riots_prevweek ln_riots ln_wave_height_prevweek asinh_wave_height_prevweek ln_wave_height ln_wave_height_next3days ln_wave_height_prev3days ln_wave_height_prev4to6days if sample==1, minmax label

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-2. Alternative types of violence and sea conditions and migrant flows to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
	
eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_viol_civ_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_battles_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek ln_viol_civ_prevweek ln_battles_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_altviolence.tex", style(tex) numbers nomtitles order(ln_riots_prevweek ln_viol_civ_prevweek ln_battles_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek  ln_viol_civ_prevweek ln_battles_prevweek   ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_viol_civ_prevweek "\textsc{Violence Against Local Civilians (ln, prior week total)}"  ln_battles_prevweek  "\textsc{Rebel-Government Violence (ln, prior week total)}"  ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-3. Alternative clustering specifications to capture potential temporal autocorrelation in migration to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek  if sample==1, absorb(month dow) cluster(edate) bw(21) 
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek  if sample==1, absorb(month dow) cluster(edate) bw(28) 
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek  if sample==1, absorb(month dow) cluster(edate) bw(56) 
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek  if sample==1, absorb(month dow) cluster(ym) 

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/robustness_StandardErrors.tex", style(tex) order(ln_riots_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$")) numbers nomtitles varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-4. Impact of sea conditions on death rates in Mediterranean Sea
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear

eststo: quietly ivreghdfe deathrate ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo: quietly ivreghdfe deathrate L.ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo: quietly ivreghdfe deathrate L(0/6).ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo: quietly ivreghdfe deathrate ln_wave_height_prev3days if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo: quietly ivreghdfe deathrate ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

// test L2.ln_wave_height L3.ln_wave_height
// ( 1)  L2.ln_wave_height = 0
// ( 2)  L3.ln_wave_height = 0
//       F(  2,   613) =    3.60
//            Prob > F =    0.0280

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/deathrates_fullmodel.tex", style(tex) numbers nomtitles order(ln_wave_height L.ln_wave_height L2.ln_wave_height L3.ln_wave_height L4.ln_wave_height L5.ln_wave_height L6.ln_wave_height ln_wave_height_prev3days ln_wave_height_prevweek) keep(ln_wave_height L.ln_wave_height L2.ln_wave_height L3.ln_wave_height L4.ln_wave_height L5.ln_wave_height L6.ln_wave_height ln_wave_height_prev3days ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_wave_height "\textsc{Wave Height (ln, current)}" L.ln_wave_height "\textsc{Wave Height (ln, lag 1)}" L2.ln_wave_height "\textsc{Wave Height (ln, lag 2)}" L3.ln_wave_height "\textsc{Wave Height (ln, lag 3)}" L4.ln_wave_height "\textsc{Wave Height (ln, lag 4)}" L5.ln_wave_height "\textsc{Wave Height (ln, lag 5)}" L6.ln_wave_height "\textsc{Wave Height (ln, lag 6)}" ln_wave_height_prev3days "\textsc{Wave Height (ln, previous 3 days)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace

eststo clear

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-5. Impact of excluding various types of potentially endogenous riot activity on migrant flows to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_riots_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_riots_noIDP_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_riots_noFUEL_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_riots_noPORTCLOSE_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_riots_noT_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

// ( 1)  ln_riots_prevweek - ln_viol_civ_prevweek = 0
//
//       F(  1,   811) =    0.11
//            Prob > F =    0.7351
	
esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_altriots.tex", style(tex) numbers nomtitles order(ln_riots_prevweek ln_riots_noIDP_prevweek ln_riots_noFUEL_prevweek ln_riots_noPORTCLOSE_prevweek ln_riots_noT_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek  ln_riots_noIDP_prevweek ln_riots_noFUEL_prevweek ln_riots_noPORTCLOSE_prevweek ln_riots_noT_prevweek   ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_riots_noIDP_prevweek "\textsc{Riots (ln, excluding IDP related events)}"  ln_riots_noFUEL_prevweek  "\textsc{Riots (ln, excluding fuel smuggling related events)}"  ln_riots_noPORTCLOSE_prevweek "\textsc{Riots (ln, excluding port closure related events)}" ln_riots_noT_prevweek "\textsc{Riots (ln, excluding economic riots in Tunisia)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-6. Correlation between sea conditions and riot activity 
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly reg ln_riots_prevweek ln_wave_height_prevweek if sample==1, vce(robust)
	eststo: quietly ivreghdfe ln_riots_prevweek ln_wave_height_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly reg ln_riots ln_wave_height if sample==1, vce(robust)
	eststo: quietly ivreghdfe ln_riots ln_wave_height  if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly reg ln_riots ln_wave_height L(1/6).ln_wave_height if sample==1, vce(robust)
	eststo: quietly ivreghdfe ln_riots ln_wave_height L(1/6).ln_wave_height  if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/riots2sea.tex", style(tex) numbers nomtitles order(ln_wave_height_prevweek ln_wave_height L.ln_wave_height L2.ln_wave_height L3.ln_wave_height L4.ln_wave_height L5.ln_wave_height L6.ln_wave_height) keep(ln_wave_height_prevweek ln_wave_height L.ln_wave_height L2.ln_wave_height L3.ln_wave_height L4.ln_wave_height L5.ln_wave_height L6.ln_wave_height) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}" ln_wave_height "\textsc{Wave Height (ln, current)}" L.ln_wave_height "\textsc{Wave Height (ln, lag 1)}" L2.ln_wave_height "\textsc{Wave Height (ln, lag 2)}" L3.ln_wave_height "\textsc{Wave Height (ln, lag 3)}" L4.ln_wave_height "\textsc{Wave Height (ln, lag 4)}" L5.ln_wave_height "\textsc{Wave Height (ln, lag 5)}" L6.ln_wave_height "\textsc{Wave Height (ln, lag 6)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace

eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-7. Evaluating the interaction of riots and sea conditions on migrant flows to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek c.ln_wave_height_prevweek#c.ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek c.ln_wave_height_prevweek#c.ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_interaction.tex", style(tex) numbers nomtitles order(ln_riots_prevweek ln_wave_height_prevweek c.ln_wave_height_prevweek#c.ln_riots_prevweek) keep(ln_riots_prevweek ln_wave_height_prevweek c.ln_wave_height_prevweek#c.ln_riots_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}" c.ln_wave_height_prevweek#c.ln_riots_prevweek "\textsc{Riots} (ln, prior week total) $\times$ \textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-8. Using alternative transformation (inverse hyperbolic sine) to evaluate relationships among riots, sea conditions and migrant flows to Italy
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
 
eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe asinh_arrival_total asinh_wave_height_prevweek asinh_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe asinh_arrivalmissing_total asinh_wave_height_prevweek asinh_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_asinh.tex", style(tex) numbers nomtitles order(ln_riots_prevweek asinh_riots_prevweek ln_wave_height_prevweek asinh_wave_height_prevweek) keep(ln_riots_prevweek asinh_riots_prevweek ln_wave_height_prevweek asinh_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" asinh_riots_prevweek "\textsc{Riots (inverse  hyperbolic sine, prior week total)}"  ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}" asinh_wave_height_prevweek "\textsc{Wave Height (inverse hyperbolic sine, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
 
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-9. Evaluating relationships among riots, sea conditions and migration using moving average of arrivals
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
 
eststo clear 

	eststo: quietly ivreghdfe ln_arrival_3dayMA ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_3dayMA ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_3dayMAarrivals.tex", style(tex) numbers nomtitles order(ln_riots_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$"))  varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-10. Alternative specifications to capture relationship between total migration (arrivals and deaths/missing) and riots and sea conditions 
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
 
eststo clear 
	
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing1L_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing2L_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing3L_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_deathvariation.tex", style(tex) order(ln_riots_prevweek ln_wave_height_prevweek) keep(ln_riots_prevweek ln_wave_height_prevweek) se stats(N r2, labels("Number of Observations" "R$^2$")) numbers nomtitles varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Table SI-11. Riots, sea conditions, and migrant flows to Italy using varying lags
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

eststo clear 

	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrival_total ln_wave_height_prev3days ln_riots_prev3days if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prevweek ln_riots_prevweek if sample==1, absorb(month dow) cluster(edate) bw(14)
	eststo: quietly ivreghdfe ln_arrivalmissing_total ln_wave_height_prev3days ln_riots_prev3days if sample==1, absorb(month dow) cluster(edate) bw(14)

esttab using "~/Dropbox/Migration_in_Med/MiM_Replication/OUTPUT/main_attenuation.tex", style(tex) numbers nomtitles order(ln_riots_prevweek  ln_riots_prev3days  ln_wave_height_prevweek ln_wave_height_prev3days) keep(ln_riots_prevweek  ln_wave_height_prev3days ln_wave_height_prevweek ln_riots_prev3days) se stats(N r2, labels("Number of Observations" "R$^2$")) varlabel(ln_riots_prevweek "\textsc{Riots (ln, prior week total)}" ln_riots_prev3days "\textsc{Riots (ln, prior 3 days total)}" ln_wave_height_prevweek "\textsc{Wave Height (ln, prior week average)}" ln_wave_height_prev3days "\textsc{Wave Height (ln, prior 3 days average)}") star(* 0.10 ** 0.05 *** 0.01) nonotes replace
eststo clear 	

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// FIGURES
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig 2. Trends in Death Rates and Sea Conditions during Sample Period
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//

tw (lowess deathrate edate if sample==1, bw(.025) lc(red) yaxis(1)) (tsline wave_height_prevweek, lp(solid) lc(blue) yaxis(2)), ///
	legend(size(vsmall) pos(12) ring(0) col(2) order(1 "Death Rate" 2 "Wave Height")) ///
	scheme(lean1) ///
	xlab(20461(60)21272 ,nogrid angle(45) format(%tdm-cy)) ylab( ,nogrid) ///
	ytitle("Death Rate", axis(1)) ytitle("Wave Height, prior week", axis(2))  xtitle("Date") 
	
gr export "~/MiM_Replication/OUTPUT/TS_deathrates.pdf", replace

***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig SI-2. Non-parametric local regression between sea conditions and riot activity
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
 
graph drop _all

tw (lpolyci ln_riots ln_wave_height if sample==1, bw(.5)), ytitle("Riots (ln)") xtitle("Wave Height (ln)") legend(order(1 "Confidence Interval" 2 "Local Polynomial Regression") pos(12) col(2)) fysize(150) yline(.1725963, lcolor(red)) name(m1)  graphregion(color(white))

hist ln_wave_height, graphregion(color(white)) yscale(reverse) frac lcolor(black) fcolor(gs12) fysize(40) xtitle("") ylabel(0(.05).1, nogrid) name(m2)

	gr combine m1 m2, col(1) imargin(0 0 0 0) ysize(10) xsize(10) scheme(lean1) graphregion(color(white))

	gr export "~/MiM_Replication/OUTPUT/riots2sea_wRugPlot.pdf", replace
	
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
// Fig SI-3. Using leads and lags of sea conditions to calibrate main effects for forecasting (anticipation) and departure delays	
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
	
eststo clear

eststo m1: ivreghdfe ln_arrival_total L(-7/7).ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m2: ivreghdfe ln_arrival_total L(1/7).ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m3: ivreghdfe ln_arrival_total L(1/2).ln_wave_height if sample==1, absorb(month dow) cluster(edate) bw(14)

coefplot (m1, label("Leads + All Lags") ciopts(recast(rcap))) (m2, label("All Lags") ciopts(recast(rcap))) (m3, label("Main Lags") ciopts(recast(rcap))), vertical yline(0) order(L7.ln_wave_height L6.ln_wave_height L5.ln_wave_height L4.ln_wave_height L3.ln_wave_height L2.ln_wave_height L.ln_wave_height ln_wave_height F.ln_wave_height F2.ln_wave_height F3.ln_wave_height F4.ln_wave_height F5.ln_wave_height F6.ln_wave_height F7.ln_wave_height) coeflabels(F7.ln_wave_height = "+7" F6.ln_wave_height = "+6" F5.ln_wave_height = "+5" F4.ln_wave_height = "+4" F3.ln_wave_height = "+3" F2.ln_wave_height = "+2" F.ln_wave_height = "+1" ln_wave_height = "0" L.ln_wave_height = "-1" L2.ln_wave_height = "-2" L3.ln_wave_height = "-3" L4.ln_wave_height = "-4" L5.ln_wave_height = "-5" L6.ln_wave_height = "-6" L7.ln_wave_height = "-7") legend(ring(0) pos(7) col(1)) title("Impact of Sea Conditions on Total Arrivals (ln)") xtitle("Days to Arrival")

	gr export "~/MiM_Replication/OUTPUT/leadslags_daily.pdf", replace

eststo clear

eststo m1: ivreghdfe ln_arrival_total ln_wave_height_next3days ln_wave_height_prev3days ln_wave_height_prev4to6days  if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m2: ivreghdfe ln_arrival_total ln_wave_height_prev3days ln_wave_height_prev4to6days  if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m3: ivreghdfe ln_arrival_total ln_wave_height_prev3days  if sample==1, absorb(month dow) cluster(edate) bw(14)

coefplot (m1, label("Leads + All Lags") ciopts(recast(rcap))) (m2, label("All Lags") ciopts(recast(rcap))) (m3, label("Main Lags") ciopts(recast(rcap))), vertical yline(0) order(ln_wave_height_prev4to6days ln_wave_height_prev3days ln_wave_height_next3days) coeflabels(ln_wave_height_prev4to6days = "4-6 days (before)" ln_wave_height_prev3days = "1-3 days (before)" ln_wave_height_next3days = "1-3 days (after)") legend(ring(0) pos(7) col(1)) title("Impact of Sea Conditions on Total Arrivals (ln)") xtitle("Time to Arrival")

	gr export "~/MiM_Replication/OUTPUT/leadslags_daybins.pdf", replace

eststo clear

eststo m1: ivreghdfe ln_arrival_total ln_wave_height_next7days ln_wave_height_prevweek  ln_wave_height_prev8to14days if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m2: ivreghdfe ln_arrival_total ln_wave_height_prevweek  ln_wave_height_prev8to14days if sample==1, absorb(month dow) cluster(edate) bw(14)
eststo m3: ivreghdfe ln_arrival_total ln_wave_height_prevweek  if sample==1, absorb(month dow) cluster(edate) bw(14)

coefplot (m1, label("Leads + All Lags") ciopts(recast(rcap))) (m2, label("All Lags") ciopts(recast(rcap))) (m3, label("Main Lags") ciopts(recast(rcap))), vertical yline(0) order(ln_wave_height_prev8to14days ln_wave_height_prevweek ln_wave_height_next7days) coeflabels(ln_wave_height_prev8to14days = "Two weeks (before)" ln_wave_height_prevweek = "One week (before)" ln_wave_height_next7days = "One week (after)") legend(ring(0) pos(7) col(1)) title("Impact of Sea Conditions on Total Arrivals (ln)") xtitle("Time to Arrival")

	gr export "~/MiM_Replication/OUTPUT/leadslags_weekbins.pdf", replace
	
***********//***********//***********//***********//***********//***********//***********//***********//***********//***********//
