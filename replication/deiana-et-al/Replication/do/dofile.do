*-------------------------------------------------------------------------------

* Setting work directory
*=============================================================================
*
clear all
set scheme s1mono
set matsize 11000
set maxiter 16000
set maxvar 32767
set more off
*

/*
The users will need to install the following packages 

ssc install reghdfe
* Install ftools (remove program if it existed previously)
cap ado uninstall ftools
net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/")

* Install reghdfe 
cap ado uninstall reghdfe
net install reghdfe, from("https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/src/")

* Install parallel, if using the parallel() option; don't install from SSC
cap ado uninstall parallel
net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace
mata mata mlib index

cap ado uninstall ivreghdfe
cap ssc install ivreg2 // Install ivreg2, the core package
net install ivreghdfe, from(https://raw.githubusercontent.com/sergiocorreia/ivreghdfe/master/src/)

* This is to check if it works correctly
cap ado uninstall ftools
cap ado uninstall reghdfe
cap ado uninstall ppmlhdfe

ssc install ftools
ssc install reghdfe
ssc install ppmlhdfe

clear all
ftools, compile
reghdfe, compile

* Test program
sysuse auto, clear
reghdfe price weight, a(turn)
ppmlhdfe price weight, a(turn)


****************
* Packages to create map  
ssc install spmap
ssc install shp2dta
ssc install mif2dta
*/ 

global whitebox 		style.editstyle boxstyle(shadestyle(color(white)) linestyle(color(white))) ///
						inner_boxstyle(shadestyle(color(white)) linestyle(color(white))) editcopy
cd "~/Replication" 


*-------------------------------------------------------------------------------
* 								MAIN TEXT  		
*-------------------------------------------------------------------------------


*=============================================================================
* Figure 1 
*=============================================================================
use "data/data_fig1.dta"

sum y, d
local min `r(min)'
local max `r(max)'
twoway 	///
		(bar cross		  y,		 ///
		lc(dknavy*.80) fc(cranberry*.80) lpattern(solid) lwidth(small)) ///
		///
		(connected dead 		  y,		 ///
		lc(dknavy*.80) mc(dknavy*.80) msymbol(o) lpattern(solid) lwidth(medthick) yaxis(2)) ///
		, ///				
		///
		ylab(0(50000)200000, axis(1) tstyle(textstyle(size(small)))) ///
		ylab(0(1000)4000, axis(2) tstyle(textstyle(size(small)))) ///
		xlab(`min'(1)`max', axis(1) tstyle(textstyle(size(small)))) ///
		xtitle("", size(small)) ///
		ytitle("Total Crossings", axis(1) size(small)) ///
		ytitle("Total Deaths", axis(2) size(small)) ///
		legend(symys(*.5) symxs(*.5) order(1 "Crossings" 2 "Deaths") ///
		size(small) region(lwidth(none)) row(1) col(2) keygap(*.2) colgap(*.2))
		gr_edit $whitebox
graph export "graphs/fig1.png", replace


*=============================================================================
* Table 1 
*=============================================================================
* No data analysis

*=============================================================================
* Figure 2 
*=============================================================================
use "data/data_fig2.dta"

twoway 	///
(bar a ym,		 ///
fcolor(dknavy*2) lcolor(dknavy*2)  lwidth(thin) yaxis(1)) ///
(rbar a B ym,		 ///
fcolor(dknavy*.75) lcolor(dknavy*.75)  lwidth(thin) yaxis(1)) ///
(rbar B C ym,		 ///
fcolor(dknavy*.25) lcolor(black)  lwidth(thin) yaxis(1)) ///
(rbar C D ym,		 ///
fcolor(dknavy*.25) lcolor(dknavy*.25)  lwidth(thin) yaxis(1)) ///
, ///
xlab( ///
1 "Hermes" 2 " " 3 " " 4 " " 5 " " 6 "Mare Nostrum" 7 " " 8 " " 9 "2014m1" 10 " " ///
11 " " 12 " " 13 " " 14 "2014m6" 15 " " 16 " " 17 " " 18 " " 19 "Triton I" 20 " " /// 
21 "2015m1" 22 " " 23 " " 24 " " 25 "Triton II" 26 " " 27 " " 28 " " 29 " " 30 " " /// 
31 " " 32 " " 33 "2016m1" 34 " " 35 " " 36 " " 37 " " 38 "2016m6" 39 " " 40 " " /// 
41 " " 42 " " 43 " " 44 " " 45 "2017m1" 46 " " 47 " " 48 " " 49 " " 50 " " /// 
51 "Minniti Code" 52 " " 53 " " 54 " " 55 " "  56 " "  ///
57 "Themis" 58 " " 59 " " 60 " " 61 " " 62 "2018m6" 63 " " 64 " " 65 " " 66 " " 67 " "  68 " " ///
69 "2019m1" 70 " " 71 " " 72 " " 73 " " 74 "2019m6" 75 " " 76 " " 77 " " 78 " " 79 " "  80 " " ///
81 "2020m1" 82 " " 83 " " 84 " " 85 " " 86 "2020m6" 87 " " 88 " " 89 " " 90 " " 91 " "  92 "2020m12" ///
) ///
					xlab(, /*format(%tdmonCCYY)*/  angle(90) axis(1) tstyle(textstyle(size(small)))) ///
					ylab(, /*format(%tdmonCCYY)*/  angle(90) axis(1) tstyle(textstyle(size(small)))) ///
					ytitle("Crossings Attempts", size(small) axis(1)) ///
					xtitle("", size(small)) ///
					xline(6 19 25 51 58, lwidth(midthin) lpattern(dash)  lcolor(black*.8) ) ///
					legend(symys(*.7) symxs(*.2) order(1 "Inflatable boat" 2 "Sturdy boat" 3 "Other boat" 4 "Unknown")  ///
					/*ring(0) pos(1)*/ col(4)  size(small) region(lwidth(none)))					
					gr_edit $whitebox
graph export "graphs/fig2.png", replace


*=============================================================================
* Figure 3 
*=============================================================================
use "data/data_fig3.dta", replace 

sum my, d
local min `r(min)'
local max `r(max)'
twoway 	///
		(line atcross 		  my if treat == 1,		 ///
		lc(cranberry) fc(none) lpattern(none) lwidth(medium)) ///
		(line atcross_inf 	  my if treat == 1,		 ///
		lc(cranberry) fc(none) lpattern(dash) lwidth(medium)) ///
		///
		(line atcross 		  my if treat == 0,		 ///
		lc(dknavy) fc(none) lpattern(none) lwidth(medium) yaxis(1)) ///
		(line atcross_noinf   my if treat == 0,		 ///
		lc(dknavy) fc(none) lpattern(dash) lwidth(medium) yaxis(1)) ///
		, ///				
		ttext(800 695 "Minniti Code", ///
		size(small) color(black)) ///
		///
		xline(690, lwidth(medium) lpattern(dash) lcolor(black*.7) ) ///
		ylab(, axis(1) tstyle(textstyle(size(small)))) ///
		xlab(`min'(6)`max', axis(1) tstyle(textstyle(size(small)))) ///
		xtitle("", size(small)) ///
		ytitle("Average Attempted Crossings", axis(1) size(small)) ///
		legend(symys(*.5) symxs(*.5) order(1 "Total Libya" 2 "Inflatable Libya" ///
										   3 "Total Tunisia" 4 "Not Inflatable Tunisia" ///
											) ///
		size(small) region(lwidth(none)) row(2) col(2) keygap(*.2) colgap(*.2))
		gr_edit $whitebox
graph export "graphs/fig3.png", replace



*=============================================================================
* Figure 4 
*=============================================================================
use "data/data_fig4.dta", clear 

forvalues y = 2010/2017 { // week
cap drop sdg_`y'
sum sdg if year == `y',d 	   
cap gen sdg_`y' = `r(mean)'
}

reg angle i.year
forvalues t = 2010/2017 { // week
	   
cap gen b_`t' = _b[`t'.year]
cap gen se_`t' = _se[`t'.year]
}

keep if _n == 1
keep b_* se_* sdg_*
g seq = 1
reshape long se_ b_ sdg_, i(seq) j(k) string
rename *_ *

cap drop ci_u
gen ci_u = b + 1.64*se

cap drop ci_l
gen ci_l = b - 1.64*se

cap drop ci_u2
gen ci_u2 = b + 1.96*se

cap drop ci_l2
gen ci_l2 = b - 1.96*se


destring k, replace
rename k year

#delimit ;
twoway
	bar sdg year, fcolor(dknavy%35*.75) fintensity(100) barw(0.95) bargap(-50) lcolor(black) lwidth(medium) yaxis(2)	 ||

	rcap ci_u ci_l year, lcolor(cranberry*0.95) yaxis(1)	  ||
	rcap ci_u2 ci_l2 year, lcolor(cranberry*0.75) yaxis(1)	  ||
	scatter b year, mcolor(cranberry) msize(small) msymbol(O) yaxis(1)	  
	
	yline(0, lcolor(black) lpattern(shortdash) lwidth(vthin))
	xlabel(2011(1)2017, labsize(small) notick angle(0) valuelabel)

	ylabel(, labsize(small) axis(1))
	ytitle("St. dev. angle of departure", size(small) axis(2))
	ylabel(, labsize(small) axis(2))
	ytitle("Angle of departure", size(small) axis(1))

	title(" ", size(small))
	text(5 2014.5  "{it:Libya}", size(small))
	text(-5 2014.5 "{it:Tunisia}", size(small))
	legend(order(1 "St. dev." 4 "Angle of departure" /*9 "Mean Angle: 41"*/))
	legend(size(small))
	legend(r(1) c(2))
	legen(region(color(white))) 
			
	xtitle(" ", size(small) )
	xscale(titlegap(5))  
	yscale(titlegap(5))  

	graphregion(color(white))
	xsize(20) ysize(12)

;
#delimit cr
graph export "graphs/fig4.png", replace



*=============================================================================
* Table 2 
*=============================================================================
* Note that you need to copy-paste the output to the latex file (tab2.tex)
use "data/data_tab2a.dta", clear
sum totacross swh_lyb_3 max_wave swh_lyb_1 swh_tun_1 swh_tun_2 swh_tun_3 swh_alg_1 fr_across3 fr_across1_3 fr_across_123


use "data/data_tab2b.dta", clear 
sum atcross onda fr_across_infl if treat==1
sum atcross onda fr_across_infl if treat==0



*=============================================================================
* Figure 5 - 6 - 7 
*=============================================================================
* No data analysis

*=============================================================================
* Table 3
*=============================================================================
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac ///
			  i.weekanno, family(poisson) vce(hac nwest 28) t(data)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tab3.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table 4
*=============================================================================
use "data/data_tab4.dta", clear
global toestout ""
label var atcross " "

ppmlhdfe atcross _treat pre_treat, absorb(w##y preMin) cluster(w##y)
estimates store reg_1
estadd scalar nobs e(N)
global toestout "${toestout} reg_1"
		
ppmlhdfe atcross c.onda##c.fr_across_infl##i.treat##i.preMin, absorb(w##y##treat) cluster(w##y)
estimates store reg_2
estadd scalar nobs e(N)
global toestout "${toestout} reg_2"

*** Export table
#delimit ;
estout ${toestout}
using "tables/tab4.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.0000000000001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{11.5cm}C{3cm}C{3cm}}
\toprule\toprule
	& (1)   & (2)      \\
	& \multicolumn{2}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
&            &                                \\
		)
prefoot(
	\midrule
Week-Year FE & \checkmark &    \\
Week-Year-From Lybia FE & & \checkmark   \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(_treat pre_treat onda c.onda#c.fr_across_infl 1.preMin#c.onda 1.preMin#c.onda#c.fr_across_infl 1.treat#c.onda 1.treat#c.onda#c.fr_across_infl  1.treat#1.preMin#c.onda 1.treat#1.preMin#c.onda#c.fr_across_infl)
order(pre_treat _treat 1.treat onda c.onda#c.fr_across_infl 1.preMin#c.onda 1.preMin#c.onda#c.fr_across_infl 1.treat#c.onda 1.treat#c.onda#c.fr_across_infl  1.treat#1.preMin#c.onda 1.treat#1.preMin#c.onda#c.fr_across_infl)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************

*=============================================================================
* Figure 8
*=============================================================================
use "data/data_fig8.dta", replace 

twoway ///
(rarea lci uci prezzi, color(cranberry%50) lw(medthin) lp(solid)) ///
(scatter beta prezzi, ///
mlw(medthin) msize(small) msy(o) mlc(cranberry*.9)  mfc(cranberry*.9)) ///	
, ///
xlab(.5(0.25)5, ///
labcolor(black) labsize(medsmall) tstyle(textstyle(size(medsmall))) ) ///
ylab(, labsize(medsmall) labcolor(black) tstyle(textstyle(size(vlarge))))  ///
xtit("(P{sub:S} - P{sub:U}) / P{sub:U}", size(medsmall)) xlab(, labgap(1) tstyle(textstyle(size(medsmall)))) ///
ytit("Simulated {&theta}", size(medsmall)) xlab(, labgap(2) tstyle(textstyle(size(medsmall))) )  ///
legend(off)
gr_edit $whitebox
graph export "graphs/fig8.png", replace



*=============================================================================
* Table 5
*=============================================================================
use "data/data_tab_main2.dta", clear 
****
global toestout ""

***

foreach outcome of varlist fr_across3 fr_across_1_3 fr_across_1_2_3 fr_across4 fr_across5   {

*********
reghdfe `outcome' policy2 policy3 policy4, ///
				  absorb(week) vce(cluster week) 
estimates store reg_`outcome'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`outcome'"

}

   su fr_across3 if policy1 == 1
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across_1_3 if policy1 == 1
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_1_2_3 if policy1 == 1
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

   su fr_across4 if policy1 == 1
local avg_fr_4 = r(mean)
local avg_fr_4 = round(`avg_fr_4',.01)

   su fr_across5 if policy1 == 1
local avg_fr_5 = r(mean)
local avg_fr_5 = round(`avg_fr_5',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tab5.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{5cm}C{3cm}C{3cm}C{3cm}C{3cm}C{3cm}}
\toprule\toprule
  & (1)         & (2)       & (3)     & (4)      & (5)       \\
Fraction of Attempted Crossings       & Inflatable   & Inflatable + & Inflatable + & Fishing & Motor                     \\
				           &    & Unknown & Unknown + &  &                      \\
				           &    &  & Other &  &                      \\
\hline
                         &        &       &          &     &                           \\
)
prefoot(
	\midrule
Pre MN Mean Outcome  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' & `avg_fr_4' & `avg_fr_5' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(policy2 policy3 policy4 ) 
order(policy2 policy3 policy4 )
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************

*=============================================================================
* Figure 9
*=============================================================================
use "data/data_fig9.dta", replace 
global y1 gom_su_tot 
global peso atcross
global funn lfit   


	sum b, d
	local mmin `r(min)'
	local mmax `r(max)'

	#delimit ;
	tw ///
	scatter  $y1 b if b<=0 [aw=$peso], mlcolor(cranberry*.85%70) mfcolor(cranberry*.85%70) msize(medsmall) msy(O) mlwidth(medsmall) || ///
	$funn    $y1 b if b<=0, color(cranberry*.85%70) lpattern(shortdash) || ///
	scatter  $y1 b if b>0  [aw=$peso],  mcolor(cranberry*.9%80) msize(medsmall) msy(Oh) mlwidth(medsmall) || /// 
	$funn    $y1 b if b>0, color(cranberry*.9%80) lpattern(longdash) ///
	///
	xline(0, lcolor(black) lpattern(dash) lwidth(medsmall))
	xlabel(`mmin'(10)`mmax', labsize(medsmall) labgap(3) notick angle(0) valuelabel)
	ylab(, labsize(medsmall) axis(1)) ///
	title(, size(medsmall)) ///
	legend(region(color(white))) ///		
	xtitle("Days relative to August 7, 2017", size(medsmall) placement(center)) ///
	graphregion(color(white)) ///
	xsize(20) ysize(12)	///
	ytitle("Average Fraction of Attempted in Inflatable", size(medsmall) placement(left)) ///
	legend(symys(*.5) symxs(*.5) order(1 " " 3 " Fraction of Inflatable") ///
	size(medsmall) region(lwidth(none)) row(1) col(4) keygap(*.2) colgap(*.2))
	; 		
	#delimit cr 
graph export "graphs/fig9.png", replace



*-------------------------------------------------------------------------------
* 								APPENDIX A  		
*-------------------------------------------------------------------------------

*=============================================================================
* Figure A1
*=============================================================================
use "data/data_figA1.dta", clear  

twoway 	///
(bar total mdate, 	 ///
			color(navy*.3) yaxis(1)) ///
(bar f_tot_cmed mdate, 	 ///
			color(cranberry*.9) yaxis(1)) ///
			, xlab(, angle(0) ///
			tstyle(textstyle(size(small)))) xtit("", margin(medsmall)) graphregion(margin(medlarge)) ///
			ytitle("Percentage of People Rescued", axis(1) size(small) margin(small)) ///
			ylab(0(.2)1, axis(1) tstyle(textstyle(size(small)))) ///
			legend(order(1 "Libya" 2 "Italy") ///
			size(small) cols(9) symxsize(6) region(lstyle(none)))
			gr_edit $whitebox
graph export "graphs/figA1.png", replace



*-------------------------------------------------------------------------------
* 								APPENDIX C 		
*-------------------------------------------------------------------------------



*=============================================================================
* Figure C1
*=============================================================================
use "data/data_figC1.dta", clear 

#delimit ;
graph bar (sum) pc_arrivals,
over(nationality,  label(angle(90) tstyle(textstyle(size(small))))) 
over(year,  label(angle(90) tstyle(textstyle(size(small))))) 
stack   asyvars 
legend(cols(6) size(small)) 
ytitle("", size(small))
ylab(, axis(1) tstyle(textstyle(size(small)))) 
legend(size(small) cols(5) symxsize(4) region(lstyle(none))) 
bar(1, color(ebblue*0.5))
bar(2, color(yellow*0.7)) 
bar(3, color(blue*.7)) 
bar(4, color(mint*.7)) 
bar(5, color(orange*0.7)) 
bar(6, color(pink*0.3)) 
bar(7, color(red*0.7)) 
bar(8, color(black*0.9)) 
bar(9, color(dkorange*0.5)) 
bar(10, color(green)) 
bar(11, color(midgreen*.5)) 
bar(12, color(emerald)) 
bar(13, color(cyan*.8))
bar(14, color(gold*.6)) 
bar(15, color(edkblue*.7)) 
bar(16, color(navy*.8)) 
bar(17, color(lime*.8)) 
bar(18, color(cranberry*.9))
bar(19, color(magenta*.5)) 
bar(20, color(gray*.8)) 
subtitle("");
#delimit cr
gr_edit $whitebox
graph export "graphs/figC1.png", replace



*=============================================================================
* Figure C2
*=============================================================================
use "data/data_figC2.dta", clear 

cap drop modate
gen modate = ym(year, month) 
format modate %tm 

tsset modate, monthly
*-------------------------------------------------------------------------------

twoway 	///
(connected atot_cmed month if year == 2009,		 ///
lcolor(grey) lpattern(solid) msymbol(O) mcolor(grey) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2010,		 ///
lcolor(cranberry*.75) lpattern(shortdash) msymbol(t) mcolor(cranberry*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2011,		 ///
lcolor(green*.75) lpattern(longdash) msymbol(D) mcolor(green*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2012,		 ///
lcolor(olive*.75) lpattern(dot) msymbol(T) mcolor(olive*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2013,		 ///
lcolor(black*.75) lpattern(longdash) msymbol(o) mcolor(black*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2014,		 ///
lcolor(blue*.75) lpattern(longdash) msymbol(X) mcolor(blue*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2015,		 ///
lcolor(ebblue*.75) lpattern(dash) msymbol(S) mcolor(ebblue*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2016,		 ///
lcolor(red*.75) lpattern(dash_dot) msymbol(d) mcolor(red*.75) msize(small) lwidth(small) yaxis(1)) ///
(connected atot_cmed month if year == 2017,		 ///
lcolor(magenta*.5) lpattern(dash) msymbol(s) mcolor(magenta*.5) msize(small) lwidth(small) yaxis(1)) ///
, xlab(1 "Jan" 2 "Feb" 3 "Mar" 4 "Apr" 5 "May" 6 "Jun" 7 "Jul" 8 "Aug" 9 "Sep" 10 "Oct" 11 "Nov" 12 "Dec", tstyle(textstyle(size(small)))) xtit("", margin(medsmall)) graphregion(margin(medlarge))  ///
ytitle("Attempted Crossings", axis(1) size(small) margin(medium) ) ///
ylab(, axis(1) tstyle(textstyle(size(small)))) ///
legend(order(1 "2009" 2 "2010" 3 "2011" 4 "2012" 5 "2013" 6 "2014" 7 "2015" 8 "2016" 9 "2017") ///
size(small) cols(9) symxsize(4) region(lstyle(none)))
gr_edit $whitebox
graph export "graphs/figC2.png", replace



*=============================================================================
* Figure C3
*=============================================================================
use "data/data_figC3.dta", clear 

sum year
local min `r(min)'
local max `r(max)'

**********
twoway  ///
(connected netimport2010 year if type==1,		 ///
		lc(dknavy*.80) mc(dknavy*.80) msymbol(o) lpattern(shortdash) lwidth(medsmall)) ///
(connected netimport2010 year if type==-1,		 ///
		lc(dknavy*.80) mc(dknavy*.80) msymbol(t) lpattern(solid) lwidth(medsmall)) ///
		, xtit("") xlab(, tstyle(textstyle(size(small))) )  ///
		ytit("") ylab(, tstyle(textstyle(size(small))) )  ///
		xlabel(`min'(1)`max',tstyle(textstyle(size(small)))) ///
		legend(label(1 "Rubber Boats and Similar Boats") ///
		label(2 "Ferries and Similar Vessels") size(small) cols(1) pos(11) ring(0) region(lwidth(none))) ///
		plotregion(margin(small)) 
		gr_edit $whitebox
graph export "graphs/figC3a.png", replace


**********
sum year
local min `r(min)'
local max `r(max)'

twoway  ///
(connected netimport2010 year if type==0,		 ///
		lc(dknavy*.80) mc(dknavy*.80) msymbol(o) lpattern(shortdash) lwidth(medsmall)) ///
		, xtit("") xlab(, tstyle(textstyle(size(small))) )  ///
		ytit("") ylab(, tstyle(textstyle(size(small))) )  ///
		xlabel(`min'(1)`max',tstyle(textstyle(size(small)))) ///
		legend(label(1 "Net Import of Lifejackets") size(small) cols(2) region(lwidth(none))) ///
		plotregion(margin(small))
		gr_edit $whitebox
graph export "graphs/figC3b.png", replace



*=============================================================================
* Figure C4
*=============================================================================
use "data/maps/attr", clear
spmap using "data/maps/coord", id(stid)  fcolor(bluishgray)  ocolor(bluishgray)  ///
		label(data("data/data_figC4.dta") xcoord(x_lon) ycoord(y_lat) /// 
		label(angle1) color(red blue) by(positivo) size(*.5 *.5) position(9) length(26)) 
graph export "graphs/figC4.png", replace



*=============================================================================
* Table C1
*=============================================================================
use "data/data_tabC1.dta", clear 
egen minDist = rowmin(distT distB distA)

global EUpolicyNM newp4 newp5 newp6 newp7 newp8 newp9 newp10 newp11   

global toestout ""

			
ppmlhdfe atot_cmed $EUpolicyNM, absorb(week) vce(cluster meseanno)
estimates store reg_1
estadd scalar nobs e(N)
global toestout "${toestout} reg_1"

reghdfe index_s_nomiss $EUpolicyNM, absorb(week) vce(cluster meseanno)	
estimates store reg_2
estadd scalar nobs e(N)
global toestout "${toestout} reg_2"

reghdfe distT $EUpolicyNM, absorb(week) vce(cluster meseanno) 
estimates store reg_3
estadd scalar nobs e(N)
global toestout "${toestout} reg_3"

reghdfe distB $EUpolicyNM, absorb(week) vce(cluster meseanno) 
estimates store reg_4
estadd scalar nobs e(N)
global toestout "${toestout} reg_4"

reghdfe distA $EUpolicyNM, absorb(week) vce(cluster meseanno) 
estimates store reg_5
estadd scalar nobs e(N)
global toestout "${toestout} reg_5"

reghdfe minDist $EUpolicyNM, absorb(week) vce(cluster meseanno) 
estimates store reg_6
estadd scalar nobs e(N)
global toestout "${toestout} reg_6"

reghdfe distL $EUpolicyNM, absorb(week) vce(cluster meseanno) 
estimates store reg_7
estadd scalar nobs e(N)
global toestout "${toestout} reg_7"


su atot_cmed if preper==1, det
local a = r(mean)
local a = round(`a',1)

su index_s_nomiss if preper==1, det	
local b = r(mean)
local b = round(`b',.01)

su distT if preper==1, det	
local c = r(mean)
local c = round(`c',1)

su distB if preper==1, det	
local d = r(mean)
local d = round(`d',1)

su distA if preper==1, det	
local f = r(mean)
local f = round(`f',1)

su minDist if preper==1, det	
local g = r(mean)
local g = round(`g',1)

su distL if preper==1, det	
local h = r(mean)
local h = round(`h',1)

*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC1.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{5cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}}
\toprule\toprule 
 & (1)         & (2)    & (3)    & (4)    & (5)    & (6) & (7)         \\
 & \textbf{Crossing Attempts} & Crossing Risk & \multicolumn{5}{c}{Distance (in km) to:} \\    \cline{4-8}
&       &     &   Tripoli    &  Bengazi    & Al Huwariyah & Min (Tripoli & Lampedusa \\
&       			  &     &       &    & &  Bengazi \& &  \\
 & &     &       &    & &  Al Huwariyah) &  \\
\hline          
&       &     &       &    &  &&  \\
		)
prefoot(
	\midrule
Week FE & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark \\
{\it Pre SAR Period Statistics} &  &  & &  &  & &  \\
Pre Mean Outcome & `a' & `b' & `c' & `d' & `f' & `g' & `h' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
drop(_cons) 
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C2
*=============================================================================
* No data analysis


*=============================================================================
* Table C3
*=============================================================================
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
ivreg2 totacross onda onda_frac onda_post onda_post_frac i.weekanno, ///
				 partial(i.weekanno) bw(28) kernel(bar) robust 
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC3.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
 &  &  &   \\
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C4
*=============================================================================
* Note that you need to manually copy-paste the different standard errors to reproduce the latex file (tabC4.tex)

********************************************************************************
* TAB C4.A - cluster month-year
********************************************************************************
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
ppmlhdfe totacross  onda onda_frac onda_post onda_post_frac, ///
				  absorb(week##anno) vce(cluster meseanno)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC4a.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


********************************************************************************
* TAB C4.B - cluster week-year
********************************************************************************

global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
ppmlhdfe totacross  onda onda_frac onda_post onda_post_frac, ///
				  absorb(week##anno) vce(cluster weekanno)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC4b.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************

********************************************************************************
* TAB C4.C - changing bandwidth 21
********************************************************************************
global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac ///
			  i.weekanno, family(poisson) vce(hac nwest 21) t(data)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC4c.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************

********************************************************************************
* TAB C4.C - changing bandwidth 21
********************************************************************************
global toestout ""
label var totacross " "

***

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac ///
			  i.weekanno, family(poisson) vce(hac nwest 14) t(data)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC4d.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C5
*=============================================================================
use "data/data_tab_main.dta", clear 

tsset data 
gen Lswh_lyb_3 = L.swh_lyb_3
cap drop max_wave
egen max_wave = rowmax(Lswh_lyb_3 swh_lyb_3)
drop swh_lyb_3

global toestout ""
label var totacross " "

***

foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = Lswh_lyb_3
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = Lswh_lyb_3 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = Lswh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = Lswh_lyb_3 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac i.weekanno, ///
			  family(poisson) vce(hac nwest 28) t(data)

estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}


foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = max_wave
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = max_wave * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = max_wave * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = max_wave * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac i.weekanno, ///
			  family(poisson) vce(hac nwest 28) t(data)

estimates store _reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} _reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC5.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{2cm}C{2cm}C{2cm}C{2cm}C{2cm}C{2cm}}
\toprule\toprule
 & (1)   & (2)   & (3)      & (4)   & (5)   & (6) \\
& \multicolumn{6}{c}{\textbf{Crossing Attempts}}    \\
\hline
      &       &       &       &       &       &             \\
& \multicolumn{6}{c}{\textbf{Definition of Unsafe Boat}} \\
\cline{2-7}
& \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\

 & \textit{} & \textit{Unknown} & \textit{Unknown +}  & \textit{} & \textit{Unknown} & \textit{Unknown +}  \\

 & \textit{} & \textit{} & \textit{Other}  \textit{} & \textit{} & \textit{Other}  \\
\hline
		)
prefoot(
	\midrule
&        &         &        &        &          &          \\
& \multicolumn{3}{c}{\textbf{Wave Height in Tripoli (t-1)}} & \multicolumn{3}{c}{\textbf{Max Wave Height  in Tripoli (t and t-1)}} \\
Week-Year FE & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark \\
{\it Pre SAR Period Statistics} &  &  & &   &  &  \\
Mean Total Attempt  & `c' & `c' & `c'  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f'  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' &  `avg_fr_3' & `avg_fr_13' & `avg_fr_123'  \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************

*=============================================================================
* Table C6
*=============================================================================
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "

***
gen swh_lyb_32 = swh_lyb_3^2

local append replace 
			foreach barca of varlist fr_across3 fr_across1_3 fr_across_123 {

cap drop onda
gen onda = swh_lyb_32
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = swh_lyb_32 * `barca'
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = swh_lyb_32 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = swh_lyb_32 * postM * `barca'
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac ///
			  i.weekanno, family(poisson) vce(hac nwest 28) t(data)
estimates store reg_`barca'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`barca'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across1_3 if postM == 0
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su onda if postM == 0
local f = r(mean)
local f = round(`f',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC6.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	  & \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	\cline{2-4}      & \textit{Inflatable} & \textit{Inflatable +} & \textit{Inflatable +} \\
	& \textit{} & \textit{Unknown} & \textit{Unknown +}  \\
	& \textit{} & \textit{} & \textit{Other}  \\
	
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
{\it Pre SAR Period Statistics} &  &  &   \\
Mean Total Attempt  & `c' & `c' & `c' \\
Mean Wave Height  & `f' & `f' & `f' \\
Mean Frac. Unsafe Boat  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C7
*=============================================================================
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "


local append replace 
	foreach onda in swh_lyb_1 swh_tun_1 swh_tun_2 swh_tun_3 swh_alg_1   {

cap drop onda
gen onda = `onda'
label var onda 	"Wave Height"
cap drop onda_frac
gen onda_frac = `onda' * fr_across_123
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post
gen onda_post = `onda' * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = `onda' * postM * fr_across_123
label var onda_post_frac 	"Wave Height * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac i.weekanno, ///
			  family(poisson) vce(hac nwest 28) t(data)
			  
estimates store reg_`onda'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`onda'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across_123 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

su swh_lyb_1 if postM == 0
local f1 = r(mean)
local f1 = round(`f1',.01)

su swh_tun_1 if postM == 0
local f2 = r(mean)
local f2 = round(`f2',.01)

su swh_tun_2 if postM == 0
local f3 = r(mean)
local f3 = round(`f3',.01)

su swh_tun_3 if postM == 0
local f4 = r(mean)
local f4 = round(`f4',.01)

su swh_alg_1 if postM == 0
local f5 = r(mean)
local f5 = round(`f5',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC7.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}C{2.5cm}}
\toprule\toprule
	& (1)   & (2)   & (3)  & (4)  & (5) \\
	& \multicolumn{5}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	&       &       &     &  &  \\
	  & \multicolumn{5}{c}{\textbf{Definition of Unsafe Boat}} \\
	  & \multicolumn{5}{c}{\textit{Inflatable + Unknown + Other}} \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark & \checkmark & \checkmark \\
Wave measured in  & \textbf{Zuwara}   & \textbf{Monastir}  & \textbf{Al Huwariyah} & \textbf{Djerba}  & \textbf{Annaba}  \\
& Libya    & \multicolumn{3}{c}{Tunisia}        & Algeria \\
{\it Pre SAR Period Statistics} &  &  &  &  & \\
Mean Total Attempt  & `c' & `c' & `c' & `c' & `c' \\
Mean Wave Height  & `f1' & `f2' & `f3' & `f4' & `f5'  \\
Mean Frac. Unsafe Boat  & `avg_fr_123' & `avg_fr_123' & `avg_fr_123' & `avg_fr_123' & `avg_fr_123' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C8
*=============================================================================
use "data/data_tab_main.dta", clear 
global toestout ""
label var totacross " "

egen x=xtile(swh_lyb_3), n(100) 

gen hw50 = x>=50
gen hw75 = x>=75
gen hw90 = x>=90

local append replace 
	foreach onda in  hw50 hw75 hw90    {

cap drop onda
gen onda = `onda'
label var onda 	"\textbf{1}[Bad weather]"
cap drop onda_frac
gen onda_frac = `onda' * fr_across3
label var onda_frac 	"\textbf{1}[Bad weather] * Fr. Boat"
cap drop onda_post
gen onda_post = `onda' * postM
label var onda_post 	"\textbf{1}[Bad weather] * Post SAR"
cap drop onda_post_frac
gen onda_post_frac = `onda' * postM * fr_across3
label var onda_post_frac 	"\textbf{1}[Bad weather] * Post SAR * Fr. Boat"

	
*********
glm totacross onda onda_frac onda_post onda_post_frac i.weekanno, ///
			  family(poisson) vce(hac nwest 28) t(data)
			  
estimates store reg_`onda'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`onda'"

}

su totacross if postM == 0
local c = r(mean)
local c = round(`c',1)

   su fr_across3 if postM == 0
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC8.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3.5cm}C{3.5cm}C{3.5cm}}
	\toprule\toprule
	& (1)   & (2)   & (3)   \\
	& \multicolumn{3}{c}{\textbf{Crossing Attempts}}                                          \\
	\hline
	&       &       &         \\
	& \multicolumn{3}{c}{\textbf{Definition of Unsafe Boat}} \\
	& \multicolumn{3}{c}{\textit{Inflatable}} \\
	\hline
		)
prefoot(
	\midrule
Week-Year FE & \checkmark & \checkmark & \checkmark  \\
Percentile of Wave Height & \textbf{$>$50 }    & \textbf{$>$75 }   & \textbf{$>$90}  \\
\textit{Pre SAR Period Statistics} &&& \\
Mean Total Attempt  & `c' & `c' & `c'  \\
Mean Frac. Unsafe Boat  & `avg_fr_123' & `avg_fr_123' & `avg_fr_123'  \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(onda_post_frac onda onda_frac onda_post) order(onda_post_frac onda onda_frac onda_post)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


*=============================================================================
* Table C9
*=============================================================================
* Note that you need to copy-paste the different models to the reproduce the latex file (tabC9.tex)


********************************************************************************
* FRAC MODEL  
********************************************************************************
use "data/data_tab_main2.dta", clear 
****
global toestout ""
global policies policy2 policy3 policy4
***

foreach outcome in fr_across3 fr_across_1_3 fr_across_1_2_3 fr_across4 fr_across5   {

*********
fracreg probit `outcome' policy2 policy3 policy4 i.week, nolog vce(cluster week)   
margins, dydx($policies) post 
estimates store reg_`outcome'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`outcome'"

}

   su fr_across3 if policy1 == 1
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su fr_across_1_3 if policy1 == 1
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su fr_across_1_2_3 if policy1 == 1
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

   su fr_across4 if policy1 == 1
local avg_fr_4 = r(mean)
local avg_fr_4 = round(`avg_fr_4',.01)

   su fr_across5 if policy1 == 1
local avg_fr_5 = r(mean)
local avg_fr_5 = round(`avg_fr_5',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC9a.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3cm}C{3cm}C{3cm}C{3cm}C{3cm}}
\toprule\toprule
  & (1)         & (2)       & (3)     & (4)      & (5)       \\
Panel A: Fraction of Attempted Crossings       & Inflatable   & Inflatable + & Inflatable + & Fishing & Motor                     \\
				           &    & Unknown & Unknown + &  &                      \\
				           &    &  & Other &  &                      \\
\hline
                         &        &       &          &     &                           \\
)
prefoot(
	\midrule
Pre MN Mean Outcome  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' & `avg_fr_4' & `avg_fr_5' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(policy2 policy3 policy4 ) 
order(policy2 policy3 policy4 )
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr
*************


********************************************************************************
* PPML 
********************************************************************************
use "data/data_tab_main2.dta", clear 
****
global toestout ""
global policies policy2 policy3 policy4
***

foreach outcome of varlist across3 across_1_3 across_1_2_3 across4 across5   {

*********
ppmlhdfe `outcome' $policies, ///
				  absorb(week) vce(cluster week) 
estimates store reg_`outcome'
estadd scalar nobs e(N)
global toestout "${toestout} reg_`outcome'"

}

   su across3 if policy1 == 1
local avg_fr_3 = r(mean)
local avg_fr_3 = round(`avg_fr_3',.01)

   su across_1_3 if policy1 == 1
local avg_fr_13 = r(mean)
local avg_fr_13 = round(`avg_fr_13',.01)

   su across_1_2_3 if policy1 == 1
local avg_fr_123 = r(mean)
local avg_fr_123 = round(`avg_fr_123',.01)

   su across4 if policy1 == 1
local avg_fr_4 = r(mean)
local avg_fr_4 = round(`avg_fr_4',.01)

   su across5 if policy1 == 1
local avg_fr_5 = r(mean)
local avg_fr_5 = round(`avg_fr_5',.01)


*** Export table
#delimit ;
estout ${toestout}
using "tables/tabC9b.tex",  
cells(b(fmt(3) star) se(par fmt(3))) style(tex) 
starlevels("" 0.00001) 
coll(,none) ml(,none)  
prehead(
\begin{tabular}{p{7cm}C{3cm}C{3cm}C{3cm}C{3cm}C{3cm}}
\toprule\toprule
  & (1)         & (2)       & (3)     & (4)      & (5)       \\
Panel B: Count of Attempted Crossings       & Inflatable   & Inflatable + & Inflatable + & Fishing & Motor                     \\
				           &    & Unknown & Unknown + &  &                      \\
				           &    &  & Other &  &                      \\
\hline
                         &        &       &          &     &                           \\
)
prefoot(
	\midrule
Pre MN Mean Outcome  & `avg_fr_3' & `avg_fr_13' & `avg_fr_123' & `avg_fr_4' & `avg_fr_5' \\
	)
postfoot(
	\bottomrule
	\end{tabular}
	)
keep(policy2 policy3 policy4) 
order(policy2 policy3 policy4)
stats(nobs, labels("Observations") fmt(0))
label replace
;
		
#delimit cr


*=============================================================================
* Figure C5
*=============================================================================
use "data/data_figC5.dta", clear

cap drop seasons
gen seasons = .
replace seasons = 1 if mese >=1 & mese <=3		// win
replace seasons = 2 if mese >=4 & mese <=6		// spr
replace seasons = 3 if mese >=7 & mese <=9 		// sum
replace seasons = 4 if mese >=10 & mese <=12 	// aut 

twoway ///
(kdensity swh_lyb_3 if seasons==1, lcolor(navy) lwidth(large) lpattern(solid)) ///
(kdensity swh_lyb_3 if seasons==2, lcolor(orange*0.5) lwidth(large) lpattern(longdash))	///
(kdensity swh_lyb_3 if seasons==3, lcolor(cranberry*.8) lwidth(large) lpattern(longdash_dot))	///
(kdensity swh_lyb_3 if seasons==4, lcolor(ebblue*0.5) lwidth(large) lpattern(dash))	///
					, ///
					ylab(, axis(1) tstyle(textstyle(size(small)))) ///
					xlab(, axis(1) tstyle(textstyle(size(small)))) ///
					xtitle("", size(small)) ///
					ytitle("Wave Height Density", size(small)) ///
legend(order(1 "Winter" 2 "Spring" 3 "Summer" 4 "Autumn") size(small) region(lwidth(none)) col(4) )
		gr_edit $whitebox
graph export "graphs/figC5.png", replace


*=============================================================================
* Figure C6
*=============================================================================
* No data analysis

*=============================================================================
* Figure C7
*=============================================================================
* No data analysis

*=============================================================================
* Figure C8
*=============================================================================
use "data/data_figC8.dta", clear

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_post_frac < -6.55 
cou if t == 1
local num `r(N)'
local pvalue = `num'/645
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------					
					
	sum be_onda_post_frac,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_post_frac,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_post_frac,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_post_frac,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))

	histogram be_onda_post_frac, ///
			  frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Post SAR * Frac. Unsafe  Boat) ///
			  xtitle(, size(small)) xlabel(-7(1)7, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(-6.55, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
	gr_edit $whitebox
	graph export "graphs/figC8b.png", replace
	
	drop p1m p1p p5m p5p
	********************************************************************************

	
*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda < -0.89
cou if t == 1
local num `r(N)'
local pvalue = `num'/645
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------	
	
	sum be_onda,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))

	histogram be_onda, ///
			  frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height) ///
			  xtitle(, size(small)) xlabel(-1(0.5)1, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(-0.89, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC8a.png", replace
	
	
		drop p1m p1p p5m p5p
	********************************************************************************

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_frac < 2.13
cou if t == 1
local num `r(N)'
local pvalue = `num'/645
local pv = round(`pvalue',.001)
*-------------------------------------------------------------------------------	


	sum be_onda_frac,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_frac,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_frac,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_frac,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))
	
	histogram be_onda_frac, ///
			  frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Frac. Unsafe  Boat) ///
			  xtitle(, size(small)) xlabel(-7(1)7, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(2.13, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC8d.png", replace
	


			drop p1m p1p p5m p5p
	********************************************************************************

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_post < 0.21
cou if t == 1
local num `r(N)'
local pvalue = `num'/645
local pv = round(`pvalue',.001)	
				
*-------------------------------------------------------------------------------	

	sum be_onda_post,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_post,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_post,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_post,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))
	
	local p5m .821
	histogram be_onda_post, ///
			  frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Post SAR) ///
			  xtitle(, size(small)) xlabel(-1(0.5)1, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(0.21, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC8c.png", replace

	
*-------------------------------------------------------------------------------



*=============================================================================
* Figure C9
*=============================================================================
use "data/data_figC9.dta", clear

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_post_frac < -6.55 
cou if t == 1
local num `r(N)'
local pvalue = `num'/78
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------
					
	sum be_onda_post_frac,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_post_frac,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_post_frac,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_post_frac,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))

	histogram be_onda_post_frac, ///
			  bin(14) frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Post SAR * Frac. Unsafe  Boat) ///
			  xtitle(, size(small)) xlabel(-7(1)7, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'.000", size(small)) ///
			  xline(-6.55, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC9b.png", replace
	
	drop p1m p1p p5m p5p
	********************************************************************************

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda < -0.89
cou if t == 1
local num `r(N)'
local pvalue = `num'/78
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------

	sum be_onda,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))
	
	local p5m .686
	histogram be_onda, ///
			  bin(14) frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height) ///
			  xtitle(, size(small)) xlabel(-1(0.5)1, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(-0.89, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC9a.png", replace
	
	
		drop p1m p1p p5m p5p
	********************************************************************************

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_frac < 2.13
cou if t == 1
local num `r(N)'
local pvalue = `num'/78
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------


	sum be_onda_frac,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_frac,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_frac,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_frac,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))
	
	histogram be_onda_frac, ///
			  bin(14) frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Frac. Unsafe  Boat) ///
			  xtitle(, size(small)) xlabel(-7(1)7, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(2.13, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC9d.png", replace
	


			drop p1m p1p p5m p5p
	********************************************************************************

*-------------------------------------------------------------------------------	
cap drop t	
gen t = be_onda_post < 0.21
cou if t == 1
local num `r(N)'
local pvalue = `num'/78
local pv = round(`pvalue',.001)					
*-------------------------------------------------------------------------------


	sum be_onda_post,det
	gen p1m  = 0-2.58*(`r(sd)')
	sum be_onda_post,det
	gen p1p  = 0+2.58*(`r(sd)')
	
	sum p1m,d
	local p1m = abs(round(p1m,.001))
	sum p1p,d
	local p1p = abs(round(p1p,.001))

	sum be_onda_post,det
	gen p5m  = 0-1.96*(`r(sd)')
	sum be_onda_post,det
	gen p5p  = 0+1.96*(`r(sd)')
	
	sum p5m,d
	local p5m = abs(round(p5m,.001))
	sum p5p,d
	local p5p = abs(round(p5p,.001))
	
	histogram be_onda_post, ///
			  bin(14) frequency fcolor(dknavy%40) lcolor(none%50) lpattern(solid) ///
			  ytitle(, size(small)) ///
			  ylabel(, labsize(small)) xtitle({&omega}{sub:s} Wave Height * Post SAR) ///
			  xtitle(, size(small)) xlabel(-1(0.5)1, labsize(small)) ///
			  note("Critical values at 1%, 5%: {&plusmn}`p1m', {&plusmn}`p5m'; Pvalue < `pv'", size(small)) ///
			  xline(0.21, lpattern(solid) lcolor(cranberry)) ///
			  xline(-`p1m' `p1p', lpattern(shortdash) lcolor(gray*.8)) ///
			  xline(-`p5m' `p5p', lpattern(longdash) lcolor(gray*.8)) ///
			  title(Title, color(white))
			  gr_edit $whitebox
	graph export "graphs/figC9c.png", replace



	
*-------------------------------------------------------------------------------



*=============================================================================
* Figure C10
*=============================================================================
use "data/data_figC10.dta", clear 

cap drop x 
gen x = fr_across3 //* 100
cap drop xx 
gen xx = swh_lyb_3 * postM 

cap drop onda 
gen onda 	  		= swh_lyb_3 
label var onda 	"Wave Height"
cap drop onda_frac 
gen onda_frac 		= swh_lyb_3 * fr_across3
label var onda_frac 	"Wave Height * Fr. Boat"
cap drop onda_post 
gen onda_post 		= swh_lyb_3 * postM
label var onda_post 	"Wave Height * Post SAR"
cap drop onda_post_frac 

xi: glm totacross onda onda_frac onda_post c.xx#c.x ///
			  i.weekanno, family(poisson) vce(hac nwest 28) t(data)
margins, eydx(xx) at(c.x=(0.1(.05)1)) 

forvalues i = 1/19 {
cap drop b_`i' 
g b_`i' = .
cap drop se_`i' 
g se_`i' = .	
}

forvalues i = 1/19 {
replace b_`i' = r(table)[1,`i']
replace se_`i' = r(table)[2,`i']	
}

keep if _n == 1
keep b_* se_* 
g seq = 1
reshape long se_ b_, i(seq) j(tipo) string
rename *_ *

cap drop ci_u
gen ci_u = b + 1.64*se

cap drop ci_l
gen ci_l = b - 1.64*se

cap drop ci_u2
gen ci_u2 = b + 1.96*se

cap drop ci_l2
gen ci_l2 = b - 1.96*se
drop seq 
destring tipo, replace 

gen tipo2 = .
replace tipo2 =  .1 if tipo == 1 
replace tipo2 = .15 if tipo == 2 
replace tipo2 =  .2 if tipo == 3 
replace tipo2 = .25 if tipo == 4 
replace tipo2 =  .3 if tipo == 5 
replace tipo2 = .35 if tipo == 6 
replace tipo2 =  .4 if tipo == 7 
replace tipo2 = .45 if tipo == 8 
replace tipo2 =  .5 if tipo == 9 
replace tipo2 = .55 if tipo == 10
replace tipo2 =  .6 if tipo == 11
replace tipo2 = .65 if tipo == 12
replace tipo2 =  .7 if tipo == 13
replace tipo2 = .75 if tipo == 14
replace tipo2 =  .8 if tipo == 15
replace tipo2 = .85 if tipo == 16
replace tipo2 =  .9 if tipo == 17
replace tipo2 = .95 if tipo == 18
replace tipo2 =   1 if tipo == 19


#delimit ;
twoway
	rcap ci_u2 ci_l2 tipo2, lcolor(dknavy*0.95) yaxis(1)	  ||
	scatter b tipo2, mcolor(dknavy) msize(small) msymbol(O) yaxis(1)	  
	
	yline(0, lcolor(black) lpattern(shortdash) lwidth(vthin))
	xlabel(0.1(0.05)1, labsize(small) notick angle(0) valuelabel)

	ylabel(, labsize(small) axis(1))
	ytitle(" ", size(small) axis(1))

	title(" ", size(small))
	legend(order(9 " " /*9 "Mean Angle: 41"*/))
	legend(size(small))
	legend(r(1) c(2))
	legen(region(color(white))) 
			
	xtitle("Fractions of Unsafe Boat", size(small) )
	xscale(titlegap(5))  
	yscale(titlegap(5))  

	graphregion(color(white))
	xsize(20) ysize(12)

;
#delimit cr
graph export "graphs/figC10.png", replace




*=============================================================================
* Figure C11
*=============================================================================
use "data/data_figC11.dta", clear
cap drop ci_u
gen ci_u = b + 1.96*se

cap drop ci_l
gen ci_l = b - 1.96*se

#delimit ;
twoway

	rcap ci_u ci_l my, lcolor(dknavy*0.75) ||
	scatter b my, mcolor(dknavy) msize(small) msymbol(O) 	

	yline(0, lcolor(black) lpattern(shortdash) lwidth(vthin))
	xline(20, lcolor(black) lpattern(shortdash) lwidth(vthin))
	xlabel(1(3)60, labsize(small) notick angle(0) valuelabel)

	xlab(1  "01/2016" ///
		  6 "06/2016" 12 "12/2016" ///
		 18 "06/2017" 24 "12/2017" ///
		 30 "06/2018" 36 "12/2018" ///
		 42 "06/2019" 48 "12/2019" ///
		 54 "06/2020" 60 "12/2020" ///
		 , ///
	tstyle(textstyle(size(small))))	///
	ylabel(, labsize(small) axis(1))
	
	title(, size(small))

	ttext(5 24 "Minniti Code", ///
	size(small) color(black)) ///

	legend(order(2))
	legend(label(2 "Month * {bf:1}[From Libya]") size(small))
	legend(r(1))
	legen(region(color(white))) 
			
	xtitle(" ", size(small) )
	xscale(titlegap(5))  
	yscale(titlegap(5))  

	graphregion(color(white))
	xsize(20) ysize(12)

;
#delimit cr
graph export "graphs/figC11.png", replace




*-------------------------------------------------------------------------------
* 								APPENDIX D 		
*-------------------------------------------------------------------------------




*=============================================================================
* Table D1
*=============================================================================
* No data analysis

*=============================================================================
* Figure D1
*=============================================================================
use "data/data_figD1.dta", clear 

twoway 	(connected pc_attempt1 modate,		 ///
lcolor(navy*.75) lpattern(solid) msymbol(O) mcolor(navy) msize(small) lwidth(small) yaxis(1)) ///
		(connected pc_attempt2 modate,		 ///
lcolor(orange*.75) lpattern(longdash) msymbol(Th) mcolor(orange) msize(small) lwidth(small) yaxis(1)) ///
		(connected pc_attempt3 modate,		 ///
lcolor(cranberry*.75) lpattern(dash) msymbol(D) mcolor(cranberry) msize(small) lwidth(small) yaxis(1)) ///
		(connected pc_attempt4 modate,		 ///
lcolor(midgreen*.75) lpattern(shortdash) msymbol(Sh) mcolor(midgreen) msize(small) lwidth(small) yaxis(1)) ///
					, ///
					ylab(, axis(1) tstyle(textstyle(size(small)))) ///
					xlab(, axis(1) tstyle(textstyle(size(small)))) ///
					xtitle("", size(small)) ///
legend(order(1 "EU Coast Patrol" 2 "Maritime Force" 3 "NGO" 4 "Commercial") size(small) region(lwidth(none)) col(4) )
					gr_edit $whitebox 
graph export "graphs/figD1.png", replace

