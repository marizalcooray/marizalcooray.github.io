/*==============================================================================
  Returns to schooling in Papua New Guinea
  ------------------------------------------------------------------------------
  Author   : Mariza Cooray

  Purpose  : Build a person-level file from the Household Income and Expenditure
             Survey (HIES) and estimate Mincer earnings equations.

             1. Parameters
             2. Build the person-level file from five HIES modules
             3. Construct variables
             4. Analysis
                  a. Mincer equation
                  b. earnings by qualification
                  c. occupation skill level by qualification
                  d. gender earnings gap, overall and by qualification
                  e. consumption relative to gross earnings
                  f. gender earnings gap by age group

  Inputs   : HIES modules in data/ (the microdata are not distributed here)
               HIES_HA1     roster: sex and age
               HIES_HB1     education: highest grade or qualification
               HIES_PC1A1   income, by pay period
               HIES_PC1B    gross earnings and consumption
               HIES_PC2     occupation
             A person is identified by x_psu x_hnum x_pnum.

  Outputs  : output/mincer_png.log
             output/fig1 to fig5 (.png)
             data/mincer_analysis_file.dta

  Requires : Stata 16 or later, and coefplot (installed below if missing)
  To run   : from the project folder,  do Mincer_analysis_PNG.do
==============================================================================*/

version 16
clear all
set more off

capture mkdir "output"
capture log close
log using "output/mincer_png.log", replace text

capture which coefplot
if _rc ssc install coefplot


**# 1. Parameters ---------------------------------------------------------------

local id          x_psu x_hnum x_pnum    // person identifier
local days_month  30                     // days per month, to convert pay periods
local school_age  6                      // age at school entry
local ref         2                      // reference qualification: Grade 12


**# 2. Build the person-level file ----------------------------------------------

tempfile roster income schooling occupation business

* Roster: sex and age
use "data/HIES_HA1", clear
rename a1_pid x_pnum
drop if missing(x_pnum)
keep `id' a1_03 a1_04
rename (a1_03 a1_04) (sex_code age)
save `roster'

* Income: convert each record to a monthly amount, then sum within person
use "data/HIES_PC1A1", clear
generate double monthly_inc = c1a15              if c1a16 == 1    // paid monthly
replace monthly_inc = c1a15 * `days_month' / 14  if c1a16 == 2    // fortnightly
replace monthly_inc = c1a15 * `days_month' / 7   if c1a16 == 3    // weekly
replace monthly_inc = c1a15 * `days_month'       if c1a16 == 4    // daily
collapse (sum) total_inc = monthly_inc, by(`id')
save `income'

* Education: highest grade or qualification
use "data/HIES_HB1", clear
rename b1_pid x_pnum
drop if missing(x_pnum)
keep `id' b1_09
rename b1_09 edu_code
save `schooling'

* Occupation
use "data/HIES_PC2", clear
keep `id' c2_06
rename c2_06 occ_code
save `occupation'

* Gross earnings and consumption, summed within person
use "data/HIES_PC1B", clear
drop if missing(c1b17)
collapse (sum) gross_earn = c1b17 consumption = c1b16, by(`id')
save `business'

* One row per person
use `roster', clear
foreach f in income schooling occupation business {
    merge 1:1 `id' using ``f'', nogenerate
}


**# 3. Construct variables ------------------------------------------------------

drop if age == 0

* Years of education. Grades count one year each; post-school qualifications
* map to 13 (trade or other certificate), 14 (diploma), 15 (bachelor's degree)
* and 16 (postgraduate). Code 98 is set to missing.
recode edu_code (14 = 13) (15 = 14) (16 = 15) (17 = 16) (98 = .), generate(years_edu)

* Highest qualification, for those with Grade 10 or above
generate byte education = edu_code - 10 if inrange(edu_code, 10, 17)
label define edulabel 0 "Grade 10" 1 "Grade 11" 2 "Grade 12" 3 "Trade Certificate" ///
    4 "Other Certificate" 5 "Diploma" 6 "Bachelor's degree" 7 "Postgraduate"
label values education edulabel

* Potential experience and income
generate exper      = age - years_edu - `school_age' if !missing(total_inc)
generate exper2     = exper^2
generate ln_inc     = ln(total_inc) if total_inc > 0
generate annual_inc = 12 * total_inc

* Sex
generate byte female = sex_code == 2 if inlist(sex_code, 1, 2)
label define female 0 "Male" 1 "Female"
label values female female

* Occupation skill level
generate byte occupation = .
replace occupation = 2 if inlist(occ_code, 11, 12, 13, 21, 22, 23, 24, 31, 32, 33, 34)
replace occupation = 1 if inlist(occ_code, 41, 42, 51, 52, 61, 71, 72, 73, 74, 81, 82)
replace occupation = 0 if inlist(occ_code, 62, 83, 91, 92, 93)
label define skill 0 "Low skilled" 1 "Semi-skilled" 2 "Highly skilled"
label values occupation skill

* Consumption relative to gross earnings, per cent
generate rel_cons = 100 * consumption / gross_earn

* Ten-year age groups
recode age (20/29 = 0 "20-29") (30/39 = 1 "30-39") (40/49 = 2 "40-49") ///
           (50/59 = 3 "50-59") (60/69 = 4 "60-69") (else = .), generate(agegroup)

label variable years_edu  "Years of education"
label variable education  "Highest qualification"
label variable exper      "Potential experience (years)"
label variable exper2     "Potential experience squared"
label variable ln_inc     "Log monthly income (kina)"
label variable annual_inc "Annual income (kina)"
label variable female     "Female"
label variable occupation "Occupation skill level"
label variable rel_cons   "Consumption relative to gross earnings (%)"
label variable agegroup   "Age group"

compress
save "data/mincer_analysis_file.dta", replace


**# 4. Analysis -----------------------------------------------------------------

**# 4a. Mincer equation
regress ln_inc years_edu exper exper2

coefplot, drop(_cons)                                                           ///
    coeflabel(years_edu = "Years of education" exper = "Experience"             ///
              exper2 = "Experience squared")                                    ///
    xline(0, lcolor(red)) xlabel(-0.15(0.05)0.15) mcolor(green)                 ///
    xtitle("{bf:Coefficient}") ytitle("{bf:Characteristics}")                   ///
    title("{bf:Wage returns}" "{bf:for additional years of education}")         ///
    graphregion(color(white))
graph export "output/fig1_mincer.png", replace width(2000)

**# 4b. Earnings by qualification, relative to the reference qualification
regress ln_inc ib`ref'.education exper exper2

* Experience-adjusted income by qualification: mean income of the reference
* group, scaled by each qualification's estimated premium exp(b)
summarize annual_inc if education == `ref' & e(sample), meanonly
local ref_mean = r(mean)
generate double adj_inc = `ref_mean' if education == `ref'
levelsof education if e(sample) & education != `ref', local(levels)
foreach k of local levels {
    replace adj_inc = `ref_mean' * exp(_b[`k'.education]) if education == `k'
}
label variable adj_inc "Experience-adjusted annual income (kina)"

graph hbar (mean) adj_inc, over(education) blabel(bar, format(%9.0fc))          ///
    ytitle("{bf:Annual income in kina}")                                        ///
    title("{bf:Average earnings by qualification}") graphregion(color(white))
graph export "output/fig2_earnings_by_qualification.png", replace width(2000)

**# 4c. Occupation skill level by qualification
mlogit occupation ib`ref'.education, baseoutcome(0)
margins education, predict(outcome(0)) predict(outcome(1)) predict(outcome(2))

**# 4d. Gender earnings gap, overall and within each qualification
regress ln_inc female ib`ref'.education exper exper2

local symbols square plus X circle triangle diamond smcircle smdiamond
local plots
forvalues k = 0/7 {
    quietly regress ln_inc female exper exper2 if education == `k'
    estimates store edu`k'
    local s : word `=`k' + 1' of `symbols'
    local plots `plots' (edu`k', msymbol(`s'))
}
coefplot `plots', asequation swapnames keep(female)                             ///
    coeflabel(edu0 = "Grade 10" edu1 = "Grade 11" edu2 = "Grade 12"             ///
              edu3 = "Trade Certificate" edu4 = "Other Certificate"             ///
              edu5 = "Diploma" edu6 = "Bachelor's degree"                       ///
              edu7 = "Postgraduate", notick)                                    ///
    xline(0, lcolor(red) lpattern(solid)) levels(95) legend(off)                ///
    xtitle("{bf:Coefficient on female}") ytitle("{bf:Qualification}")           ///
    title("{bf:The Influence of Sex on Wages}" "{bf:By Education}")             ///
    graphregion(color(white))
graph export "output/fig3_gender_gap_by_qualification.png", replace width(2000)

**# 4e. Consumption relative to gross earnings
regress rel_cons ib`ref'.education

local ref_label : label edulabel `ref'
coefplot, drop(_cons) baselevels                                                ///
    coeflabel(`ref'.education = "`ref_label' {it:(Reference)}")                 ///
    xline(0, lcolor(red)) xlabel(-300(50)300) mcolor(green)                     ///
    xtitle("{bf:Coefficient}") ytitle("{bf:Type of educational qualification}") ///
    title("{bf:The Influence of Education}" "{bf:on Relative Consumption}")     ///
    graphregion(color(white))
graph export "output/fig4_relative_consumption.png", replace width(2000)

**# 4f. Gender earnings gap by age group
regress ln_inc i.female##i.agegroup
margins female#agegroup
marginsplot, xdimension(agegroup) title("{bf:Wage differential}")
graph export "output/fig5_gender_gap_by_age.png", replace width(2000)

log close
