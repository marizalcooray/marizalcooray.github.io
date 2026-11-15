/*==============================================================================
  Oil price pass-through to consumer prices: Tonga, Fiji and Samoa
  ------------------------------------------------------------------------------
  Author   : Mariza Cooray
  Date     : November 2026

  Purpose  : Estimate how movements in an oil price index pass through to the
             CPI divisions of three Pacific island economies. The estimates
             parameterise a simulation model of the 2026 oil price shock.

             1. Parameters
             2. Programs
             3. Model specifications
             4. Estimation
                  a. short-run (log-log) equations for the fuel-exposed
                     divisions, each with an out-of-sample check
                  b. transport: test for asymmetric pass-through
                  c. transport: bounds test for a long-run relationship
                  d. long-run pass-through for all 12 divisions
             5. Results table and figure

  Inputs   : Stata files in data/ (names in section 1). Monthly, one row per
             month in chronological order, no date variable. Variables:
               food ... misc   the 12 CPI divisions
               oilind          oil price index
               firstr          aggregate first-round effect
               xrate           exchange rate (Fiji and Samoa)

  Outputs  : output/oil_passthrough.log
             output/passthrough_estimates.dta, .csv   every estimate, long form
             output/lr_passthrough_table.csv          division by country
             output/lr_passthrough.png

  Requires : Stata 16 or later, and ardl (installed below if missing)
  To run   : from the project folder,  do oil_passthrough.do
             Run the file in one go: the sections share local macros.
==============================================================================*/

version 16
clear all
set more off

capture mkdir "output"
capture log close
log using "output/oil_passthrough.log", replace text

capture which ardl
if _rc ssc install ardl


**# 1. Parameters ---------------------------------------------------------------

local countries Tonga Fiji Samoa
local divs      food alc cloth hous fur health trans comm recr edu rest misc
local fuel_divs trans hous food rest    // divisions with a direct fuel or freight channel

local nwlag     12          // Newey-West truncation lag (months)
local holdout   2023m12     // last month of the estimation window in the backtest
local nardl_p   4           // own lags in the asymmetry test
local nardl_q   4           // oil price lags in the asymmetry test

* Data files and their first month: sr = short-run block, lr = long-run block
local file_sr_Tonga  oil_price_ind_CPI
local t0_sr_Tonga    2021m1
local file_lr_Tonga  oil_price_ind_inflation
local t0_lr_Tonga    2022m1

local file_sr_Fiji   oil_price_ind_CPI_Fiji
local t0_sr_Fiji     2019m1
local file_lr_Fiji   oil_price_ind_CPI_Fiji
local t0_lr_Fiji     2019m1

local file_sr_Samoa  oil_price_ind_CPI_Samoa
local t0_sr_Samoa    2020m1
local file_lr_Samoa  oil_price_ind_CPI_Samoa
local t0_lr_Samoa    2020m1


**# 2. Programs -----------------------------------------------------------------

* Load one country file and declare the monthly time series
capture program drop load_data
program define load_data
    args file first_month
    use "data/`file'.dta", clear
    generate int mdate = tm(`first_month') + _n - 1
    format mdate %tm
    tsset mdate
end

* Estimate
*     y_t = a + sum_j b_j x_{t-j} + sum_k c_k y_{t-k} [+ controls] + e_t
* by OLS with Newey-West standard errors and return
*     impact   = sum_j b_j
*     longrun  = sum_j b_j / (1 - sum_k c_k), with a delta-method standard error
*     rho      = sum_k c_k (the long-run multiplier is finite only if rho < 1)
*     rmse_oos = one-step-ahead RMSE after holdout(), if requested
capture program drop passthrough
program define passthrough, rclass
    syntax varname(numeric), xvar(varname numeric) xlags(numlist integer >=0) ///
        ylags(numlist integer >0)                                             ///
        [controls(string) noCONStant nwlag(integer 12) holdout(string)]

    local y `varlist'
    local xterms
    local yterms
    local num 0
    local den 1
    foreach j of numlist `xlags' {
        local xterms `xterms' L`j'.`xvar'
        local num    `num' + _b[L`j'.`xvar']
    }
    foreach k of numlist `ylags' {
        local yterms `yterms' L`k'.`y'
        local den    `den' - _b[L`k'.`y']
    }
    local rhs `xterms' `yterms' `controls'

    * Backtest: estimate up to the holdout month, predict one step ahead after it
    if "`holdout'" != "" {
        quietly tsset
        local t `r(timevar)'
        tempvar yhat err2
        quietly newey `y' `rhs' if `t' <= `holdout', lag(`nwlag') `constant'
        quietly predict double `yhat' if `t' > `holdout'
        quietly generate double `err2' = (`y' - `yhat')^2
        quietly summarize `err2'
        return scalar rmse_oos = sqrt(r(mean))
    }

    * Full-sample estimate, impact effect and long-run multiplier
    quietly newey `y' `rhs', lag(`nwlag') `constant'
    return scalar N   = e(N)
    return scalar rho = 1 - (`den')
    quietly nlcom (impact: `num') (longrun: (`num') / (`den'))
    tempname b V
    matrix `b' = r(b)
    matrix `V' = r(V)
    return scalar impact  = `b'[1,1]
    return scalar longrun = `b'[1,2]
    return scalar se_long = sqrt(`V'[2,2])
end


**# 3. Model specifications -----------------------------------------------------

* One specification per block (sr, lr), country and division:
*     x_*   lags of the oil price
*     y_*   lags of the dependent variable
*     o_*   further options: controls() and noconstant
* The default is the contemporaneous oil price and one own lag. The departures
* below are the lags retained after general-to-specific reduction of a model
* with oil lags 0-12 and own lags 1-6 (1-12 in the short-run block).

foreach b in sr lr {
    foreach c of local countries {
        foreach v of local divs {
            local x_`b'_`c'_`v' 0
            local y_`b'_`c'_`v' 1
            local o_`b'_`c'_`v'
        }
    }
}

* Short-run block (logs). Food and restaurants also respond to the aggregate
* first-round effect; transport in Fiji and Samoa to the exchange rate.
local y_sr_Tonga_food   1 6
local o_sr_Tonga_food   controls(L6.lfirstr)
local x_sr_Tonga_rest   6
local o_sr_Tonga_rest   controls(L1.lfirstr) noconstant

local x_sr_Fiji_trans   0 1
local o_sr_Fiji_trans   controls(L5.lxrate)
local y_sr_Fiji_food    1 6
local o_sr_Fiji_food    controls(L6.lfirstr)
local x_sr_Fiji_rest    6
local o_sr_Fiji_rest    controls(L1.lfirstr) noconstant

local x_sr_Samoa_trans  1
local o_sr_Samoa_trans  controls(L5.lxrate)
local o_sr_Samoa_food   controls(L2.lfirstr)
local o_sr_Samoa_rest   controls(L0.lfirstr)

* Long-run block
local x_lr_Tonga_food   0 1 4
local y_lr_Tonga_food   1 3 4
local x_lr_Tonga_alc    1
local y_lr_Tonga_cloth  1 2
local x_lr_Tonga_hous   0 4
local y_lr_Tonga_fur    1 2
local x_lr_Tonga_health 0 1 2
local y_lr_Tonga_trans  1 2
local x_lr_Tonga_rest   2 6
local x_lr_Tonga_misc   0 5

local y_lr_Fiji_food    1 3
local x_lr_Fiji_alc     0 1
local o_lr_Fiji_alc     noconstant
local y_lr_Fiji_cloth   1 2
local x_lr_Fiji_comm    1
local x_lr_Fiji_edu     0 1 2 4 5
local o_lr_Fiji_edu     noconstant
local x_lr_Fiji_rest    1 2
local y_lr_Fiji_rest    1 2
local o_lr_Fiji_rest    noconstant

local x_lr_Samoa_food   0 1 2
local o_lr_Samoa_food   noconstant
local x_lr_Samoa_alc    1
local o_lr_Samoa_alc    noconstant
local y_lr_Samoa_cloth  1 2
local x_lr_Samoa_health 1
local x_lr_Samoa_comm   1
local x_lr_Samoa_edu    0 1 2 4 5
local o_lr_Samoa_edu    noconstant
local x_lr_Samoa_rest   1 2
local y_lr_Samoa_rest   1 2
local o_lr_Samoa_rest   noconstant
local x_lr_Samoa_misc   1
local o_lr_Samoa_misc   noconstant


**# 4. Estimation ---------------------------------------------------------------

tempname est
postfile `est' str2 block str5 country byte div str6 code       ///
    double(impact longrun se_long rho rmse_oos) int nobs        ///
    using "output/passthrough_estimates.dta", replace

local ncountries : word count `countries'
matrix ASYM = J(`ncountries', 2, .)
matrix rownames ASYM = `countries'
matrix colnames ASYM = F p
local row 0

foreach c of local countries {
    local ++row
    display as text _n "{hline 78}" _n "`c'" _n "{hline 78}"

    **# 4a. Short-run elasticities, fuel-exposed divisions
    load_data "`file_sr_`c''" `t0_sr_`c''
    confirm numeric variable `fuel_divs' oilind firstr
    foreach v of varlist `fuel_divs' oilind firstr {
        assert `v' > 0 if !missing(`v')
        generate double l`v' = ln(`v')
    }
    capture confirm numeric variable xrate
    if !_rc generate double lxrate = ln(xrate)

    foreach v of local fuel_divs {
        local i : list posof "`v'" in divs
        passthrough l`v', xvar(loilind) xlags(`x_sr_`c'_`v'')       ///
            ylags(`y_sr_`c'_`v'') nwlag(`nwlag')                    ///
            holdout(tm(`holdout')) `o_sr_`c'_`v''
        post `est' ("sr") ("`c'") (`i') ("`v'") (r(impact)) (r(longrun)) ///
            (r(se_long)) (r(rho)) (r(rmse_oos)) (r(N))
    }

    **# 4b. Transport: is pass-through symmetric?
    * NARDL decomposition: cumulate oil price increases and decreases separately
    generate double d_oil   = D.loilind
    generate double oil_pos = sum(max(d_oil, 0))
    generate double oil_neg = sum(min(d_oil, 0))
    quietly newey ltrans L(0/`nardl_q').oil_pos L(0/`nardl_q').oil_neg ///
        L(1/`nardl_p').ltrans, lag(`nwlag')

    * Both long-run multipliers share a denominator, so they are equal if and
    * only if the summed coefficients are equal
    local pos _b[oil_pos]
    local neg _b[oil_neg]
    forvalues j = 1/`nardl_q' {
        local pos `pos' + _b[L`j'.oil_pos]
        local neg `neg' + _b[L`j'.oil_neg]
    }
    test `pos' = `neg'
    matrix ASYM[`row', 1] = r(F)
    matrix ASYM[`row', 2] = r(p)

    **# 4c. Transport: bounds test for a long-run (levels) relationship
    ardl ltrans loilind, maxlags(4) bic ec
    estat ectest

    **# 4d. Long-run pass-through, all divisions
    load_data "`file_lr_`c''" `t0_lr_`c''
    confirm numeric variable `divs' oilind
    foreach v of local divs {
        local i : list posof "`v'" in divs
        passthrough `v', xvar(oilind) xlags(`x_lr_`c'_`v'')         ///
            ylags(`y_lr_`c'_`v'') nwlag(`nwlag') `o_lr_`c'_`v''
        local rho = r(rho)
        post `est' ("lr") ("`c'") (`i') ("`v'") (r(impact)) (r(longrun)) ///
            (r(se_long)) (r(rho)) (.) (r(N))
        if `rho' >= 1 {
            display as error "`c' `v': own-lag coefficients sum to " ///
                %5.3f `rho' ", so the long-run multiplier is not finite"
        }
    }
}
postclose `est'


**# 5. Results table and figure -------------------------------------------------

use "output/passthrough_estimates.dta", clear

label define div  1 "Food and non-alcoholic beverages"  2 "Alcohol, tobacco and kava" ///
                  3 "Clothing and footwear"             4 "Housing and utilities"     ///
                  5 "Furnishings and maintenance"       6 "Health"                    ///
                  7 "Transport"                         8 "Communication"             ///
                  9 "Recreation and culture"           10 "Education"                 ///
                 11 "Restaurants and hotels"           12 "Miscellaneous"
label values div div

label variable block    "sr = short-run (logs), lr = long-run"
label variable impact   "Sum of oil price coefficients"
label variable longrun  "Long-run pass-through"
label variable se_long  "Delta-method standard error of longrun"
label variable rho      "Sum of own-lag coefficients"
label variable rmse_oos "One-step-ahead RMSE after `holdout'"
label variable nobs     "Observations"

generate double lr_lo = longrun - invnormal(0.975) * se_long
generate double lr_hi = longrun + invnormal(0.975) * se_long
format impact longrun se_long lr_lo lr_hi rho rmse_oos %9.3f

list block country div impact longrun se_long rho rmse_oos nobs, ///
    sepby(block country) noobs abbreviate(12)
matrix list ASYM, format(%9.3f) ///
    title("Transport: Wald test of symmetric long-run pass-through")

save "output/passthrough_estimates.dta", replace
export delimited using "output/passthrough_estimates.csv", replace

* Figure: long-run pass-through by division, one panel per country
graph hbar (asis) longrun if block == "lr", over(div, label(labsize(small))) ///
    by(country, rows(1) yrescale note("")                                     ///
       title("Long-run oil price pass-through by CPI division"))              ///
    ytitle("Long-run pass-through coefficient") yline(0, lcolor(gs8))
graph export "output/lr_passthrough.png", replace width(2400)

* Table: divisions in rows, countries in columns
keep if block == "lr"
keep country div longrun
reshape wide longrun, i(div) j(country) string
rename longrun* *
export delimited using "output/lr_passthrough_table.csv", replace

log close
