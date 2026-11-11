/*=============================================================================
       Pass-through estimates for oil shocks' simulation model 
					     F I J I
					 v. July 4, 2026
=============================================================================*/

use "oil_price_ind_CPI_Fiji.dta", clear
gen mdate = tm(2019m1) + _n - 1
format mdate %tm
tsset mdate

* Log variables
gen lfood   = ln(food)
gen lalc    = ln(alc)
gen lcloth  = ln(cloth)
gen lhous   = ln(hous)
gen lfur	= ln(fur)
gen lhealth = ln(health)
gen ltrans  = ln(trans)
gen lcomm   = ln(comm)
gen lrecr  	= ln(recr)
gen ledu    = ln(edu)
gen lrest   = ln(rest)
gen lmisc   = ln(misc)
gen loilind = ln(oilind)
gen lfao    = ln(fao)
gen lxrate    = ln(xrate)

* First round aggregated effects: 
gen lfirstr = ln(firstr)

* Decomposing oil price changes into positive and negative
gen doilind = D.oilind
gen doilind_pos = max(doilind, 0)   // price increases
gen doilind_neg = min(doilind, 0)   // price decreases

* Constructing cumulative partial sums (NARDL approach)
* positive and negative cumulative oil price changes
gen oilind_pos = 0
gen oilind_neg = 0
replace oilind_pos = sum(max(doilind, 0))
replace oilind_neg = sum(min(doilind, 0))


**# (1) Pass-through estimates: oil prices to consumer prices, short-run effects
**# ------------- Transport -----------------------------
	varsoc ltrans loilind lxrate, maxlag(12) 
	*  optimal lags: 3
	newey ltrans L(0/6).loilind L(1/6).ltrans L(0/6).lxrate, lag(12)
	newey ltrans L0.loilind L1.loilind ///
				 L1.ltrans ///
				 L5.lxrate, lag(12)	
	matrix b = e(b)'
	matrix list b
	
	
			* ~~~~~~~ Final model selection ~~~~~~~
			quiet reg ltrans L0.loilind L1.loilind L1.ltrans L5.xrate				 
			estat ic
			quiet reg ltrans L0.loilind L1.ltrans L5.xrate				 
			estat ic
	
			* ~~~~~~~ Backtesting forecasting accuracy ~~~~~~~
			* Estimation using the 2021-2023 sample only (using the best model)
			newey ltrans L0.loilind L1.loilind ///
				 L1.ltrans ///
				 L5.xrate  ///
				 if mdate <= tm(2023m12), lag(12)
			* Predicted values for 2024-2025
			predict ltrans_hat if mdate > tm(2023m12)
			* Comparing to actuals
			gen forecast_error_trans = ltrans - ltrans_hat
			sum forecast_error_trans
			
			* ~~~~~~~ Testing asymmetric pass-through effects ~~~~~~~	
			* Approach 1: using ARDL
			newey ltrans L(0/3).doilind_pos L(0/3).doilind_neg, lag(12)
			test _b[doilind_pos] = _b[doilind_neg]
			* Null hypothesis: no asymmetric effect of shocks (no asymmetric pass-through effects)
	
			* ~~~~~~~ Cointegration Testing ~~~~~~~
			varsoc ltrans loilind, maxlag(12) 
			ardl ltrans loilind lxrate, lags(3) maxlags(4) ec bic
			estat ectest
			estat btest    // Bounds test for cointegration
			* Null: no cointegration
	
**# ------------- Electricity, housing
	newey lhous L(0/12).loilind L(1/12).lhous, lag(12)
	newey lhous L0.loilind L1.lhous, lag(12)
	matrix b = e(b)'
	matrix list b

**# ------------- Food	
	newey lfood L(0/12).loilind L(1/12).lfood L(0/12).lfirstr L(0/12).lfao, lag(12)
	newey lfood L0.loilind L4.loilind L1.lfood L6.lfood L6.lfirstr, lag(12)
	newey lfood L0.loilind L1.lfood L6.lfood L6.lfirstr, lag(12)	
	matrix b = e(b)'
	matrix list b
	
		* ~~~~~~~ Backtesting forecasting accuracy ~~~~~~~
		* Estimation using the 2021-2023 sample only
		newey lfood L0.loilind L1.lfood L6.lfood L6.lfirstr ///
		  if mdate <= tm(2023m12), lag(12)
		* Predicted values for 2024-2025
		predict lfood_hat if mdate > tm(2023m12)
		* Comparing to actuals
		gen forecast_error_food = lfood - lfood_hat
		sum forecast_error_food
		
		* ~~~~~~~ Testing asymmetric pass-through effects ~~~~~~~	
		* Approach 1: using ARDL
		newey lfood L(0/4).doilind_pos L(0/4).doilind_neg L(1/4).food, lag(12)
		test _b[doilind_pos] = _b[doilind_neg]
		* Null hypothesis: no asymmetric effect of shocks (no asymmetric pass-through effects)
		* Approach 2: NARDL
		* Estimation in levels (this preserves the long-run relationship)
		newey lfood L(0/4).oilind_pos L(0/4).oilind_neg ///
					L(1/4).lfood L(0/2).lfirstr, lag(12)
		* Long-run asymmetry test:
		* LR+ = sum of pos coefficients / (1 - sum of AR coefficients)
		* LR- = sum of neg coefficients / (1 - sum of AR coefficients)
		nlcom (_b[oilind_pos] + _b[L1.oilind_pos] + _b[L2.oilind_pos]) / ///
			  (1 - _b[L1.lfood] - _b[L2.lfood] - _b[L3.lfood] - _b[L4.lfood])
		nlcom (_b[oilind_neg] + _b[L1.oilind_neg] + _b[L2.oilind_neg]) / ///
			  (1 - _b[L1.lfood] - _b[L2.lfood] - _b[L3.lfood] - _b[L4.lfood])
			  
**# ------------- Restaurants
	newey lrest L(0/12).loilind L(1/12).lrest L(0/12).lfirstr, lag(12)
	newey lrest L6.loilind L1.lrest L1.lfirstr, lag(12) nocons
	matrix b = e(b)'
	matrix list b

* -----------------------------------------------------------------------------
**# (2) Pass-through estimates: oil prices to inflation, long-run effects	
use "oil_price_ind_CPI_Fiji.dta", clear
gen mdate = tm(2021m1) + _n - 1
format mdate %tm
tsset mdate

* Matrix to store long-run estimates
matrix LR = J(12, 1, .)
matrix colnames LR = b
matrix rownames LR = food alc cloth hous fur health trans comm recr edu rest misc
local r = 1
	
**# Long-rung pass-through estimates: channels	
**# ------------- Food and non-alcoholic beverages
newey food L(0/6).oilind L(1/6).food, lag(12)
newey food L(0/4).oilind L(1/4).food, lag(12)
newey food L0.oilind L1.food L3.food, lag(12)
	nlcom (_b[oilind]) / ///
      (1 - _b[L1.food] - _b[L3.food])
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	
**# ------------- Alcoholic beverages, tobacco and kava
newey alc L(0/12).oilind L(1/6).alc, lag(12) 
newey alc L0.oilind L1.oilind L1.alc, lag(12) nocons
	nlcom (_b[L0.oilind] + _b[L1.oilind]) / ///
      (1 - _b[L1.alc])	  
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Clothing and footwear
newey cloth L(0/12).oilind L(1/6).cloth, lag(12)
newey cloth L0.oilind L(1/2).cloth, lag(12)
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.cloth] + _b[L2.cloth])
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Housing, water, electricity and gas
newey hous L(0/12).oilind L(1/6).hous, lag(12)
newey hous L0.oilind L1.hous, lag(12)
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.hous])
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Furnishings, household equipment and household maintenance
newey fur L(0/12).oilind L(1/6).fur, lag(12)	  
newey fur L0.oilind L1.fur, lag(12)
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.fur])
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	
**# ------------- Health	  
newey health L(0/12).oilind L(1/6).health, lag(12)	
newey health L0.oilind L1.health, lag(12) 
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.health])	  
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Transport
newey trans L(0/12).oilind L(1/6).trans, lag(12)
newey trans L0.oilind L1.trans, lag(12)
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.trans])
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	 
**# ------------- Communication
newey comm L(0/12).oilind L(1/6).comm, lag(12)
newey comm L1.oilind L1.comm, lag(12)	  
	nlcom (_b[L1.oilind]) / ///
      (1 - _b[L1.comm])	  
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	
**# ------------- Recreation and culture
newey recr L(0/12).oilind L(1/6).recr, lag(12)
newey recr L0.oilind L1.recr, lag(12)	  
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.recr])	  
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  	  
**# ------------- Education
newey edu L(0/12).oilind L(1/6).edu, lag(12)
newey edu L0.oilind L1.oilind L2.oilind L4.oilind L5.oilind L1.edu, lag(12) nocons
	nlcom (_b[L0.oilind]+_b[L1.oilind] +_b[L2.oilind] +_b[L4.oilind] +_b[L5.oilind] ) / ///
      (1 - _b[L1.edu])	  	 
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Restaurants and hotels
newey rest L(0/12).oilind L(1/6).rest, lag(12)
newey rest L1.oilind L2.oilind L1.rest L2.rest, lag(12)	nocons  
	nlcom (_b[L1.oilind] + _b[L2.oilind]) / ///
      (1 - _b[L1.rest] - _b[L2.rest])	
	matrix LR[`r',1] = r(b)[1,1]
	local ++r
	  
**# ------------- Miscellaneous goods and services
newey misc L(0/12).oilind L(1/6).misc, lag(12)
newey edu L0.oilind L1.misc, lag(12)	  
	nlcom (_b[L0.oilind]) / ///
      (1 - _b[L1.misc])	  	 
	matrix LR[`r',1] = r(b)[1,1]
	local ++r

matrix list LR

* -----------------------------------------------------------------------------	  
**# (3) Trend component (weights of the Bates-Granger combination)
gen trend = _n
reg firstr trend

* ---------------------- end of script ---------------------------------------