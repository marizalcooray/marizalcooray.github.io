/*******************************************************************************
DFM for Lao PDR
Version: October, 2024
Working code
*******************************************************************************/
clear all
* Set the working directory to the folder containing this file and the data

*-------------------------------------------------------------------------------
**# Declaring globals
	glo path_in "."
	glo var_list FDI EDLT_CH EDLT DS RPD_CH RPD XRUSD XRCNY XRTHB INF CAB GDPG_CH GDPG CAB_CH CAB_TH covid19
	glo var_lbop CA	CAPI CASI KA FA DI FADIL FAPI NEO OB
	glo var_lqrt FA CAGE CAGI CASH INF INFCORE DEPCAP DEPCLEAR DEPENP DEPOSITS DEPRES FX FXRESUS LSXVAL LSXVOL M1 M2 MB MKKR PRATE QM CA CAPI CASI KA DI FADIL FAPI NEO OB


**# ---------------------------- XR analysis -----------------------------------
* Loading dataset
	import excel $path_in\dby.xlsx, sheet("Sheet1") firstrow
	tsset year

	lab var FDI 	"Foreign direct investment (net, BoP, US$)"
	lab var INF 	"Inflation (%)"
	lab var EDLT	"External debt stock (long term, DOD, US$)"
	* DOD: disbursed and outstanding debt
	lab var DS		"Public and publicly debt service (% GNI)"
	lab var RPD		"Total reserves (% of total debt)"
	lab var GDPG    "GDP growth"
	lab var CAB    	"Current account balance (BoP, US$)"
	lab var XRUSD	"Exchange rate (Kip per US$, period average)"
	lab var XRCNY	"Exchange rate (Kip per CNY, period average)"
	lab var XRTHB	"Exchange rate (Kip per THB, period average)"

	gen 	warukr = 0
	replace warukr = 1 if year > 2021

	gen 	covid19 = 0
	replace covid19 = 41 	if year == 2020
	replace covid19 = 106190 if year == 2021
	replace covid19 = 111411 if year == 2022
	replace covid19 = 1264 	if year == 2023
	replace covid19 = 154 	if year == 2024
	lab var covid19 "Cases of COVID-19"

* Standarization
	foreach vari of global var_list {
		egen z`vari' = std(`vari')
	}

* Dynamic factor models
	dfactor (zXRTHB zXRCNY, noconstant) (f = zINF zFDI zDS zRPD zCAB zGDPG zcovid19, ar(1)), method(dejong) diff
	estat ic
	predict factor, factor
	gen xriskfactor = factor
	lab var xriskfactor "Latent exchange rate risk factor"
	tsline xriskfactor, tline(2005)
	
	tsline FDI, tline(2005)
	tsline RPD, tline(2010)
	tsline CAB, tline(2010)
	tsline DS, tline(2015)
	
**# ------------------------ Quarterly XR analysis -----------------------------
* Loading dataset
	import excel $path_in\dbq.xlsx, sheet("Sheet1") firstrow clear
	generate time = tq(2012q3) + _n-1
	format time %tq
	tsset time, quarterly
	
	lab var M1 	"Narrow Money (M1, LAK)"
	lab var OB 	"Overall balance (BoP)"	
	
* Standarization
	foreach vari of global var_lqrt {
		egen z`vari' = std(`vari')
	}

* Dynamic factor model
	dfactor (zFX zFXRESUS, noconstant) (f = zFA zINFCORE zDEPCAP zM1 zPRATE zCA zCAPI zDI zFADIL zFAPI zOB, ar(1)), method(dejong) diff
	estat ic
	predict factor, factor
	gen XRQfactor = factor
	lab var XRQfactor "Latent exchange rate risk factor (quarterly)"
	tsline XRQfactor, tline(2020q3)
	
	tsline M1, tline(2019q3)
	tsline OB, tline(2020q3)
	
**# ------------------------ Quarterly BoP analysis -----------------------------
* Loading dataset
	import excel $path_in\dbbop.xlsx, sheet("Sheet1") firstrow clear
	generate time = tq(2012q1) + _n-1
	format time %tq
	tsset time, quarterly

	lab var CA 		"Current account"
	lab var KA 		"Capital account (USD)"
	lab var FA 		"Financial account"
	lab var CAPI 	"Current account primary income (USD)"
	lab var CASI 	"Current account secondary income"
	lab var DI 		"Direct investment (USD)"
	lab var FADIL 	"Liabilities (USD)"
	lab var FAPI 	"Portfolio investment"

* Standarization
	foreach vari of global var_lbop {
		egen z`vari' = std(`vari')
	}

* Dynamic factor models
	dfactor (zCA zKA zFA zCAPI zCASI zDI zFADIL zFAPI, noconstant) (f =, ar(1)), method(dejong) diff
	estat ic
	dfactor (zCA zKA zFA zCAPI zFAPI, noconstant) (f =, ar(1)), method(dejong) diff
	estat ic
	predict factor, factor
	gen BoPfactor = factor
	lab var BoPfactor "Latent factor of BoP risk"
	tsline BoPfactor, tline(2016q1)
	tsline KA, tline(2015q3)
	tsline CAPI, tline(2016q1)

* ------------------------- e n d  o f  s c r i p t ----------------------------
