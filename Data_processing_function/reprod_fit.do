* fit reproduction ratio
import delimited using "D:\data\COVID-19\ready2fit\trimmed\combined_5_1.csv"
*import delimited using "D:\data\COVID-19\ready2fit\trimmed\combined_R_weighted_gamma20.csv"
*import delimited using "D:\data\COVID-19\ready2fit\trimmed\AU_CIR_5_1.csv"
*import delimited using "D:\data\COVID-19\ready2fit\trimmed\CAN_CIR_5_1.csv"
replace trans_ratio=. if trans_ratio==-9999
* if a very large R suggests infecion data error
*replace trans_ratio=. if trans_ratio>4
*replace trans_ratio=. if trans_ratio>5
* remove dates with no community spreading
replace date =. if date==-9999
*1 case senario
replace region=11 if region==8 & obj==7 
replace region=11 if region==8 & obj==16 
replace region=11 if region==8 & obj==17
replace obj=1 if region==11 & obj==7
replace obj=2 if region==11 & obj==16
replace obj=3 if region==11 & obj==17
replace obj=obj+30 if region==6 
replace region=5 if region==6

*2 case senario


**** NO2 ***********
*replace no2 =. if no2==0
generate logNo2=log(no2)
generate scaledNo2=no2*1e4
*****mete****************
generate logPrec=log(p)
generate scaledPrec=p*1e5
generate scaledUV=uv/1e5
generate logUV=log(uv/3600)
generate logSolar=log(sr)
generate ctemp=t-273.15
generate logsh=log(sh)
generate scaledSH=sh*1e3

*******standardize**************
*egen sLogNo2=std(logNo2)
*egen sT=std(ctemp)
*egen sLogUV=std(logUV)
*egen sWind=std(ws)
*egen srh=std(rh)
*******OLS regression***************
* regression with all global records
*reghdfe trans_ratio sLogNo2 sLogUV sT sWind srh, absorb( region#date,savefe)

*gen panelname=string(region, "%02.0f") +string(obj, "%02.0f")
*encode panelname, gen(pn)
*sort pn date
*duplicates report pn date
*duplicates list pn date


reghdfe trans_ratio logNo2 logUV ctemp ws rh, absorb (region#date,savefe)
*********regression with developed contries and China****************
*reghdfe trans_ratio logNo2 ctemp if region <=3 | region==7 | region==10 , absorb( region#date,savefe)

* regression with controlled temperature
*reghdfe trans_ratio rh ctemp ws if ctemp>=25, absorb( region#date,savefe)
*reghdfe trans_ratio logNo2 ctemp logUV rh ws if ctemp>=25, absorb( region#date,savefe)


predict R_pred
generate pred_trans_ratio=R_pred+__hdfe1__
scatter trans_ratio pred_trans_ratio 
export delimited date obj region __hdfe1__ using "D:\data\COVID-19\ready2fit\out\coef.csv"
clear
outreg2 using x.xls

