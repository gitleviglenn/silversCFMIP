;-------------------------------------------------------------------------
; radflux_am4CRE_panA.ncl 
;-------------------------------------------------------------------------
;
; compute CRE as a function of latitude for particular CFMIP experiments
;
; feedbacks are also computed for these experiments using matlab
; glbmn_toa_cfmip_fb_driver.m which also uses
; glbmn_toa_cfmip_fb.m
;
; biases can be computed relative to either version of the CERES data
; or relative to the piControl experiment of CM4.0
;
; The cloud radiative effects are normalized for different experiments according
; to differing time periods.  I still need to determine an appropriate normalization
; for the 1pct CO2 ramp experiment and the two ceres profiles.  
;
; the temperature change of the air near the surface is not exactly 4K but:
; delT_amip_fut = 4.6816
; delT_amip_p4K = 4.4938
; delT_amip_m4K = 4.5217 ; verified in this script
; delT_aqua     = 4.1391 ; verified in this script;  why isn't this exactly 4K?  
;
; levi silvers                                            june 2020
;-------------------------------------------------------------------------

begin

path="/Users/silvers/data/cfmip_toafluxes/"

path_obs="/Users/silvers/data/SatData/CERES/"

; everything = 0; plotbias = 0 --> Panel on left, showing the raw CRE for CERES-EBAF, coupled experiments, and amip
; everything = 2; plotbias = 0 --> Panel on far right showing amip type feedback experiments

; as of November 2020
; everything = 0; plotbias = 1 --> middle panel of figure 7, with CRE feedback relative to piControl

everything = 0 ; switch defining which profiles are plotted
; 0 plots the CRE from: ceres, amip, historical, abrupt, and 1pct
; 1 plots the shortwave, longwave, and net CRE for amip and aqua
; 2 plots CRE Feedback figure showing aqua and amip-based warming experiments

plotbias = 0 ; default is 0; for plotbias=1, a bias is plotted 

; years for the abrupt and ramp warming experiments
early="000101-010012"
late="010101-015012"

; years for the historical experiment
firsthalf="185001-194912"
seconhalf="195001-201412"

hist_yrs=seconhalf
warm_yrs=late
;
;  David has recommended that we use the years 1996-2015 for the historical 
;  and amip experiments to exclude some of the volcanic activity from the 
;  early 90's.  He doesn't have confidence in the GCMs ability to deal with
;  that period.  
;  The coupled runs from CMIP6 only extend through december of 2014, so 
;  we are likely looking at jan 1996 through dec 2014.  19 years, or 228 months.
;
;  as originally written, the historical files over the later periods cover 65 
;  years, which is 780 months. i want to grab the last 228 of them.  
;
lastmonth=779
timeper1=228 ;228 = 19 years
;timeper1=777
firstmonth=lastmonth-timeper1-1  
;
; since amip starts in 1979, the range is different...
; 432 months is 36 years
lastmonth_amip = 431
timeper2=timeper1
;timeper2=430 ;228 = 19 years
firstmonth_amip=lastmonth_amip-timeper2-1  
;-------------------------------------------------------------------------
; define variables borrowed from the matlab scripts: 
delT_amip_fut = 4.6816
delT_amip_p4K = 4.4938
delT_amip_m4K = 4.5217
;delT_aqua     = 4.1391
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
; for amip 
;-------------------------------------------------------------------------

; define filenames
filename="atmos_cmip.197901-201412.rsut.nc"
infile1=path+filename
filename="atmos_cmip.197901-201412.rsutcs.nc"
infile2=path+filename
filename="atmos_cmip.197901-201412.rlut.nc"
infile3=path+filename
filename="atmos_cmip.197901-201412.rlutcs.nc"
infile4=path+filename
filename="atmos_cmip.197901-201412.rsdt.nc"
infile5=path+filename
filename="tas_Amon_GFDL-AM4_amip_r1i1p1f1_gr1_198001-201412.nc"
;filename="atmos_cmip.197901-201412.tas.nc"
infile6=path+filename

  print("infile1 is: "+infile1)

; ceres data paths: 
filename="CERES_EBAF-TOA_Ed2.8_Subset_200003-201607.nc"
filename_4p1="CERES_EBAF_Ed4.1_Subset_200003-201809.nc"
infile_obs=path_obs+filename
infile_obs_4p1=path_obs+filename_4p1
  print("infile_obs ceres v2.8 is: "+infile_obs)
  print("infile_obs_4p1 ceres v4.1 is: "+infile_obs_4p1)

crs = addfile(infile_obs,"r")

; read data
add_rsut=addfile(infile1,"r")
amip_rsut=add_rsut->rsut(firstmonth_amip:lastmonth_amip,:,:)
add_rsutcs=addfile(infile2,"r")
amip_rsutcs=add_rsutcs->rsutcs(firstmonth_amip:lastmonth_amip,:,:)
add_rlut=addfile(infile3,"r")
amip_rlut=add_rlut->rlut(firstmonth_amip:lastmonth_amip,:,:)
add_rlutcs=addfile(infile4,"r")
amip_rlutcs=add_rlutcs->rlutcs(firstmonth_amip:lastmonth_amip,:,:)
add_rsdt=addfile(infile5,"r")
amip_rsdt=add_rsdt->rsdt(firstmonth_amip:lastmonth_amip,:,:)
add_tas=addfile(infile6,"r")
amip_tas=add_tas->tas(:,:,:)

; compute time average of tas
amip_tas_tmn=dim_avg_n(amip_tas,0)

; compute weighted global mean of tas
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=
latam      = add_tas->lat
rad      = 4.0*atan(1.0)/180.0
clatam     = cos(latam*rad)

amip_tas_wave=wgt_areaave(amip_tas_tmn,clatam,1.0,0)

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("weighted mean of amip over particular years is: "+amip_tas_wave)
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=

;-------------------------------------------------------------------------
; piControl
;-------------------------------------------------------------------------

; define filenames
filename="rsut_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_1=path+filename
filename="rsutcs_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_2=path+filename
filename="rlut_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_3=path+filename
filename="rlutcs_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_4=path+filename
filename="rsdt_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_5=path+filename
filename="tas_Amon_GFDL-CM4_piControl_r1i1p1f1_gr1_025101-035012.nc"
infile_pic_6=path+filename

; read data
add_rsut=addfile(infile_pic_1,"r")
pic_rsut=add_rsut->rsut(:,:,:)
add_rsutcs=addfile(infile_pic_2,"r")
pic_rsutcs=add_rsutcs->rsutcs(:,:,:)
add_rlut=addfile(infile_pic_3,"r")
pic_rlut=add_rlut->rlut(:,:,:)
add_rlutcs=addfile(infile_pic_4,"r")
pic_rlutcs=add_rlutcs->rlutcs(:,:,:)
add_rsdt=addfile(infile_pic_5,"r")
pic_rsdt=add_rsdt->rsdt(:,:,:)
add_tas=addfile(infile_pic_6,"r")
pic_tas=add_tas->tas(:,:,:)

;-------------------------------------------------------------------------

; preparation for cosine weighted global mean averages...
latitude=add_rsut->lat(:)

lat      = crs->lat
rad      = 4.0*atan(1.0)/180.0
clat     = cos(lat*rad)
clat_am4 = cos(latitude*rad)


;----------- average data   
amip_rsut_tmn=dim_avg_n(amip_rsut,0)
amip_rsut_tgmn=wgt_areaave(amip_rsut_tmn,clat_am4,1.0,0)
amip_rsut_tzmn=dim_avg_n(amip_rsut_tmn,1)

amip_rsutcs_tmn=dim_avg_n(amip_rsutcs,0)
amip_rsutcs_tgmn=wgt_areaave(amip_rsutcs_tmn,clat_am4,1.0,0)
amip_rsutcs_tzmn=dim_avg_n(amip_rsutcs_tmn,1)

amip_rlut_tmn=dim_avg_n(amip_rlut,0)
amip_rlut_tgmn=wgt_areaave(amip_rlut_tmn,clat_am4,1.0,0)
amip_rlut_tzmn=dim_avg_n(amip_rlut_tmn,1)

amip_rlutcs_tmn=dim_avg_n(amip_rlutcs,0)
amip_rlutcs_tgmn=wgt_areaave(amip_rlutcs_tmn,clat_am4,1.0,0)
amip_rlutcs_tzmn=dim_avg_n(amip_rlutcs_tmn,1)

amip_rsdt_tmn=dim_avg_n(amip_rsdt,0)
amip_rsdt_tgmn=wgt_areaave(amip_rsdt_tmn,clat_am4,1.0,0)
amip_rsdt_tzmn=dim_avg_n(amip_rsdt_tmn,1)

pic_rsut_tmn=dim_avg_n(pic_rsut,0)
pic_rsut_tgmn=wgt_areaave(pic_rsut_tmn,clat_am4,1.0,0)
pic_rsut_tzmn=dim_avg_n(pic_rsut_tmn,1)

pic_rsutcs_tmn=dim_avg_n(pic_rsutcs,0)
pic_rsutcs_tgmn=wgt_areaave(pic_rsutcs_tmn,clat_am4,1.0,0)
pic_rsutcs_tzmn=dim_avg_n(pic_rsutcs_tmn,1)

pic_rlut_tmn=dim_avg_n(pic_rlut,0)
pic_rlut_tgmn=wgt_areaave(pic_rlut_tmn,clat_am4,1.0,0)
pic_rlut_tzmn=dim_avg_n(pic_rlut_tmn,1)

pic_rlutcs_tmn=dim_avg_n(pic_rlutcs,0)
pic_rlutcs_tgmn=wgt_areaave(pic_rlutcs_tmn,clat_am4,1.0,0)
pic_rlutcs_tzmn=dim_avg_n(pic_rlutcs_tmn,1)

pic_rsdt_tmn=dim_avg_n(pic_rsdt,0)
pic_rsdt_tgmn=wgt_areaave(pic_rsdt_tmn,clat_am4,1.0,0)
pic_rsdt_tzmn=dim_avg_n(pic_rsdt_tmn,1)

pic_tas_tmn=dim_avg_n(pic_tas,0)
pic_tas_tgmn=wgt_areaave(pic_tas_tmn,clat_am4,1.0,0)
;  print("global time mean of tas of PICONTROL is: "+pic_tas_tgmn)
;  print("dimensions of amip_rsut are: "+dimsizes(amip_rsut))
;  print("dimensions of pic_rsut are:  "+dimsizes(pic_rsut))

;------------ calculate CRE: cloud radiative effect
; for amip
cre_lw_amip=amip_rlut_tzmn-amip_rlutcs_tzmn
cre_sw_amip=amip_rsut_tzmn-amip_rsutcs_tzmn
cre_net_amip=cre_sw_amip+cre_lw_amip

; for cosine weighted global time mean values...
cre_lw_amip_gmn=amip_rlut_tgmn-amip_rlutcs_tgmn
cre_sw_amip_gmn=amip_rsut_tgmn-amip_rsutcs_tgmn
cre_net_amip_tgmn=cre_sw_amip_gmn+cre_lw_amip_gmn

;  print("dimsizes of amip_rlut_tzmn are: "+dimsizes(amip_rlut_tzmn))
;  print("dimsizes of cre_net_amip are: "+dimsizes(cre_net_amip))

cre_lw_amip_tmn =amip_rlut_tmn-amip_rlutcs_tmn
cre_sw_amip_tmn  =amip_rsut_tmn-amip_rsutcs_tmn
cre_net_amip_tmn =cre_sw_amip_tmn+cre_lw_amip_tmn

cre_lw_amip!0="lat"
cre_lw_amip&lat=amip_rsut&lat

cre_sw_amip!0="lat"
cre_sw_amip&lat=amip_rsut&lat

cre_net_amip!0="lat"
cre_net_amip&lat=amip_rsut&lat

cre_lw_amip_tmn!0="lat"
cre_lw_amip_tmn&lat=amip_rsut&lat
cre_lw_amip_tmn!1="lon"
cre_lw_amip_tmn&lon=amip_rsut&lon

cre_sw_amip_tmn!0="lat"
cre_sw_amip_tmn&lat=amip_rsut&lat
cre_sw_amip_tmn!1="lon"
cre_sw_amip_tmn&lon=amip_rsut&lon

cre_net_amip_tmn!0="lat"
cre_net_amip_tmn&lat=amip_rsut&lat
cre_net_amip_tmn!1="lon"
cre_net_amip_tmn&lon=amip_rsut&lon

; for pic
cre_lw_pic=pic_rlut_tzmn-pic_rlutcs_tzmn
cre_sw_pic=pic_rsut_tzmn-pic_rsutcs_tzmn
cre_net_pic=cre_sw_pic+cre_lw_pic

cre_lw_pic_tmn =pic_rlut_tmn-pic_rlutcs_tmn
cre_sw_pic_tmn  =pic_rsut_tmn-pic_rsutcs_tmn
cre_net_pic_tmn =cre_sw_pic_tmn+cre_lw_pic_tmn

cre_lw_pic!0="lat"
cre_lw_pic&lat=pic_rsut&lat

cre_sw_pic!0="lat"
cre_sw_pic&lat=pic_rsut&lat

cre_net_pic!0="lat"
cre_net_pic&lat=pic_rsut&lat

cre_lw_pic_tmn!0="lat"
cre_lw_pic_tmn&lat=pic_rsut&lat
cre_lw_pic_tmn!1="lon"
cre_lw_pic_tmn&lon=pic_rsut&lon

cre_sw_pic_tmn!0="lat"
cre_sw_pic_tmn&lat=pic_rsut&lat
cre_sw_pic_tmn!1="lon"
cre_sw_pic_tmn&lon=pic_rsut&lon

cre_net_pic_tmn!0="lat"
cre_net_pic_tmn&lat=pic_rsut&lat
cre_net_pic_tmn!1="lon"
cre_net_pic_tmn&lon=pic_rsut&lon

;-------------------------------------------------------------------------
;---------------------------------------------------------------------
; ceres data read below
;---------------------------------------------------------------------
;  crs = addfile(infile_obs,"r")
varname1="toa_cre_sw_mon"
varname2="toa_cre_lw_mon"
varname3="toa_cre_net_mon"
ceres_swcre = crs->$varname1$(:,:,:)
ceres_lwcre = crs->$varname2$(:,:,:)
ceres_netcre = crs->$varname3$(:,:,:)

print("dims before averaging "+dimsizes(ceres_swcre))
ceres_swcre_tmn  =dim_avg_n(ceres_swcre,0)
ceres_lwcre_tmn  =dim_avg_n(ceres_lwcre,0)
ceres_netcre_tmn =dim_avg_n(ceres_netcre,0)
print("dims after averaging "+dimsizes(ceres_swcre_tmn))

;
  crs_4p1 = addfile(infile_obs_4p1,"r")
varname1="toa_sw_all_mon"
varname2="toa_sw_clr_t_mon"; what is the difference between the clr_t and clr_c vars?
;varname2="toa_sw_clr_c_mon"; what is the difference between the clr_t and clr_c vars?
varname3="toa_lw_all_mon"
varname4="toa_lw_clr_t_mon"
;varname4="toa_lw_clr_c_mon"
ceres_sw_all = crs_4p1->$varname1$(:,:,:)
ceres_sw_clr = crs_4p1->$varname2$(:,:,:)
ceres_lw_all = crs_4p1->$varname3$(:,:,:)
ceres_lw_clr = crs_4p1->$varname4$(:,:,:)

ceres_lw_cre=-ceres_lw_all+ceres_lw_clr
ceres_lw_cre_tmn=dim_avg_n(ceres_lw_cre,0)
ceres_sw_cre=-ceres_sw_all+ceres_sw_clr
ceres_sw_cre_tmn=dim_avg_n(ceres_sw_cre,0)

ceres_2p8_cre=ceres_netcre_tmn
;ceres_2p8_cre=ceres_lwcre_tmn+ceres_swcre_tmn
ceres_4p1_cre=ceres_lw_cre_tmn+ceres_sw_cre_tmn

; ceres data is at a higher resolution (180x360 than) so we need to interpolate 
; CERES to the AM4 grid if we want to plot global maps

ceres_2p8_cre_ztmn=dim_avg_n(ceres_2p8_cre,1)
ceres_4p1_cre_ztmn=dim_avg_n(ceres_4p1_cre,1)

;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile5=path+filename
filename="tas_Amon_GFDL-CM4_historical_r1i1p1f1_gr1_"+hist_yrs+".nc"
infile6=path+filename
;
;; read data
add_rsut_historical=addfile(infile1,"r")
h2_rsut_historical=add_rsut_historical->rsut(firstmonth:lastmonth,:,:)
add_rsutcs_historical=addfile(infile2,"r")
h2_rsutcs_historical=add_rsutcs_historical->rsutcs(firstmonth:lastmonth,:,:)
add_rlut_historical=addfile(infile3,"r")
h2_rlut_historical=add_rlut_historical->rlut(firstmonth:lastmonth,:,:)
add_rlutcs_historical=addfile(infile4,"r")
h2_rlutcs_historical=add_rlutcs_historical->rlutcs(firstmonth:lastmonth,:,:)
add_rsdt_historical=addfile(infile5,"r")
h2_rsdt_historical=add_rsdt_historical->rsdt(firstmonth:lastmonth,:,:)
add_tas_historical=addfile(infile6,"r")
h2_tas_historical=add_tas_historical->tas(firstmonth:lastmonth,:,:)
;
print("!@#$Q$#^#$^%&$^*")
print("Years? one of the historical files is: "+infile2)
print("!@#$Q$#^#$^%&$^*")
;
;; average data
h2_tas_historical_tmn=dim_avg_n(h2_tas_historical,0)
h2_rsut_historical_tmn=dim_avg_n(h2_rsut_historical,0)
h2_rsut_historical_tzmn=dim_avg_n(h2_rsut_historical_tmn,1)
h2_rsutcs_historical_tmn=dim_avg_n(h2_rsutcs_historical,0)
h2_rsutcs_historical_tzmn=dim_avg_n(h2_rsutcs_historical_tmn,1)
h2_rlut_historical_tmn=dim_avg_n(h2_rlut_historical,0)
h2_rlut_historical_tzmn=dim_avg_n(h2_rlut_historical_tmn,1)
h2_rlutcs_historical_tmn=dim_avg_n(h2_rlutcs_historical,0)
h2_rlutcs_historical_tzmn=dim_avg_n(h2_rlutcs_historical_tmn,1)
h2_rsdt_historical_tmn=dim_avg_n(h2_rsdt_historical,0)
h2_rsdt_historical_tzmn=dim_avg_n(h2_rsdt_historical_tmn,1)
;
h2_tas_historical_tgmn=wgt_areaave(h2_tas_historical_tmn,clat_am4,1.0,0)
del_tas_hist=h2_tas_historical_tgmn-pic_tas_tgmn
;
;; calculate CRE: cloud radiative effect
cre_lw_h2_historical=h2_rlut_historical_tzmn-h2_rlutcs_historical_tzmn
cre_sw_h2_historical=h2_rsut_historical_tzmn-h2_rsutcs_historical_tzmn
cre_net_h2_historical=cre_sw_h2_historical+cre_lw_h2_historical
;

cre_lw_h2_historical!0="lat"
cre_lw_h2_historical&lat=h2_rsut_historical&lat
cre_sw_h2_historical!0="lat"
cre_sw_h2_historical&lat=h2_rsut_historical&lat
cre_net_h2_historical!0="lat"
cre_net_h2_historical&lat=h2_rsut_historical&lat
;
;;-------------------------------------------------------------------------
; for 1pctCO2 
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile5=path+filename
filename="tas_Amon_GFDL-CM4_1pctCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile6=path+filename
  print("Years?  One of the 1pctCO2 files is: "+filename)
;
;; read data
add_rsut_1pctCO2=addfile(infile1,"r")
abr2_rsut_1pctCO2=add_rsut_1pctCO2->rsut(:,:,:)
add_rsutcs_1pctCO2=addfile(infile2,"r")
abr2_rsutcs_1pctCO2=add_rsutcs_1pctCO2->rsutcs(:,:,:)
add_rlut_1pctCO2=addfile(infile3,"r")
abr2_rlut_1pctCO2=add_rlut_1pctCO2->rlut(:,:,:)
add_rlutcs_1pctCO2=addfile(infile4,"r")
abr2_rlutcs_1pctCO2=add_rlutcs_1pctCO2->rlutcs(:,:,:)
add_rsdt_1pctCO2=addfile(infile5,"r")
abr2_rsdt_1pctCO2=add_rsdt_1pctCO2->rsdt(:,:,:)
add_tas_1pctCO2=addfile(infile6,"r")
abr2_tas_1pctCO2=add_tas_1pctCO2->tas(:,:,:)

;; average data
abr2_tas_1pctCO2_tmn=dim_avg_n(abr2_tas_1pctCO2,0)
abr2_rsut_1pctCO2_tmn=dim_avg_n(abr2_rsut_1pctCO2,0)
abr2_rsut_1pctCO2_tzmn=dim_avg_n(abr2_rsut_1pctCO2_tmn,1)
abr2_rsutcs_1pctCO2_tmn=dim_avg_n(abr2_rsutcs_1pctCO2,0)
abr2_rsutcs_1pctCO2_tzmn=dim_avg_n(abr2_rsutcs_1pctCO2_tmn,1)
abr2_rlut_1pctCO2_tmn=dim_avg_n(abr2_rlut_1pctCO2,0)
abr2_rlut_1pctCO2_tzmn=dim_avg_n(abr2_rlut_1pctCO2_tmn,1)
abr2_rlutcs_1pctCO2_tmn=dim_avg_n(abr2_rlutcs_1pctCO2,0)
abr2_rlutcs_1pctCO2_tzmn=dim_avg_n(abr2_rlutcs_1pctCO2_tmn,1)
abr2_rsdt_1pctCO2_tmn=dim_avg_n(abr2_rsdt_1pctCO2,0)
abr2_rsdt_1pctCO2_tzmn=dim_avg_n(abr2_rsdt_1pctCO2_tmn,1)
;
; calculate cosine weighted global mean tas
abr2_tas_1pctCO2_tgmn=wgt_areaave(abr2_tas_1pctCO2_tmn,clat_am4,1.0,0)

;; calculate CRE: cloud radiative effect
cre_lw_1pct=abr2_rlut_1pctCO2_tzmn-abr2_rlutcs_1pctCO2_tzmn
cre_sw_1pct=abr2_rsut_1pctCO2_tzmn-abr2_rsutcs_1pctCO2_tzmn
cre_net_1pct=cre_sw_1pct+cre_lw_1pct
;
cre_lw_1pct!0="lat"
cre_lw_1pct&lat=abr2_rsut_1pctCO2&lat
cre_sw_1pct!0="lat"
cre_sw_1pct&lat=abr2_rsut_1pctCO2&lat
cre_net_1pct!0="lat"
cre_net_1pct&lat=abr2_rsut_1pctCO2&lat
;
;;-------------------------------------------------------------------------
; for abrupt 4xCO2 
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile5=path+filename
filename="tas_Amon_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gr1_"+warm_yrs+".nc"
infile6=path+filename
;
  print("Years?  One of the 1pctCO2 files is: "+filename)
;
;; read data
add_rsut_abrupt4xCO2=addfile(infile1,"r")
abr2_rsut_abrupt4xCO2=add_rsut_abrupt4xCO2->rsut(:,:,:)
add_rsutcs_abrupt4xCO2=addfile(infile2,"r")
abr2_rsutcs_abrupt4xCO2=add_rsutcs_abrupt4xCO2->rsutcs(:,:,:)
add_rlut_abrupt4xCO2=addfile(infile3,"r")
abr2_rlut_abrupt4xCO2=add_rlut_abrupt4xCO2->rlut(:,:,:)
add_rlutcs_abrupt4xCO2=addfile(infile4,"r")
abr2_rlutcs_abrupt4xCO2=add_rlutcs_abrupt4xCO2->rlutcs(:,:,:)
add_rsdt_abrupt4xCO2=addfile(infile5,"r")
abr2_rsdt_abrupt4xCO2=add_rsdt_abrupt4xCO2->rsdt(:,:,:)
add_tas_abrupt4xCO2=addfile(infile6,"r")
abr2_tas_abrupt4xCO2=add_tas_abrupt4xCO2->tas(:,:,:)
;
;; average data
abr2_tas_abrupt4xCO2_tmn=dim_avg_n(abr2_tas_abrupt4xCO2,0)
abr2_rsut_abrupt4xCO2_tmn=dim_avg_n(abr2_rsut_abrupt4xCO2,0)
abr2_rsut_abrupt4xCO2_tzmn=dim_avg_n(abr2_rsut_abrupt4xCO2_tmn,1)
abr2_rsutcs_abrupt4xCO2_tmn=dim_avg_n(abr2_rsutcs_abrupt4xCO2,0)
abr2_rsutcs_abrupt4xCO2_tzmn=dim_avg_n(abr2_rsutcs_abrupt4xCO2_tmn,1)
abr2_rlut_abrupt4xCO2_tmn=dim_avg_n(abr2_rlut_abrupt4xCO2,0)
abr2_rlut_abrupt4xCO2_tzmn=dim_avg_n(abr2_rlut_abrupt4xCO2_tmn,1)
abr2_rlutcs_abrupt4xCO2_tmn=dim_avg_n(abr2_rlutcs_abrupt4xCO2,0)
abr2_rlutcs_abrupt4xCO2_tzmn=dim_avg_n(abr2_rlutcs_abrupt4xCO2_tmn,1)
abr2_rsdt_abrupt4xCO2_tmn=dim_avg_n(abr2_rsdt_abrupt4xCO2,0)
abr2_rsdt_abrupt4xCO2_tzmn=dim_avg_n(abr2_rsdt_abrupt4xCO2_tmn,1)

; calculate cosine weighted global mean tas
abr2_tas_abrupt4xCO2_tgmn=wgt_areaave(abr2_tas_abrupt4xCO2_tmn,clat_am4,1.0,0)
;
;; calculate CRE: cloud radiative effect
cre_lw_abr2=abr2_rlut_abrupt4xCO2_tzmn-abr2_rlutcs_abrupt4xCO2_tzmn
cre_sw_abr2=abr2_rsut_abrupt4xCO2_tzmn-abr2_rsutcs_abrupt4xCO2_tzmn
cre_net_abr2=cre_sw_abr2+cre_lw_abr2
;
cre_lw_abr2!0="lat"
cre_lw_abr2&lat=abr2_rsut_abrupt4xCO2&lat
cre_sw_abr2!0="lat"
cre_sw_abr2&lat=abr2_rsut_abrupt4xCO2&lat
cre_net_abr2!0="lat"
cre_net_abr2&lat=abr2_rsut_abrupt4xCO2&lat
;
;--------------------------------------------------------------------------
;
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=

;-------------------------------------------------------------------------
; compute bias of cre relative to a particular control exp.
;-------------------------------------------------------------------------

; control CRE:
cont_cre=-cre_net_pic 
;-------------------------------------------------------------------------

;-------------------------------------------------------------------------
; compute global mean values of data from ceres and am4
;-------------------------------------------------------------------------
axis_coslat=new((/180/),float)

ceres_2p8_ave=wgt_areaave(ceres_2p8_cre,clat,1.0,0)
ceres_4p1_ave=wgt_areaave(ceres_4p1_cre,clat,1.0,0)

  print("mean of ceres 2p8: "+dim_avg_n(ceres_2p8_cre_ztmn,0))
  print("w mean of ceres 2p8: "+ceres_2p8_ave)
  print("mean of ceres 4p1: "+dim_avg_n(ceres_4p1_cre_ztmn,0))
  print("w mean of ceres 4p1: "+ceres_4p1_ave)

; amip p4K 
cre_amip=wgt_areaave(cre_net_amip_tmn,clat_am4,1.0,0)

  print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  print(" The following values indicate that the sign of the CRE and Feedbacks should be switched...")
  print(" Clouds should cool the earth overall, by roughly 20 W/m2 ")
  print("w mean of amip time series net:          "+cre_amip)
  print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")

;-------------------------------------------------------------------------
; make plot
  ;wks   = gsn_open_wks ("newPDF","testplot")
  wks = gsn_open_wks("eps","testplot")          ; send graphics to ps file
  plot = new(1,graphic)

; resources for the panels, and the legend
  res                        = True
  lgres                      = True
 
  fontsize    = 0.03
  lthick      = 2.0
  lthin       = 1.0
  ldash       = 16
  ldash_solid = 0
 
  res@tiMainFont = 21
  res@tiMainFontHeightF = fontsize
  ;res@tiYAxisString  = "[W/m2K]" 
  ;res@tiYAxisString  = "[W/m2]" 
  res@tiYAxisFont    = 21
  res@tiYAxisFontHeightF = fontsize 
  res@tiXAxisString  = "Latitude"
  ;res@tiXAxisString  = ""
  res@tiXAxisFont    = 21
  res@tiXAxisFontHeightF = fontsize
 
  res@vpWidthF          = 1.2 
  res@vpHeightF         = 0.5 
 
  res@tmYLLabelFontHeightF = fontsize 
  res@tmXBLabelFontHeightF = fontsize
 
  lthick=3.0
  res@xyLineThicknessF     = lthick
  res@tmXBMode = "explicit"
  res@tmXBValues = (/-60,-30,0,30,60/)
  res@tmXBLabels = (/-60,-30,0,30,60/)
  res@tmXBMinorValues =(/-80.,-70.,-50.,-40.,-20.,-10.,10.,20.,40.,50.,70.,80./)

; here is where the colors, line styles, and string names are defined 
 ; colors = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid4","darkslategray4","gold"/) 
    colors = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid3","darkslategray4","cyan3","red"/) 
  title_1="net CRE"
  title_2="sw CRE"
  title_3="lw CRE"
  labels = (/title_1,title_2,title_3/) ; bottom to top
  ; the plotting order is: observations,r,s,u ; top to bottom
 
 res@gsnDraw         = False
 res@gsnFrame         = False

if (everything .ge. 2) then ; plot amip based feedback of CRE
 res@tiMainString = "SW CRE Feedback"
 res@tiYAxisString  = "[W/m2K]" 
 res@trYMinF = -4
 res@trYMaxF = 1 
else if (everything .ge. 1) then ; plot the whole shebang
 res@tiMainString = "CRE "
 res@tiYAxisString  = "[W/m2]" 
 res@trYMinF = -70 
 res@trYMaxF = 110 
else
 res@tiMainString = "Cloud Radiative Effect (CRE)"
; res@tiYAxisString  = "[W/m2]" 
 res@tiYAxisString  = "CRE (W "+"m~S~-2"+"~N~"+")" 
 res@trYMinF = -50 
 res@trYMaxF = 10 
end if
end if
if (plotbias .gt. 0) then ; plot the cre feedback
end if

 res@xyLineColors=colors(3)
 res@xyDashPattern="0"

 res@xyLineColors=colors(1)
 res@xyDashPattern="0"
 plot1a=gsn_csm_xy(wks,latitude(:),cre_net_amip(:),res)
 res@xyLineColors=colors(1)
 plot1b=gsn_csm_xy(wks,latitude(:),cre_sw_amip(:),res)
 res@xyLineColors=colors(1)
 plot1c=gsn_csm_xy(wks,latitude(:),cre_lw_amip(:),res)


 res@xyDashPattern="0"
 res@xyLineColors=colors(4)
 plot5a=gsn_csm_xy(wks,latitude(:),cre_net_h2_historical(:),res)
 res@xyDashPattern="16"
 res@xyLineColors=colors(4)
 plot5b=gsn_csm_xy(wks,latitude(:),cre_sw_h2_historical(:),res)
 res@xyLineColors=colors(4)
 plot5c=gsn_csm_xy(wks,latitude(:),cre_lw_h2_historical(:),res)

 res@xyDashPattern="0"
 res@xyLineColors=colors(5)
 plot6a=gsn_csm_xy(wks,latitude(:),cre_net_abr2(:),res)
 res@xyDashPattern="16"
 res@xyLineColors=colors(5)
 plot6b=gsn_csm_xy(wks,latitude(:),cre_sw_abr2(:),res)
 res@xyLineColors=colors(5)
 plot6c=gsn_csm_xy(wks,latitude(:),cre_lw_abr2(:),res)

;;; compute bias
;
 res@xyLineColors=colors(1)
 res@xyDashPattern="0"
 plot_alt_1=gsn_csm_xy(wks,latitude(:),-cre_net_amip(:),res)
 res@xyDashPattern="0"
 res@xyLineColors=colors(4)
 plot_alt_2=gsn_csm_xy(wks,latitude(:),-cre_net_h2_historical(:),res)
 res@xyDashPattern="0"
 res@xyLineColors=colors(5)
 plot_alt_3=gsn_csm_xy(wks,latitude(:),-cre_net_abr2(:),res)
 res@xyLineColors=colors(6)
 plot_alt_4=gsn_csm_xy(wks,latitude(:),-cre_net_1pct(:),res)
 res@xyDashPattern="0"
 res@xyLineColors="grey70"
 plot_alt_5=gsn_csm_xy(wks,latitude(:),ceres_2p8_cre_ztmn(:),res)
 res@xyDashPattern="0"
 res@xyLineColors=colors(0)
 plot_alt_6=gsn_csm_xy(wks,latitude(:),ceres_4p1_cre_ztmn(:),res)

print("plotted ceres 2p8 mean is: "+dim_avg_n(ceres_2p8_cre_ztmn,0))
print("plotted ceres 4p1 mean is: "+dim_avg_n(ceres_4p1_cre_ztmn,0))
print("plotted cre net amip mean is: "+dim_avg_n(cre_net_amip,0))
print("plotted cre net h2 historical mean is: "+dim_avg_n(cre_net_h2_historical,0))
print("plotted cre net abr2 mean is: "+dim_avg_n(cre_net_abr2,0))
print("plotted cre net 1pct mean is: "+dim_avg_n(cre_net_1pct,0))

if (everything .ge. 2) then ; plot the cre feedback
  colors_cre = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid3","darkslategray4","cyan3","red","steelblue1"/) 
  labels_cre = (/"amip +future4K","amip +4K","amip -4K","aqua +4K"/) ; bottom to top

testaxis=cos(latitude*rad)
testaxis=where(latitude.lt.0,-testaxis,testaxis)

  plot(0)=plot0a
else
  overlay(plot_alt_1,plot_alt_2)
  overlay(plot_alt_1,plot_alt_3)
  overlay(plot_alt_1,plot_alt_4)
  overlay(plot_alt_1,plot_alt_5)
  overlay(plot_alt_1,plot_alt_6)
  plot(0)=plot_alt_1
end if

;-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=--=-=-=-=-=-=-=-=-=-=-=-=-
;       add a legend

 lgres                      = True
 lgres@xyLabelMode            = "Custom"
 lgres@xyLineLabelFontHeightF = 0.020                ; font height
 lgres@lgItemType           = "Lines"
 lgres@lgLabelFontHeightF   = .07
 lgres@vpWidthF             = 0.3        ; width of legend
 lgres@vpHeightF            = 0.20        ; height of legend
 lgres@lgLineThicknessF     = lthick
 lgres@lgPerimThicknessF    = 2.0
 lgres@lgMonoDashIndex      = False 
 lgres@lgPerimOn            = False

  resP                 = True
; draw panel with white space added
if (plotbias .le. 0) then ; plot the cre feedback
  if (everything .ge. 2) then
   print("everything is greather than 2")
   legend  = gsn_create_legend (wks, 4, labels_cre, lgres)
   lgres@lgLineColors         = (/colors_cre(7),colors_cre(2),colors_cre(3),colors_cre(8)/)
   lgres@lgDashIndexes = (/"0","0","0","0"/)
   legend  = gsn_create_legend (wks, 4, labels_cre, lgres)
   ;resP                 = True
   resP@amJust = "BottomRight"
   ;resP@amParallelPosF   = 0.42    ; Move legend to right
   ;resP@amOrthogonalPosF = -0.1     ; Move legend down more negative moves higher
   resP@amParallelPosF   = 0.32    ; Move legend to right
   resP@amOrthogonalPosF = 0.35     ; Move legend down more negative moves higher
   annoid = gsn_add_annotation(plot0a,legend,resP) ; add legend to plot
  else if(everything .ge. 1) then ; plot the whole shebang
   legend  = gsn_create_legend (wks, 3, labels, lgres)
   lgres@lgDashIndexes = (/"0","0","0","0"/)
   ;resP                 = True
   resP@amJust = "BottomRight"
   ;resP@amParallelPosF   = 0.42    ; Move legend to right
   resP@amParallelPosF   = -0.12    ; Move legend to right
   resP@amOrthogonalPosF = -0.85     ; Move legend down more negative moves higher
   annoid = gsn_add_annotation(plot0a,legend,resP) ; add legend to plot
  else
   gsres               = True
   txres               = True
   gsres@gsMarkerIndex = 16
   txres@txFontHeightF = 0.015
   ;colors_leg = (/"darkgoldenrod","mediumorchid3","cyan3","darkslategray4","grey70","black"/) 
   colors_leg = (/"darkgoldenrod","mediumorchid3","darkslategray4","cyan3","grey70","black"/) 
   labels_b = (/"amip","historical","abrupt 4xCO2","1% CO2 Ramp","CERES-EBAF v2.8","CERES-EBAF v4.1"/)
   yb=0.37
   yt=0.41
   xleg = (/0.40,0.40,0.54,0.54,0.74,0.74/)   ; Location of
   xtxt = (/0.46,0.48,0.65,0.65,0.87,0.87/)   ; legend markers
   yleg = (/yt,yb,yt,yb,yt,yb/)   ; and text
   ytxt = (/yt,yb,yt,yb,yt,yb/)   ; strings.
   do i = 0,5
     gsres@gsLineColor              = colors_leg(i)
     gsres@gsLineThicknessF         = 4.5
     gsn_polyline_ndc(wks,            (/xleg(i),xleg(i)+0.03/),(/yleg(i),yleg(i)/),gsres)
     gsn_text_ndc    (wks,labels_b(i),xtxt(i),ytxt(i),txres) 
   end do
  end if 
  end if
  ;resP                 = True
  ;frame(wks)
;  annoid = gsn_add_annotation(plot_alt_1,legend,resP) ; add legend to plot
end if

 gsn_panel(wks,plot,(/1,1/),resP)

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("global time mean of tas of PICONTROL is: "+pic_tas_tgmn)
  print("global time mean of tas of historical is: "+h2_tas_historical_tgmn)
  print("weighted mean of amip over particular years is: "+amip_tas_wave)
  print("global time mean of tas of 1pctCO2 is: "+abr2_tas_1pctCO2_tgmn)
  print("global time mean of tas of abrupt4xCO2 is: "+abr2_tas_abrupt4xCO2_tgmn)
  print("change of tas for selected historical period: "+del_tas_hist)
del_tas_1pct=abr2_tas_1pctCO2_tgmn-pic_tas_tgmn
  print("change of tas for selected 1pctCO2 period: "+del_tas_1pct)
del_tas_abrupt=abr2_tas_abrupt4xCO2_tgmn-pic_tas_tgmn
  print("change of tas for selected abrupt4xCO2 period: "+del_tas_abrupt)
del_tas_amip=amip_tas_wave-pic_tas_tgmn
  print("change of tas for selected amip period: "+del_tas_amip)
  print("BLUEFACE")
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")

;;-------------------------------------------------------------------------
end

