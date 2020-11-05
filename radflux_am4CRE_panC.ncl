;-------------------------------------------------------------------------
; radflux_am4CRE_panC.ncl 
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
; everything = 0; plotbias =0; Corresponds to the top panel in which the ceres CRE is plotted
;                              along with amip historical, abrupt4xCO2, and the 1pctCO2 experiments.
; everything = 0; plotbias =1; Corresponds to the middle panel
; everything = 2; plotbias =1; Corresponds to the bottom panel of the CRE figure in the CFMIP paper
;                              showing amip p4k, m4k, future4K, and aqua p4K
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

everything = 2 ; switch defining which profiles are plotted
; 0 plots the CRE from: ceres, amip, historical, abrupt, and 1pct
; 1 plots the shortwave, longwave, and net CRE for amip and aqua
; 2 plots CRE Feedback figure showing aqua and amip-based warming experiments

plotbias = 1 ; default is 0; for plotbias=1, a bias is plotted 

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
delT_aqua     = 4.1391
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
; for amip m4K
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile5=path+filename
filename="tas_Amon_GFDL-CM4_amip-m4K_r1i1p1f1_gr1_197901-201412.nc"
infile6=path+filename

print("infile1 is: "+infile1)

; read data
add_rsut_m4K=addfile(infile1,"r")
amip_rsut_m4K=add_rsut_m4K->rsut(:,:,:)
add_rsutcs_m4K=addfile(infile2,"r")
amip_rsutcs_m4K=add_rsutcs_m4K->rsutcs(:,:,:)
add_rlut_m4K=addfile(infile3,"r")
amip_rlut_m4K=add_rlut_m4K->rlut(:,:,:)
add_rlutcs_m4K=addfile(infile4,"r")
amip_rlutcs_m4K=add_rlutcs_m4K->rlutcs(:,:,:)
add_rsdt_m4K=addfile(infile5,"r")
amip_rsdt_m4K=add_rsdt_m4K->rsdt(:,:,:)
add_tas_m4K=addfile(infile6,"r")
amip_tas_m4K=add_tas_m4K->tas(:,:,:)

; average data
amip_rsut_m4K_tmn=dim_avg_n(amip_rsut_m4K,0)
amip_rsut_m4K_tzmn=dim_avg_n(amip_rsut_m4K_tmn,1)
amip_rsutcs_m4K_tmn=dim_avg_n(amip_rsutcs_m4K,0)
amip_rsutcs_m4K_tzmn=dim_avg_n(amip_rsutcs_m4K_tmn,1)
amip_rlut_m4K_tmn=dim_avg_n(amip_rlut_m4K,0)
amip_rlut_m4K_tzmn=dim_avg_n(amip_rlut_m4K_tmn,1)
amip_rlutcs_m4K_tmn=dim_avg_n(amip_rlutcs_m4K,0)
amip_rlutcs_m4K_tzmn=dim_avg_n(amip_rlutcs_m4K_tmn,1)
amip_rsdt_m4K_tmn=dim_avg_n(amip_rsdt_m4K,0)
amip_rsdt_m4K_tzmn=dim_avg_n(amip_rsdt_m4K_tmn,1)
amip_tas_m4K_tmn=dim_avg_n(amip_tas_m4K,0)

; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=
amip_tas_m4K_wave=wgt_areaave(amip_tas_m4K_tmn,clatam,1.0,0)

delT_m4K=amip_tas_wave-amip_tas_m4K_wave

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("weighted mean of amip m4K 1979-2014 is: "+amip_tas_m4K_wave)
  print("del T for amip m4K is: "+delT_m4K)
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=

;print("dimensions of amip_rsut_m4k are: "+dimsizes(amip_rsut_m4K))
;print("dimensions of amip_rsut_m4k_tmn are: "+dimsizes(amip_rsut_m4K_tmn))
;print("dimensions of amip_rsut_m4k_tzmn are: "+dimsizes(amip_rsut_m4K_tzmn))

; calculate CRE: cloud radiative effect
cre_lw_amip_m4K=amip_rlut_m4K_tzmn-amip_rlutcs_m4K_tzmn
cre_sw_amip_m4K=amip_rsut_m4K_tzmn-amip_rsutcs_m4K_tzmn
cre_net_amip_m4K=cre_sw_amip_m4K+cre_lw_amip_m4K

cre_lw_amip_m4K_tmn =amip_rlut_m4K_tmn-amip_rlutcs_m4K_tmn
cre_sw_amip_m4K_tmn  =amip_rsut_m4K_tmn-amip_rsutcs_m4K_tmn
cre_net_amip_m4K_tmn =cre_sw_amip_m4K_tmn+cre_lw_amip_m4K_tmn

cre_lw_amip_m4K!0="lat"
cre_lw_amip_m4K&lat=amip_rsut_m4K&lat

cre_sw_amip_m4K!0="lat"
cre_sw_amip_m4K&lat=amip_rsut_m4K&lat

cre_net_amip_m4K!0="lat"
cre_net_amip_m4K&lat=amip_rsut_m4K&lat

cre_lw_amip_m4K_tmn!0="lat"
cre_lw_amip_m4K_tmn&lat=amip_rsut_m4K&lat
cre_lw_amip_m4K_tmn!1="lon"
cre_lw_amip_m4K_tmn&lon=amip_rsut_m4K&lon

cre_sw_amip_m4K_tmn!0="lat"
cre_sw_amip_m4K_tmn&lat=amip_rsut_m4K&lat
cre_sw_amip_m4K_tmn!1="lon"
cre_sw_amip_m4K_tmn&lon=amip_rsut_m4K&lon

cre_net_amip_m4K_tmn!0="lat"
cre_net_amip_m4K_tmn&lat=amip_rsut_m4K&lat
cre_net_amip_m4K_tmn!1="lon"
cre_net_amip_m4K_tmn&lon=amip_rsut_m4K&lon
;-------------------------------------------------------------------------

;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
; for amip p4K
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_amip-p4K_r1i1p1f1_gr1_197901-201412.nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_amip-p4K_r1i1p1f1_gr1_197901-201412.nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_amip-p4K_r1i1p1f1_gr1_197901-201412.nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_amip-p4K_r1i1p1f1_gr1_197901-201412.nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_amip-p4K_r1i1p1f1_gr1_197901-201412.nc"
infile5=path+filename

; read data
add_rsut_p4K=addfile(infile1,"r")
amip_rsut_p4K=add_rsut_p4K->rsut(:,:,:)
add_rsutcs_p4K=addfile(infile2,"r")
amip_rsutcs_p4K=add_rsutcs_p4K->rsutcs(:,:,:)
add_rlut_p4K=addfile(infile3,"r")
amip_rlut_p4K=add_rlut_p4K->rlut(:,:,:)
add_rlutcs_p4K=addfile(infile4,"r")
amip_rlutcs_p4K=add_rlutcs_p4K->rlutcs(:,:,:)
add_rsdt_p4K=addfile(infile5,"r")
amip_rsdt_p4K=add_rsdt_p4K->rsdt(:,:,:)

; average data
amip_rsut_p4K_tmn=dim_avg_n(amip_rsut_p4K,0)
amip_rsut_p4K_tzmn=dim_avg_n(amip_rsut_p4K_tmn,1)
amip_rsutcs_p4K_tmn=dim_avg_n(amip_rsutcs_p4K,0)
amip_rsutcs_p4K_tzmn=dim_avg_n(amip_rsutcs_p4K_tmn,1)
amip_rlut_p4K_tmn=dim_avg_n(amip_rlut_p4K,0)
amip_rlut_p4K_tzmn=dim_avg_n(amip_rlut_p4K_tmn,1)
amip_rlutcs_p4K_tmn=dim_avg_n(amip_rlutcs_p4K,0)
amip_rlutcs_p4K_tzmn=dim_avg_n(amip_rlutcs_p4K_tmn,1)
amip_rsdt_p4K_tmn=dim_avg_n(amip_rsdt_p4K,0)
amip_rsdt_p4K_tzmn=dim_avg_n(amip_rsdt_p4K_tmn,1)

; calculate CRE: cloud radiative effect
cre_lw_amip_p4K=amip_rlut_p4K_tzmn-amip_rlutcs_p4K_tzmn
cre_sw_amip_p4K=amip_rsut_p4K_tzmn-amip_rsutcs_p4K_tzmn
cre_net_amip_p4K=cre_sw_amip_p4K+cre_lw_amip_p4K

cre_lw_amip_p4K_tmn=amip_rlut_p4K_tmn-amip_rlutcs_p4K_tmn
cre_sw_amip_p4K_tmn=amip_rsut_p4K_tmn-amip_rsutcs_p4K_tmn
cre_net_amip_p4K_tmn=cre_sw_amip_p4K_tmn+cre_lw_amip_p4K_tmn

cre_lw_amip_p4K!0="lat"
cre_lw_amip_p4K&lat=amip_rsut_p4K&lat

cre_sw_amip_p4K!0="lat"
cre_sw_amip_p4K&lat=amip_rsut_p4K&lat

cre_net_amip_p4K!0="lat"
cre_net_amip_p4K&lat=amip_rsut_p4K&lat

cre_lw_amip_p4K_tmn!0="lat"
cre_lw_amip_p4K_tmn&lat=amip_rsut_p4K&lat
cre_lw_amip_p4K_tmn!1="lon"
cre_lw_amip_p4K_tmn&lon=amip_rsut_p4K&lon

cre_sw_amip_p4K_tmn!0="lat"
cre_sw_amip_p4K_tmn&lat=amip_rsut_p4K&lat
cre_sw_amip_p4K_tmn!1="lon"
cre_sw_amip_p4K_tmn&lon=amip_rsut_p4K&lon

cre_net_amip_p4K_tmn!0="lat"
cre_net_amip_p4K_tmn&lat=amip_rsut_p4K&lat
cre_net_amip_p4K_tmn!1="lon"
cre_net_amip_p4K_tmn&lon=amip_rsut_p4K&lon
;-------------------------------------------------------------------------
;-------------------------------------------------------------------------
; for amip future
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_amip-future4K_r1i1p1f1_gr1_197901-201412.nc"
infile5=path+filename

; read data
add_rsut_future4K=addfile(infile1,"r")
amip_rsut_future4K=add_rsut_future4K->rsut(:,:,:)
add_rsutcs_future4K=addfile(infile2,"r")
amip_rsutcs_future4K=add_rsutcs_future4K->rsutcs(:,:,:)
add_rlut_future4K=addfile(infile3,"r")
amip_rlut_future4K=add_rlut_future4K->rlut(:,:,:)
add_rlutcs_future4K=addfile(infile4,"r")
amip_rlutcs_future4K=add_rlutcs_future4K->rlutcs(:,:,:)
add_rsdt_future4K=addfile(infile5,"r")
amip_rsdt_future4K=add_rsdt_future4K->rsdt(:,:,:)

; average data
amip_rsut_future4K_tmn=dim_avg_n(amip_rsut_future4K,0)
amip_rsut_future4K_tgmn=wgt_areaave(amip_rsut_future4K_tmn,clat_am4,1.0,0)
amip_rsut_future4K_tzmn=dim_avg_n(amip_rsut_future4K_tmn,1)

amip_rsutcs_future4K_tmn=dim_avg_n(amip_rsutcs_future4K,0)
amip_rsutcs_future4K_tgmn=wgt_areaave(amip_rsutcs_future4K_tmn,clat_am4,1.0,0)
amip_rsutcs_future4K_tzmn=dim_avg_n(amip_rsutcs_future4K_tmn,1)

amip_rlut_future4K_tmn=dim_avg_n(amip_rlut_future4K,0)
amip_rlut_future4K_tgmn=wgt_areaave(amip_rlut_future4K_tmn,clat_am4,1.0,0)
amip_rlut_future4K_tzmn=dim_avg_n(amip_rlut_future4K_tmn,1)

amip_rlutcs_future4K_tmn=dim_avg_n(amip_rlutcs_future4K,0)
amip_rlutcs_future4K_tgmn=wgt_areaave(amip_rlutcs_future4K_tmn,clat_am4,1.0,0)
amip_rlutcs_future4K_tzmn=dim_avg_n(amip_rlutcs_future4K_tmn,1)

amip_rsdt_future4K_tmn=dim_avg_n(amip_rsdt_future4K,0)
amip_rsdt_future4K_tgmn=wgt_areaave(amip_rsdt_future4K_tmn,clat_am4,1.0,0)
amip_rsdt_future4K_tzmn=dim_avg_n(amip_rsdt_future4K_tmn,1)
;
; calculate CRE: cloud radiative effect
cre_lw_amip_future4K=amip_rlut_future4K_tzmn-amip_rlutcs_future4K_tzmn
cre_sw_amip_future4K=amip_rsut_future4K_tzmn-amip_rsutcs_future4K_tzmn
cre_net_amip_future4K=cre_sw_amip_future4K+cre_lw_amip_future4K

cre_lw_amip_future4K_tmn =amip_rlut_future4K_tmn-amip_rlutcs_future4K_tmn
cre_sw_amip_future4K_tmn  =amip_rsut_future4K_tmn-amip_rsutcs_future4K_tmn
cre_net_amip_future4K_tmn =cre_sw_amip_future4K_tmn+cre_lw_amip_future4K_tmn

; compute global time mean cre
cre_lw_amip_future4K_tgmn =amip_rlut_future4K_tgmn-amip_rlutcs_future4K_tgmn
cre_sw_amip_future4K_tgmn =amip_rsut_future4K_tgmn-amip_rsutcs_future4K_tgmn
cre_net_amip_future4K_tgmn =cre_sw_amip_future4K_tgmn+cre_lw_amip_future4K_tgmn

cre_lw_amip_future4K!0="lat"
cre_lw_amip_future4K&lat=amip_rsut_future4K&lat
cre_sw_amip_future4K!0="lat"
cre_sw_amip_future4K&lat=amip_rsut_future4K&lat
cre_net_amip_future4K!0="lat"
cre_net_amip_future4K&lat=amip_rsut_future4K&lat

cre_lw_amip_future4K_tmn!0="lat"
cre_lw_amip_future4K_tmn&lat=amip_rsut_future4K&lat
cre_lw_amip_future4K_tmn!1="lon"
cre_lw_amip_future4K_tmn&lon=amip_rsut_future4K&lon

cre_sw_amip_future4K_tmn!0="lat"
cre_sw_amip_future4K_tmn&lat=amip_rsut_future4K&lat
cre_sw_amip_future4K_tmn!1="lon"
cre_sw_amip_future4K_tmn&lon=amip_rsut_future4K&lon

cre_net_amip_future4K_tmn!0="lat"
cre_net_amip_future4K_tmn&lat=amip_rsut_future4K&lat
cre_net_amip_future4K_tmn!1="lon"
cre_net_amip_future4K_tmn&lon=amip_rsut_future4K&lon
;;;-------------------------------------------------------------------------
; for aqua control 
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile5=path+filename

; get near surface air temperature...
filename="tas_Amon_GFDL-CM4_aqua-control_r1i1p1f1_gr1_198001-198912.nc"
infile6=path+filename

; read data
add_rsut_aqua=addfile(infile1,"r")
amip_rsut_aqua=add_rsut_aqua->rsut(:,:,:)
add_rsutcs_aqua=addfile(infile2,"r")
amip_rsutcs_aqua=add_rsutcs_aqua->rsutcs(:,:,:)
add_rlut_aqua=addfile(infile3,"r")
amip_rlut_aqua=add_rlut_aqua->rlut(:,:,:)
add_rlutcs_aqua=addfile(infile4,"r")
amip_rlutcs_aqua=add_rlutcs_aqua->rlutcs(:,:,:)
add_rsdt_aqua=addfile(infile5,"r")
amip_rsdt_aqua=add_rsdt_aqua->rsdt(:,:,:)
add_tas_aqua=addfile(infile6,"r")
amip_tas_aqua=add_tas_aqua->tas(:,:,:)

; average data
amip_rsut_aqua_tmn=dim_avg_n(amip_rsut_aqua,0)
amip_rsut_aqua_tzmn=dim_avg_n(amip_rsut_aqua_tmn,1)
amip_rsutcs_aqua_tmn=dim_avg_n(amip_rsutcs_aqua,0)
amip_rsutcs_aqua_tzmn=dim_avg_n(amip_rsutcs_aqua_tmn,1)
amip_rlut_aqua_tmn=dim_avg_n(amip_rlut_aqua,0)
amip_rlut_aqua_tzmn=dim_avg_n(amip_rlut_aqua_tmn,1)
amip_rlutcs_aqua_tmn=dim_avg_n(amip_rlutcs_aqua,0)
amip_rlutcs_aqua_tzmn=dim_avg_n(amip_rlutcs_aqua_tmn,1)
amip_rsdt_aqua_tmn=dim_avg_n(amip_rsdt_aqua,0)
amip_rsdt_aqua_tzmn=dim_avg_n(amip_rsdt_aqua_tmn,1)
; tas (check cosine weighting...)
amip_tas_aqua_tmn=dim_avg_n(amip_tas_aqua,0)
amip_tas_aqua_tzmn=dim_avg_n(amip_tas_aqua_tmn,1)

; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=
latap      = add_tas_aqua->lat
rad      = 4.0*atan(1.0)/180.0
clatap     = cos(latap*rad)

aqua_tas_wave=wgt_areaave(amip_tas_aqua_tmn,clatap,1.0,0)

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("weighted mean of aqua control is: "+aqua_tas_wave)
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=

; calculate CRE: cloud radiative effect
cre_lw_aqua=amip_rlut_aqua_tzmn-amip_rlutcs_aqua_tzmn
cre_sw_aqua=amip_rsut_aqua_tzmn-amip_rsutcs_aqua_tzmn
cre_net_aqua=cre_sw_aqua+cre_lw_aqua

cre_lw_aqua_tmn=amip_rlut_aqua_tmn-amip_rlutcs_aqua_tmn
cre_sw_aqua_tmn=amip_rsut_aqua_tmn-amip_rsutcs_aqua_tmn
cre_net_aqua_tmn=cre_sw_aqua_tmn+cre_lw_aqua_tmn

cre_lw_aqua!0="lat"
cre_lw_aqua&lat=amip_rsut_aqua&lat

cre_sw_aqua!0="lat"
cre_sw_aqua&lat=amip_rsut_aqua&lat

cre_net_aqua!0="lat"
cre_net_aqua&lat=amip_rsut_aqua&lat

cre_lw_aqua_tmn!0="lat"
cre_lw_aqua_tmn&lat=amip_rsut_aqua&lat
cre_lw_aqua_tmn!1="lon"
cre_lw_aqua_tmn&lon=amip_rsut_aqua&lon

cre_sw_aqua_tmn!0="lat"
cre_sw_aqua_tmn&lat=amip_rsut_aqua&lat
cre_sw_aqua_tmn!1="lon"
cre_sw_aqua_tmn&lon=amip_rsut_aqua&lon

cre_net_aqua_tmn!0="lat"
cre_net_aqua_tmn&lat=amip_rsut_aqua&lat
cre_net_aqua_tmn!1="lon"
cre_net_aqua_tmn&lon=amip_rsut_aqua&lon
;
;--------------------------------------------------------------------------
; for aqua p4K 
;-------------------------------------------------------------------------
filename="rsut_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile1=path+filename
filename="rsutcs_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile2=path+filename
filename="rlut_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile3=path+filename
filename="rlutcs_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile4=path+filename
filename="rsdt_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile5=path+filename
; get near surface air temperature...
filename="tas_Amon_GFDL-CM4_aqua-p4K_r1i1p1f1_gr1_198001-198912.nc"
infile6=path+filename


; read data
add_rsut_aqua_p4K=addfile(infile1,"r")
amip_rsut_aqua_p4K=add_rsut_aqua_p4K->rsut(:,:,:)
add_rsutcs_aqua_p4K=addfile(infile2,"r")
amip_rsutcs_aqua_p4K=add_rsutcs_aqua_p4K->rsutcs(:,:,:)
add_rlut_aqua_p4K=addfile(infile3,"r")
amip_rlut_aqua_p4K=add_rlut_aqua_p4K->rlut(:,:,:)
add_rlutcs_aqua_p4K=addfile(infile4,"r")
amip_rlutcs_aqua_p4K=add_rlutcs_aqua_p4K->rlutcs(:,:,:)
add_rsdt_aqua_p4K=addfile(infile5,"r")
amip_rsdt_aqua_p4K=add_rsdt_aqua_p4K->rsdt(:,:,:)
add_tas_aqua_p4K=addfile(infile6,"r")
amip_tas_aqua_p4K=add_tas_aqua_p4K->tas(:,:,:)

;--------------------------------------------------------------------
; average data      e.g. amip_rsut_aqua_p4K(time|lat|lon)
;  printVarSummary(amip_rsut_aqua_p4K)
amip_rsut_aqua_p4K_tmn=dim_avg_n(amip_rsut_aqua_p4K,0)
;  printVarSummary(amip_rsut_aqua_p4K_tmn)
;   this is the point at which global mean should be taken...
amip_rsut_aqua_p4K_tzmn=dim_avg_n(amip_rsut_aqua_p4K_tmn,1)
;  printVarSummary(amip_rsut_aqua_p4K_tzmn)

amip_rsutcs_aqua_p4K_tmn=dim_avg_n(amip_rsutcs_aqua_p4K,0)
amip_rsutcs_aqua_p4K_tzmn=dim_avg_n(amip_rsutcs_aqua_p4K_tmn,1)

amip_rlut_aqua_p4K_tmn=dim_avg_n(amip_rlut_aqua_p4K,0)
amip_rlut_aqua_p4K_tzmn=dim_avg_n(amip_rlut_aqua_p4K_tmn,1)

amip_rlutcs_aqua_p4K_tmn=dim_avg_n(amip_rlutcs_aqua_p4K,0)
amip_rlutcs_aqua_p4K_tzmn=dim_avg_n(amip_rlutcs_aqua_p4K_tmn,1)

amip_rsdt_aqua_p4K_tmn=dim_avg_n(amip_rsdt_aqua_p4K,0)
amip_rsdt_aqua_p4K_tzmn=dim_avg_n(amip_rsdt_aqua_p4K_tmn,1)

amip_tas_aqua_p4K_tmn=dim_avg_n(amip_tas_aqua_p4K,0)
amip_tas_aqua_p4K_tzmn=dim_avg_n(amip_tas_aqua_p4K_tmn,1)


; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=
latap      = add_tas_aqua_p4K->lat
rad      = 4.0*atan(1.0)/180.0
clatap     = cos(latap*rad)

aqua_tas_p4K_wave=wgt_areaave(amip_tas_aqua_p4K_tmn,clatap,1.0,0)

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("weighted mean of aqua p4K is: "+aqua_tas_p4K_wave)
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  snake=aqua_tas_p4K_wave-aqua_tas_wave
  print("del T for aqua +4K is: "+snake)
; =0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=

;--------------------------------------------------------------------
; calculate CRE: cloud radiative effect
cre_lw_aqua_p4K=amip_rlut_aqua_p4K_tzmn-amip_rlutcs_aqua_p4K_tzmn
cre_sw_aqua_p4K=amip_rsut_aqua_p4K_tzmn-amip_rsutcs_aqua_p4K_tzmn
cre_net_aqua_p4K=cre_sw_aqua_p4K+cre_lw_aqua_p4K

cre_lw_aqua_p4K_tmn=amip_rlut_aqua_p4K_tmn-amip_rlutcs_aqua_p4K_tmn
cre_sw_aqua_p4K_tmn=amip_rsut_aqua_p4K_tmn-amip_rsutcs_aqua_p4K_tmn
;  printVarSummary(cre_lw_aqua_p4K_tmn)
cre_net_aqua_p4K_tmn=cre_sw_aqua_p4K_tmn+cre_lw_aqua_p4K_tmn

cre_lw_aqua_p4K!0="lat"
cre_lw_aqua_p4K&lat=amip_rsut_aqua_p4K&lat

cre_sw_aqua_p4K!0="lat"
cre_sw_aqua_p4K&lat=amip_rsut_aqua_p4K&lat

cre_net_aqua_p4K!0="lat"
cre_net_aqua_p4K&lat=amip_rsut_aqua_p4K&lat

cre_net_aqua_p4K_tmn!0="lat"
cre_net_aqua_p4K_tmn&lat=amip_rsut_aqua&lat
cre_net_aqua_p4K_tmn!1="lon"
cre_net_aqua_p4K_tmn&lon=amip_rsut_aqua&lon

;
;-------------------------------------------------------------------------
; compute bias of cre relative to a particular control exp.
;-------------------------------------------------------------------------

; control CRE:
;cont_cre=ceres_2p8_cre_ztmn
;cont_cre=ceres_4p1_cre_ztmn
cont_cre=-cre_net_pic 

;-------------------------------------------------------------------------
; compute CRE feedbacks....
;-------------------------------------------------------------------------

cre_aqua_p4K_fb     =cre_net_aqua_p4K-cre_net_aqua
cresw_aqua_p4K_fb   =cre_sw_aqua_p4K-cre_sw_aqua
crelw_aqua_p4K_fb   =cre_lw_aqua_p4K-cre_lw_aqua

cre_amip_base_fb    =cre_net_amip-cre_net_amip
crelw_amip_base_fb  =cre_lw_amip-cre_lw_amip
cresw_amip_base_fb  =cre_sw_amip-cre_sw_amip

cre_amip_fut_fb     =cre_net_amip_future4K-cre_net_amip
cresw_amip_fut_fb   =cre_sw_amip_future4K-cre_sw_amip
crelw_amip_fut_fb   =cre_lw_amip_future4K-cre_lw_amip

cre_amip_p4K_fb     =cre_net_amip_p4K-cre_net_amip
cresw_amip_p4K_fb   =cre_sw_amip_p4K-cre_sw_amip
crelw_amip_p4K_fb   =cre_lw_amip_p4K-cre_lw_amip

;cre_amip_m4K_fb    =cre_net_amip_m4K-cre_net_amip
cre_amip_m4K_fb     =cre_net_amip-cre_net_amip_m4K
cresw_amip_m4K_fb   =cre_sw_amip-cre_sw_amip_m4K
crelw_amip_m4K_fb   =cre_lw_amip-cre_lw_amip_m4K

;cre_amip_m4K_fb=cre_net_amip_m4K-cre_net_amip

cresw_aqua_p4K_fb!0="lat"
cresw_aqua_p4K_fb&lat=amip_rsut_future4K&lat
crelw_aqua_p4K_fb!0="lat"
crelw_aqua_p4K_fb&lat=amip_rsut_future4K&lat

cresw_amip_base_fb!0="lat"
cresw_amip_base_fb&lat=amip_rsut_future4K&lat
crelw_amip_base_fb!0="lat"
crelw_amip_base_fb&lat=amip_rsut_future4K&lat

cresw_amip_fut_fb!0="lat"
cresw_amip_fut_fb&lat=amip_rsut_future4K&lat
crelw_amip_fut_fb!0="lat"
crelw_amip_fut_fb&lat=amip_rsut_future4K&lat

cresw_amip_p4K_fb!0="lat"
cresw_amip_p4K_fb&lat=amip_rsut_future4K&lat
crelw_amip_p4K_fb!0="lat"
crelw_amip_p4K_fb&lat=amip_rsut_future4K&lat

cresw_amip_m4K_fb!0="lat"
cresw_amip_m4K_fb&lat=amip_rsut_future4K&lat
crelw_amip_m4K_fb!0="lat"
crelw_amip_m4K_fb&lat=amip_rsut_future4K&lat

cre_aqua_p4K_fb!0="lat"
cre_aqua_p4K_fb&lat=amip_rsut_future4K&lat
cre_amip_base_fb!0="lat"
cre_amip_base_fb&lat=amip_rsut_future4K&lat
cre_amip_fut_fb!0="lat"
cre_amip_fut_fb&lat=amip_rsut_future4K&lat
cre_amip_p4K_fb!0="lat"
cre_amip_p4K_fb&lat=amip_rsut_future4K&lat
cre_amip_m4K_fb!0="lat"
cre_amip_m4K_fb&lat=amip_rsut_future4K&lat

;-------------------------------------------------------------------------
; compute global mean values of data from ceres and am4
;-------------------------------------------------------------------------

axis_coslat=new((/180/),float)

; amip p4K 
cre_p4K=wgt_areaave(cre_net_amip_p4K_tmn,clat_am4,1.0,0)
cre_amip=wgt_areaave(cre_net_amip_tmn,clat_am4,1.0,0)
amip_p4K_a=(cre_p4K-cre_amip)/delT_amip_p4K

; amip m4K
cre_m4K=wgt_areaave(cre_net_amip_m4K_tmn,clat_am4,1.0,0)
amip_m4K_a=(cre_m4K-cre_amip)/delT_amip_m4K
  ;print("dimsizes of cre_net_amip_m4K_tmn are: "+dimsizes(cre_net_amip_m4K_tmn))
  ;print("dimsizes of cre_m4K are: "+dimsizes(cre_m4K))
  ;print("dimsizes of amip_m4K_a are: "+dimsizes(amip_m4K_a))
  print("delT_amip_m4K is: "+delT_amip_m4K)

; amip future4K
cre_future4K=wgt_areaave(cre_net_amip_future4K_tmn,clat_am4,1.0,0)
amip_future4K_a=(cre_future4K-cre_amip)/delT_amip_fut
amip_future4K_b=(cre_net_amip_future4K_tgmn-cre_net_amip_tgmn)/delT_amip_fut

; aqua p4K 
cre_aqua=wgt_areaave(cre_net_aqua_tmn,clat_am4,1.0,0)
cre_aqua4K=wgt_areaave(cre_net_aqua_p4K_tmn,clat_am4,1.0,0)
aqua_p4K_a=(cre_aqua4K-cre_aqua)/delT_aqua

  print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  print(" The following values indicate that the sign of the CRE and Feedbacks should be switched...")
  print(" Clouds should cool the earth overall, by roughly 20 W/m2 ")
  print("w mean of CRE fb for aqua p 4K :      "+aqua_p4K_a)
  print("w mean of amip time series net:          "+cre_amip)
  print("w mean of amip p 4K net:     "+cre_p4K)
  print("w mean of amip p 4K fb:      "+amip_p4K_a)
  print("w mean of amip m 4K fb:      "+amip_m4K_a)
  print("old weighted mean of amip future 4K fb: "+amip_future4K_a)
  print("new weighted mean of amip future 4K fb: "+amip_future4K_b)
  print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")

amip_cre_fb_p4K =cre_amip_p4K_fb/delT_amip_p4K
amip_b=amip_cre_fb_p4K*clat_am4/180

  print("w jerry-rigged mean of amip p4K fb: "+sum(amip_b))

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
 res@tiYAxisString  = "CRE Feedback (W "+"m~S~-2~NN~K~S~-1"+"~N~"+")" 
 res@tiMainString = "SW Feedback relative to amip"
 res@trYMinF = -4 
 res@trYMaxF = 1 
end if

 res@xyLineColors=colors(3)
 res@xyDashPattern="0"
 plot0a=gsn_csm_xy(wks,latitude(:),cre_net_amip_m4K(:),res)
 res@xyDashPattern="2"
 res@xyLineColors=colors(3)
 plot0b=gsn_csm_xy(wks,latitude(:),cre_sw_amip_m4K(:),res)
 res@xyLineColors=colors(3)
 plot0c=gsn_csm_xy(wks,latitude(:),cre_lw_amip_m4K(:),res)

 res@xyLineColors=colors(1)
 res@xyDashPattern="0"
 plot1a=gsn_csm_xy(wks,latitude(:),cre_net_amip(:),res)
 res@xyLineColors=colors(1)
 plot1b=gsn_csm_xy(wks,latitude(:),cre_sw_amip(:),res)
 res@xyLineColors=colors(1)
 plot1c=gsn_csm_xy(wks,latitude(:),cre_lw_amip(:),res)

 res@xyLineColors=colors(2)
 res@xyDashPattern="0"
 plot2a=gsn_csm_xy(wks,latitude(:),cre_net_amip_p4K(:),res)
 res@xyDashPattern="16"
 res@xyLineColors=colors(2)
 plot2b=gsn_csm_xy(wks,latitude(:),cre_sw_amip_p4K(:),res)
 res@xyLineColors=colors(2)
 plot2c=gsn_csm_xy(wks,latitude(:),cre_lw_amip_p4K(:),res)

 res@xyDashPattern="0"
 res@xyLineColors=colors(0)
 plot3a=gsn_csm_xy(wks,latitude(:),cre_net_aqua(:),res)
 res@xyLineColors=colors(0)
 plot3b=gsn_csm_xy(wks,latitude(:),cre_sw_aqua(:),res)
 res@xyLineColors=colors(0)
 plot3c=gsn_csm_xy(wks,latitude(:),cre_lw_aqua(:),res)
 
 res@xyDashPattern="0"
 res@xyLineColors=colors(0)
 plot4a=gsn_csm_xy(wks,latitude(:),cre_net_aqua_p4K(:),res)
 res@xyDashPattern="16"
 res@xyLineColors=colors(0)
 plot4b=gsn_csm_xy(wks,latitude(:),cre_sw_aqua_p4K(:),res)
 res@xyLineColors=colors(0)
 plot4c=gsn_csm_xy(wks,latitude(:),cre_lw_aqua_p4K(:),res)


;;; compute bias
if (plotbias .gt. 0) then ; plot the cre feedback
 print("using a bias.  what bias? ")
 cre_net_amip          =(cre_net_amip+cont_cre)/0.76  ; why not normalized by 1.3K??
end if

  colors_cre = (/"black","darkgoldenrod","darkolivegreen3","steelblue","mediumorchid3","darkslategray4","cyan3","red","steelblue1"/) 
  labels_cre = (/"amip +future4K","amip +4K","amip -4K","aqua +4K"/) ; bottom to top

testaxis=cos(latitude*rad)
testaxis=where(latitude.lt.0,-testaxis,testaxis)

  res@xyLineColors=colors_cre(0)
  res@xyLineThicknessF  = lthin
  res@xyDashPattern = ldash_solid
  plot0a=gsn_csm_xy(wks,latitude(:),-cre_amip_base_fb(:),res)
  res@xyDashPattern = ldash
  plot0ac=gsn_csm_xy(wks,latitude(:),-crelw_amip_base_fb(:),res)
  res@xyDashPattern = ldash_solid
  res@xyLineThicknessF  = lthick
  plot0ab=gsn_csm_xy(wks,latitude(:),-cresw_amip_base_fb(:),res)

  res@xyLineColors=colors_cre(2)
  res@xyLineThicknessF  = lthin
  res@xyDashPattern = ldash_solid
  plot0b=gsn_csm_xy(wks,latitude(:),-cre_amip_p4K_fb(:)/delT_amip_p4K,res)
  res@xyDashPattern = ldash
  plot0bc=gsn_csm_xy(wks,latitude(:),-crelw_amip_p4K_fb(:)/delT_amip_p4K,res)
  res@xyDashPattern = ldash_solid
  res@xyLineThicknessF  = lthick
  plot0bb=gsn_csm_xy(wks,latitude(:),-cresw_amip_p4K_fb(:)/delT_amip_p4K,res)

  res@xyLineColors=colors_cre(3)
  res@xyLineThicknessF  = lthin
  res@xyDashPattern = ldash_solid
  plot0c=gsn_csm_xy(wks,latitude(:),-cre_amip_m4K_fb(:)/delT_amip_m4K,res)
  res@xyDashPattern = ldash
  plot0cc=gsn_csm_xy(wks,latitude(:),-crelw_amip_m4K_fb(:)/delT_amip_m4K,res)
  res@xyDashPattern = ldash_solid
  res@xyLineThicknessF  = lthick
  plot0cb=gsn_csm_xy(wks,latitude(:),-cresw_amip_m4K_fb(:)/delT_amip_m4K,res)

  res@xyLineColors=colors_cre(7)
  res@xyLineThicknessF  = lthin
  res@xyDashPattern = ldash_solid
  plot0d=gsn_csm_xy(wks,latitude(:),-cre_amip_fut_fb(:)/delT_amip_fut,res)
  res@xyDashPattern = ldash
  plot0dc=gsn_csm_xy(wks,latitude(:),-crelw_amip_fut_fb(:)/delT_amip_fut,res)
  res@xyDashPattern = ldash_solid
  res@xyLineThicknessF  = lthick
  plot0db=gsn_csm_xy(wks,latitude(:),-cresw_amip_fut_fb(:)/delT_amip_fut,res)

  res@xyLineColors=colors_cre(8)
  res@xyLineThicknessF  = lthin
  res@xyDashPattern = ldash_solid
  plot0e=gsn_csm_xy(wks,latitude(:),-cre_aqua_p4K_fb(:)/delT_aqua,res)
  res@xyDashPattern = ldash
  plot0ec=gsn_csm_xy(wks,latitude(:),-crelw_aqua_p4K_fb(:)/delT_aqua,res)
  res@xyDashPattern = ldash_solid
  res@xyLineThicknessF  = lthick
  plot0eb=gsn_csm_xy(wks,latitude(:),-cresw_aqua_p4K_fb(:)/delT_aqua,res)

; plot sw cre
;  overlay(plot0a,plot0bb) ; amip p4k
;  overlay(plot0bb,plot0b) ; amip p4k total
  overlay(plot0bb,plot0cb) ; amip m4k
;  overlay(plot0bb,plot0c) ; amip m4k
  overlay(plot0bb,plot0db) ; amip fut
;  overlay(plot0bb,plot0d) ; amip fut
  overlay(plot0bb,plot0eb) ; aqua p4k
;  overlay(plot0bb,plot0e) ; aqua p4k

  plot(0)=plot0bb

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
  gsres               = True
  txres               = True
  gsres@gsMarkerIndex = 16
  txres@txFontHeightF = 0.015
    colors_leg = (/"darkolivegreen3","steelblue","red","steelblue1"/) 
    labels_b = (/"amip p4K","amip m4K","amip Future4K","Aqua p4K"/)
  ;colors_leg = (/"darkgoldenrod","mediumorchid3","darkolivegreen3","darkslategray4","grey49","black"/) 
  ;labels_b = (/"amip","historical","abrupt 4xCO2","1% CO2 Ramp","CERES-EBAF v2.8","CERES-EBAF v4.1"/)
  yb=0.37
  yt=0.41
  xleg = (/0.40,0.40,0.54,0.54/)   ; Location of
  xtxt = (/0.48,0.48,0.65,0.62/)   ; legend markers
  yleg = (/yt,yb,yt,yb,yt,yb/)   ; and text
  ytxt = (/yt,yb,yt,yb,yt,yb/)   ; strings.
  do i = 0,3
    gsres@gsLineColor              = colors_leg(i)
    gsres@gsLineThicknessF         = 4.5
    gsn_polyline_ndc(wks,            (/xleg(i),xleg(i)+0.03/),(/yleg(i),yleg(i)/),gsres)
    gsn_text_ndc    (wks,labels_b(i),xtxt(i),ytxt(i),txres) 
  end do

 gsn_panel(wks,plot,(/1,1/),resP)

  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")
  print("global time mean of tas of PICONTROL is: "+pic_tas_tgmn)
  print("weighted mean of amip over particular years is: "+amip_tas_wave)
del_tas_amip=amip_tas_wave-pic_tas_tgmn
  print("change of tas for selected amip period: "+del_tas_amip)
  print("BLUEFACE")
  print("=0=-0=-0=-0=-0=0=-0=-0=-0=-0=-0=-0=-0-=0=")

;; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
;; compute average cre values and print to screen
; net_cre_amip             = dim_avg_n(cre_net_amip,0)
; print("net cre amip      = "+net_cre_amip)
; net_cre_amip_p4K         = dim_avg_n(cre_net_amip_p4K,0)
; print("net cre amip p4K  = "+net_cre_amip_p4K)
; net_cre_amip_m4K         = dim_avg_n(cre_net_amip_m4K,0)
; print("net cre amip m4K  = "+net_cre_amip_m4K)
;
;print("dimensions of amip future fb are: "+dimsizes(cre_amip_fut_fb)) 
;net_amip_fut_fb            = dim_avg_n(cre_amip_fut_fb,0)
;print("net cre fut fb      = "+net_amip_fut_fb)
;net_amip_p4K_fb            = dim_avg_n(cre_amip_p4K_fb,0)
;print("net cre p4K fb      = "+net_amip_p4K_fb)
;net_amip_m4K_fb            = dim_avg_n(cre_amip_m4K_fb,0)
;print("net cre m4K fb      = "+net_amip_m4K_fb)
;
; net_cre_aqua   =dim_avg_n(cre_net_aqua,0)
; print("net cre aqua = "+net_cre_aqua)
; net_cre_aquap4K=dim_avg_n(cre_net_aqua_p4K,0)
; print("net cre aqua+4K = "+net_cre_aquap4K)
;
;; compute the cre feedback values and print to screen
; amip_cre_fb_p4K =(net_cre_amip_p4K-net_cre_amip)/delT_amip_p4K
; print("amip cre fb p4K: "+amip_cre_fb_p4K)
; amip_cre_fb_m4K =(net_cre_amip-net_cre_amip_m4K)/delT_amip_m4K
; print("amip cre fb m4K: "+amip_cre_fb_m4K)
; aqua_cre_fb=(net_cre_aquap4K-net_cre_aqua)/delT_aqua
; print("net cre aqua: "+aqua_cre_fb)
;
;;-------------------------------------------------------------------------
end

