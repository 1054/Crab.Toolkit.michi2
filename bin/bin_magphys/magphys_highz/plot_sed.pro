; This code was made to interpret grafically the Spectral Energy Distribution
; for a galaxy, fitted with the code of E. da Cunha (2008).
; It takes the imputs from the DATA file of the code and creates a .ps file for 
; each galaxy that includes the theoretical SED of this galaxy, oveploted the data 
; points and also the probability distribution functions for some important parameters.
; T.Bitsakis(17-9-2010)
; 
; version Oct. 2011
; bugs fixed by Kate Rownlands, Univ. Nottigham
; 
; largely modified by dzliu
; 2016
; 20180618 added the output of converted wavelength (um) and flux (mJy).
; 
; Usage: idl84 -e "plot_sed" -args 10001 ; added by dzliu
; Usage: idl84 -e "plot_sed" -args "best-fit_SED" ; added by dzliu
; 






; Dependent functions

forward_function TeXtoIDL ; added by dzliu -- see www.idlcoyote.com/tips/func_var.html

function dzliu_lumtoflux, input_z, input_lambda_um, input_lum_Lsun
  dzliu_lumdist = lumdist(input_z)
  dzliu_lumdisq = 4*!PI*dzliu_lumdist^2
  ;dzliu_flux_mJy = (1.+input_z)*input_lum_Lsun/dzliu_lumdisq*40.31970d/(2.99792458e5/input_lambda_um) ; 1 Lsun Mpc-2 = 40.31970 mJy GHz
  dzliu_flux_mJy = input_lum_Lsun/dzliu_lumdisq*40.31970d/(2.99792458e5/input_lambda_um) ; 1 Lsun Mpc-2 = 40.31970 mJy GHz
  return, dzliu_flux_mJy
end

function dzliu_logtickformat, axis, index, number
    IF number LT 0 THEN RETURN, '-Inf'
    IF number GE 1E-2 - 1E-6 AND number LT 1E-1 - 1E-6 THEN RETURN, STRING(float(number), FORMAT='(F0.2)') ; - 1E-6 is for precision problem
    IF number GE 1E-1 - 1E-6 AND number LT 1E-0 - 1E-6 THEN RETURN, STRING(float(number), FORMAT='(F0.1)') ; - 1E-6 is for precision problem
    IF number GE 1E-0 - 1E-6 AND number LE 1E+2 - 1E-6 THEN RETURN, STRING(float(number), FORMAT='(I)')    ; - 1E-6 is for precision problem
    ex = STRING(number, FORMAT='(e8.0)')
    pt = STRPOS(ex, '.')
    first = STRTRIM(STRMID(ex, 0, pt),2)
    sign = STRMID(ex, pt+2, 1)
    exp = STRMID(ex, pt+3)
    WHILE STRMID(exp, 0, 1) EQ '0' DO exp = STRMID(exp, 1) ; shave off leading zero in exponent
    IF (LONG(exp) EQ 0) THEN BEGIN & sign = '' & exp = '0' & ENDIF
    IF sign EQ '+' THEN sign = ''
    IF first EQ '1' THEN first = '' ELSE first = first + 'x'
    RETURN, first + '10!U' + sign + exp + '!N'
end






; MAIN PROGRAM

pro plot_sed, galaxy
  
  ;================================================INPUTS===============================================
 
  ;galaxy='' ; commmented by dzliu, now read from input argument
  ;readcol,'DATA.dat',f=(A,F),galaxy,nn  ;Reads the data file to find the number of the elements as well the names of galaxies
  ;n=n_elements(nn)
  
  ; Read command line arguments
  if n_elements((command_line_args())) gt 0 then begin ; added by dzliu
      galaxy=(command_line_args())[0] ; added by dzliu
  endif ; added by dzliu
  if n_elements(galaxy) eq 0 then begin
      read,galaxy,prompt='Give the name of the source: ' 
  endif
  
  ; Check system variables
  if ~strlen(GETENV('magphys')) then message, 'Error! GETENV(magphys) was not defined!' ; added by dzliu
  if size(galaxy,/tname) ne 'STRING' then galaxy=STRING(format='(I0)',galaxy) ; added by dzliu
  if file_test(GETENV('magphys')+path_sep()+'textoidl',/dir) then !PATH=expand_path('+'+GETENV('magphys')+path_sep()+'textoidl')+':'+!PATH ; added by dzliu
  if file_test(GETENV('magphys')+path_sep()+'coyote',/dir) then !PATH=expand_path('+'+GETENV('magphys')+path_sep()+'coyote')+':'+!PATH ; added by dzliu
  resolve_all, /quiet ; added by dzliu
  print, 'Output to '+galaxy+'.ps'
  
  ; Read USER_FILTERS
  fmt='A,F,I,I'
  filter_name=''
  if ~strlen(GETENV('USER_FILTERS')) then message, 'Error! GETENV(USER_FILTERS) was not defined!'               ; added by dzliu
  if ~file_test(GETENV('USER_FILTERS'),/read) then message, 'Error! '+GETENV('USER_FILTERS')+' was not found!'  ; added by dzliu
  readcol,GETENV('USER_FILTERS'),f=fmt,filter_name,lambda_eff,filter_id,fit                                            ; Reads the filter file, lambda_eff are in units of um
  m = n_elements(lambda_eff)
  
  ; Read MAGPHYS *.fit data file
  name_ps=galaxy+'.ps'
  name_fit=galaxy+'.fit'
  name_sed=galaxy+'.sed'
  name_sed_out=galaxy+'.sed.um.mJy.txt'
  name_xy=galaxy
  if not file_test(name_fit,/read) then message, 'Error! '+name_fit+' was not found!' ; added by dzliu
  if not file_test(name_sed,/read) then message, 'Error! '+name_sed+' was not found!' ; added by dzliu
  oti1=''
  oti2=''
  oti3=''
  Lv_obs = fltarr(m)
  Lv_err = fltarr(m)
  Lv_fit = fltarr(m)
  openr,1,name_fit
  for i=1,2 do readf,1,oti1
  readf,1,Lv_obs                     ; notes by dzliu: observed flux values in *.fit are L_{\nu} in units of L_{\odot} Hz^{-1}
  readf,1,Lv_err                     ; notes by dzliu: observed flux errors in *.fit are L_{\nu} in units of L_{\odot} Hz^{-1}
  for i=1,4 do readf,1,oti2
  readf,1,i_sfh,i_ir,chi2,z
  for i=1,3 do readf,1,oti3
  readf,1,Lv_fit                     ; notes by dzliu: bset fit fluxes in *.fit are L_{\nu} in units of L_{\odot} Hz^{-1}
  
  ; Read MAGPHYS *.sed data file
  readcol,name_sed,lambda,lum_at,lum_un,SKIPLINE =10      ; notes by dzliu: lambda are in log and in units of {\AA}; lum_at are the L_{\lambda} in log and in units of L_{\odot} {\AA}^{-1}. 
  
  ;==============================================CALCULATIONS===========================================
  
  ; Convert the wavelengths
  w_obs = lambda_eff                   ; obs-frame wavelength in um. Note that "lambda_eff" are the observations data points, in linear and in units of um.
  w_sed = 10^lambda / 1D4              ; obs-frame wavelength in um. Note that "lambda" are the SED library data points, in log and in units of {\AA}.
  lambda_err = lambda_eff*0
  
  ; Convert the best fit SED fluxes to wLw ({\lambda}L_{\lambda}) or vLv ({\nu}L_{\nu})
  Lw_sed_at  = 10^lum_at                  ; L_{\lambda} in units of L_{\odot} {\AA}^{-1}, total attenuated SED
  Lw_sed_un  = 10^lum_un                  ; L_{\lambda} in units of L_{\odot} {\AA}^{-1}, unattenuated SED
  wLw_sed_at = 10^lambda*Lw_sed_at        ; {\lambda}L_{\lambda} in units of L_{\odot}, total attenuated SED
  wLw_sed_un = 10^lambda*Lw_sed_un        ; {\lambda}L_{\lambda} in units of L_{\odot}, unattenuated SED
  vLv_sed_at = wLw_sed_at                 ; {\nu}L_{\nu} = {\lambda}L_{\lambda} in units of L_{\odot}, total attenuated SED
  vLv_sed_un = wLw_sed_un                 ; {\nu}L_{\nu} = {\lambda}L_{\lambda} in units of L_{\odot}, unattenuated SED
  
  ;Lv_sed_at  = 10^lum_at                  ; L_{\nu} in units of L_{\odot} {Hz}^{-1}, total attenuated SED
  ;Lv_sed_un  = 10^lum_un                  ; L_{\nu} in units of L_{\odot} {Hz}^{-1}, unattenuated SED
  ;vLv_sed_at = (1.+z)*(Lv_sed_at)*2.99792458e8/(10^lambda/1e10) ; {\nu}L_{\nu} in units of L_{\odot}, total attenuated SED
  ;vLv_sed_un = (1.+z)*(Lv_sed_un)*2.99792458e8/(10^lambda/1e10) ; {\nu}L_{\nu} in units of L_{\odot}, unattenuated SED
  
  ; Convert the observed fluxes and best fit fluxes to vLv ({\nu}L_{\nu})
  vLv_obs     = (1.+z)*(Lv_obs       )*2.99792458e8/(w_obs/1e6)  ; notes by dzliu: vLv_obs are in units of L_{\odot}; (1.+z) is needed according to the original code of plot_sed.pro and the fit_sed_highz.f line 275.
  vLv_fit     = (1.+z)*(Lv_fit       )*2.99792458e8/(w_obs/1e6)
  vLv_obs_err = (1.+z)*(Lv_err       )*2.99792458e8/(w_obs/1e6)
  vLv_obs_lo  = (1.+z)*(Lv_obs-Lv_err)*2.99792458e8/(w_obs/1e6)
  vLv_obs_hi  = (1.+z)*(Lv_obs+Lv_err)*2.99792458e8/(w_obs/1e6)
  
  ; For error bars which go down to infinity (i.e. error is bigger than flux)
  IF N_ELEMENTS(where(vLv_obs_lo LE 1D-30, /NULL)) GT 0 THEN BEGIN
      vLv_obs_lo[where(vLv_obs_lo LE 1D-30, /NULL)] = 1D-30
  ENDIF
  
  ; Compute the difference/residual between observed vLv and best fit vLv
  vLv_res = alog10(vLv_obs) - alog10(vLv_fit)
  
  ; Compute flux density from vLv
  f_sed_at  = dzliu_lumtoflux(z,w_sed,vLv_sed_at)  ; -- the best fit dust attenuated SED flux
  f_sed_un  = dzliu_lumtoflux(z,w_sed,vLv_sed_un)  ; -- the best fit dust unattenuated SED flux
  f_fit     = dzliu_lumtoflux(z,w_obs,vLv_fit)     ; -- the best fit flux
  f_obs     = dzliu_lumtoflux(z,w_obs,vLv_obs)     ; -- the observed flux
  f_obs_err = dzliu_lumtoflux(z,w_obs,vLv_obs_err) ; -- the observed flux error
  f_obs_lo  = dzliu_lumtoflux(z,w_obs,vLv_obs_lo)  ; -- the observed flux lowest value for 1-sigma confidence
  f_obs_hi  = dzliu_lumtoflux(z,w_obs,vLv_obs_hi)  ; -- the observed flux highest value for 1-sigma confidence
  
  ; For values below 0, we set them to 1D-30
  IF N_ELEMENTS(where(f_fit LE 1D-30, /NULL)) GT 0 THEN BEGIN
      f_fit[where(f_fit LE 1D-30, /NULL)] = 1D-30
  ENDIF
  
  ; Compute the difference/residual between observed fluxes and best fit fluxes
  f_res    = alog10(f_obs   ) - alog10(f_fit)
  f_res_lo = alog10(f_obs_lo) - alog10(f_fit)
  f_res_hi = alog10(f_obs_hi) - alog10(f_fit)
  print, 'f_obs', f_obs
  print, 'f_fit', f_fit
  print, 'f_res', f_res
  
    
  ;================================================HISTOGRAM============================================
  
  f_muSFH=dblarr(2,20)
  f_muIR=dblarr(2,20)
  mu=dblarr(2,20)
  tau_V=dblarr(2,160) ; highz is 246-86=160 instead of 48
  sSFR=dblarr(2,70)
  Mstars=dblarr(2,60)
  Ldust=dblarr(2,60)
  Tc_ISM=dblarr(2,40) ;<dzliu modified> Tc_ISM=dblarr(2,10) -> Tc_ISM=dblarr(2,40), see also "dzliu modified" in "fit_sed_highz.f"
  Tw_BC=dblarr(2,60) ;<dzliu modified> Tw_BC=dblarr(2,30) -> Tw_BC=dblarr(2,60), see also "dzliu modified" in "fit_sed_highz.f"
  xi_C_tot=dblarr(2,20) ; highz added this
  xi_PAH_tot=dblarr(2,20) ; highz added this
  xi_MIR_tot=dblarr(2,20) ; highz added this
  xi_W_tot=dblarr(2,20) ; highz added this
  tau_V_ISM=dblarr(2,80)
  Mdust=dblarr(2,60)
  SFR=dblarr(2,60)
  age_M=dblarr(2,50) ; highz added this
  Mstar_LH_ratio=dblarr(2,100) ; highz added this
  Mstar_LK_ratio=dblarr(2,100) ; highz added this
  A_V=dblarr(2,80) ; highz added this
  Tdust=dblarr(2,14) ; highz added this
  ;;xi_C_tot=dblarr(2,20)
  ;;xi_W_tot=dblarr(2,20)
  skips=''
  
  name=''
  openr,20,name_fit
  
  for j=1,16 do readf,20,skips & readf,20,f_muSFH
  for j=1,3 do readf,20,skips & readf,20,f_muIR   & f_mu=0.5*(f_muIR+f_muSFH)
  for j=1,3 do readf,20,skips & readf,20,mu
  for j=1,3 do readf,20,skips & readf,20,tau_V
  for j=1,3 do readf,20,skips & readf,20,sSFR
  for j=1,3 do readf,20,skips & readf,20,Mstars
  for j=1,3 do readf,20,skips & readf,20,Ldust
  for j=1,3 do readf,20,skips & readf,20,Tc_ISM
  for j=1,3 do readf,20,skips & readf,20,Tw_BC
  for j=1,3 do readf,20,skips & readf,20,xi_C_tot
  for j=1,3 do readf,20,skips & readf,20,xi_PAH_tot
  for j=1,3 do readf,20,skips & readf,20,xi_MIR_tot
  for j=1,3 do readf,20,skips & readf,20,xi_W_tot
  for j=1,3 do readf,20,skips & readf,20,tau_V_ISM
  for j=1,3 do readf,20,skips & readf,20,Mdust
  for j=1,3 do readf,20,skips & readf,20,SFR
  for j=1,3 do readf,20,skips & readf,20,age_M
  for j=1,3 do readf,20,skips & readf,20,Mstar_LH_ratio
  for j=1,3 do readf,20,skips & readf,20,Mstar_LK_ratio
  for j=1,3 do readf,20,skips & readf,20,A_V
  for j=1,3 do readf,20,skips & readf,20,Tdust
  for j=1,2 do readf,20,skips
  
  xtitle01=TeXtoIDL('f_\mu') & xrange01=[0.01,0.98]
  ;xtitle03=TeXtoIDL('\mu')
  xtitle04=TeXtoIDL('\tau_V') & xrange04=[0.01,5.8]
  xtitle05=TeXtoIDL('log(M_{stars}/M_{o})') & xrange05=[8.1,12.9]
  xtitle06=TeXtoIDL('log(L_{dust}/L_{o})') & xrange06=[8.01,13.99]
  xtitle07=TeXtoIDL('T_{C}^{ISM}/K') & xrange07=[15.1,55.99]
  xtitle08=TeXtoIDL('T_{W}^{BC}/K') & xrange08=[30.5,99.5]
  ;<dzliu>;xtitle08=TeXtoIDL('T_{dust}/K')
  xtitle09=TeXtoIDL('\mu\tau_V') & xrange09=[0.01,2.95]
  xtitle10=TeXtoIDL('log(M_{dust}/M_{o})') & xrange10=[4.1,10.5]
  xtitle11=TeXtoIDL('log(sSFR) yr^{-1}') & xrange11=[-12.9,-7.1]
  xtitle12=TeXtoIDL('log(SFR/M_{o} yr^{-1})') & xrange12=[-3.9,5.0]
  xtitle13=TeXtoIDL('\xi_C^{tot}') & xrange13=[0.01,0.99]
  xtitle14=TeXtoIDL('\xi_W^{tot}') & xrange14=[0.01,0.99]
  
  xrange01=[]
  xrange04=[]
  xrange05=[]
  xrange06=[]
  xrange07=[]
  xrange08=[]
  xrange10=[]
  xrange11=[]
  xrange12=[]
  xrange13=[]
  xrange14=[]
  
  
  
  
  
  
  
  
  ;================================================SAVING=============================================
  save, FILENAME=name_sav, w_sed, f_sed_at, f_sed_un, filter_name, w_obs, f_obs, f_obs_err, f_obs_lo, f_obs_hi, f_fit, f_res, f_res_lo, f_res_hi, Mstars, SFR, sSFR, Ldust, Mdust, Tdust, age_M, A_V, tau_V, tau_V_ISM, mu, Tc_ISM, Tw_BC
  
  openw, lun, name_sed_out, /get_lun
  printf, lun, "wave_um", "f_attenu_mJy", "f_unattenu_mJy", "vLv_attenu_Lsun", "vLv_unattenu_Lsun", format='("# ",A-14," ",A16," ",A16," ",A16," ",A16)'
  for i=0,n_elements(w_sed)-1 do begin
    printf, lun, w_sed[i], f_sed_at[i], f_sed_un[i], vLv_sed_at[i], vLv_sed_un[i], format='("  ",E-14.6," ",E16.6," ",E16.6," ",E16.6," ",E16.6)'
  endfor
  close, lun
  free_lun, lun
  
  
  
  
  ;================================================PLOTTING=============================================
  
  device,decomposed=0
  
  set_plot, 'ps'
  !P.MULTI=[0,0,20]
  TVLCT, 0, 0, 0 , 0
  TVLCT, 200,55, 70 , 1
  TVLCT,  0,255, 0 , 2
  TVLCT, 0, 150, 255 , 3
  TVLCT, 100, 100, 100 , 4
  
  ; Define Positions
  ; upper panel
  positionU1=[0.120,0.55,0.920,0.99]
  positionU2=[0.120,0.45,0.920,0.55]
  ; 1st row
  position01=[0.150,0.25,0.275,0.35]
  position02=[0.275,0.25,0.400,0.35]
  position03=[0.400,0.25,0.525,0.35]
  position04=[0.525,0.25,0.650,0.35]
  position05=[0.650,0.25,0.775,0.35]
  position06=[0.775,0.25,0.900,0.35]
  ; 2nd row
  position07=[0.150,0.07,0.275,0.17]
  position08=[0.275,0.07,0.400,0.17]
  position09=[0.400,0.07,0.525,0.17]
  position10=[0.525,0.07,0.650,0.17]
  position11=[0.650,0.07,0.775,0.17]
  position12=[0.775,0.07,0.900,0.17]
  
  ; Open PS file to plot
  device,filename=name_ps,/color,XSIZE=21,YSIZE=21,/ENCAPSULATED ; size unit in centimeters, US letter 8.5 x 11.0 inches
  
  ; Define Axis Parameters
  ;ytitle = TeXtoIDL("log(\lambdaL_{\lambda}/L"+sunsymbol()+")")
  ;ytitle = "flux [mJy]" ; TeXtoIDL("flux [mJy]")
  xtitle = TextoIDL("\lambda/\mum [observed-frame]")
  xrange = [0.07,3e5]
  ;yrange = [7.1,14.]
  yrange = [1e-6,1e4]
  
  ; Plot frame
  plot,POSITION=positionU1,[1D-30],[1D-30],/NoData,/xlog,/ylog,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,ytickformat='dzliu_logtickformat',xminor=9,thick=3,xthick=3,ythick=3,charthick=3,charsize=1.75,XTICKFORMAT="(A1)"
  xyouts,0.05,(0.65+0.99)/2.0,/NORMAL,TeXtoIDL("flux [mJy]"),charthick=3,charsize=1.0,align=0.5,orient=90 ; ytitle
  
  ; Plot best fit SED
  w_flag = WHERE(f_sed_at GT 0 AND f_sed_un GT 0, /NULL)
  IF N_ELEMENTS(w_flag) GT 0 THEN BEGIN
      oplot,w_sed[w_flag],f_sed_at[w_flag]
      oplot,w_sed[w_flag],f_sed_un[w_flag],color=3
  ENDIF
  
  ; Plot observed fluxes
  cid_detected = []
  cid_undetect = []
  cid_detected_ALMA = []
  cid_undetect_ALMA = []
  FOR k=0,N_ELEMENTS(f_obs)-1 DO BEGIN
    IF f_obs[k]/f_obs_err[k] GE 3.0 THEN BEGIN
      IF STRMATCH(filter_name[k],'ALMA*') THEN BEGIN
        cid_detected_ALMA = [cid_detected_ALMA, k]
      ENDIF ELSE BEGIN
        cid_detected = [cid_detected, k]
      ENDELSE
    ENDIF ELSE BEGIN
      IF STRMATCH(filter_name[k],'ALMA*') THEN BEGIN
        cid_undetect_ALMA = [cid_undetect_ALMA, k]
      ENDIF ELSE BEGIN
        cid_undetect = [cid_undetect, k]
      ENDELSE
    ENDELSE
  ENDFOR
  IF N_ELEMENTS(cid_detected) GT 0 THEN BEGIN
    plotsym,8,0.9,/fill
    oploterror,w_obs[cid_detected],f_obs[cid_detected],lambda_err[cid_detected],(f_obs[cid_detected]-f_obs_lo[cid_detected]),psym=8,symsize=0.5,/LOBAR
    oploterror,w_obs[cid_detected],f_obs[cid_detected],lambda_err[cid_detected],(f_obs_hi[cid_detected]-f_obs[cid_detected]),psym=8,symsize=0.5,/HIBAR
  ENDIF
  IF N_ELEMENTS(cid_undetect) GT 0 THEN BEGIN
    plotsym,1
    oplot,w_obs[cid_undetect],3.0*f_obs_err[cid_undetect],psym=8,symsize=0.8
  ENDIF
  IF N_ELEMENTS(cid_detected_ALMA) GT 0 THEN BEGIN
    plotsym,8,0.9,/fill
    oploterror,w_obs[cid_detected_ALMA],f_obs[cid_detected_ALMA],lambda_err[cid_detected_ALMA],(f_obs[cid_detected_ALMA]-f_obs_lo[cid_detected_ALMA]),psym=8,symsize=0.5,color=1,errcolor=1,/LOBAR
    oploterror,w_obs[cid_detected_ALMA],f_obs[cid_detected_ALMA],lambda_err[cid_detected_ALMA],(f_obs_hi[cid_detected_ALMA]-f_obs[cid_detected_ALMA]),psym=8,symsize=0.5,color=1,errcolor=1,/HIBAR
  ENDIF
  IF N_ELEMENTS(cid_undetect_ALMA) GT 0 THEN BEGIN
    plotsym,1
    oplot,w_obs[cid_undetect_ALMA],3.0*f_obs_err[cid_undetect_ALMA],psym=8,symsize=0.8,color=1
  ENDIF
  xyouts,xrange[1]-0.17*(xrange[1]-xrange[0]),yrange[1]/10^(1.0),name_xy,charthick=3,charsize=0.8, align=1
  xyouts,xrange[1]-0.17*(xrange[1]-xrange[0]),yrange[1]/10^(1.0+0.9),TeXtoIDL("z=")+STRTRIM(STRING(z,FORMAT='(F0.4)'),2),charthick=3,charsize=0.8, align=1
  xyouts,xrange[1]-0.17*(xrange[1]-xrange[0]),yrange[1]/10^(1.0+0.9+0.9),TeXtoIDL("\chi^{2}=")+STRTRIM(STRING(chi2,FORMAT='(G10)'),2),charthick=3,charsize=0.8, align=1
  
  ; Plot the flux residual panel
  yrange = [-1.0,1.0]
  plot,POSITION=positionU2,[0.07,2500],[0,0],lines=1,/xlog,xrange=xrange,yrange=yrange,xtitle=xtitle,xstyle=1,ystyle=1,xtickformat='dzliu_logtickformat',xminor=9,thick=3,xthick=3,ythick=3,xticklen=0.1,yticks=2,yminor=5,charthick=3,charsize=1.75
  w_flag = WHERE(f_res GE yrange[0] AND f_res LE yrange[1], /NULL)
  IF N_ELEMENTS(w_flag) GT 0 THEN BEGIN
    plotsym,8,0.9,/fill
    oploterror,w_obs[w_flag],f_res[w_flag],lambda_err,(f_res-f_res_lo),psym=8,symsize=0.5,color=0,/LOBAR
    oploterror,w_obs[w_flag],f_res[w_flag],lambda_err,(f_res_hi-f_res),psym=8,symsize=0.5,color=0,/HIBAR
  ENDIF
  w_flag = WHERE(f_res LT yrange[0], /NULL)
  IF N_ELEMENTS(w_flag) GT 0 THEN BEGIN
    plotsym,2,0.9,/fill ; upward arrow
    oplot,w_obs[w_flag],f_res[w_flag]*0+yrange[0],psym=8,symsize=1.5,thick=3,color=0
  ENDIF
  w_flag = WHERE(f_res GT yrange[1], /NULL)
  IF N_ELEMENTS(w_flag) GT 0 THEN BEGIN
    plotsym,1,0.9,/fill ; downward arrow
    oplot,w_obs[w_flag],f_res[w_flag]*0+yrange[1],psym=8,symsize=1.5,thick=3,color=0
  ENDIF
  xyouts,0.05,(0.55+0.65)/2.0,/NORMAL,TeXtoIDL("log(f_{OBS}/f_{SED})"),charthick=3,charsize=0.8,align=0.5,orient=90 ; ytitle
  
  ; Histograms 
  plot,f_mu[0,*],f_mu[1,*],psym=10,xrange=xrange01,xtitle=xtitle01,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,ytitle='Likehood Distr.',POSITION=position01
  plot,tau_V_ISM[0,*],tau_V_ISM[1,*],psym=10,xrange=xrange09,xtitle=xtitle09,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position03,YTICKFORMAT="(A1)"
  plot,tau_V[0,*],tau_V[1,*],psym=10,xrange=xrange04,xtitle=xtitle04,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position02,YTICKFORMAT="(A1)"
  plot,Ldust[0,*],Ldust[1,*],psym=10,xrange=xrange06,xtitle=xtitle06,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position07,ytitle='Likehood Distr.'
  ;<dzliu>;plot,Tc_ISM[0,*],Tc_ISM[1,*],psym=10,xrange=[15.1,28.99],xtitle=xtitle07,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position09,YTICKFORMAT="(A1)"
  ;<dzliu>;plot,Tw_BC[0,*],Tw_BC[1,*],psym=10,xrange=[30.5,64.5],xtitle=xtitle08,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position10,YTICKFORMAT="(A1)"
  ;<dzliu>;plot,Tdust[0,*],Tdust[1,*],psym=10,xrange=[15.1,79.9],xtitle=xtitle08,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position10,YTICKFORMAT="(A1)"
  plot,Tc_ISM[0,*],Tc_ISM[1,*],psym=10,xrange=xrange07,xtitle=xtitle07,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position09,YTICKFORMAT="(A1)"
  plot,Tw_BC[0,*],Tw_BC[1,*],psym=10,xrange=xrange08,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position10,YTICKFORMAT="(A1)"
  plot,Mdust[0,*],Mdust[1,*],psym=10,xrange=xrange10,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position08,YTICKFORMAT="(A1)"
  plot,Mstars[0,*],Mstars[1,*],psym=10,xrange=xrange05,xtitle=xtitle05,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position04,YTICKFORMAT="(A1)"
  plot,sSFR[0,*],sSFR[1,*],psym=10,xrange=xrange11,xtitle=xtitle11,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,xminor=5,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position05,YTICKFORMAT="(A1)"
  plot,SFR[0,*],SFR[1,*],psym=10,xrange=xrange12,xtitle=xtitle12,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position06,YTICKFORMAT="(A1)"
  plot,xi_C_tot[0,*],xi_C_tot[1,*],psym=10,xrange=xrange13,xtitle=xtitle13,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position11,YTICKFORMAT="(A1)"
  plot,xi_W_tot[0,*],xi_W_tot[1,*],psym=10,xrange=xrange14,xtitle=xtitle14,xstyle=1,thick=5,xthick=3,ythick=3,xticklen=0.1,ystyle=1,yrange=[0,1.18],charthick=3,charsize=1.5,xcharsize=0.9,POSITION=position12,YTICKFORMAT="(A1)"
  
  xyouts, (position08[0]+position08[3])/2.0, (position08[1]-0.01), xtitle08, /NORMAL, ALIGNMENT=0.5, CHARSIZE=1.0
  xyouts, (position10[0]+position10[3])/2.0, (position10[1]-0.01), xtitle10, /NORMAL, ALIGNMENT=0.5, CHARSIZE=1.0
  
  
  device,/close
  set_plot,'x'

close,1,20
;endfor
end
 
