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
; 20180618 copied "plot_sed.pro" to this code "convert_sed.pro" to convert the units and save as text file.
; 
; Usage: idl84 -e "convert_sed" -args "best-fit_SED" ; added by dzliu
; 






; Dependent functions

;forward_function TeXtoIDL ; added by dzliu -- see www.idlcoyote.com/tips/func_var.html

function dzliu_lumtoflux, input_z, input_lambda_um, input_lum_Lsun
  dzliu_lumdist = lumdist(input_z)
  dzliu_lumdisq = 4.D * !PI*dzliu_lumdist^2
  ;dzliu_flux_mJy = (1.+input_z)*input_lum_Lsun/dzliu_lumdisq*40.31970d/(2.99792458e5/input_lambda_um) ; 1 Lsun Mpc-2 = 40.31970 mJy GHz
  dzliu_flux_mJy = input_lum_Lsun/dzliu_lumdisq*40.31970D/(2.99792458D5/input_lambda_um) ; 1 Lsun Mpc-2 = 40.31970 mJy GHz
  return, dzliu_flux_mJy
end






; MAIN PROGRAM

pro convert_sed, galaxy
  
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
  if size(galaxy,/tname) ne 'STRING' then galaxy=STRING(format='(I0)',galaxy) ; added by dzliu
  ;if ~strlen(GETENV('magphys')) then message, 'Error! GETENV(magphys) was not defined!' ; added by dzliu
  ;if file_test(GETENV('magphys')+path_sep()+'textoidl',/dir) then !PATH=expand_path('+'+GETENV('magphys')+path_sep()+'textoidl')+':'+!PATH ; added by dzliu
  ;if file_test(GETENV('magphys')+path_sep()+'coyote',/dir) then !PATH=expand_path('+'+GETENV('magphys')+path_sep()+'coyote')+':'+!PATH ; added by dzliu
  resolve_all, /quiet ; added by dzliu
  print, 'Output to '+galaxy+'.sed.um.mJy.txt'
  
  ; Read USER_FILTERS
  fmt='A,F,I,I'
  filt=''
  USER_FILTERS='magphys_input_filters.dat'
  if ~strlen(USER_FILTERS) then message, 'Error! USER_FILTERS was not defined!'               ; added by dzliu
  if ~file_test(USER_FILTERS,/read) then message, 'Error! '+USER_FILTERS+' was not found!'  ; added by dzliu
  readcol,USER_FILTERS,f=fmt,filt,lambda_eff,filter_id,fit                                            ; Reads the filter file, lambda_eff are in units of um
  m = n_elements(lambda_eff)
  
  ; Read MAGPHYS *.fit data file
  name_sav=galaxy+'.sav'
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
  vLv_obs    = (1.D + z)*(Lv_obs       )*2.99792458D8/(w_obs/1D6)  ; notes by dzliu: vLv_obs are in units of L_{\odot}; (1.+z) is needed according to the original code of plot_sed.pro and the fit_sed_highz.f line 275.
  vLv_fit    = (1.D + z)*(Lv_fit       )*2.99792458D8/(w_obs/1D6)
  vLv_obs_lo = (1.D + z)*(Lv_obs-Lv_err)*2.99792458D8/(w_obs/1D6)
  vLv_obs_hi = (1.D + z)*(Lv_obs+Lv_err)*2.99792458D8/(w_obs/1D6)
  
  ; For error bars which go down to infinity (i.e. error is bigger than flux)
  IF N_ELEMENTS(where(vLv_obs_lo LE 1D-99, /NULL)) GT 0 THEN BEGIN
      vLv_obs_lo[where(vLv_obs_lo LE 1D-99, /NULL)] = 1D-99
  ENDIF
  
  ; Compute the difference/residual between observed vLv and best fit vLv
  vLv_res = alog10(vLv_obs) - alog10(vLv_fit)
  
  ; Compute flux density from vLv
  f_sed_at = dzliu_lumtoflux(z,w_sed,vLv_sed_at) ; -- the best fit dust attenuated SED flux
  f_sed_un = dzliu_lumtoflux(z,w_sed,vLv_sed_un) ; -- the best fit dust unattenuated SED flux
  f_fit    = dzliu_lumtoflux(z,w_obs,vLv_fit)    ; -- the best fit flux
  f_obs    = dzliu_lumtoflux(z,w_obs,vLv_obs)    ; -- the observed flux
  f_obs_lo = dzliu_lumtoflux(z,w_obs,vLv_obs_lo) ; -- the observed flux lowest value for 1-sigma confidence
  f_obs_hi = dzliu_lumtoflux(z,w_obs,vLv_obs_hi) ; -- the observed flux highest value for 1-sigma confidence
  
  ; For values below 0, we set them to 1D-30
  IF N_ELEMENTS(where(f_fit LE 1D-99, /NULL)) GT 0 THEN BEGIN
      f_fit[where(f_fit LE 1D-99, /NULL)] = 1D-99
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
  
  
  
  
  
  
  
  save, FILENAME=name_sav, w_sed, f_sed_at, f_sed_un, f_obs, f_obs_lo, f_obs_hi, f_fit, f_res, f_res_lo, f_res_hi, Mstars, SFR, sSFR, Ldust, Mdust, Tdust, age_M, A_V, tau_V, tau_V_ISM, mu, Tc_ISM, Tw_BC
  
  openw, lun, name_sed_out, /get_lun
  printf, lun, "wave_um", "f_attenu_mJy", "f_unattenu_mJy", "vLv_attenu_Lsun", "vLv_unattenu_Lsun", format='("# ",A-14," ",A16," ",A16," ",A16," ",A16)'
  for i=0,n_elements(w_sed)-1 do begin
    if i gt 0 then begin
      ;print, 'Checking non-monochromatic wavelengths:', string(format='(E-14.6)', w_sed[i]), string(format='(E-14.6)', w_sed[i-1])
      if string(format='(E-14.6)', w_sed[i]) eq string(format='(E-14.6)', w_sed[i-1]) then begin
        continue
      endif
    endif
    printf, lun, w_sed[i], f_sed_at[i], f_sed_un[i], vLv_sed_at[i], vLv_sed_un[i], format='("  ",E-14.6," ",E16.6," ",E16.6," ",E16.6," ",E16.6)'
  endfor
  close, lun
  free_lun, lun
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  close,1,20
end
 
