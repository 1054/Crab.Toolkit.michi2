forward_function TeXtoIDL ; added by dzliu -- see www.idlcoyote.com/tips/func_var.html

function dzliu_lumtoflux, input_z, input_lambda_um, input_lum_Lsun
  dzliu_lumdist = lumdist(input_z)
  dzliu_lumsang = 4*!PI*dzliu_lumdist^2
  dzliu_flux_mJy = (1.+input_z)*input_lum_Lsun/dzliu_lumsang*40.31970/(2.99792458e5/input_lambda_um) ; 1 Lsun Mpc-2 = 40.31970 mJy GHz
  return, dzliu_flux_mJy
end

; Usage: idl84 -e "plot_sed" -args 10001 ; added by dzliu

pro plot_sed, galaxy

;This code was made to interpret grafically the Spectral Energy Distribution
;for a galaxy, fitted with the code of E. da Cunha (2008).
;It takes the imputs from the DATA file of the code and creates a .ps file for 
;each galaxy that includes the theoretical SED of this galaxy, oveploted the data 
;points and also the probability distribution functions for some important parameters.
;T.Bitsakis(17-9-2010)

; version Oct. 2011
; bugs fixed by Kate Rownlands, Univ. Nottigham

;================================================INPUTS===============================================
 
;  galaxy='' ; commmented by dzliu, now read from input argument
;  readcol,'DATA.dat',f=(A,F),galaxy,nn  ;Reads the data file to find the number of the elements as well the names of galaxies
;  n=n_elements(nn)
 if n_elements((command_line_args())) gt 0 then begin ; added by dzliu
    galaxy=(command_line_args())[0] ; added by dzliu
 endif ; added by dzliu
 if n_elements(galaxy) eq 0 then begin
    read,galaxy,prompt='Give the name of the source: ' 
 endif
 
 if size(galaxy,/tname) ne 'STRING' then galaxy=STRING(format='(I0)',galaxy) ; added by dzliu
 if file_test('textoidl',/dir) then !PATH=expand_path('+'+(file_search('textoidl',/full))[0])+':'+!PATH ; added by dzliu
 if file_test('coyote',/dir) then !PATH=expand_path('+'+(file_search('coyote',/full))[0])+':'+!PATH ; added by dzliu
 resolve_all, /quiet ; added by dzliu
 print, 'Output to '+galaxy+'.ps'

 fmt='A,F,I,I'
  filt=''
  ;readcol,'/data/kate/MAGPHYS/magphys/eg_user_files/filters_HATLAS.dat',f=fmt,filt,lambda_eff,filter_id,fit  ;Reads the filter file
  if ~strlen(GETENV('USER_FILTERS')) then message, 'Error! GETENV(USER_FILTERS) was not defined!' ; added by dzliu
  if ~file_test(GETENV('USER_FILTERS'),/read) then message, 'Error! '+GETENV('USER_FILTERS')+' was not found!' ; added by dzliu
  readcol,GETENV('USER_FILTERS'),f=fmt,filt,lambda_eff,filter_id,fit  ;Reads the filter file, modified by dzliu 
  m=n_elements(lambda_eff)

;for jj=1,n do begin
   name_ps=galaxy+'.ps'
   name_fit=galaxy+'.fit'
   name_sed=galaxy+'.sed'
   name_xy=galaxy
   
   if not file_test(name_fit,/read) then message, 'Error! '+name_fit+' was not found!' ; added by dzliu
   if not file_test(name_sed,/read) then message, 'Error! '+name_sed+' was not found!' ; added by dzliu
   
   oti1=''
   oti2=''
   oti3=''
   flux=fltarr(m)
   e_flux=fltarr(m)
   p_flux=fltarr(m)

   openr,1,name_fit
   for i=1,2 do readf,1,oti1
   readf,1,flux
   readf,1,e_flux
   for i=1,4 do readf,1,oti2
   readf,1,i_sfh,i_ir,chi2,z
   for i=1,3 do readf,1,oti3
   readf,1,p_flux

   readcol,name_sed,lambda,lum_at,lum_un,SKIPLINE =10

;==============================================CALCULATIONS===========================================
  x=10^lambda ;wavelength
  y_at=alog10(x*10^lum_at) ;total attenuated SED
  y_un=alog10(x*10^lum_un) ;unattenuated SED
  xx=x/1.e+4 ;lambda in microns
  lambda_err=lambda_eff*0
  device,decomposed=0

; Convert the fluxes in the different photometric bands into luminosities lambda*L
; OBSERVED FLUXES
  L_flux=alog10((1.+z)*flux*3e+14/lambda_eff)
  
; ERRORS
  ;L_eflux=alog10((1.+z)*flux*3e+14/lambda_eff)-alog10((1.+z)*flux*3e+14/lambda_eff-(1.+z)*e_flux*3e+14/lambda_eff)

  L_eflux_lo=alog10((1.+z)*flux*3e+14/lambda_eff)-alog10((1.+z)*flux*3e+14/lambda_eff-e_flux*(1.+z)*3e+14/lambda_eff)

  ;for error bars which go down to infinity (i.e. error is bigger than flux)
  loc = [where(~finite(L_eflux_lo))]
  if (loc[0] gt -0.5) then L_eflux_lo[where(~finite(L_eflux_lo))] = L_flux[where(~finite(L_eflux_lo))]

  L_eflux_hi=-alog10((1.+z)*flux*3e+14/lambda_eff)+alog10((1.+z)*flux*3e+14/lambda_eff+e_flux*(1.+z)*3e+14/lambda_eff)

; PREDICTED LUMINOSITIES
  L_pflux=alog10((1.+z)*flux*3e+14/lambda_eff)-alog10((1.+z)*p_flux*3e+14/lambda_eff) 

;================================================HISTOGRAM============================================
f_muSFH=dblarr(2,20)
f_muIR=dblarr(2,20)
mu=dblarr(2,20)
tau_V=dblarr(2,48)
sSFR=dblarr(2,70)
Mstars=dblarr(2,60)
Ldust=dblarr(2,60)
Tc_ISM=dblarr(2,10)
Tw_BC=dblarr(2,30)
tau_V_ISM=dblarr(2,80)
Mdust=dblarr(2,60)
SFR=dblarr(2,60)
xi_C_tot=dblarr(2,20)
xi_W_tot=dblarr(2,20)
skip1=''
skip2=''
skip3=''
skip4=''
skip5=''
skip6=''
skip7=''
skip8=''
skip9=''
skip10=''
skip11=''
skip12=''
skip13=''
skip14=''
skip15=''

name=''
openr,20,name_fit

for j=1,16 do begin
readf,20,skip1
endfor
readf,20,f_muSFH
for j=1,3 do begin
readf,20,skip2
endfor
readf,20,f_muIR
f_mu=0.5*(f_muIR+f_muSFH)
for j=1,3 do begin
readf,20,skip3
endfor
readf,20,mu
for j=1,3 do begin
readf,20,skip4
endfor
readf,20,tau_V
for j=1,3 do begin
readf,20,skip5
endfor
readf,20,sSFR
for j=1,3 do begin
readf,20,skip6
endfor
readf,20,Mstars
for j=1,3 do begin
readf,20,skip7
endfor
readf,20,Ldust
for j=1,3 do begin
readf,20,skip8
endfor
readf,20,Tc_ISM
for j=1,3 do begin
readf,20,skip9
endfor
readf,20,Tw_BC
for j=1,3 do begin
readf,20,skip10
endfor
readf,20,xi_C_tot
for j=1,49 do begin
readf,20,skip11
endfor
readf,20,xi_W_tot
for j=1,3 do begin
readf,20, skip12
endfor
readf,20,tau_V_ISM
for j=1,3 do begin
readf,20,skip13
endfor
readf,20,Mdust
for j=1,3 do begin
readf,20,skip14
endfor
readf,20,SFR
for j=1,2 do begin
readf,20,skip15
endfor

t1=TeXtoIDL('f_\mu')
t3=TeXtoIDL('\mu')
t4=TeXtoIDL('\tau_V')
t5=TeXtoIDL('log(M_{stars}/M_{o})')
t6=TeXtoIDL('log(L_{dust}/L_{o})')
t7=TeXtoIDL('T_{C}^{ISM}/K')
t8=TeXtoIDL('T_{W}^{BC}/K')
t9=TeXtoIDL('\mu\tau_V')
t10=TeXtoIDL('log(M_{dust}/M_{o})')
t11=TeXtoIDL('log(sSFR) yr^{-1}')
t12=TeXtoIDL('log(SFR/M_{o} yr^{-1})')
t13=TeXtoIDL('\xi_C^{tot}')
t14=TeXtoIDL('\xi_W^{tot}')


;================================================PLOTTING=============================================
  yt=TeXtoIDL("log(\lambdaL_{\lambda}/L"+sunsymbol()+")")
  xt=TextoIDL("\lambda/\mum [observed-frame]")
  chi_t=TeXtoIDL("\chi^{2}=")
  plotsym,8,0.9,/fill

  set_plot, 'ps'
  !P.MULTI=[0,0,20]
  TVLCT, 0, 0, 0 , 0
  TVLCT, 200,55, 70 , 1
  TVLCT,  0,255, 0 , 2
  TVLCT, 0, 150, 255 , 3
  TVLCT, 100, 100, 100 , 4

; positions
; 1st row
position1=[0.15,0.3,0.275,0.45]
position2=[0.275,0.3,0.40,0.45]
position3=[0.40,0.3,0.525,0.45]
position4=[0.525,0.3,0.65,0.45]
position5=[0.65,0.3,0.775,0.45]
position6=[0.775,0.3,0.9,0.45]
; 2nd row
position7=[0.15,0.07,0.275,0.22]
position8=[0.275,0.07,0.40,0.22]
position9=[0.40,0.07,0.525,0.22]
position10=[0.525,0.07,0.65,0.22]
position11=[0.65,0.07,0.775,0.22]
position12=[0.775,0.07,0.90,0.22]

  device,filename=name_ps,/color

  ; dzliu Y axis flux mJy
  xr = [0.07,1500]
  yr = [7.1,14.]
  if 1 then begin ; -- dzliu -- TODO -- plot Y axis as flux in mJy
  y_at = alog10(dzliu_lumtoflux(z,lambda,10^y_at/lambda)) ; -- the dust attenuated SED
  y_un = alog10(dzliu_lumtoflux(z,lambda,10^y_un/lambda)) ; -- the dust unattenuated SED
  L_pflux = alog10(dzliu_lumtoflux(z,lambda_eff,10^(L_flux-L_pflux)/lambda_eff)) ; must put before rewritting L_flux
  L_flux = alog10(dzliu_lumtoflux(z,lambda_eff,10^L_flux/lambda_eff))
  L_pflux = L_flux - L_pflux
  yt = TeXtoIDL("log_{10} flux [mJy]")
  yr = [-6,4]
  endif
  
  ; dzliu
  ;plot,xx,y_at,/xlog,xr=[0.07,1500],yr=[7.1,12.],ytitle=yt,xstyle=1,ystyle=1,thick=2,POSITION=[0.12,0.65,0.92,0.99],charsize=1.75,XTICKFORMAT="(A1)"
  plot,xx,y_at,/xlog,xr=xr,yr=yr,ytitle=yt,xstyle=1,ystyle=1,thick=2,POSITION=[0.12,0.65,0.92,0.99],charsize=1.75,XTICKFORMAT="(A1)"
  oplot,xx,y_un,color=3
  oploterror,lambda_eff,L_flux,lambda_err,L_eflux_lo,psym=8,symsize=0.5,color=1,errcolor=1, /LOBAR
  oploterror,lambda_eff,L_flux,lambda_err,L_eflux_hi,psym=8,symsize=0.5,color=1,errcolor=1, /HIBAR
   xyouts,500,yr[1]-1.0,name_xy,charthick=2,charsize=0.8
   xyouts,500,yr[1]-1.7,chi_t,charsize=0.8
   xyouts,700,yr[1]-1.7,STRTRIM(STRING(chi2, format='(F5.2)'),1),charsize=0.8
  
  ; dzliu
  ;plot,[0.07,1500],[0,0],lines=1,xr=[0.07,1500],yr=[-1.0,1.0],xtitle=xt,xstyle=1,ystyle=1,/xlog, $
  plot,[0.07,1500],[0,0],lines=1,xr=xr,yr=[-1.0,1.0],xtitle=xt,xstyle=1,ystyle=1,/xlog, $  
       POSITION=[0.12,0.55,0.92,0.65],yticks=2,yminor=5,xticklen=0.1,charsize=1.75
  oploterror,lambda_eff,L_pflux,lambda_err,L_eflux_lo,psym=8,symsize=0.5,color=0, /LOBAR
  oploterror,lambda_eff,L_pflux,lambda_err,L_eflux_hi,psym=8,symsize=0.5,color=0, /HIBAR

;Histograms 
  
   plot,f_mu[0,*],f_mu[1,*],psym=10,xtitle=t1,thick=5,ystyle=1,yr=[0,1],xr=[0.01,0.98],xstyle=1,charsize=1.3,ytitle='Likehood Distr.',POSITION=position1
   plot,tau_V_ISM[0,*],tau_V_ISM[1,*],xr=[0.01,1.95],psym=10,xtitle=t9,xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position3,YTICKFORMAT="(A1)"
   plot,tau_V[0,*],tau_V[1,*],psym=10,xtitle=t4,thick=5,ystyle=1,yr=[0,1],charsize=1.3,xr=[0.01,4.8],xstyle=1,POSITION=position2,YTICKFORMAT="(A1)"
   plot,Ldust[0,*],Ldust[1,*],xr=[8.01,11.99],psym=10,xtitle=t6,xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position7,ytitle='Likehood Distr.'
   plot,Tc_ISM[0,*],Tc_ISM[1,*],xr=[15.1,24.99],psym=10,xtitle=t7,xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position9,YTICKFORMAT="(A1)"

   plot,Tw_BC[0,*],Tw_BC[1,*],psym=10,xtitle=t8,thick=5,ystyle=1,yr=[0,1],xr=[30.5,59.5],xstyle=1,charsize=1.3,POSITION=position10,YTICKFORMAT="(A1)"

   plot,Mdust[0,*],Mdust[1,*],psym=10,xtitle=t10,xr=[4.1,8.9],xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position8,YTICKFORMAT="(A1)"
   plot,Mstars[0,*],Mstars[1,*],psym=10,xtitle=t5,xr=[8.1,11.9],xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position4,YTICKFORMAT="(A1)"
   plot,sSFR[0,*],sSFR[1,*],psym=10,xtitle=t11,xr=[-12.9,-8.1],xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position5,YTICKFORMAT="(A1)"
   plot,SFR[0,*],SFR[1,*],psym=10,xtitle=t12,thick=5,ystyle=1,yr=[0,1],xr=[-3.9,4],xstyle=1,charsize=1.3,POSITION=position6,YTICKFORMAT="(A1)"

   plot,xi_C_tot[0,*],xi_C_tot[1,*],xr=[0.01,0.99],psym=10,xtitle=t13,xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position11,YTICKFORMAT="(A1)"
   plot,xi_W_tot[0,*],xi_W_tot[1,*],xr=[0.01,0.99],psym=10,xtitle=t14,xstyle=1,thick=5,ystyle=1,yr=[0,1],charsize=1.3,POSITION=position12,YTICKFORMAT="(A1)"


  device,/close
  set_plot,'x'

close,1,20
;endfor
end
 
