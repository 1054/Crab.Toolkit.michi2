;+
; NAME:
;       READ_SPEC 
;
; PURPOSE:
;      Function to convert a *.spec file into an IDL structure.  
;      Input may be an array of filenames.
;
; CALLING SEQUENCE:
;   res = read_spec(file)
;
; KEYWORD PARAMETERS:
;   miles - set this keyword when reading in spectra based on Miles 
;   pickles - set this keyword when reading in spectra based on Pickles
;   NB: you should not need to set these keywords if you are running
;       the latest version of fsps
;
; MODIFICATION HISTORY: 
;   ? - created by CFC
;   09/09/11 - Updated to read in files that have arbitrary wavelength
;              arrays and where the first full line is the wavelength array
;
;-
;-----------------------------------------------------------------;
PRO PLOT_SPEC
    ;spec = read_spec1('/Users/dliu/Programming/sps/fsps/fsps_ssp/SSP_Padova_BaSeL_Chabrier_allZ/SSP_Padova_BaSeL_Chabrier_Z0.0002.out.spec')
    spec = read_spec1('/Users/dliu/Programming/sps/fsps/fsps_ssp/SSP_Padova_BaSeL_Chabrier_allZ/SSP_Padova_BaSeL_Chabrier_Z0.0025.out.spec')
    Nlambda = (size((spec.lambda),/dim))[0]
    NageGyr = (size((spec.lambda),/dim))[1]
    print, FORMAT='(A,I0)', "Nlambda = ", Nlambda
    print, FORMAT='(A,I0)', "NageGyr = ", NageGyr
    xrange=[3e1,1e8]*1e-4 & yrange=[1e-21,2e-12]
    xtitle="wave [um]" & ytitle="flux density [Lsun/Hz]" ; xtitle="wave [angstroms]"
    iage = 0
    p = plot((spec.lambda)[*,iage]*1e-4, (spec.spec)[*,iage], /xlog, /ylog, /nodata, xtitle=xtitle, ytitle=ytitle, font_name='NGC', xtickint=1, xstyle=1, ystyle=1, xrange=xrange, yrange=yrange)
    FOR iage=7,NageGyr-1,12 DO BEGIN
        vcolor = 10+iage*(250.0/(NageGyr+10.0))
        print, FORMAT='(A,I-5,A,F-8.5,A,F0.5,A,E0.5,A,I0)', "index = ", iage, " age = ", (spec.agegyr)[iage], " Gyr     mass = ", 10^(spec.logmass)[iage], "    sfr = ", 10^(spec.logsfr)[iage], "    color = ", vcolor
        p = plot((spec.lambda)[*,iage]*1e-4, (spec.spec)[*,iage], linestyle='solid', rgb_table=25, vert_colors=vcolor, /overplot, xrange=xrange, yrange=yrange)
    ENDFOR
    
END



FUNCTION READ_SPEC1, file

  openr,lun,file,/get_lun

  ;burn the header
  char = '#'
  WHILE strmid(char,0,1) EQ '#' DO BEGIN
     readf,lun,char
  ENDWHILE 

  ;check if the spec file is of the "new" type, where both 
  ;the number of age steps and the number of spectral elements are included
  char = strsplit(char,' ',/extr)
  IF n_elements(char) GT 1 THEN BEGIN
     nt = long(char[0])
     nl = long(char[1]) 
  ENDIF ELSE BEGIN
     print,'ERROR: you are not passing a properly formatted *spec file'
     STOP
  ENDELSE

  str  = {agegyr:0.0,logmass:0.0,loglbol:0.0,logsfr:0.0,spec:fltarr(nl),$
          lambda:fltarr(nl)}
  str   = replicate(str,nt)
  tspec = fltarr(nl)
  t = 0.
  m = 0.
  l = 0.
  s = 0.

  ;if the number of spectral elements is passed, then the first
  ;line here is the wavelength array.
  IF n_elements(char) GT 1 THEN BEGIN
     readf,lun,tspec
     lam = tspec
  ENDIF

  FOR i=0,nt-1 DO BEGIN
     
     readf,lun,t,m,l,s
     str[i].agegyr  = 10^t/1E9
     str[i].logmass = m
     str[i].loglbol = l
     str[i].logsfr  = s

     readf,lun,tspec
     str[i].spec   = tspec
     str[i].lambda = lam
     
  ENDFOR
  
  close,lun
  free_lun,lun

  RETURN,str

END

;------------------------------------------------------------;
;------------------------------------------------------------;

FUNCTION READ_SPEC, file

  ct = n_elements(file)

  ; spsdir = getenv('SPS_HOME')
  spsdir = ''
  IF spsdir EQ '' THEN BEGIN
     print,'READ_SPEC ERROR: spsdir environment '+$
           'variable not set, returning...'
     return,0
  ENDIF
  spsdir = spsdir+'/SPECTRA/'

     
  ff = findfile(file[0],count=ctf)
  IF ctf EQ 0 THEN BEGIN
     print,'READ_SPEC ERROR: file not found: ',file
     return,0
  ENDIF

  spec = strpos(file[0],'spec')
  IF spec EQ -1 THEN BEGIN
     print,'READ_SPEC ERROR: you did not pass a .spec file, returning'
     return,0
  ENDIF

  astr = read_spec1(file[0])

  IF ct GT 1 THEN BEGIN
     str = replicate(astr[0],n_elements(astr),ct)
     str[*,0] = astr
     FOR i=1,ct-1 DO str[*,i] = read_spec1(file[i])
  ENDIF ELSE str = astr
  
  RETURN, str

END

