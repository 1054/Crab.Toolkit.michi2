
PRO makelibSED_FSPS
    ; 
    ; speclist = FILE_SEARCH('/Users/dliu/Programming/sps/fsps/fsps_ssp/SSP_Padova_BaSeL_Chabrier_allZ/SSP_Padova_BaSeL_Chabrier_Z*.out.spec')
    ; 
    ; outfile = '/Users/dliu/Working/SpireLines/Tool/Level_3_SciData/07_SED_Synthesis/02_Make_Lib_SED/FSPS.Padova.BaSeL.lib.SED'
    ;
    ;speclist = FILE_SEARCH('/Users/dliu/Programming/sps/fsps/fsps_ssp/SSP_Padova_BaSeL_Chabrier_allZ/SSP_Padova_BaSeL_Chabrier_Z0.0190.out.spec')
    speclist = FILE_SEARCH('FSPSspec/SSP_Padova_BaSeL_Chabrier_Z0.0190.out.spec')
    ;
    ;outfile = '/Users/dliu/Working/SpireLines/Tool/Level_3_SciData/07_SED_Synthesis/02_Make_Lib_SED/FSPS.Padova.BaSeL.Z0.0190.lib.SED'
    outfile = 'outputs/FSPS.Padova.BaSeL.Z0.0190.lib.SED'
    ; 
    spawn, 'mkdir -p outputs 2>/dev/null'
    ; 
    OPENW, fp, outfile, /GET_LUN
    PRINTF, fp, FORMAT='(A,A)', "# ", CrabStringCurrentTime()
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR  = 2 # Wave & Flux"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR1 = band wavelength (rest-frame) (lambda) [um]"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR2 = flux density in unit of solar luminosity per Hz"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR1 = 1 # the colomn number of VAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR2 = 2 # the colomn number of VAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR1 = 1963 # NAXIS1 -- the number of wavelength in one spectrum"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR2 = 1    # NAXIS2 -- the number of spectra -- should equal to NPAR1*NPAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NPAR  = 3 # Metallicity, Age & Mass"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR1 = Absolute Metallicity (note that solar metallicity is 0.0190)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR2 = Age in unit of Gyr"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR3 = Mass in unit of solar mass (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR1 = 3 # the colomn number of PAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR2 = 4 # the colomn number of PAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR3 = 5 # the colomn number of PAR2"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR1 = ", N_ELEMENTS(speclist), " # Number of different PAR1"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR2 = ", ((188-1-0)/11)+1, " # Number of different PAR2"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR3 = ", ((188-1-0)/11)+1, " # Number of different PAR3"
    PRINTF, fp, FORMAT='(A)', "# "
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A12)', "# ", "lambda", "flux", "metallicity", "age", "mass"
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A12)', "# ", "[um]", "[Lsun/Hz]", " ", "[Gyr]", "[Msun]"
    
    tmpVar = 0
    
    ; 
    FOREACH specfile, speclist DO BEGIN
        
        TmpPos1 = STRPOS(specfile,'_Chabrier_Z')+STRLEN('_Chabrier_Z')
        TmpPos2 = STRPOS(specfile,'.out.spec')
        ZMeta = STRMID(specfile,TmpPos1,TmpPos2-TmpPos1)
        
        PRINT, specfile, " ", ZMeta
        
        libSED = read_spec1(specfile)
        Nlambda = (SIZE((libSED.lambda),/dim))[0]
        NageGyr = (SIZE((libSED.lambda),/dim))[1]
        ; PRINT, FORMAT='(A,I0)', "Nlambda = ", Nlambda
        ; PRINT, FORMAT='(A,I0)', "NageGyr = ", NageGyr ; <TODO> SHOULD BE 188 !!! THIS DEFINES NPAR2 AND NPAR3 !!!
        
        ; match wave unit to [um]
        libSED.lambda = libSED.lambda * 1e-4
        
        FOR iage=0,NageGyr-1,11 DO BEGIN
            
            ; PRINT, FORMAT='(A,I-5,A,F-12.5,A,F-8.5,A,F0.5,A,E0.5)', "index = ", iage, " Z = ", ZMeta, " age = ", (libSED.ageGyr)[iage], " Gyr     Mass = ", 10^(libSED.logmass)[iage], " Msun    SFR = ", 10^(libSED.logsfr)[iage]
            
            wave = (libSED.lambda)[*,iage] ; [um]
            flux = (libSED.spec)[*,iage] ; [L_solar/Hz]
            
            ; PRINT, N_ELEMENTS(wave), N_ELEMENTS(flux)
            
            
            FOR i=0,N_ELEMENTS(flux)-1 DO BEGIN
                PRINTF, fp, FORMAT='(E15.6,E15.6,F15.5,F12.5,F12.5)', wave[i], flux[i], DOUBLE(ZMeta), (libSED.ageGyr)[iage], 10^(libSED.logmass)[iage]
            ENDFOR
            
            
            
        ENDFOR
        
    ENDFOREACH
    
    CLOSE, fp & FREE_LUN, fp
    
    PRINT, 'Done!'
    
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

