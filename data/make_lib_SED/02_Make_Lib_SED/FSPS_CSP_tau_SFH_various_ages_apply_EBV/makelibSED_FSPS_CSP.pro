
PRO makelibSED_FSPS_CSP
    ; 
    ; 
    speclist = FILE_SEARCH('../../01_SED_Models/FSPS/fsps_csp/Z_EQ_1.0_Zsun_tau_EQ_1.0_fburst_EQ_0_tburst_EQ_99/CSP*.out.spec')
    ;
    outfile = 'outputs/FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.lib.SED'
    outdatafile = 'outputs/FSPS.CSP.tau.1Gyr.Padova.BaSeL.Z0.0190.lib.DAT'
    ; 
    spawn, 'mkdir -p outputs 2>/dev/null'
    ; 
    OPENW, fpdatafile, outdatafile, /GET_LUN
    ; 
    Nallpar = 0
    ; 
    ; 
    FOREACH specfile, speclist DO BEGIN
        
        ZMeta = ''
        IF STRPOS(specfile,'_Chabrier_Z') GT 0 THEN BEGIN
            TmpPos1 = STRPOS(specfile,'_Chabrier_Z')+STRLEN('_Chabrier_Z')
            TmpPos2 = STRPOS(specfile,'.out.spec')
            ZMeta = STRMID(specfile,TmpPos1,TmpPos2-TmpPos1)
        ENDIF
        if STRTRIM(ZMeta,2) EQ '' THEN ZMeta = '0.0190'
        
        PRINT, specfile + " (ZMeta = "+ZMeta+")"
        ;break
        
        libSED = read_spec1(specfile)
        Nlambda = (SIZE((libSED.lambda),/dim))[0]
        NageGyr = (SIZE((libSED.lambda),/dim))[1]
        PRINT, FORMAT='(A,I0)', "Nlambda = ", Nlambda
        PRINT, FORMAT='(A,I0)', "NageGyr = ", NageGyr ; <TODO> SHOULD BE 94 !!! THIS DEFINES NPAR2 AND NPAR3 !!!
        
        ; match wave unit to [um]
        libSED.lambda = libSED.lambda * 1e-4
        
        list_iage = []
        FOR iage = 0,NageGyr-1,11 DO BEGIN
            list_iage = [list_iage, iage]
        ENDFOR
        IF list_iage[N_ELEMENTS(list_iage)-1] NE iage THEN BEGIN
            list_iage = [list_iage, NageGyr-1]
        ENDIF
        PRINT, 'Looping ages: ', list_iage
        
        FOREACH iage, list_iage DO BEGIN
            
            ; PRINT, FORMAT='(A,I-5,A,F-12.5,A,F-8.5,A,F0.5,A,E0.5)', "index = ", iage, " Z = ", ZMeta, " age = ", (libSED.ageGyr)[iage], " Gyr     Mass = ", 10^(libSED.logmass)[iage], " Msun    SFR = ", 10^(libSED.logsfr)[iage]
            
            wave = (libSED.lambda)[*,iage] ; [um]
            flux = (libSED.spec)[*,iage] ; [L_solar/Hz]
            
            ; PRINT, N_ELEMENTS(wave), N_ELEMENTS(flux)
            
            
            FOR i=0,N_ELEMENTS(flux)-1 DO BEGIN
                PRINTF, fpdatafile, FORMAT='(E15.6,E15.6,F15.5,E13.5,E13.5,E13.5,E13.5)', wave[i], flux[i], DOUBLE(ZMeta), (libSED.ageGyr)[iage], 10^(libSED.logmass)[iage], 10^(libSED.logsfr)[iage], 10^(libSED.loglbol)[iage]
            ENDFOR
            
            Nallpar = Nallpar + 1
            
        ENDFOREACH
        
    ENDFOREACH
    
    CLOSE, fpdatafile & FREE_LUN, fpdatafile
    
    
    
    
    
    OPENR, fpdatafile, outdatafile, /GET_LUN
    
    POINT_LUN, fpdatafile, 0
    
    OPENW, fp, outfile, /GET_LUN
    
    POINT_LUN, fp, 0
    
    PRINTF, fp, FORMAT='(A,A)', "# ", CrabStringCurrentTime()
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR  = 2 # Wave & Flux"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR1 = band wavelength (rest-frame) (lambda) [um]"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR2 = flux density in unit of solar luminosity per Hz"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR1 = 1 # the colomn number of VAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR2 = 2 # the colomn number of VAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR1 = 1963 # NAXIS1 -- the number of wavelength in one spectrum"
    PRINTF, fp, FORMAT='(A,A,I-4,A)', "# ", "NVAR2 = ", N_ELEMENTS(speclist) * Nallpar, " # NAXIS2 -- the number of spectra -- should equal to NPAR1*NPAR2 (only consider non-independent PAR)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NPAR  = 5    # Metallicity, Age, Mass, SFR and Lbol"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR1 = Absolute Metallicity (note that solar metallicity is 0.0190)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR2 = Age in unit of Gyr"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR3 = Mass in unit of solar mass, Charbier2003 IMF (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR4 = SFR in unit of solar mass per year, Charbier2003 IMF (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR5 = Lbol in unit of solar luminosity (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR1 = 3 # the colomn number of PAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR2 = 4 # the colomn number of PAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR3 = 5 # the colomn number of PAR3"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR4 = 6 # the colomn number of PAR4"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR5 = 7 # the colomn number of PAR5"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR1 = ", N_ELEMENTS(speclist), " # Number of different PAR1"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR2 = ", Nallpar, " # Number of different PAR2"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR3 = ",  1, " # Number of different PAR3 (non-independent)"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR4 = ",  1, " # Number of different PAR4 (non-independent)"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR5 = ",  1, " # Number of different PAR5 (non-independent)"
    PRINTF, fp, FORMAT='(A)', "# "
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A13,A13,A13,A13)', "# ", "lambda",    "flux",  "metallicity",   "age",   "mass",       "SFR",   "Lbol"
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A13,A13,A13,A13)', "# ", "[um]", "[Lsun/Hz]", "[12+lg(O/H)]", "[Gyr]", "[Msun]", "[Msun/yr]", "[Lsun]"
    
    WHILE NOT EOF(fpdatafile) DO BEGIN
        READU, fpdatafile, dataline
        WRITEU, fp, dataline
    ENDWHILE
    
    CLOSE, fp & FREE_LUN, fp
    
    CLOSE, fpdatafile & FREE_LUN, fpdatafile
    
    
    
    PRINT, 'Done!'
    
END






