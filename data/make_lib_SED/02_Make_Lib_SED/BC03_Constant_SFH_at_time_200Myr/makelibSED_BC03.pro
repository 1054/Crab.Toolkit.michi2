
PRO makelibSED_BC03
    ; 
    speclist = ['datatable_wave_AA_flux_Lsun_per_AA.txt']
    ;
    outfile = 'lib.BC03.Padova1994.BaSeL.Z0.0190.SED'
    ; 
    OPENW, fp, outfile, /GET_LUN
    PRINTF, fp, FORMAT='(A,A)', "# ", CrabStringCurrentTime()
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR  = 2 # Wave & Flux"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR1 = band wavelength (rest-frame) (lambda) [um]"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR2 = flux density in unit of solar luminosity per Hz"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR1 = 1 # the colomn number of VAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR2 = 2 # the colomn number of VAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR1 = 2023 # NAXIS1 -- the number of wavelength in one spectrum"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR2 = 1    # NAXIS2 -- the number of spectra -- should equal to NPAR1*NPAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NPAR  = 3    # Metallicity, Time, Mass"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR1 = Absolute Metallicity (note that solar metallicity is 0.0190)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR2 = Time in unit of Gyr"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR3 = Mass in unit of solar mass (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR1 = 3 # the colomn number of PAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR2 = 4 # the colomn number of PAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR3 = 5 # the colomn number of PAR3"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR1 = ", 1, " # Number of different PAR1"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR2 = ", 1, " # Number of different PAR2"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR3 = ", 1, " # Number of different PAR3"
    PRINTF, fp, FORMAT='(A)', "# "
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A15)', "# ", "lambda", "flux", "metallicity", "age", "mass"
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A15)', "# ", "[um]", "[Lsun/Hz]", " ", "[Gyr]", "[Msun]"
    
    tmpVar = 0
    
    ; 
    FOREACH specfile, speclist DO BEGIN
        
        ;TmpPos1 = STRPOS(specfile,'_Chabrier_Z')+STRLEN('_Chabrier_Z')
        ;TmpPos2 = STRPOS(specfile,'.out.spec')
        ;ZMeta = STRMID(specfile,TmpPos1,TmpPos2-TmpPos1)
        ZMeta = '0.0190'
        
        PRINT, specfile, " ", ZMeta
        
        readcol, specfile, wave_AA, flux_LsunAA
        
        lambdaLlambda = wave_AA * flux_LsunAA
        vLv = lambdaLlambda
        wave = wave_AA / 1e4 ; [um]
        freq = 2.99792458d5 / wave * 1e9 ; Hz
        flux = vLv / freq ; [L_solar/Hz]
        
        age = 2e6 / 1e9 ; Gyr
        
        ; the BC03 CSP is constant SFH of 1 Msun/yr, age is 2Myr, mass loss is 
        sfr = 1.0 ; Msun/yr, see '/Users/dzliu/Softwares/BC03/a_dzliu_note_20171004_make_constant_SFH.txt'
        mass = sfr * age*1e9 ; <TODO> mass loss?
        
        FOR i=0,N_ELEMENTS(flux)-1 DO BEGIN
            PRINTF, fp, FORMAT='(E15.6,E15.6,F15.5,F12.5,E15.6)', wave[i], flux[i], DOUBLE(ZMeta), age, mass
        ENDFOR
        
        
        
        
    ENDFOREACH
    
    CLOSE, fp & FREE_LUN, fp
    
    PRINT, 'Done!'
    
END
