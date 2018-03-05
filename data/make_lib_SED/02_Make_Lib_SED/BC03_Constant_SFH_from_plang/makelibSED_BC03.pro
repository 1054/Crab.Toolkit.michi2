
PRO makelibSED_BC03
    ; 
    restore, 'Zsolar_cst_SF1_spec_MIRI.save', /verbose
    SPEC_SFR = SFR
    SPEC_MASS = M_TOT
    SPEC_AGE = 10^LOGAGE
    ;
    outfile = 'lib.BC03.CstSFH.Z0.0190.SED'
    ; 
    OPENW, fp, outfile, /GET_LUN
    PRINTF, fp, FORMAT='(A,A)', "# ", CrabStringCurrentTime()
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR  = 2 # Wave & Flux"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR1 = band wavelength (rest-frame) (lambda) [um]"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TVAR2 = flux density in unit of solar luminosity per Hz"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR1 = 1 # the colomn number of VAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CVAR2 = 2 # the colomn number of VAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR1 = "+STRING(FORMAT='(I0)',(SIZE(SPEC_FLUX,/DIM))[1])+" # NAXIS1 -- the number of wavelength in one spectrum"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NVAR2 = "+STRING(FORMAT='(I0)',(SIZE(SPEC_FLUX,/DIM))[0])+"    # NAXIS2 -- the number of templates -- should equal to NPAR1*NPAR2*NPAR3*NPAR4"
    PRINTF, fp, FORMAT='(A,A)', "# ", "NPAR  = 4           # Metallicity, Age, Mass, SFR"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR1 = Metallicity # (12+log(O/H), note that solar metallicity is 0.0190)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR2 = Age         # in unit of Gyr"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR3 = Mass        # in unit of solar mass (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "TPAR4 = SFR         # in unit of solar mass per year (associated with Age, not independent)"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR1 = 3           # the colomn number of PAR1"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR2 = 4           # the colomn number of PAR2"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR3 = 5           # the colomn number of PAR3"
    PRINTF, fp, FORMAT='(A,A)', "# ", "CPAR4 = 6           # the colomn number of PAR4"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR1 = ", 1, " # Number of different PAR1"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR2 = ", (SIZE(SPEC_FLUX,/DIM))[0], " # Number of different PAR2"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR3 = ", 1, " # Number of different PAR3"
    PRINTF, fp, FORMAT='(A,A,I0,A)', "# ", "NPAR4 = ", 1, " # Number of different PAR4"
    PRINTF, fp, FORMAT='(A)', "# "
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A15,A15)', "# ", "lambda", "flux", "metallicity", "age", "mass", "SFR"
    PRINTF, fp, FORMAT='(A,A13,A15,A15,A12,A15,A15)', "# ", "[um]", "[Lsun/Hz]", "[---]", "[Gyr]", "[Msun]", "[Msun/yr]"
    
    tmpVar = 0
    
    ; 
    FOR specindex = 0, (SIZE(SPEC_FLUX,/DIM))[0]-1 DO BEGIN
        
        ;TmpPos1 = STRPOS(specfile,'_Chabrier_Z')+STRLEN('_Chabrier_Z')
        ;TmpPos2 = STRPOS(specfile,'.out.spec')
        ;ZMeta = STRMID(specfile,TmpPos1,TmpPos2-TmpPos1)
        ZMeta = '0.0190'
        
        PRINT, 'Processing spectrum ' + STRING(FORMAT='(I0)',specindex) + ' Z=' + ZMeta + ' age='+STRING(FORMAT='(E0.3)',SPEC_AGE[specindex])
        
        wave_AA = SPEC_LAMBDA[specindex,*]
        flux_LsunAA = SPEC_FLUX[specindex,*]
        
        lambdaLlambda = wave_AA * flux_LsunAA
        vLv = lambdaLlambda
        wave = wave_AA / 1e4 ; [um]
        freq = 2.99792458d5 / wave * 1e9 ; Hz
        flux = vLv / freq ; [L_solar/Hz]
        
        age = SPEC_AGE[specindex] / 1e9 ; Gyr
        
        sfr = SPEC_SFR[specindex] ; Msun/yr
        mass = SPEC_MASS[specindex] ; Msun, cumulated
        
        FOR i=0,N_ELEMENTS(flux)-1 DO BEGIN
            PRINTF, fp, FORMAT='(E15.6,E15.6,F15.5,F12.5,E15.6,E15.6)', wave[i], flux[i], DOUBLE(ZMeta), age, mass, sfr
        ENDFOR
        
        
        
        
    ENDFOR
    
    CLOSE, fp & FREE_LUN, fp
    
    PRINT, 'Done!'
    
END
