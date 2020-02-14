PRO makeChi2ContourTkVnH2
    
    
    ;CD, '/Users/dliu/Working/SpireLines/Tool/Level_3_SciData/06_LVG_Synthesis/61_Fit_CO_Ladder/go_fit_coir_normalization_sled'
    
    SPAWN, 'rm xyz_Tk_nH2_*_excitation*.csv'
    
    ; Extract the jpeak1 jpeak2 param-space
    IF NOT FILE_TEST("xyz_Tk_nH2_low_excitation.csv") THEN BEGIN
 
;       readcol, "fit_double.out", format='(A,A,A,A,A,A,A,A,A,A,A,A)', skipline=2, id, Chi, i1, a1, i2, a2, Tk1, nH1, MH1, Tk2, nH2, MH2
        readcol, "fit_double.out", format='(A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A,A)', skipline=2, id, Chi, i1, a1, i2, a2, Tk1, nH1, NULL_1, NULL_2, NULL_3, NULL_4, NULL_5, NULL_6, NULL_7, NULL_8, NULL_9, NULL_10, NULL_11, NULL_12, NULL_13, Tk2, nH2, MH2
        id  = STRING(FORMAT='(I0)',INDGEN(N_ELEMENTS(id),/LONG))
        nH1 = STRING(FORMAT='(F0.2)',ALOG10(DOUBLE(nH1)))
        nH2 = STRING(FORMAT='(F0.2)',ALOG10(DOUBLE(nH2)))
        
;       Where1IsLow2IsHigh = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 LE jpeak2, /NULL)
        Where1IsLow2IsHigh = WHERE(FLOAT(a1) GT 0 AND FLOAT(a2) GT 0 AND FLOAT(nH1) LE FLOAT(nH2), /NULL)
        ColdArrI =  id[Where1IsLow2IsHigh]
        ColdArrX = nH1[Where1IsLow2IsHigh]
        ColdArrY = Tk1[Where1IsLow2IsHigh]
        ColdArrZ = Chi[Where1IsLow2IsHigh]
        CrabTablePrintC, "xyz_Tk_nH2_low_excitation_full.csv", ColdArrI, ColdArrX, ColdArrY, ColdArrZ;, FORMAT='I15,F15.5,F15.5,G15.5'
        makeChi2MinArray, ColdArrX, ColdArrY, ColdArrZ, IndexArray=ColdArrI
        CrabTablePrintC, "xyz_Tk_nH2_low_excitation.csv", ColdArrI, ColdArrX, ColdArrY, ColdArrZ;, FORMAT='I15,F15.5,F15.5,G15.5'
        
        WarmArrI =  id[Where1IsLow2IsHigh]
        WarmArrX = nH2[Where1IsLow2IsHigh]
        WarmArrY = Tk2[Where1IsLow2IsHigh]
        WarmArrZ = Chi[Where1IsLow2IsHigh]
        CrabTablePrintC, "xyz_Tk_nH2_high_excitation_full.csv", WarmArrI, WarmArrX, WarmArrY, WarmArrZ;, FORMAT='I15,F15.5,F15.5,G15.5'
        makeChi2MinArray, WarmArrX, WarmArrY, WarmArrZ, IndexArray=WarmArrI
        CrabTablePrintC, "xyz_Tk_nH2_high_excitation.csv", WarmArrI, WarmArrX, WarmArrY, WarmArrZ;, FORMAT='I15,F15.5,F15.5,G15.5'
        
    ENDIF
    
    
    
;    cgLoadCT, 25, /Brewer, /Reverse
;    cgLoadCT, 70, /Reverse
;    cgLoadCT, 72, /Reverse
;    cgLoadCT, 73, /Reverse  ; my fav
    cgLoadCT, 25
    
;   levels = [0.0, 0.1, 0.2, 0.5, 2.0, 10.0, 50.0, 100.0]
;   levels = [0.0, 1.0, 1.5, 2.0, 10.0, 50.0, 100.0, 200.0]
    levels = [1.5, 2.0, 10.0, 50.0, 100.0, 200.0]
    
    FigYRange = [10.,295.]
    FigXRange = [2.0,6.1]
    
    FigTitle = cgSymbol("chi")+"!E2!N"
    
    FigXTitle = "n!DH!L2!D (low excitation)!N"
    FigYTitle = "T!Dkin (low excitation)!N"
    
    readcol, "xyz_Tk_nH2_low_excitation.csv", ArrayI, ArrayX, ArrayY, ArrayZ
    
    makeChi2Contour, ArrayX, ArrayY, ArrayZ, Title=FigTitle, XTitle=FigXTitle, YTitle=FigYTitle, $
                     ColorbarTitle=FigTitle, ColorbarColor='black', $
                     XRange=FigXRange, YRange=FigYRange, $
                     Levels=levels, SaveEPS='fit_chi2_distribution_low_excitation.eps'
    
    FigXTitle = "n!DH!L2!D (high excitation)!N"
    FigYTitle = "T!Dkin (high excitation)!N"
    
    readcol, "xyz_Tk_nH2_high_excitation.csv", ArrayI, ArrayX, ArrayY, ArrayZ
    
    makeChi2Contour, ArrayX, ArrayY, ArrayZ, Title=FigTitle, XTitle=FigXTitle, YTitle=FigYTitle, $
                     ColorbarTitle=FigTitle, ColorbarColor='black', ColorbarLocation = 'left', $
                     XRange=FigXRange, YRange=FigYRange, $
                     Levels=levels, SaveEPS='fit_chi2_distribution_high_excitation.eps'
    
    
END
