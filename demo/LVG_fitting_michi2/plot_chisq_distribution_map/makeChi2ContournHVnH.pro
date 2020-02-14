PRO makeChi2ContournHVnH
    
    
    CD, '/Users/dliu/Working/SpireLines/Tool/Level_3_SciData/06_LVG_Synthesis/61_Fit_CO_Ladder/go_fit_coir_normalization_sled'
    
    SPAWN, 'rm xyz_nH_nH.csv'
    
    ; Extract the jpeak1 jpeak2 param-space
    IF NOT FILE_TEST("xyz_nH_nH.csv") THEN BEGIN
        
        readcol, "fit-double.csv", format='(A,A,A,A,A,A,A,A,A,A)', skipline=2, id, Chi, i1, a1, i2, a2, Tk1, nH1, Tk2, nH2
;       readcol, "fit-double_jpeak.csv", idj, jpeak1, jpeak2
        id  = STRING(FORMAT='(I0)',INDGEN(N_ELEMENTS(id),/LONG))
        nH1 = STRING(FORMAT='(F0.2)',ALOG10(DOUBLE(nH1)))
        nH2 = STRING(FORMAT='(F0.2)',ALOG10(DOUBLE(nH2)))
        
;       Where1IsLow = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 LE jpeak2, /NULL)
;       Where2IsLow = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 GT jpeak2, /NULL)
        Where1IsLow = WHERE(FLOAT(a1) GT 0 AND FLOAT(a2) GT 0 AND FLOAT(nH1) LE FLOAT(nH2), /NULL)
        Where2IsLow = WHERE(FLOAT(a1) GT 0 AND FLOAT(a2) GT 0 AND FLOAT(nH1) GT FLOAT(nH2), /NULL)
        ArrayI = [ id[Where1IsLow], id[Where2IsLow]]
        ArrayX = [nH1[Where1IsLow],nH2[Where2IsLow]]
        ArrayY = [nH2[Where1IsLow],nH1[Where2IsLow]]
        ArrayZ = [Chi[Where1IsLow],Chi[Where2IsLow]]
        CrabTablePrintC, "xyz_nH_nH_full.csv", ArrayI, ArrayX, ArrayY, ArrayZ ;FORMAT='I15,F15.5,F15.5,G15.5'
        makeChi2MinArray, ArrayX, ArrayY, ArrayZ, IndexArray=ArrayI
        CrabTablePrintC, "xyz_nH_nH.csv", ArrayI, ArrayX, ArrayY, ArrayZ ;, FORMAT='I15,F15.5,F15.5,G15.5'
        
    ENDIF
    
    
    
    LoadCT, 25
;   levels = [0.0, 0.1, 0.2, 0.5, 2.0, 10.0, 50.0, 100.0]
    levels = [0.0, 1.0, 1.5, 2.0, 10.0, 50.0, 100.0, 200.0]
    
    FigXRange = [2.0,6.1]
    FigYRange = [2.0,6.1]
    
    FigTitle = cgSymbol("chi")+"!E2!N"
    
    ; Axis
    FigXTitle = "n!DH!L2!D (low excitation)!N"
    FigYTitle = "n!DH!L2!D (high excitation)!N"
    
    ; Readcol
    readcol, "xyz_nH_nH.csv", ArrayI, ArrayX, ArrayY, ArrayZ
    
    ; Contour
    makeChi2Contour, ArrayX, ArrayY, ArrayZ, Title=FigTitle, XTitle=FigXTitle, YTitle=FigYTitle, $
                     ColorbarTitle=FigTitle, ColorbarColor='black', $
                     XRange=FigXRange, YRange=FigYRange, $
                     Levels=levels, SaveEPS='fit_5_chi2_nH_nH.eps'
    
    
    
END
