PRO makeChi2ContourTkVTk
    
;    CD, "/Users/dliu/Working/2014-CEA/Tool/Level_3_SciData/LVG/"
    CD, "/Users/dliu/Working/2014-CEA/Data/Level_3_SciData/BzK_CO_SED/BzK-All"
    
    
    ; Extract the jpeak1 jpeak2 param-space
    IF NOT FILE_TEST("xyz_Tk_Tk.csv") THEN BEGIN
        readcol, "dl_out_jpeak.csv", id, Tk1, nH1, Tk2, nH2, a1, a2, i1, i2, Chi, ChiNum, jpeak1, jpeak2
        readcol, "dl_jf1.csv", id1s, f1co1
        readcol, "dl_jf2.csv", id2s, f2co1
        
        
        Where1IsLow = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 LE jpeak2, /NULL)
        Where2IsLow = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 GT jpeak2, /NULL)
        Where2IsHigh = Where1IsLow
        Where1IsHigh = Where2IsLow
        ArrayI = [ id[Where1IsLow], id[Where2IsLow]]
        ArrayX = [Tk1[Where1IsLow],Tk2[Where2IsLow]]
        ArrayY = [Tk2[Where1IsLow],Tk1[Where2IsLow]]
        ArrayZ = [Chi[Where1IsLow],Chi[Where2IsLow]]
        makeChi2MinArray, ArrayX, ArrayY, ArrayZ, IndexArray=ArrayI
        
;        Filters = WHERE(jpeak1 GT 0 AND jpeak2 GT 0 AND jpeak1 LE jpeak2, /NULL)
;        ArrayI = [] & ArrayX = [] & ArrayY = [] & ArrayZ = []
;        FOR i=0,N_ELEMENTS(Filters)-1 DO BEGIN
;            Filter = Filters[i]
;            IF Filter GE 0 THEN BEGIN
;                TempArrX = Tk1[Filters] & TempVarX = Tk1[Filter]
;                TempArrY = Tk2[Filters] & TempVarY = Tk2[Filter]
;                TempFlag = WHERE(TempArrX EQ TempVarX AND TempArrY EQ TempVarY)
;                TempArrZ = (ChiSq[Filters])[TempFlag]
;                ArrayI = [ArrayI,id[Filter]]
;                ArrayX = [ArrayX,Tk1[Filter]]
;                ArrayY = [ArrayY,Tk2[Filter]]
;                ArrayZ = [ArrayZ,MIN(TempArrZ)]
;                Filters[TempFlag] = -9
;            ENDIF
;        ENDFOR
        CrabTablePrintC, "xyz_Tk_Tk.csv", ArrayI, ArrayX, ArrayY, ArrayZ
        
        RETURN
    ENDIF
    
    
    
    
    
    
    
    LoadCT, 25
    levels = [0.0, 0.1, 0.2, 0.5, 2.0, 10.0, 50.0, 100.0]
    FigFile = 'fit_6_chi2_Tk_Tk.eps'
    FigTitle = cgSymbol("chi")+"!E2!N"
    FigXTitle = "T!Dkin (low excitation)!N"
    FigYTitle = "T!Dkin (high excitation)!N"
    FigXRange = [10,100]
    FigYRange = [10,100]
    FigColorBarTitle = cgSymbol("chi")+"!E2!N"
    
    
    
    
    
    
    
    ; ReadCSV
    readcol, "xyz_Tk_Tk.csv", ArrayI, ArrayX, ArrayY, ArrayZ
    
    ; Contour
    makeChi2Contour, ArrayX, ArrayY, ArrayZ, Title=FigTitle, XTitle=FigXTitle, YTitle=FigYTitle, $
                     ColorbarTitle=FigTitle, ColorbarColor='black', $ 
                     XRange=FigXRange, YRange=FigYRange, $
                     Levels=levels, SaveEPS=FigFile
    
    
    
    
    ; Plot
    ; CrabImageQuickPlot, AxisZ, TVZoom=500
    
END
