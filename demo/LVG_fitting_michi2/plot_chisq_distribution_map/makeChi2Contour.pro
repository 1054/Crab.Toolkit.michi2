PRO makeChi2Contour, Title=Title, XTitle=XTitle, YTitle=YTitle, $
                     ColorbarTitle = ColorbarTitle, ColorbarColor = ColorbarColor, ColorbarLocation = ColorbarLocation, $
                     XArray, YArray, ZArray, Levels=Levels, $
                     XRange = XRange, YRange = YRange, $
                     XTICKINTERVAL = XTICKINTERVAL, YTICKINTERVAL = YTICKINTERVAL, $
                     SaveEPS = SaveEPS
    
    ; Check
    IF N_ELEMENTS(SaveEPS) EQ 0 THEN BEGIN
        MESSAGE, 'makeChi2Contour: Please set SaveEPS! <TODO>'
    ENDIF
    
    
    ; Contour Levels
    IF N_ELEMENTS(Levels) EQ 0 THEN BEGIN
        histos = CrabArrayFenBin(ZArray, Count=10, BinCentres=Levels, BinEdges=BinEdges) ; BinEdges=Levels
        Levels = levels[0:N_ELEMENTS(WHERE(Levels GT 0))-1]
        Levels = [0.0,Levels]
    ENDIF
    
    ; Axis
    IF N_ELEMENTS(XRange) EQ 2 THEN XStyle=1
    IF N_ELEMENTS(YRange) EQ 2 THEN YStyle=1
    
    
    ; Color Table
    ; LoadCT, 25
    TVLCT, CT_R, CT_G, CT_B, /GET
    
    
    ; Colors index
    ; Colors = BYTSCL([Levels[1:-1],MAX(ZArray)]) ; <TODO> need to be sure 
    Colors = BYTSCL(FINDGEN(N_ELEMENTS(Levels)))
    
    
    ; Contour
    PRINT, "contour levels: ", CrabStringPrintArray(Levels)
    SET_PLOT, 'PS'
    DEVICE, FILENAME=SaveEPS, /Color, /ENCAPSULATED, DECOMPOSED=0, XSIZE=24, YSIZE=20 ; cm
    
    Contour, XTITLE=XTitle, YTITLE=YTitle, TITLE=Title, $
             ZArray, XArray, YArray, /NODATA, LEVELS=Levels, /IRREGULAR, $
             COLOR=cgColor("black"), $
             XRange=XRange, YRange=YRange, $
             XStyle=XStyle, YStyle=YStyle, $
             XTICKINTERVAL=XTICKINTERVAL, YTICKINTERVAL=YTICKINTERVAL, $
             XThick=1.5, YThick=1.5, $
             Thick=1.5, CharSize=2.2, CharThick=4.2
    
    
    Contour, $
             ZArray, XArray, YArray, /OVERPLOT, /FILL, LEVELS=Levels, /IRREGULAR
    
    
    FOR ci=0,N_ELEMENTS(Levels)-1 DO BEGIN
        
;        ; plot some dashed shadow for high Chi2 regions <TODO>
;        
;        IF Levels[ci] LT 2.0 THEN CONTINUE
;        IF N_ELEMENTS(WHERE(ZArray GE Levels[ci],/NULL)) EQ 0 THEN CONTINUE
;        
;        Contour, $
;             ZArray, XArray, YArray, /OVERPLOT, /FILL, LEVELS=[Levels[ci],MAX([Levels,ZArray])], /IRREGULAR, $
;             COLOR=cgColor("black"), BACKGROUND=Colors[ci], $
;             C_ORIENTATION=[90], C_LINESTYLE=1, C_COLORS=[cgColor("black"),cgColor("black")]
        
        ; BREAK
        
    ENDFOR
    
    Contour, $
             ZArray, XArray, YArray, /OVERPLOT, LEVELS=Levels, /IRREGULAR, $
             Color = cgColor("black"), $
             C_Colors = [REPLICATE(cgColor("black"),N_ELEMENTS(Levels))], $
             C_Labels =  REPLICATE(1,N_ELEMENTS(Levels)), $
             C_CharSize=1.0, C_CharThick=2.0
    
    
    ; finally plot the axis ticket -- they are usually covered up by contour
    PLOT, [0], [0], /NoErase, /NoData, $
             COLOR=cgColor("black"), $
             XRange=XRange, YRange=YRange, $
             XStyle=XStyle, YStyle=YStyle, $
             XTICKFORMAT='(A1)', YTICKFORMAT='(A1)', $
             XTICKINTERVAL=XTICKINTERVAL, YTICKINTERVAL=YTICKINTERVAL, $
             XThick=1.5, YThick=1.5, $
             Thick=1.5, CharSize=2.2, CharThick=4.2
    
    
    
    
    IF !D.NAME EQ 'PS' THEN BEGIN
        dl_levels = Levels
        dl_barBytes = bytscl(FINDGEN(N_ELEMENTS(dl_levels)),TOP=255)
        dl_barPixel = FINDGEN(256)/255.d * (N_ELEMENTS(dl_levels))
        dl_barValue = dl_barPixel*0.0+255.d
        dl_barTicks = []
        dl_barPosit = [0.83, 0.65, 0.86, 0.90] ; x0 y0 x1 y1
        IF N_ELEMENTS(ColorbarLocation) EQ 1 AND SIZE(ColorbarLocation,/TNAME) EQ 'STRING' THEN BEGIN
            IF STRLOWCASE(ColorbarLocation) EQ 'left' THEN BEGIN
                dl_barPosit = [0.23, 0.65, 0.26, 0.90] ; x0 y0 x1 y1
            ENDIF
        ENDIF
        FOR dl_i=0,N_ELEMENTS(dl_levels)-1 DO BEGIN
            dl_barFlag = WHERE(dl_barPixel GE dl_i AND dl_barPixel LT dl_i+1)
            IF N_ELEMENTS(dl_barFlag) GT 0 THEN dl_barValue[dl_barFlag] = dl_barBytes[dl_i]
            dl_barChar = ''
            IF ABS(dl_levels[dl_i]) LT 10. THEN dl_barChar = STRING(FORMAT='(F0.1)',dl_levels[dl_i]) ELSE dl_barChar = STRING(FORMAT='(G0.3)',dl_levels[dl_i])
            dl_barTicks = [dl_barTicks, dl_barChar] 
        ENDFOR
        dl_barImage = REPLICATE(1,20) # dl_barValue
        dl_barTicks = [dl_barTicks, " "] 
        dl_barPosX0 = dl_barPosit[0]*!D.X_SIZE
        dl_barPosY0 = dl_barPosit[1]*!D.Y_SIZE
        dl_barPosX1 = dl_barPosit[2]*!D.X_SIZE
        dl_barPosY1 = dl_barPosit[3]*!D.Y_SIZE
        TV, dl_barImage, dl_barPosX0, dl_barPosY0, XSIZE=dl_barPosX1-dl_barPosX0, YSIZE=dl_barPosY1-dl_barPosY0
        
        dl_barColor = cgColor('black')
        dl_barColor = 255
        IF N_ELEMENTS(ColorbarColor) EQ 1 AND SIZE(ColorbarColor,/TNAME) EQ 'INT' THEN dl_barColor = ColorbarColor
        IF N_ELEMENTS(ColorbarColor) EQ 1 AND SIZE(ColorbarColor,/TNAME) EQ 'LONG' THEN dl_barColor = ColorbarColor
        IF N_ELEMENTS(ColorbarColor) EQ 1 AND SIZE(ColorbarColor,/TNAME) EQ 'STRING' THEN dl_barColor = cgColor(ColorbarColor)
        
        Plot, [0], [0], /NoData, /NOERASE, $
              Color=dl_barColor, $
              CharSize=1.2, CharThick=3.0, $
              XRange=[1,20], YRange=[1,N_ELEMENTS(Levels)], $
              XTicks=1, YTicks=N_ELEMENTS(Levels), $
              XTICKFORMAT='(A1)', YTICKFORMAT='(A1)', $
              XMinor=1, YMinor=1, $
              Position=dl_barPosit
        
        Axis, XTitle="", YTitle="", /YAXIS, $
              Color=dl_barColor, $
              CharSize=1.2, CharThick=3.0, $
              YRange=[1,N_ELEMENTS(Levels)], $
              YTicks=N_ELEMENTS(Levels), $
              YTickName=dl_barTicks, $
              YTickLen=0.25, $
              YMinor=1, $
              YStyle=1
        
        IF N_ELEMENTS(ColorbarTitle) EQ 1 THEN BEGIN
              XYOUTS, dl_barPosX1+250, dl_barPosY1-100, ColorbarTitle, /DEVICE, $
                      Color=dl_barColor, $
                      CharSize=1.2, CharThick=3.0, $
                      Alignment=0.0
        ENDIF
        
    ENDIF ELSE BEGIN
    
        cgColorBar, POSITION=[0.85, 0.65, 0.89, 0.90], /VERTICAL, /RIGHT, $
                    TITLE=ColorBarTitle, CharThick=3.0, $
                    RANGE=CrabMinMax(Levels), dl_levels=Levels
    
    ENDELSE
    
    DEVICE, /CLOSE
    SET_PLOT, 'X'
   ;ConvertPS2PDF, SaveEPS
    
END
