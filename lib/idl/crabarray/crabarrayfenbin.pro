; 
; Divide bins according to given levels or count
; 
FUNCTION CrabArrayFenBin, InputArray, Levels=InputLevels, Count=InputCount, Verbose=Verbose, Debug=Debug, Continue=Continue, Iteration=Iteration, $
                                      MIN=binMIN, MAX=binMAX, NotIncludingMIN=NotIncludingMIN, NotIncludingMAX=NotIncludingMAX, $
                                      BinEdges=binEdges, BinCentres=binCents, BinWidthes=binWidth, BinVolumes=binVoles
    
    IF KEYWORD_SET(Debug) THEN BEGIN
        InputArray = RANDOMU(3.17,3000,/DOUBLE)
        InputLevels = [0.1,0.2,0.3,0.3,0.2,0.1,0.05,0.01]
        Verbose = 1
    ENDIF
    
    IF N_ELEMENTS(InputArray) THEN RETURN, !NULL
    
    fCount = 8
    IF N_ELEMENTS(InputCount) GT 0 THEN fCount = InputCount
    
    fVolume = REPLICATE(0.1,fCount)
    IF N_ELEMENTS(InputLevels) GT 0 THEN fVolume = InputLevels ; InputLevels has higher priority.
    
    fArray = InputArray
    fFilter = WHERE(FINITE(InputArray),/NULL)
    IF N_ELEMENTS(fFilter) GT 0 THEN fArray = fArray[fFilter]
    IF N_ELEMENTS(binMIN)  EQ 1 AND KEYWORD_SET(NotIncludingMIN) THEN fFilter = WHERE(fArray GT binMIN, /NULL)
    IF N_ELEMENTS(binMIN)  EQ 1 AND NOT KEYWORD_SET(NotIncludingMIN) THEN fFilter = WHERE(fArray GE binMIN, /NULL)
    IF N_ELEMENTS(fFilter) GT 0 THEN fArray = fArray[fFilter]
    IF N_ELEMENTS(binMAX)  EQ 1 AND KEYWORD_SET(NotIncludingMAX) THEN fFilter = WHERE(fArray LT binMAX, /NULL)
    IF N_ELEMENTS(binMAX)  EQ 1 AND NOT KEYWORD_SET(NotIncludingMAX) THEN fFilter = WHERE(fArray LE binMAX, /NULL)
    IF N_ELEMENTS(fFilter) GT 0 THEN fArray = fArray[fFilter]
    
    fCount = N_ELEMENTS(fVolume)
    fVolume = ROUND(fVolume/TOTAL(fVolume)*N_ELEMENTS(fArray)) ; number of data points in each bin
    binMIN  = MIN(fArray)
    binMAX  = MAX(fArray)
    binCents = REPLICATE(0.0d,fCount)
    binWidth = REPLICATE(0.0d,fCount)
    binVoles = REPLICATE(0.0d,fCount)
    binEdges = REPLICATE(0.0d,fCount+1)
    binSteps = MIN([ABS(MEAN(fArray)-binMIN),ABS(binMAX-MEAN(fArray)),binMAX-binMIN])/1000.0
    
    IF binSteps EQ 0.0 THEN MESSAGE, 'CrabArrayFenBin: InputArray is not good enough to fen bin. Would you like to set MIN and MAX?', CONTINUE=Continue
    
    IF N_ELEMENTS(Iteration) NE 1 THEN Iteration = 30
    
    FOR i=0,fCount-1 DO BEGIN
        IF i EQ 0 THEN binEdges[0]=binMIN
        tempFlag = WHERE(fArray GE binEdges[i],tempCount,/NULL) ; extract the sub array that fArray >= binEdge
        IF tempCount GT 0 THEN BEGIN
            tArray = fArray[tempFlag] ; tArray is the sub array that fArray >= binEdge
            ; if MAX(tArray) EQ MIN(tArray) then it means that all items are equal -- this sub array is single valued!
            IF MAX(tArray) EQ MIN(tArray) OR i EQ fCount-1 THEN BEGIN
                binEdges[i+1] = MAX(tArray)
                tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
            ENDIF ELSE BEGIN
                ; <Updated><20140707><DzLiu>
                tStep = MAX(tArray) - MIN(tArray)
                tEdge = binEdges[i]
                tCont = 0
                WHILE ABS(TOTAL(tArray LT tEdge)-fVolume[i]) GE 1 DO BEGIN
                    IF TOTAL(tArray LT tEdge+tStep) GT fVolume[i] THEN BEGIN
                        tStep = tStep/10.0d
                        IF tCont GE Iteration THEN BREAK ; <Added><TODO><20140716><DzLIU>
                    ENDIF ELSE BEGIN
                        tEdge = tEdge+tStep
                        tStep = MAX(tArray)-tEdge
                        tCont = tCont+1
                    ENDELSE
                ENDWHILE
                IF tEdge GT binMAX THEN tEdge = binMAX
                binEdges[i+1] = tEdge
                tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
                ; PRINT, i, tempCount, binEdges[i], binEdges[i+1], binSteps
                ; binSteps = MIN([ABS(MEDIAN(tArray)-MIN(tArray)),ABS(MAX(tArray)-MEDIAN(tArray)),$
                ;                 ABS((MEAN(tArray)-MIN(tArray))),ABS((MAX(tArray)-MEAN(tArray))),$
                ;                 MAX(tArray)-MIN(tArray)])/1000.0
                ; binSteps = MIN([ABS(MEAN(tArray)-MIN(tArray)),ABS(MAX(tArray)-MEAN(tArray)),MAX(tArray)-MIN(tArray)])/1000.0
                ; test binSteps -- in case the step is too large
                ; binEdges[i+1] = binEdges[i] + binSteps
                ; tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
                ; WHILE tempCount GT fVolume[i] AND binEdges[i+1] LT binMAX DO BEGIN
                ;     binSteps = binSteps/10.0
                ;     binEdges[i+1] = binEdges[i] + binSteps
                ;     tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
                ; ENDWHILE
                ; calc binEdge -- within current bin count the number
                ; tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
                ; tempLPCT = 10 ; loop count
                ; WHILE tempCount LT fVolume[i] AND binEdges[i+1] LT binMAX DO BEGIN
                ;     binEdges[i+1] = binEdges[i+1] + binSteps ; * ALOG10(tempLPCT)
                ;     tempFlag = WHERE(fArray GE binEdges[i] AND fArray LT binEdges[i+1], tempCount, /NULL)
                ;     tempLPCT = tempLPCT + 1
                ; ENDWHILE
            ENDELSE
            binVoles[i] = tempCount
            binCents[i] = (binEdges[i]+binEdges[i+1])/2.0
            binWidth[i] = ABS(binEdges[i]-binEdges[i+1])
            IF KEYWORD_SET(Verbose) THEN BEGIN
                IF i EQ 0 THEN PRINT, FORMAT='("#",A9,A17,A17,A17,A15,A15,A17)', "BinId", "BinLeft", "BinCentre", "BinRight", "BinVolume", "InputVolume", "BinStep"
                PRINT, FORMAT='(I10,g17.7,g17.7,g17.7,I15,I15,g17.7)', i+1, binEdges[i], binCents[i], binEdges[i+1], binVoles[i], fVolume[i], binSteps
            ENDIF
            IF binEdges[i+1] GE binMAX THEN BEGIN
                BREAK
            ENDIF
        ENDIF
    ENDFOR
    
    ; PRINT, binEdges
    RETURN, binVoles
    
END