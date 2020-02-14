; 
; This program will remove all duplicated X,Y values, only keeping minimum Chi2 value. 
; 
PRO makeChi2MinArray, XArray, YArray, Chi2Array, IndexArray=IndexArray
    
    ; Check
    IF N_ELEMENTS(Chi2Array) EQ 0 THEN BEGIN
        MESSAGE, 'makeChi2Contour: Please set SaveEPS! <TODO>'
    ENDIF
    
    Filters = WHERE(Chi2Array GT 0.0, /NULL)
    
    ; Deduplicated Array
    DdupArrX = []
    DdupArrY = []
    DdupArrZ = []
    DdupArrI = []
    
    IF N_ELEMENTS(Filters) GT 0 THEN BEGIN
        WHILE N_ELEMENTS(WHERE(Filters GE 0,/NULL)) GT 0 DO BEGIN
            TempVarX = XArray[(WHERE(Filters GE 0))[0]]
            TempVarY = YArray[(WHERE(Filters GE 0))[0]]
            TempFlag = WHERE(Filters GE 0 AND XArray EQ TempVarX AND YArray EQ TempVarY, /NULL) ; corrected 20141217 dzliu ; <TODO> exactly equal?
            IF N_ELEMENTS(TempFlag) GT 0 THEN BEGIN
                TempArrZ = DOUBLE(Chi2Array[TempFlag])
                TempVarZ = (Chi2Array[TempFlag])[(WHERE(TempArrZ EQ MIN(TempArrZ)))[0]]
                TempVarI = (IndexArray[TempFlag])[(WHERE(TempArrZ EQ MIN(TempArrZ)))[0]]
                DdupArrX = [DdupArrX,TempVarX]
                DdupArrY = [DdupArrY,TempVarY]
                DdupArrZ = [DdupArrZ,TempVarZ]
                IF N_ELEMENTS(IndexArray) EQ N_ELEMENTS(Chi2Array) THEN DdupArrI = [DdupArrI,TempVarI] ; DdupArrI = [DdupArrI,IndexArray[Filter]]
                Filters[TempFlag] = -9
            ENDIF
        ENDWHILE
        
;        FOR fi=0,N_ELEMENTS(Filters)-1 DO BEGIN
;            Filter = Filters[fi]
;            IF Filter GE 0 THEN BEGIN
;                TempFilters = Filters[WHERE(Filters GE 0)] ; TempFilters contains at least Filter
;                TempArrX = XArray[TempFilters] ; corrected 20141217 dzliu
;                TempArrY = YArray[TempFilters] ; corrected 20141217 dzliu
;                TempVarX = XArray[Filter]
;                TempVarY = YArray[Filter]
;                
;                TempFlag = WHERE(TempArrX EQ TempVarX AND TempArrY EQ TempVarY, /NULL) ; <TODO> exactly equal?
;                
;                IF N_ELEMENTS(TempFlag) GT 0 THEN BEGIN
;                    TempArrZ = (Chi2Array[TempFilters])[TempFlag]
;                    TempArrI = (IndexArray[TempFilters])[TempFlag] ; corrected 20141217 dzliu
;                    
;                    TempFilter = WHERE(TempArrZ EQ MIN(TempArrZ))
;                    TempVarI = TempArrI[TempFilter]
;                    TempVarZ = TempArrZ[TempFilter]
;                    
;                    DdupArrX = [DdupArrX,TempVarX]
;                    DdupArrY = [DdupArrY,TempVarY]
;                    DdupArrZ = [DdupArrZ,TempVarZ]
;                    IF N_ELEMENTS(IndexArray) EQ N_ELEMENTS(Chi2Array) THEN DdupArrI = [DdupArrI,TempVarI] ; DdupArrI = [DdupArrI,IndexArray[Filter]] 
;                    
;                    Filters[TempFlag] = -9
;                ENDIF
;            ENDIF
;        ENDFOR
    ENDIF
    
    ; 
    XArray = DdupArrX
    YArray = DdupArrY
    Chi2Array = DdupArrZ
    IndexArray = DdupArrI ; corrected 20141217 dzliu
    
    
END
