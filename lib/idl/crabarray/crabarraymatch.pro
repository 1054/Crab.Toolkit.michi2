; 
; Match two array
; 
FUNCTION CrabArrayMatch, ArrayToMatch, ArrayAsRefer, YES_OR_NO=YES_OR_NO
    
    IF N_ELEMENTS(ArrayToMatch) EQ 0 THEN BEGIN
        PRINT, 'CrabArrayMatch, ArrayToMatch, ArrayAsReference'
        RETURN, !NULL
    ENDIF
    
    ArrayMatched = MAKE_ARRAY(N_ELEMENTS(ArrayToMatch),/INTEGER,VALUE=0)
    
    IF N_ELEMENTS(ArrayAsRefer) EQ 0 THEN BEGIN
        
        ArrayMatched = ArrayMatched*0+1 ;<TODO><20140805><DzLIU> if reference array is empty, return all true!
        
    ENDIF ELSE BEGIN
        
        Arr1 = ArrayToMatch
        
        IF SIZE(ArrayToMatch,/TNAME) EQ 'STRING' AND SIZE(ArrayAsRefer,/TNAME) NE 'STRING' THEN Arr2 = STRING(ArrayAsRefer) ELSE Arr2 = ArrayAsRefer
            
        FOREACH Item2, Arr2 DO BEGIN
            ArrayMatched = (ArrayMatched OR (ArrayToMatch EQ Item2))
        ENDFOREACH
        
    ENDELSE
    
    IF KEYWORD_SET(YES_OR_NO) THEN ArrayMatched = (TOTAL(ArrayMatched) GT 0)
    
    RETURN, ArrayMatched
END