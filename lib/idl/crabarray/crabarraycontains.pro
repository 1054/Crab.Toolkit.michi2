; 
; check whether an array contains some values
; 
FUNCTION CrabArrayContains, InputArray, ValuesToMatch, YES_OR_NO=YES_OR_NO
    
    IF N_Params() LT 2 THEN BEGIN
        PRINT, 'Usage: CrabArrayContains, InputArray, ValuesToMatch'
    ENDIF
    
    IF N_ELEMENTS(InputArray) EQ 0 THEN BEGIN
        IF KEYWORD_SET(YES_OR_NO) THEN RETURN, 0
        RETURN, !NULL
    ENDIF
    
    NewArray = InputArray
    FOREACH ValueToMatch, ValuesToMatch DO BEGIN
        IdMatched = WHERE(NewArray EQ ValueToMatch, /NULL)
        IsMatched = (N_ELEMENTS(IdMatched) GT 0)
        IF N_ELEMENTS(IsMatchedArray) EQ 0 THEN BEGIN
            IsMatchedArray = IsMatched
        ENDIF ELSE BEGIN
            IsMatchedArray = [ IsMatchedArray, IsMatched ]
        ENDELSE
    ENDFOREACH
    
    IF KEYWORD_SET(YES_OR_NO) THEN BEGIN
        IsMatched = (TOTAL(IsMatchedArray) GT 0)
        RETURN, IsMatched
    ENDIF
    
    RETURN, IsMatchedArray
    
END