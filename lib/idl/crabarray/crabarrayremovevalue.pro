; 
; remove some values in an array
; 
FUNCTION CrabArrayRemoveValue, InputArray, ValuesToRemove
    
    NewArray = InputArray
    FOREACH ValueToRem, ValuesToRemove DO BEGIN
        IdMatched = WHERE(NewArray EQ ValueToRem, /NULL)
        IF N_ELEMENTS(IdMatched) GT 0 THEN BEGIN
            IdToKeep = WHERE(NewArray NE ValueToRem, /NULL)
            IF N_ELEMENTS(IdToKeep) GT 0 THEN BEGIN
                NewArray = NewArray[IdToKeep]
            ENDIF ELSE BEGIN
                ; all items are removed!
                RETURN, []
            ENDELSE
        ENDIF
    ENDFOREACH
    RETURN, NewArray
    
END