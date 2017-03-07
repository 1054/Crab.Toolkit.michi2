; ------------------------------------------------------------------------------------------------------
; CrabStringPrintFloat ---
;                 CrabStringPrintFloat( 1.34343, PREC=3)
;                      -->             "1.343"
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringPrintFloat, Number, FORMAT=FORMAT, PRECISION=PRECISION
    IF SIZE(Number,/TNAME) NE 'FLOAT' AND SIZE(Number,/TNAME) NE 'DOUBLE' THEN RETURN, ''
    IF N_ELEMENTS(PRECISION) EQ 1 THEN BEGIN
        IF PRECISION GT 0 THEN FORMAT = '(F0.'+STRTRIM(PRECISION,2)+')'
    ENDIF ELSE BEGIN
        FORMAT = '(F0.3)'
    ENDELSE
    PrintText = ""
    PrintText = STRTRIM(STRING(Number[0],FORMAT=FORMAT),2)
    RETURN, PrintText
END