; ------------------------------------------------------------------------------------------------------
; CrabStringPadSpace ---
;                 CrabStringPadSpace( "YES!", PREC=3)
;                      -->             "1.343"
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringPadSpace, InputTexts, PadTrailing=PadTrailing, PadLeading=PadLeading, PadLength=PadLength, TotalLength=TotalLength
    FOREACH InputText, InputTexts DO BEGIN
        InputLength = STRLEN(InputText)
        OutputText = ""
        IF KEYWORD_SET(TotalLength) THEN BEGIN
            PadLength = TotalLength - InputLength
        ENDIF
        IF KEYWORD_SET(PadLeading) OR NOT KEYWORD_SET(PadTrailing) THEN BEGIN
            OutputText = string(replicate(32b,PadLength)) + InputText
        ENDIF
        IF KEYWORD_SET(PadTrailing) THEN BEGIN
            OutputText = InputText + string(replicate(32b,PadLength))
        ENDIF
    ENDFOREACH
    RETURN, OutputText
END