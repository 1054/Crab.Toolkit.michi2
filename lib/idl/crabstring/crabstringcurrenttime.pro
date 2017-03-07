; ------------------------------------------------------------------------------------------------------
; CrabStringCurrentTime --- 
;                 PRINT, CrabStringCurrentTime()
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringCurrentTime, UTC=UTC, Location=Location, Precision=Precision, Compress=Compress
    IF N_ELEMENTS(Precision) EQ 0 THEN PrecOffset = 0 ELSE PrecOffset = Precision+1
    CurrentDateTime = TIMESTAMP(/UTC)
    CurrentDate = STRMID(CurrentDateTime,0,STRPOS(CurrentDateTime,'T'))
    CurrentTime = STRMID(CurrentDateTime,STRPOS(CurrentDateTime,'T')+1)
    CurrentTime = STRMID(CurrentTime,0,STRPOS(CurrentTime,'.')+PrecOffset) ; only need 3 digits of seconds
    IF N_ELEMENTS(Location) EQ 1 THEN BEGIN
        ; TODO Location
    ENDIF
    IF N_ELEMENTS(UTC) EQ 0 THEN BEGIN
        UTC = 0
    ENDIF ELSE BEGIN
        CurrentHour = STRMID(CurrentTime,0,STRPOS(CurrentTime,':'))
        CurrentTime = STRMID(CurrentTime,STRPOS(CurrentTime,':'))
        ; PRINT, CurrentHour
        ; PRINT, CurrentTime
        CurrentTime = STRING(FORMAT='(I0)',FIX(CurrentHour)+FIX(UTC)) + CurrentTime
    ENDELSE
    OutputDateTime = CurrentDate + ' ' + CurrentTime + ' ' + STRING(FORMAT='("UTC",I+0)',UTC)
    ; Compress the string
    IF KEYWORD_SET(Compress) THEN OutputDateTime = CrabStringClean(CurrentDate,TextsToRemove=['-']) + '.' + CrabStringClean(CurrentTime,TextsToRemove=[':']) + '.' + STRING(FORMAT='("UTC",I+0)',UTC)
    RETURN, OutputDateTime
    ; TODO
    ; Another method
    ; PRINT, STRING(SYSTIME(/JULIAN),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
END