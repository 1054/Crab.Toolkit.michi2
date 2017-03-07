; ------------------------------------------------------------------------------------------------------
; CrabStringClean ---
;                 if input '123^456^789^^ABC' then output '123456789ABC'
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringClean, Str, TextsToRemove=TextsToRemove
    DumpStr = Str
    CleanStr = Str
    IF N_ELEMENTS(TextsToRemove) EQ 0 THEN BEGIN
        BadStrs = ['^', '/', '\\',',','{','}','_']
    ENDIF ELSE BEGIN
        BadStrs = TextsToRemove
    ENDELSE
    FOR SI=0,N_ELEMENTS(CleanStr)-1 DO BEGIN
        CleanSubStr = CleanStr[SI]
        DumpSubStr = DumpStr[SI]
        FOREACH BadStr, BadStrs DO BEGIN
            IF STRPOS(DumpSubStr,BadStr) GE 0 THEN BEGIN
                SplitStr = STRSPLIT(CleanSubStr,BadStr,/EXTRACT) ; <TODO> can not FOLD_CASE
                CleanSubStr = STRJOIN(SplitStr,'')
                CleanSubStr = STRTRIM(CleanSubStr,2)
            ENDIF
        ENDFOREACH
        CleanStr[SI] = CleanSubStr
    ENDFOR
    ; a = '123^456^789^^ABC'
    ; print, CrabStringClean(a)
    ; CleanStr = STRTRIM(CleanStr,2
    RETURN, CleanStr
END
FUNCTION CrabStringCleaning, Str, TextsToRemove
    RETURN, CrabStringClean(Str, TextsToRemove=TextsToRemove)
END