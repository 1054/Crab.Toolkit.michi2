; ------------------------------------------------------------------------------------------------------
; CrabStringFindWholeWord --- 
;                 
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringFindWholeWord, Text, SearchStr, YES_OR_NO=YES_OR_NO
    ; exam input
    IF SIZE(Text,/TYPE) NE 7 OR SIZE(SearchStr,/TYPE) NE 7 THEN BEGIN
        IF KEYWORD_SET(YES_OR_NO) THEN BEGIN
            RETURN, 0
        ENDIF ELSE BEGIN
            Print, 'Usage: position = CrabStringFindWholeWord(text,searchstr)'
            Return, []
        ENDELSE
    ENDIF
    ; initial search
    FirstPos = STRPOS(Text, SearchStr)
    IF FirstPos EQ -1 THEN BEGIN ; unfortunately, we got no any single match
        Return, []
    ENDIF
    ; search all matched positions from text
    PosList = []
    FoundPos = -2
    WHILE FoundPos NE -1 AND FoundPos LT STRLEN(Text) DO BEGIN
        FoundPos = STRPOS(Text, SearchStr, FoundPos+1) ; start character position
        LeftPos  = FoundPos - 1
        RightPos = FoundPos + STRLEN(SearchStr)
        LeftCheck = 0   ; check whether the character at the left of the found position is blank space. 1 for yes.
        RightCheck = 0  ; check whether the character at the right of the found position is blank space. 1 for yes.
        IF LeftPos GE 0 THEN BEGIN
            IF STRMID(Text,LeftPos,1) EQ " " THEN BEGIN
                LeftCheck = 1
            ENDIF
        ENDIF ELSE LeftCheck = 1
        IF RightPos LT  STRLEN(Text) THEN BEGIN
            IF STRMID(Text,RightPos,1) EQ " " THEN BEGIN
                RightCheck = 1
            ENDIF
        ENDIF ELSE RightCheck = 1
        IF LeftCheck EQ 1 AND RightCheck EQ 1 THEN BEGIN ; check done, no problem. this is a whole word.
            PosList = [ PosList , FoundPos]
        ENDIF
    ENDWHILE
    IF KEYWORD_SET(YES_OR_NO) THEN BEGIN
        IF N_ELEMENTS(PosList) EQ 0 THEN RETURN, 0 ELSE RETURN, 1
    ENDIF
    
    RETURN, PosList
END