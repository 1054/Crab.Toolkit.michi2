; ----------------------------------------------------------------------------------------------------------------------------------------
; CrabStringSearch ---
;                 We can search two string array. 
; 2016-07-16      
;                 RETURN_FINAL_POS: return the final pos if not found
; 
; ----------------------------------------------------------------------------------------------------------------------------------------
FUNCTION CrabStringSearch, Text, SearchStr, StartPos, REVERSE_SEARCH=REVERSE_SEARCH, FOLD_CASE=FOLD_CASE, ONE_ELEMENT=ONE_ELEMENT,$ 
                                                      RETURN_FINAL_POS=RETURN_FINAL_POS
    ; exam input
    IF SIZE(Text,/TYPE) NE 7 OR SIZE(SearchStr,/TYPE) NE 7 THEN BEGIN
        Print, 'Usage: pos = CrabStringSearch(text,searchstr)'
        Return, []
    ENDIF
    ; exam StartPos
    IF N_ELEMENTS(StartPos) EQ 0 THEN PosStart = 0 ELSE PosStart = StartPos
    IF N_ELEMENTS(Text) GT 1 AND N_ELEMENTS(PosStart) EQ 1 THEN PosStart = REPLICATE(PosStart,N_ELEMENTS(Text))
    IF N_ELEMENTS(Text) GT N_ELEMENTS(PosStart) THEN PosStart = [ PosStart, REPLICATE(0,N_ELEMENTS(Text)-N_ELEMENTS(PosStart)) ]
    ; search each input element
    IF N_ELEMENTS(Text) EQ 1 AND N_ELEMENTS(SearchStr) EQ 1 THEN BEGIN
        PosFound = STRPOS(Text, SearchStr, PosStart, REVERSE_SEARCH=REVERSE_SEARCH)
        IF PosFound EQ -1 THEN BEGIN
            IF KEYWORD_SET(RETURN_FINAL_POS) THEN BEGIN
                IF KEYWORD_SET(REVERSE_SEARCH) THEN BEGIN
                    PosFound = 0
                ENDIF ELSE BEGIN
                    PosFound = STRLEN(Text)-1
                ENDELSE
            ENDIF
        ENDIF
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(RETURN_FINAL_POS) THEN BEGIN
            ; return the final pos if not found
            IF KEYWORD_SET(REVERSE_SEARCH) THEN BEGIN
                PosFound = MAKE_ARRAY(N_ELEMENTS(Text),/INTEGER,VALUE=0)
            ENDIF ELSE BEGIN
                PosFound = MAKE_ARRAY(N_ELEMENTS(Text),/INTEGER,VALUE=STRLEN(Text)-1)
            ENDELSE
        ENDIF ELSE BEGIN
            PosFound = MAKE_ARRAY(N_ELEMENTS(Text),/INTEGER,VALUE=-1)
        ENDELSE
        FOR IDRow=0,N_ELEMENTS(Text)-1 DO BEGIN
            IF STRTRIM(Text[IDRow],2) EQ '' THEN CONTINUE ;<Added><20160716>
            IF STRLEN(Text[IDRow]) LE PosStart[IDRow] THEN CONTINUE ;<Added><20160716>
            FOR IDCol=0,N_ELEMENTS(SearchStr)-1 DO BEGIN
                TmpFound = STRPOS(Text[IDRow], SearchStr[IDCol], PosStart[IDRow], REVERSE_SEARCH=REVERSE_SEARCH)
                IF TmpFound GE 0 AND PosFound[IDRow] EQ -1 THEN PosFound[IDRow] = TmpFound
                IF KEYWORD_SET(REVERSE_SEARCH) THEN BEGIN
                    IF TmpFound GE 0 AND TmpFound GT PosFound[IDRow] THEN PosFound[IDRow] = TmpFound
                ENDIF ELSE BEGIN
                    IF TmpFound GE 0 AND TmpFound LT PosFound[IDRow] THEN PosFound[IDRow] = TmpFound
                ENDELSE
            ENDFOR
        ENDFOR
    ENDELSE
    ; return PosFound
    IF KEYWORD_SET(ONE_ELEMENT) THEN RETURN, PosFound[0]
    Return, PosFound
END