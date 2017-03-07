; ----------------------------------------------------------------------------------------------------------------------------------------
; CrabStringMatch ---
;                 We can match two string array. 
;                 
;                 By setting /YES_OR_NO, if there has one matched string in two array, we will return 1.
;                 e.g. PRINT, CrabStringMatch(['NGC1068','NGC0253','NGC1365-NE'],$
;                                              'NGC1365','NGC5194','NGC4418'],$
;                                              /YES_OR_NO)
;                 ===> 1
;                 
;                 By setting /USE_WILDCARD, we will add a pair of '*' around each item of SearchStr
;                 e.g. PRINT, CrabStringMatch(['NGC1068','NGC0253','NGC1365-NE'],$
;                                              'NGC1365',$
;                                              /USE_WILDCARD)
;                 ===> [0,0,1]
; 
; 2014-12-11      ONE_DIMENSION
; 2016-07-16      USE_WILDCARD set 2 to only append wildcard, 3 to only prepend wildcard, other values to both append and prepend. 
; 
; ----------------------------------------------------------------------------------------------------------------------------------------
FUNCTION CrabStringMatch, Text, SearchStr, FOLD_CASE=FOLD_CASE, YES_OR_NO=YES_OR_NO, USE_WILDCARD=USE_WILDCARD, ONE_DIMENSION=ONE_DIMENSION
    ; exam input
    IF SIZE(Text,/TYPE) NE 7 OR SIZE(SearchStr,/TYPE) NE 7 THEN BEGIN
        IF KEYWORD_SET(YES_OR_NO) THEN BEGIN
            RETURN, 0
        ENDIF ELSE BEGIN
            Print, 'Usage: check = CrabStringMatch(text,searchstr)'
            Return, []
        ENDELSE
    ENDIF
    ; search each input element
    IF N_ELEMENTS(Text) EQ 1 AND N_ELEMENTS(SearchStr) EQ 1 THEN BEGIN
        OneSStr = SearchStr
        IF KEYWORD_SET(USE_WILDCARD) THEN BEGIN
            ;<Modified><20160716><DzLIU> IF KEYWORD_SET(USE_WILDCARD) THEN OneSStr='*'+OneSStr+'*'
            IF USE_WILDCARD EQ 2 THEN BEGIN
                OneSStr=OneSStr+'*' ; <Added><20160716><DzLIU> USE_WILDCARD Only Append Wild Card!
            ENDIF ELSE IF USE_WILDCARD EQ 3 THEN BEGIN
                OneSStr='*'+OneSStr ; <Added><20160716><DzLIU> USE_WILDCARD Only Prepend Wild Card!
            ENDIF ELSE BEGIN
                OneSStr='*'+OneSStr+'*' ; <Corrected><20140905><DzLIU> USE_WILDCARD!
            ENDELSE
        ENDIF
        Matched = STRMATCH(Text,OneSStr)
    ENDIF ELSE BEGIN
        Matched = MAKE_ARRAY(N_ELEMENTS(SearchStr), N_ELEMENTS(Text),/INTEGER,VALUE=0)
        FOR IDRow=0,N_ELEMENTS(Text)-1 DO BEGIN
            IF STRTRIM(Text[IDRow],2) EQ '' THEN CONTINUE ;<Added><20160716>
            FOR IDCol=0,N_ELEMENTS(SearchStr)-1 DO BEGIN
                OneBStr = Text[IDRow]
                OneSStr = SearchStr[IDCol]
                IF KEYWORD_SET(USE_WILDCARD) THEN BEGIN
                    ;<Modified><20160716><DzLIU> IF KEYWORD_SET(USE_WILDCARD) THEN OneSStr='*'+OneSStr+'*'
                    IF USE_WILDCARD EQ 2 THEN BEGIN
                        OneSStr=OneSStr+'*' ; <Added><20160716><DzLIU> USE_WILDCARD Only Append Wild Card!
                    ENDIF ELSE IF USE_WILDCARD EQ 3 THEN BEGIN
                        OneSStr='*'+OneSStr ; <Added><20160716><DzLIU> USE_WILDCARD Only Prepend Wild Card!
                    ENDIF ELSE BEGIN
                        OneSStr='*'+OneSStr+'*' ; <Corrected><20140905><DzLIU> USE_WILDCARD!
                    ENDELSE
                ENDIF
                Matched[IDCol,IDRow] = STRMATCH(OneBStr, OneSStr, FOLD_CASE=FOLD_CASE)
            ENDFOR
        ENDFOR
    ENDELSE
    IF KEYWORD_SET(ONE_DIMENSION) THEN BEGIN
        Matched = TOTAL(Matched,1,/PRESERVE_TYPE)
    ENDIF
    IF KEYWORD_SET(YES_OR_NO) THEN BEGIN
        IF TOTAL(Matched) GT 0 THEN Matched = 1 ELSE Matched = 0
    ENDIF
    ; return matched index
    Return, Matched
    ; e.g. PRINT, CrabStringMatch(['P1','P2'],['M1','P1','P2','P3'])
    ;             0   1   0   0
    ;             0   0   1   0
END