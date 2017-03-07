; ------------------------------------------------------------------------------------------------------
; CrabStringFindNumber ---
;                      find a number value inside a string
;                      PRINT, CrabStringFindNumber("adfb -30293")
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringFindNumber, Text, TextAfter=TextAfter, TextBefore=TextBefore, DOUBLE=DOUBLE
    ; exam input
    IF SIZE(Text,/TYPE) NE 7 THEN BEGIN
        Print, 'Usage: position = CrabStringFindWholeWord(text,searchstr)'
        IF N_ELEMENTS(Text) EQ 1 THEN Return, '' ELSE RETURN, []
    ENDIF
    ; exam input TextAfter
    IF N_ELEMENTS(TextAfter) EQ 1 AND N_ELEMENTS(Text) GT 1 THEN BEGIN
        TextAfter = REPLICATE(TextAfter[0],N_ELEMENTS(Text))
    ENDIF
    ; exam input TextBefore
    IF N_ELEMENTS(TextBefore) EQ 1 AND N_ELEMENTS(Text) GT 1 THEN BEGIN
        TextBefore = REPLICATE(TextBefore[0],N_ELEMENTS(Text))
    ENDIF
    ; get the SubText which is defined by the text between TextAfter and TextBefore
    PosAfter = REPLICATE(0,N_ELEMENTS(Text))
    IF SIZE(TextAfter,/TNAME) EQ 'STRING' AND N_ELEMENTS(TextAfter) EQ N_ELEMENTS(Text) THEN BEGIN
        FOR i=0,N_ELEMENTS(Text)-1 DO BEGIN
            TextUpCase = STRUPCASE(Text[i])
            TextAfterUpCase = STRUPCASE(TextAfter[i])
            PosAfter[i] = STRPOS(TextUpCase,TextAfterUpCase)
            IF PosAfter[i] EQ -1 THEN PosAfter[i] = STRLEN(Text[i])-1
        ENDFOR
    ENDIF
    PosBefore = REPLICATE(-1,N_ELEMENTS(Text))
    IF SIZE(TextBefore,/TNAME) EQ 'STRING' AND N_ELEMENTS(TextBefore) EQ N_ELEMENTS(Text) THEN BEGIN
        FOR i=0,N_ELEMENTS(Text)-1 DO BEGIN
            TextUpCase = STRUPCASE(Text[i])
            TextBeforeUpCase = STRUPCASE(TextBefore[i])
            PosBefore[i] = STRPOS(TextUpCase,TextBeforeUpCase)
        ENDFOR
    ENDIF
    SubText = REPLICATE('',N_ELEMENTS(Text))
    SubLength = REPLICATE(0,N_ELEMENTS(Text))
    FOR i=0,N_ELEMENTS(Text)-1 DO BEGIN
        IF PosBefore[i] EQ -1 THEN SubLength[i] = STRLEN(Text[i])-PosAfter[i] ELSE SubLength[i] = PosBefore[i] - PosAfter[i] + 1
        SubText[i] = STRMID(Text[i],PosAfter[i],SubLength[i])
    ENDFOR
    ; find the first meet number from the begining of SubText
    SubSign = REPLICATE('',N_ELEMENTS(Text))
    SubNumber = REPLICATE('',N_ELEMENTS(Text))
    SubNumberFound = REPLICATE(0,N_ELEMENTS(Text))
    SubNumberAsDouble = REPLICATE(!VALUES.D_NAN,N_ELEMENTS(Text))
    FOR i=0,N_ELEMENTS(Text)-1 DO BEGIN
        FOR OneCharId=0,STRLEN(SubText[i])-1 DO BEGIN
            OneChar = STRMID(SubText[i],OneCharId,1)
            IF OneChar GE '0' AND OneChar LE '9' THEN BEGIN
                IF SubNumberFound[i] EQ 0 AND OneCharId-1 GE 0 THEN BEGIN ; judge the sign of first meet number
                    IF STRMID(SubText[i],OneCharId-1,1) EQ '-' THEN BEGIN
                        SubSign[i]='-'
                        SubNumber[i]='-'
                    ENDIF
                ENDIF
                IF SubNumberFound[i] EQ 0 THEN SubNumberFound[i] = 1
                IF SubNumberFound[i] EQ 1 THEN SubNumber[i] = SubNumber[i] + OneChar
            ENDIF ELSE IF OneChar EQ 'E' OR OneChar EQ 'e' THEN BEGIN
                IF SubNumberFound[i] EQ 1 THEN SubNumber[i] = SubNumber[i] + OneChar
            ENDIF ELSE IF SubNumberFound[i] EQ 1 THEN BREAK
        ENDFOR
        ; BREAK to here
        IF SubNumber[i] NE '' THEN SubNumberAsDouble[i] = DOUBLE(SubNumber[i])
    ENDFOR
    ; if only 1 element then return scalar value
    IF N_ELEMENTS(Text) EQ 1 THEN BEGIN
        SubNumber = SubNumber[0]
        SubNumberFound = SubNumberFound[0]
        SubNumberAsDouble = SubNumberAsDouble[0]
    ENDIF
    ; if required double type then return double type value 
    IF KEYWORD_SET(DOUBLE) THEN BEGIN
        RETURN, SubNumberAsDouble
    ENDIF
    RETURN, SubNumber
END