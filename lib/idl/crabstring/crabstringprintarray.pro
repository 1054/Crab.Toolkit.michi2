; --------------------------------------------------------------------------------------------------------------
; CrabStringPrintArray --- 
;                 CrabStringPrintArray( [1.34343, 2.3434, 4.5454], PREC=2)
;                      -->             "[1.34, 2.34, 4.55]"
; --------------------------------------------------------------------------------------------------------------
FUNCTION CrabStringPrintArray, Array, NELEMENTS=NElements, FORMAT=Format, PRECISION=Precision, SPLITTER=Splitter, $
                                      NoBracket=NoBracket
    
    Dim = SIZE(Array, /DIM)
    
    IF N_ELEMENTS(Dim) NE 1 THEN MESSAGE, "Usage: CrabStringPrintArray([1,2,3,4,5],FORMAT=['F0.2','F0.3'])"
    
    IF Dim LE 0 THEN RETURN, ''
    
    Formats = MAKE_ARRAY(Dim[0],/STRING)
    
    ; switch cases of input array types
    IF SIZE(Array,/TNAME) EQ 'UNDEFINED' THEN BEGIN
        RETURN, ''
    ENDIF ELSE IF SIZE(Array,/TNAME) EQ 'INT' OR SIZE(Array,/TNAME) EQ 'LONG' THEN BEGIN
        DefaultPositiveFormat = '(I0)'
        DefaultNegativeFormat = '(I0)'
        DefaultZeroZeroFormat = '(I0)'
        FOR I=0,N_ELEMENTS(Array)-1 DO BEGIN
            IF N_ELEMENTS(Format) EQ 0 THEN BEGIN
                IF Array[I] GT 0 AND SIZE(DefaultPositiveFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultPositiveFormat
                IF Array[I] LT 0 AND SIZE(DefaultNegativeFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultNegativeFormat
                IF Array[I] EQ 0 AND SIZE(DefaultZeroZeroFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultZeroZeroFormat
            ENDIF ELSE BEGIN
                IF I LE N_ELEMENTS(Format)-1 THEN Formats[I] = Format[I] ELSE Formats[I] = Format[-1]
            ENDELSE
        ENDFOR
    ENDIF ELSE IF SIZE(Array,/TNAME) EQ 'DOUBLE' OR SIZE(Array,/TNAME) EQ 'FLOAT' THEN BEGIN
        IF N_ELEMENTS(Precision) EQ 0 THEN Precision=4
;        DefaultPositiveFormat = '(G'+STRTRIM(Precision+6,2)+'.'+STRTRIM(Precision,2)+')'
;        DefaultNegativeFormat = '(G'+STRTRIM(Precision+7,2)+'.'+STRTRIM(Precision,2)+')'
;        DefaultZeroZeroFormat = '(G'+STRTRIM(Precision+6,2)+'.'+STRTRIM(Precision,2)+')'
        DefaultPositiveFormat = '(G0.'+STRTRIM(Precision,2)+')'
        DefaultNegativeFormat = '(G0.'+STRTRIM(Precision,2)+')'
        DefaultZeroZeroFormat = '(G0.'+STRTRIM(Precision,2)+')'
        FOR I=0,N_ELEMENTS(Array)-1 DO BEGIN
            IF N_ELEMENTS(Format) EQ 0 THEN BEGIN
                IF Array[I] GT 0 AND SIZE(DefaultPositiveFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultPositiveFormat
                IF Array[I] LT 0 AND SIZE(DefaultNegativeFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultNegativeFormat
                IF Array[I] EQ 0 AND SIZE(DefaultZeroZeroFormat,/TNAME) EQ 'STRING' THEN Formats[I]=DefaultZeroZeroFormat
            ENDIF ELSE BEGIN
                IF I LE N_ELEMENTS(Format)-1 THEN Formats[I] = Format[I] ELSE Formats[I] = Format[-1]
            ENDELSE
        ENDFOR
    ENDIF ELSE IF SIZE(Array,/TNAME) EQ 'STRING' THEN BEGIN
        FOR I=0,N_ELEMENTS(Array)-1 DO BEGIN
            IF N_ELEMENTS(Format) EQ 0 THEN BEGIN
                Formats[I]='(A)'
            ENDIF ELSE BEGIN
                IF I LE N_ELEMENTS(Format)-1 THEN Formats[I] = Format[I] ELSE Formats[I] = Format[-1]
            ENDELSE
        ENDFOR
    ENDIF ELSE RETURN, ''
    
    
    ; define splitter
    IF N_ELEMENTS(Splitter) EQ 0 THEN Splitter=", "
    
    
    ; prepare output text
    PrintText = ""
    IF N_ELEMENTS(Dim) EQ 1 THEN BEGIN
        IF NOT KEYWORD_SET(NoBracket) THEN PrintText = PrintText + "["
        IF N_ELEMENTS(NELEMENTS) EQ 0 THEN NELEMENTS=Dim[0]
        IF NELEMENTS GT Dim[0] THEN NELEMENTS=Dim[0]
        FOR I=0,NELEMENTS-1 DO BEGIN
            IF I NE 0 THEN PrintText = PrintText + Splitter
            IF NOT STRMATCH(Formats[I],'\(*\)') THEN Formats[I]='('+Formats[I]+')'
           ;PrintText = PrintText + STRTRIM(STRING(Array[I],FORMAT=Formats[I]),2)
            PrintText = PrintText + STRING(Array[I],FORMAT=Formats[I])
        ENDFOR
        IF NOT KEYWORD_SET(NoBracket) THEN PrintText = PrintText + "]"
    ENDIF
    
    
    RETURN, PrintText
END