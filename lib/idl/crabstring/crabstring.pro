; Crab String
; CrabStringReplace
; CrabStringClean
; CrabStringSplit
; CrabStringFindWholeWord
; CrabStringFindNumber
;        Aim: Find out a number inside of any text
;        e.g. PRINT, CrabStringFindNumber("what if we input a -3.1415926 here?")
;             -3.1415926
; 
; CrabStringMatch
; CrabStringPrintArray
;        Aim: Print an array easily with certain format, and with a pair of bracket with it. 
;        e.g. PRINT, CrabStringPrintArray([3.14,3.14159,3.141])
;             [3.140, 3.142, 3.141]
; 
; CrabStringPrintFloat
; CrabStringPrintExpVar
; CrabStringPrintReturn
; CrabStringPrintSigma
; CrabStringPrintChar


; ------------------------------------------------------------------------------------------------------
; CrabString ---
;                      
; ------------------------------------------------------------------------------------------------------
PRO CrabString
    
    RETURN
END



; ------------------------------------------------------------------------------------------------------
; CrabStringPrintExpVar ---
;                 CrabStringPrintExpVar( 1.34343, PREC=3)
;                      -->              "1.343E+10"
;                 CrabStringPrintExpVar( 1.34343, PREC=3, /SIMPLE)
;                      -->              "1.34e10"
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringPrintExpVar, Number, FORMAT=FORMAT, PRECISION=PRECISION, SIMPLE=SIMPLE
    IF SIZE(Number,/TNAME) NE 'FLOAT' AND SIZE(Number,/TNAME) NE 'DOUBLE' THEN RETURN, ''
    IF N_ELEMENTS(PRECISION) EQ 1 THEN BEGIN
        IF PRECISION LE 0 THEN PRECISION=3
    ENDIF ELSE BEGIN
        PRECISION = 3
    ENDELSE
    IF N_ELEMENTS(SIMPLE) EQ 1 THEN BEGIN 
        FORMAT = '(F0.'+STRTRIM(PRECISION,2)+',"e",I0)'
    ENDIF ELSE BEGIN
        FORMAT = '(F0.'+STRTRIM(PRECISION,2)+',"E",I+03)'
    ENDELSE
    IF Number[0] NE 0.0 THEN BEGIN
        Numb_N = FIX(ALOG10(ABS(Number[0])))
        Numb_A = Number[0]/(10.0D^Numb_N)
    ENDIF ELSE BEGIN
        Numb_N = 0
        Numb_A = 0.0
    ENDELSE
    PrintText = ""
    PrintText = STRTRIM(STRING(Numb_A,Numb_N,FORMAT=FORMAT),2)
    RETURN, PrintText
END



; ------------------------------------------------------------------------------------------------------
; CrabStringPrintReturn ---
;                 CrabStringPrintReturn
;                      -->             "CR"
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringPrintReturn
    Str = ''
    IF STRMATCH(!Version.OS_FAMILY,'*Windows*',/FOLD_CASE) THEN Str = STRING(BYTE(['0A'xL,'0D'xL]))
    IF STRMATCH(!Version.OS_FAMILY,'*unix*',/FOLD_CASE) THEN Str = STRING(BYTE(['0A'xL]))
    RETURN, Str
END


FUNCTION CrabStringPrintSigma
    Str = ''
    IF STRMATCH(!Version.OS_FAMILY,'Windows',/FOLD_CASE) THEN Str = STRING(BYTE(['A6'xL,'D2'xL]))
    IF STRMATCH(!Version.OS_FAMILY,'unix',/FOLD_CASE) THEN Str = STRING(BYTE(['CF'xL,'83'xL]))
    RETURN, Str
END


FUNCTION CrabStringPrintPlusMinus
    Str = ''
    IF STRMATCH(!Version.OS_FAMILY,'Windows',/FOLD_CASE) THEN Str = STRING(BYTE(['A1'xL,'C0'xL]))
    IF STRMATCH(!Version.OS_FAMILY,'unix',/FOLD_CASE) THEN Str = STRING(BYTE(['C2'xL,'B1'xL]))
    RETURN, Str
END


FUNCTION CrabStringPrintChar, CharText
    IF STRUPCASE(CharText) EQ "RETURN" THEN BEGIN
        RETURN, CrabStringPrintReturn()
    ENDIF ELSE IF STRUPCASE(CharText) EQ "SIGMA" THEN BEGIN
        RETURN, CrabStringPrintSigma()
    ENDIF ELSE IF STRUPCASE(CharText) EQ "RIGHTARROW" THEN BEGIN
        RETURN, STRING(BYTE(['A1'xL,'FA'xL]))
    ENDIF ELSE IF STRUPCASE(CharText) EQ "" THEN BEGIN
        RETURN, ""
    ENDIF
END


