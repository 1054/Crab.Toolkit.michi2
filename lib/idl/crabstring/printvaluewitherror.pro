;FUNCTION print_value_error, InputValue, InputError, LOG=LOG
;    
;    IF N_ELEMENTS(InputValue) NE 1 THEN RETURN, ''
;    IF N_ELEMENTS(InputError) NE 1 THEN RETURN, ''
;    Value = Double(InputValue)
;    Error = Double(InputError)

FUNCTION PrintValueWithError, ValueWithError, LOG=LOG, EXP=EXP, ShowAsExp=ShowAsExp, Format=Format
    
    IF N_ELEMENTS(ValueWithError) NE 2 THEN RETURN, ''
    Value = ValueWithError[0]
    Error = ValueWithError[1]
    
    IF KEYWORD_SET(LOG) THEN BEGIN ;;; the value and error are log values and we need to use Value=10^Value
        Error = Error/Value
        Value = ALOG10(Value)
    ENDIF
    
    IF KEYWORD_SET(EXP) THEN BEGIN ;;; the value and error are log values and we need to use Value=10^Value
;        ValueHigh = 10^(Value+Error)
;        ValueLow = 10^(Value-Error)
        Value = 10^Value
;        Error = MEAN([ValueHigh-Value, Value-ValueLow])
        Error = Error * Value
    ENDIF
    
    
    IF KEYWORD_SET(ShowAsExp) THEN BEGIN
        ValueWithError_Up    = ValueWithError[0] + ValueWithError[1]
        ValueWithError_Low   = ValueWithError[0] - ValueWithError[1]
        ValueWithError_Err_1 = ABS( ALOG10(ValueWithError_Up) - ALOG10(ValueWithError[0])  )
        ValueWithError_Err_2 = ABS( ALOG10(ValueWithError[0]) - ALOG10(ValueWithError_Low) )
        ValueWithError_Err   = MAX([ValueWithError_Err_1,ValueWithError_Err_2]) ; <TODO> why this is smaller than below??
        ValueWithError_Err   = ValueWithError[1]/ValueWithError[0] ; error propagation!
        IF N_ELEMENTS(Format) EQ 0 THEN Format = '(F0.2)'
        RETURN, '10!U'+STRING(FORMAT=Format,ALOG10(ValueWithError[0]))+'!Z(00B1)'+$ ; cgsymbol('+-')+$
                       STRING(FORMAT=Format,       ValueWithError_Err)+'!N'
    ENDIF
    
    ; <TODO> what if very big data?
    IF Error GT Value*1E6 THEN BEGIN
        Error = !VALUES.D_NAN
    ENDIF
    
    FormCode1 = 'F0.2'
    FormCode2 = 'F0.2'
    
    IF Value LT 0.001 THEN BEGIN
        FormCode1 = 'F0.6'
    ENDIF ELSE IF Value GE 0.001 AND Value LT 0.01 THEN BEGIN ; 1 digits
        FormCode1 = 'F0.5'
    ENDIF ELSE IF Value GE 0.01 AND Value LT 0.1 THEN BEGIN ; 1 digits
        FormCode1 = 'F0.4'
    ENDIF ELSE IF Value GE 0.1 AND Value LT 1 THEN BEGIN ; 1 digits
        FormCode1 = 'F0.3'
    ENDIF ELSE IF Value GE 1 AND Value LT 10 THEN BEGIN ; 1 digits
        FormCode1 = 'F0.2'
    ENDIF ELSE IF Value GE 10 AND Value LT 100 THEN BEGIN ; 2 digits
        FormCode1 = 'F0.1'
    ENDIF ELSE IF Value GE 1E2 AND Value LT 1E3 THEN BEGIN ; 3 digits
        FormCode1 = 'I0'
    ENDIF ELSE IF Value GE 1E3 AND Value LT 1E4 THEN BEGIN ; 4 digits
        FormCode1 = 'I0'
    ENDIF
    
    IF Error LT 1 THEN BEGIN
        FormCode2 = 'F0.2'
    ENDIF ELSE IF Error GE 1 AND Error LT 10 THEN BEGIN ; 1 digits
        FormCode2 = 'F0.2'
    ENDIF ELSE IF Error GE 10 AND Error LT 100 THEN BEGIN ; 2 digits
        FormCode2 = 'F0.1'
    ENDIF ELSE IF Error GE 1E2 AND Error LT 1E3 THEN BEGIN ; 3 digits
        FormCode2 = 'I0'
    ENDIF ELSE IF Error GE 1E3 AND Error LT 1E4 THEN BEGIN ; 4 digits
        FormCode2 = 'I0'
    ENDIF
    
    IF N_ELEMENTS(Format) EQ 1 THEN BEGIN
        Format = Format
        IF NOT STRMATCH(Format,'(*,A,*)') THEN Format='('+Format+',A,'+Format+')'
        IF NOT STRMATCH(Format,'(*)') THEN Format='('+Format+')'
    ENDIF ELSE BEGIN
        Format = '('+FormCode1+',A,'+FormCode2+')'
    ENDELSE
    StrToPrint = STRING(FORMAT=Format,Value,'!Z(00B1)', Error) ; cgsymbol('+-')
    RETURN, StrToPrint
END