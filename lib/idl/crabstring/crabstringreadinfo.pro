; ------------------------------------------------------------------------------------------------------
; CrabStringReadInfo ---
;                 treat one str line as an "Key = Value # comment" format, and read the value.  
;                 e.g. input "* Column density [cm-2]:  1.000E+20" 
;                 PRINT, CrabStringReadInfo(['* Column density [cm-2]:  1.000E+20', '* T(background)     [K]:    2.730'], $
;                                            '* Column density', Equaller=':', /Double)
;                 then output "1.000E+20"
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringReadInfo,   StrArray,  KeyName,  DOUBLE=DOUBLE, $
                                               Equaller=Equaller, $
                                              Commenter=Commenter, $
                                               Splitter=Splitter, $
                                               AllValues=AllValues, $
                                               RemoveSuffix=RemoveSuffix, $
                                               KeyNamePattern=KeyNamePattern, $
                                               Verbose=Verbose
    
    ;;; Exam input
    IF KEYWORD_SET(KeyNamePattern) EQ 0 AND SIZE(KeyName,/TNAME) NE 'STRING' THEN RETURN,''
    IF KEYWORD_SET(Equaller) EQ 0 THEN Equaller = '='
    IF KEYWORD_SET(Commenter) EQ 0 THEN Commenter = '#'
    
    
    ;;; Convert KeyName to KeyPattern
    ;;; use STRMATCH to match one text line with the KeyPattern rather than the KeyName
    IF KEYWORD_SET(KeyNamePattern) THEN BEGIN
        IF SIZE(KeyNamePattern,/TNAME) EQ 'STRING' THEN KeyPattern = KeyNamePattern
    ENDIF ELSE BEGIN
        KeyPattern = KeyName
        ;;;Take care of the single char [
        IF STRMATCH(KeyPattern,'*\[*') THEN KeyPattern = CrabStringReplace(KeyPattern,'[','\[')
        ;;;Take care of the single char ]
        IF STRMATCH(KeyPattern,'*\]*') THEN KeyPattern = CrabStringReplace(KeyPattern,']','\]')
        ;;;Take care of the single char *
        IF STRMATCH(KeyPattern,'*[*]*') THEN KeyPattern = CrabStringReplace(KeyPattern,'*','[*]') ; must after [ and ]
        KeyPattern = KeyPattern + '*'
    ENDELSE
    
    
    ;;; Read each line, skip comment lines starting with '#'
    AllValues = []
    FOREACH StrLine, StrArray DO BEGIN
        TempLine = ''
        TempValue= ''
        ReadS,StrLine,TempLine
        TempLine = STRTRIM(TempLine,2)
        IF STRLEN(TempLine)                                      EQ 0 THEN CONTINUE
        IF STRMATCH(TempLine,Commenter+'*') AND STRMATCH(KeyPattern,Commenter+'*') EQ 0 THEN CONTINUE ; sometimes KeyName = "# RMS"
        IF STRMATCH(TempLine,KeyPattern,/FOLD_CASE) THEN BEGIN
            ;;;;;; Print,TempLine
            IF KEYWORD_SET(Verbose) THEN BEGIN
                PRINT, 'CrabStringReadInfo: Found line "'+TempLine+'"'
            ENDIF
            ;;;;;; Exam "KeyName  =  KeyValue  # Comment"
            foundEq  = STRPOS(TempLine,Equaller) ; <Update><20131201><DzLiu> use Equaller to replace the default '=' mark. 
            IF foundEq NE -1 THEN BEGIN
                TempValue    = STRTRIM(STRMID(TempLine,foundEq+1),2)
                foundComment = STRPOS(TempValue,Commenter)  ; If there has in-line comment marked after #
                IF foundComment NE -1 THEN BEGIN
                    TempValue = STRTRIM(STRMID(TempValue,0,foundComment),2)
                ENDIF
;                IF KeyWord_Set(DOUBLE) THEN BEGIN
;                    TempValue = TempValue+'D'
;                ENDIF
                AllValues = [AllValues,TempValue]
            ENDIF
        ENDIF
    ENDFOREACH
    ;;; If nothing found
    IF N_ELEMENTS(AllValues) EQ 0 THEN BEGIN 
        KeyValue = ""
    ;;; If found then get first value
    ENDIF ELSE BEGIN
        KeyValue = AllValues[0]
        ;;; If set removesuffix then remove suffix, e.g. KeyValue = "300 [Km/s]", RemoveSuffix="[Km/s]". 
        IF SIZE(RemoveSuffix,/TNAME) EQ 'STRING' THEN BEGIN
            KeyValue = CrabStringReplace(KeyValue, RemoveSuffix, '')
        ENDIF
        ;;; If set output as an array
        IF KEYWORD_SET(Splitter) THEN BEGIN
            IF SIZE(Splitter,/TNAME) EQ 'STRING' THEN BEGIN
                KeyValue = CrabStringSplit(KeyValue,Splitter=Splitter,/RemoveBrackets,/Compress)
            ENDIF ELSE BEGIN
                KeyValue = CrabStringSplit(KeyValue,Splitter=',',/RemoveBrackets,/Compress)
            ENDELSE
        ENDIF
    ENDELSE
    
    IF KeyWord_Set(DOUBLE) THEN BEGIN
        AllValues = DOUBLE(AllValues)
        RETURN, DOUBLE(KeyValue)
    ENDIF
    RETURN, KeyValue
    
END