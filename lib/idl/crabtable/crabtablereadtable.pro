FUNCTION CrabTableReadTable, FilePath,  Formats=Formats, $
                                        Commenter=Commenter, $
                                        Splitter=Splitter, $
                                        AllValues=AllValues, $
                                        Verbose=Verbose, $
                                        ContinuOnError=ContinuOnError
    
    ;;; Check input
    IF N_ELEMENTS(FilePath) LE 0 THEN MESSAGE, "Usage: CrabTableReadTable, FilePath"
    
    ;;; Exam whether file exists or not, then open it
    IF FILE_TEST(FilePath,/Read) LE 0 THEN MESSAGE, "Error! FilePath "+FilePath+" is not found!"
    OPENR, FileUnit, FilePath, /GET_LUN
    
    ;;; Exam input
    IF KEYWORD_SET(KeyNamePattern) EQ 0 AND SIZE(KeyName,/TNAME) NE 'STRING' THEN RETURN,''
    IF KEYWORD_SET(Equaller) EQ 0 THEN Equaller = '='
    
    ;;; Read each line, set first comment lines starting with '#'
    
    ;;; Read each line, skip comment lines starting with '#'
    AllValues = []
    WHILE (~EOF(FileUnit)) DO BEGIN
        TempLine = ''
        TempValue= ''
        ReadF,FileUnit,TempLine
        TempLine = STRTRIM(TempLine,2)
        IF STRLEN(TempLine)                                      EQ 0 THEN CONTINUE
        IF STRMATCH(TempLine,'#*') AND STRMATCH(KeyPattern,'#*') EQ 0 THEN CONTINUE ; sometimes KeyName = "# RMS"
        IF STRMATCH(TempLine,KeyPattern,/FOLD_CASE) THEN BEGIN
            ;;;;;; Print,TempLine
            IF KEYWORD_SET(Verbose) THEN BEGIN
                PRINT, 'CrabTableReadInfo: Found line "'+TempLine+'"'
            ENDIF
            ;;;;;; Exam "KeyName  =  KeyValue  # Comment"
            foundEq  = STRPOS(TempLine,Equaller) ; <Update><20131201><DzLiu> use Equaller to replace the default '=' mark. 
            IF foundEq NE -1 THEN BEGIN
               ;TempValue    = STRTRIM(STRMID(TempLine,foundEq+1),2)
                TempValue    = STRTRIM(STRMID(TempLine,foundEq+STRLEN(Equaller)),2) ; <Corrected><20140716><DzLIU>
                foundComment = STRPOS(TempValue,'#')  ; If there has in-line comment marked after #
                IF foundComment NE -1 THEN BEGIN
                    TempValue = STRTRIM(STRMID(TempValue,0,foundComment),2)
                ENDIF
                IF STRMATCH(TempValue,'"*"') THEN BEGIN
                    TempValue = STRMID(TempValue,1,STRLEN(TempValue)-2) ; remove quotes pair ""
                ENDIF
                IF STRMATCH(TempValue,"'*'") THEN BEGIN
                    TempValue = STRMID(TempValue,1,STRLEN(TempValue)-2) ; remove quotes pair ""
                ENDIF
;                IF KeyWord_Set(DOUBLE) THEN BEGIN
;                    TempValue = TempValue+'D'
;                ENDIF
                AllValues = [AllValues,TempValue]
            ENDIF
        ENDIF
    ENDWHILE
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
    
    
    ;;;Close
    CLOSE, FileUnit
    FREE_LUN, FileUnit
    
    ;;;Warn if empty
    IF KeyWord_Set(GiveWarning) AND N_ELEMENTS(KeyValue) EQ 0 THEN BEGIN
        IF KeyValue EQ "" THEN BEGIN
            MESSAGE, "CrabTableReadInfo: Warning: Key "+KeyName+" was not found or is empty!"
        ENDIF
    ENDIF
    
    ;;;Return
    IF KeyWord_Set(DOUBLE) THEN BEGIN
        IF N_ELEMENTS(AllValues) GT 0 THEN AllValues = DOUBLE(AllValues) ELSE AllValues=[0.0d]
        IF KeyValue NE "" THEN KeyValue = DOUBLE(KeyValue) ELSE KeyValue=0.0d
    ENDIF
    IF KeyWord_Set(INT) THEN BEGIN
        IF N_ELEMENTS(AllValues) GT 0 THEN AllValues = FIX(AllValues)
        KeyValue = FIX(KeyValue)
    ENDIF
    RETURN, KeyValue
END