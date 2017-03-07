PRO CrabTableWriteInfo, FilePath, KeyName, KeyValue, KeyComment = KeyComment, LineChanged = LineChanged
    ; 2013-10-11
    ; Open an ini file (information file), find the key by KeyName and write its KeyValue. 
    ; e.g. "Source  = NGC6946"
    ;      "RA      = 123123123"
    ;      "Dec     = 123123123"
    ;
    
    ;;;Exam whether file exists or not
    IF FILE_TEST(FilePath,/Read) LE 0 THEN RETURN
    OPENU, FileUnit, FilePath, /GET_LUN
    
    ;;;Exam input
    IF STRLEN(STRTRIM(KeyName,2)) LE 0 THEN RETURN
    
    ;;;Exam input
    IF STRLEN(STRTRIM(KeyValue,2)) LE 0 THEN RETURN
    
    ;;;Read each line, skip comment lines starting with '#'
    NewContents = []
    LineChanged = 0
    WhereInsert =-1L
    WHILE (~EOF(FileUnit)) DO BEGIN
        TempLine = ''
        POINT_LUN,-FileUnit,TempPtrPos
        ReadF,FileUnit,TempLine
        
        IF LineChanged NE 0 THEN $
            NewContents=[NewContents,TempLine] ; record those lines after the inserted line.
        
        TempLine = STRTRIM(TempLine,2)
        IF STRLEN(TempLine)                                   EQ 0 THEN CONTINUE
        IF STRMATCH(TempLine,'#*') AND STRMATCH(KeyName,'#*') EQ 0 THEN CONTINUE
        IF STRMATCH(TempLine,KeyName+'*',/FOLD_CASE)          NE 1 THEN CONTINUE
        ;Print,TempLine
        ;;;;;; Exam "KeyName  =  KeyValue  # Comment"
        foundEq  = STRPOS(TempLine,'=')
        IF foundEq NE -1 AND LineChanged EQ 0 THEN BEGIN ; <TODO> we will only proceed with one line.
            NewLine = STRMID(TempLine,0,foundEq+1) + ' ' + KeyValue
            IF SIZE(KeyComment,/TNAME) EQ 'STRING' THEN NewLine = NewLine + ' # ' + KeyComment
            WhereInsert=TempPtrPos                                 ; record the start position of the found will-be-modified line
            NewContents=[NewContents,NewLine] ; recored current modified line content
            LineChanged++
        ENDIF
    ENDWHILE
    
    POINT_LUN,FileUnit,WhereInsert
    TRUNCATE_LUN,FileUnit
    FOREACH NewLine,NewContents DO BEGIN
        PRINTF,FileUnit,NewLine
    ENDFOREACH
    
    
    ;;;Close
    CLOSE, FileUnit
    FREE_LUN, FileUnit
    
    ;;;Return
    RETURN
END