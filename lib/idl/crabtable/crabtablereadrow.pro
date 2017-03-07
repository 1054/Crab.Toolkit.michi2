;
; Updated: 2013-03-27
; Updated: 2013-06-29 : FREE_LUN !!!
; 
FUNCTION CrabTableReadRow, FilePath, RowText, RowID=RowID
    
    
    ;;;Exam whether file exists or not
    IF FILE_TEST(FilePath,/Read) LE 0 THEN RETURN,[]
    OPENR, FileUnit, FilePath, /GET_LUN
    
    ;;;Exam input
    IF N_ELEMENTS(RowText) LT 1 OR SIZE(RowText,/TYPE) NE 7 THEN BEGIN
        RowID = -1
        CLOSE, FileUnit, /FORCE
        FREE_LUN, FileUnit   ;;; <Corrected> 2013-06-29 FREE_LUN !
        RETURN,[]
    ENDIF
    
    ;RowContentSearchPattern should be a pattern for strmatch, so if the input is "SOMETEXT", we will make it "*SOMETEXT*"
    IF STRPOS(RowText,'*') EQ -1 THEN BEGIN
        RowContentSearchPattern = '*'+RowText+'*'
    ENDIF ELSE RowContentSearchPattern = RowText
    
    ;;;Read each line, skip comment lines starting with '#'
    MatchedRows = []
    MatchedRowIDs = []
    TempID = -2 ; note that CrabTableReadColumn will skip the first ColumnHeader line, so here TempID starts with -2. 
    WHILE (~EOF(FileUnit)) DO BEGIN
        TempLine = ''
        ReadF,FileUnit,TempLine
        IF STRLEN(STRTRIM(TempLine,2)) EQ 0 THEN CONTINUE
        IF STRMATCH(RowContentSearchPattern,'#*') EQ 0 AND STRMATCH(TempLine,'#*') EQ 1 THEN CONTINUE
        IF STRMATCH(RowContentSearchPattern,'#*') EQ 0 AND STRMATCH(STRTRIM(TempLine,2),'#*') EQ 1 THEN CONTINUE
        IF STRMATCH(TempLine,'======*') EQ 1 THEN CONTINUE
        IF STRMATCH(TempLine,'------*') EQ 1 THEN CONTINUE
        TempID = TempID + 1 ; after skipping all invalid lines, here we got valid lines. First valid line is the column header line, with ID=-1. 
        ;;;;;; Match Row Header
        IF STRMATCH(TempLine,RowContentSearchPattern) EQ 0 THEN CONTINUE
        MatchedRows = [MatchedRows,TempLine]
        MatchedRowIDs = [MatchedRowIDs,TempID] ; 
    ENDWHILE
    
    ;;;Close
    CLOSE, FileUnit
    CLOSE, FileUnit, /FORCE
    FREE_LUN, FileUnit
    FREE_LUN, FileUnit, /FORCE
    
    ;;;
    RowID = MatchedRowIDs
    RETURN, MatchedRows
    
END