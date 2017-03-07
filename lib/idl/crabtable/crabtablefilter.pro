; Filter the content of a text file. 
; 
PRO CrabTableFilter, FilePath, Filters, REGEX=REGEX, FOLD_CASE=FOLD_CASE, OutputFilePath=OutputFilePath
    
    ; check input file
    IF NOT CrabFileCheck(FilePath,/NOWARNING) THEN BEGIN
        MESSAGE, 'CrabTableFilter: FilePath is invalid! '+FilePath
    ENDIF
    IF N_ELEMENTS(Filters) EQ 0 THEN BEGIN
        MESSAGE, 'CrabTableFilter: Filter is empty!'
    ENDIF
    
    ; prepare output file
    IF NOT CrabFileCheck(OutputFilePath,/NOEXIST) THEN BEGIN
        OutputFilePath=FILE_DIRNAME(FilePath,/MARK_DIR)+FILE_BASENAME(FilePath)+'.filtered'
    ENDIF
    
    ; open output file to write
    OPENW, OutputFileUnit, OutputFilePath, /GET_LUN
    
    ; open input file to read
    OPENR, FileUnit, FilePath, /GET_LUN
    TempLine = ''
    WHILE(~EOF(FileUnit)) DO BEGIN
        ; read each line
        ReadF,FileUnit,TempLine
        ; match the filters
        IF CrabStringMatch(TempLine,Filters,/USE_WILDCARD,/YES_OR_NO,FOLD_CASE=FOLD_CASE) THEN BEGIN
            PRINTF,OutputFileUnit,TempLine
        ENDIF
    ENDWHILE
    
    ;;;Close
    CLOSE, FileUnit
    FREE_LUN, FileUnit
    CLOSE, OutputFileUnit
    FREE_LUN, OutputFileUnit
    
END