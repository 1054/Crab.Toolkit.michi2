; Print the content of a table to a file. 
; 
PRO CrabTablePrintF, FilePathOrFileUnit, TableDataArray, Format=Format, APPEND=APPEND, $
                                                         PrependColumn1=PrependColumn1, PrependFormat1=PrependFormat1, $
                                                         PrependColumn2=PrependColumn2, PrependFormat2=PrependFormat2, $
                                                         PrependColumn3=PrependColumn3, PrependFormat3=PrependFormat3
    
    ; check table data
    IF SIZE(TableDataArray,/N_DIM) NE 2 THEN BEGIN
      MESSAGE, 'CrabTablePrintF: TableDataArray is not two-dimension! '
    ENDIF
    
    ; check output file
    IF SIZE(FilePathOrFileUnit,/TNAME) EQ 'UNDEFINED' THEN BEGIN
        MESSAGE, 'CrabTablePrintF: FilePathOrFileUnit is not set! '
    ENDIF ELSE IF SIZE(FilePathOrFileUnit,/TNAME) EQ 'LONG' THEN BEGIN
        OutputFileUnit = FilePathOrFileUnit
    ENDIF ELSE IF SIZE(FilePathOrFileUnit,/TNAME) EQ 'STRING' THEN BEGIN
        OutputFilePath = FilePathOrFileUnit
        OPENW, OutputFileUnit, OutputFilePath, /GET_LUN, APPEND=APPEND
    ENDIF ELSE BEGIN
        MESSAGE, 'CrabTablePrintF: FilePathOrFileUnit is invalid! '+FilePathOrFileUnit
    ENDELSE
    
    ; check table dimension
    TableDim = SIZE(TableDataArray,/DIM)
    
    ; check format
    IF SIZE(Format,/TNAME) EQ 'STRING' THEN BEGIN
        OutputFormat = Format
    ENDIF
    
    ; loop each item
    FOR YId = 0, TableDim[1]-1 DO BEGIN
        OutputLine = ''
        IF N_ELEMENTS(PrependColumn1) EQ TableDim[1] THEN BEGIN
            OutputLine = OutputLine + STRING(FORMAT=PrependFormat1,PrependColumn1[YId])
        ENDIF
        IF N_ELEMENTS(PrependColumn2) EQ TableDim[1] THEN BEGIN
            OutputLine = OutputLine + STRING(FORMAT=PrependFormat2,PrependColumn2[YId])
        ENDIF
        IF N_ELEMENTS(PrependColumn3) EQ TableDim[1] THEN BEGIN
            OutputLine = OutputLine + STRING(FORMAT=PrependFormat3,PrependColumn3[YId])
        ENDIF
        FOR XId = 0, TableDim[0]-1 DO BEGIN
            OutputText = ''
            OutputText = STRING(FORMAT=OutputFormat,TableDataArray[XId,YId])
            OutputLine = OutputLine + OutputText
        ENDFOR
        PRINTF, OutputFileUnit, OutputLine
    ENDFOR
    
    ;;;Close
    CLOSE, OutputFileUnit
    FREE_LUN, OutputFileUnit
    
END