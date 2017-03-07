;
; Updated: 2013-03-27
; Updated: 2013-04-06
; Updated: 2013-06-06 : the first character of ColumnHeader is a greek letter
; 
; 
FUNCTION CrabTableCheckValidLine, TempLine, Commente, TempLineId, SkipLine
    IsValidContent = 1
    ; check empty line
    IF IsValidContent EQ 1 THEN BEGIN
        IF STRLEN(TempLine) EQ 0 THEN BEGIN
            IsValidContent=0
        ENDIF
    ENDIF
    ; check commented or not
    IF IsValidContent EQ 1 THEN BEGIN
        IF CrabStringMatch(TempLine,Commente,USE_WILDCARD=2,/YES_OR_NO) THEN BEGIN
            IsValidContent=0
        ENDIF
    ENDIF
    ; Skip lines by the given SkipLine, which can be an array of string patterns, or an array of integer indicating which line number to be skipped, starting from 1. 
    IF IsValidContent EQ 1 AND N_ELEMENTS(SkipLine) GT 0 THEN BEGIN
        IF SIZE(SkipLine,/TNAME) EQ 'STRING' THEN BEGIN
            IsValidContent = CrabStringMatch(TempLine,SkipLine,/YES_OR_NO)
            ;<20160716>FOREACH SkipLinePattern,SkipLine DO BEGIN
            ;<20160716>    IF STRMATCH(TempLine,SkipLinePattern) EQ 1 THEN BEGIN
            ;<20160716>        IsValidContent=0
            ;<20160716>    ENDIF
            ;<20160716>ENDFOREACH
        ENDIF ELSE IF SIZE(SkipLine,/TNAME) EQ 'INT' THEN BEGIN
            FOREACH SkipLineIndex,SkipLine DO BEGIN ; <TODO> SkipLineIndex starts from 0.
                IF TempLineId+1 EQ SkipLineIndex THEN BEGIN
                    IsValidContent=0
                ENDIF
            ENDFOREACH
        ENDIF
    ENDIF
    ; Return
    RETURN, IsValidContent
END
; 
; 
FUNCTION CrabTableReadColumn, FilePath, ColumnHeader, ColumnHeaderBefore = ColumnHeaderBefore, $
                                                      ColumnPositionFound = ColumnPositionFound, $
                                                      ColumnPositionBefore = ColumnPositionBefore, $
                                                      Splitter = Splitter, $
                                                      Commente = Commente, $
                                                      SkipLine = SkipLine, $
                                                      SkipValue = SkipValue, $
                                                      SkipValidLines = SkipValidLines, $
                                                      ReadValidLines = ReadValidLines, $
                                                      NoHeaderLine = NoHeaderLine, $
                                                      FastReadDataBlock = FastReadDataBlock, $
                                                      FirstReadLineContent = FirstReadLineContent, $ ; the first read valid line as an output, might be the header line or first data line
                                                      IncludeHeader = IncludeHeader, $
                                                      Verbose = Verbose
    ; 
    ; Read a column from a well formatted text table contents.
    ; e.g.            # some comments
    ;                 # some comments
    ;                 Source      RA     DEC     Flux     Luminosity
    ;                 NGC6946   123123  123123  123123    123123123123
    ;                 NGC1234   12312   1231       123          123123
    ; This program can read a whole column like Column "Flux" = ['123123','123'].
    ;  
    ;    ;;;<TEST>
    ;    FilePath = 'E:\Working\SpireLines\ListBySource2\Arp220\sfit\Arp220_FromHIPE_SLWC3.CSV'
    ;    ColumnHeader = 'LineName'
    ;
    ; 2013-10-28 added keyword Splitter
    ; 2016-07-14 FastReadDataBlock: set to make it faster by skipping checking comment line in data block lines
    ; 
    
    ;;; Check input
    IF N_ELEMENTS(FilePath) EQ 0 THEN RETURN,[]
    
    ;;; Check whether file exists or not
    IF FILE_TEST(FilePath,/Read) LE 0 THEN BEGIN
        IF KEYWORD_SET(Verbose) THEN PRINT, 'CrabTableReadColumn: '+FilePath+' not found! Exit!'
        RETURN,[]
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(Verbose) THEN PRINT, 'CrabTableReadColumn: '+FilePath+' ready to read!'
    ENDELSE
    
    ;;; Check input column header
    IF SIZE(ColumnHeader,/TNAME) EQ 'STRING' THEN BEGIN
        IF STRLEN(STRTRIM(ColumnHeader,2)) LE 0 THEN BEGIN
            IF KEYWORD_SET(Verbose) THEN PRINT, 'CrabTableReadColumn: Column Header is empty! Exit!'
            CLOSE, FileUnit, /FORCE
            FREE_LUN, FileUnit   ;;; <Corrected> 2013-06-08 FREE_LUN !
            RETURN,[]
        ENDIF
        ColumnIndex = [-1, -1] ; if ColumnIndex is a two-element array, then it's actually a byte range. 
    ENDIF ELSE BEGIN
        ; ColumnHeader is actually ColumnIndex, starting from 1. 
        ColumnIndex = FIX(ColumnHeader) ; if ColumnIndex is a single-element array, then it's the column number index starting from 1. 
    ENDELSE
    
    ;;; Check 
    IF KEYWORD_SET(FastReadDataBlock) THEN BEGIN
        ; TODO: Only works for high version IDL > 8.5.1 ?
        IF KEYWORD_SET(Verbose) THEN PRINT, 'CrabTableReadColumn: FastReadDataBlock set! Assuming that all data lines are continuous and no comment lines in between!'
    ENDIF
    
    ;;; Define Commente
    IF NOT KEYWORD_SET(Commente) THEN BEGIN
        Commente = ['#','====','----']
    ENDIF
    ;<20160716>ELSE BEGIN
    ;<20160716>    IF STRPOS(Commente,'*') GT -1 AND STRMATCH(Commente,'\*') EQ -1 THEN BEGIN
    ;<20160716>        Commente = CrabStringReplace(Commente,'*','\*')
    ;<20160716>    ENDIF
    ;<20160716>ENDELSE
    IF N_ELEMENTS(Commente) EQ 1 THEN Commente = [Commente] ; Commente i.e. the comment line leading mark can be an array
    
    
    
    ;;; Read each line, read leading comment lines starting with '#' as TableHeadLine
    OPENR, FileUnit, FilePath, /GET_LUN
    TableHeadLine  = ''
    TableDataLine = ''
    TableLineId   = 0L
    WHILE (~EOF(FileUnit)) DO BEGIN
        TempLine = ''
        ReadF,FileUnit,TempLine
        ;<20160716>IsValidContent = 1
        ;<20160716>IF IsValidContent EQ 1 AND STRLEN(STRTRIM(TempLine,2))                EQ 0 THEN IsValidContent=0
        ;<20160716>IF IsValidContent EQ 1 AND STRMATCH(TempLine,Commente+'*')            EQ 1 THEN IsValidContent=0
        ;<20160716>IF IsValidContent EQ 1 AND STRMATCH(TempLine,'======*')               EQ 1 THEN IsValidContent=0
        ;<20160716>IF IsValidContent EQ 1 AND STRMATCH(TempLine,'------*')               EQ 1 THEN IsValidContent=0
        ;<20160716>IF IsValidContent EQ 1 AND STRMATCH(STRTRIM(TempLine,2),Commente+'*') EQ 1 THEN IsValidContent=0
        IsValidContent = CrabTableCheckValidLine(TempLine,Commente,TableLineId,SkipLine) ;<20160716> added new function
        ;IF KEYWORD_SET(Verbose) THEN PRINT, TempLine
        ; The first line we meet is the TableHeader <TODO><20131201><DzLiu> this may have some confusions though~~
        IF IsValidContent EQ 1 AND STRLEN(TableHeadLine) EQ 0 THEN BEGIN
            TableHeadLine = TempLine
        ENDIF ELSE BEGIN
          ; Then we skip all commented lines until we met the first valid line ; Then we take all valid contents as TableContens, and prepend the header line
          IF IsValidContent EQ 1 THEN BEGIN
            TableDataLine = TempLine & BREAK
          ENDIF
        ENDELSE
        TableLineId = TableLineId + 1L
    ENDWHILE
    
    
    
    ;;; ColumnHeader是字符串还是整型数
    IF N_ELEMENTS(ColumnIndex) EQ 2 THEN BEGIN
        ;;; Correct Greek Letters In TableHeadLine
        GreekPos = 0
        TestTabHeader = STRMID(TableHeadLine,69,8)
        WHILE GreekPos NE -1 DO BEGIN
            GreekPos = STRPOS(TableHeadLine,'σ') ; <Corrected><2013-03-21> One greek letter has a strlen of 2. 
            IF GreekPos GE 1 THEN $
                TableHeadLine = STRMID(TableHeadLine,0,GreekPos)+'s'+STRMID(TableHeadLine,GreekPos+2,STRLEN(TableHeadLine)-GreekPos-1)
            IF GreekPos EQ 0 THEN $ ; if the first character is a greek letter
                TableHeadLine = 's'+STRMID(TableHeadLine,GreekPos+2,STRLEN(TableHeadLine)-GreekPos-1)
        ENDWHILE
        ;;; Correct Greek Letters In TableHeadLine
        GreekPos = 0
        TestTabData = STRMID(TableDataLine,69,8)
        WHILE GreekPos NE -1 DO BEGIN
            GreekPos = STRPOS(TableDataLine,'σ') ; <Corrected><2013-03-21> One greek letter has a strlen of 2.
            IF GreekPos GE 1 THEN $
                TableDataLine = STRMID(TableDataLine,0,GreekPos)+'s'+STRMID(TableDataLine,GreekPos+2,STRLEN(TableDataLine)-GreekPos-1)
            IF GreekPos EQ 0 THEN $ ; if the first character is a greek letter
                TableDataLine = 's'+STRMID(TableDataLine,GreekPos+2,STRLEN(TableDataLine)-GreekPos-1)
        ENDWHILE
        ;;; Correct Greek Letters In ColumnHeader
        GreekPos = 0
        WHILE GreekPos NE -1 DO BEGIN
            GreekPos = STRPOS(ColumnHeader,'σ') ; <Corrected><2013-03-21> One greek letter has a strlen of 2. 
            IF GreekPos GE 1 THEN $
                ColumnHeader = STRMID(ColumnHeader,0,GreekPos)+'s'+STRMID(ColumnHeader,GreekPos+2,STRLEN(ColumnHeader)-GreekPos-1)
            IF GreekPos EQ 0 THEN $ ; if the first character is a greek letter
                ColumnHeader = 's'+STRMID(ColumnHeader,GreekPos+2,STRLEN(ColumnHeader)-GreekPos-1)
        ENDWHILE
        
        ;;; 在TableHeadLine行中找到ColumnHeader字符.  (20160716:Obsolete:注意TableContents第一行就是TableHeadLine.)
        ColumnFound = CrabStringFindWholeWord(TableHeadLine, ColumnHeader) ; Match the whole word
        IF N_ELEMENTS(ColumnFound) EQ 0 THEN BEGIN
            IF KEYWORD_SET(Verbose) THEN PRINT, 'Looking for "'+ColumnHeader+'" in "'+TableHeadLine+'"'
            PRINT, 'Error! Column Header '+ColumnHeader+' not found! Exit!'
            CLOSE, FileUnit, /FORCE ; column header not found!
            FREE_LUN, FileUnit ; column header not found!
            RETURN, [] ; column header not found!
        ENDIF
        HPosL = -1
        IF N_ELEMENTS(ColumnHeaderBefore) EQ 1 AND SIZE(ColumnHeaderBefore,/TYPE) EQ 7 THEN $
            ColumnBefore = CrabStringFindWholeWord(TableHeadLine, ColumnHeaderBefore)
        IF N_ELEMENTS(ColumnPositionBefore) GE 1 THEN $
            ColumnBefore = LONG(ColumnPositionBefore[0])
        IF N_ELEMENTS(ColumnBefore) GE 1 THEN BEGIN
            ; find the ColumnFound right after the ColumnBefore
            HPosL = ColumnFound[WHERE(ColumnFound-ColumnBefore[0] GE 0)] 
            HPosL = HPosL[0]
        ENDIF ELSE BEGIN
            ; just take the first found position
            HPosL = ColumnFound[0]
        ENDELSE
        HPosR = HPosL + STRLEN(ColumnHeader)-1
        ColumnPositionFound = [HPosL,HPosR] ; output the found position
        ColumnIndex = ColumnPositionFound
    
    
        ;;; IncludeHeader
        IF NOT KEYWORD_SET(NoHeaderLine) THEN BEGIN
            IF NOT KEYWORD_SET(IncludeHeader) THEN BEGIN ; if not set NoHeaderLine and not set IncludeHeader, then drop the first valid line, i.e. TableHeadLine. 
                TableHeadLine = ''
            ENDIF
        ENDIF
        
        
    ENDIF
    
    ;;; Splitter
    IF N_ELEMENTS(Splitter) EQ 0 THEN Splitter = [' '] ; default Splitter is one single space
    IF N_ELEMENTS(Splitter) EQ 1 THEN Splitter = [Splitter]
    ;<20160716>IF N_ELEMENTS(Splitter) GT 1 THEN Splitter = Splitter[0]
    ;<20160716>IF SIZE(Splitter,/TNAME) NE 'STRING' THEN Splitter = STRING(Splitter)
    
    
    
    ;;; 循环TableContents.  (20160716:Obsolete:注意跳掉第一行标题行.)
    ColumnContents = []
    TepmValidLineId = 0L ; SkipValidLines
    TempValidLineCount = 0L ; ReadValidLines
    WHILE (~EOF(FileUnit)) DO BEGIN
        IF NOT KEYWORD_SET(FastReadDataBlock) THEN BEGIN
            IF TableHeadLine NE '' THEN BEGIN ; if there are some undealed data line in memory
                TempLine = TableHeadLine
                TableHeadLine = ''
                TepmValidLineId = TepmValidLineId + 1L
            ENDIF ELSE IF TableDataLine NE '' THEN BEGIN ; if there are some undealed data line in memory
                TempLine = TableDataLine
                TableDataLine = ''
                TepmValidLineId = TepmValidLineId + 1L
            ENDIF ELSE BEGIN
                TempLine = ''
                ReadF,FileUnit,TempLine
                TableLineId = TableLineId + 1L
                TepmValidLineId = TepmValidLineId + 1L
                ; Skip comment line (if not set FastReadDataBlock)
                IsValidContent = CrabTableCheckValidLine(TempLine,Commente,TableLineId,SkipLine) ;<20160716> added new function
                IF IsValidContent EQ 0 THEN CONTINUE
            ENDELSE
            ; Skip valid lines <added><20160730><dzliu>
            IF N_ELEMENTS(SkipValidLines) EQ 1 THEN BEGIN
                ;PRINT, "<TODO><DEBUG><20160731> TepmValidLineId SkipValidLines", TepmValidLineId, SkipValidLines
                IF TepmValidLineId LE SkipValidLines THEN CONTINUE
            ENDIF
            ; Count valid lines <added><20160730><dzliu>
            TempValidLineCount = TempValidLineCount + 1L
            ; Get First valid line content <added><20160730><dzliu>
            IF TempValidLineCount EQ 1L THEN BEGIN
                FirstReadLineContent = TempLine
            ENDIF
            ; Read specified number of lines <added><20160730><dzliu>
            IF N_ELEMENTS(ReadValidLines) EQ 1 THEN BEGIN
                ;PRINT, "<TODO><DEBUG><20160731> TempValidLineCount ReadValidLines", TempValidLineCount, ReadValidLines
                IF TempValidLineCount GT ReadValidLines THEN BREAK
            ENDIF
        ENDIF ELSE BEGIN
            ; Read all rest lines as data block
            TempFileLineNumb = FILE_LINES(FilePath) - TableLineId - 1L
            TempLine = MAKE_ARRAY(TempFileLineNumb,/STRING)
            ReadF,FileUnit,TempLine
            TableLineId = TableLineId + TempFileLineNumb
            IF TableDataLine NE '' THEN BEGIN ; if there are some undealed data line in memory
                TempLine = [TableDataLine, TempLine]
                IF KEYWORD_SET(Verbose) THEN PRINT, "TableDataLine = ", TableDataLine
                TableDataLine = ''
            ENDIF
            IF TableHeadLine NE '' THEN BEGIN ; if there are some undealed data line in memory
                TempLine = [TableHeadLine, TempLine]
                IF KEYWORD_SET(Verbose) THEN PRINT, "TableHeadLine = ", TableHeadLine
                TableHeadLine = ''
            ENDIF
            ; Skip valid lines <added><20160730><dzliu>
            IF N_ELEMENTS(SkipValidLines) EQ 1 THEN BEGIN
                IF SkipValidLines GT N_ELEMENTS(TempLine) THEN MESSAGE, "No line to read after skipping "+STRING(FORMAT='(I0)',SkipValidLines)+" lines!"
                TempLine = TempLine[SkipValidLines:N_ELEMENTS(TempLine)-1L]
            ENDIF
            ; Get First valid line content <added><20160730><dzliu>
            IF TempValidLineCount EQ 0L THEN BEGIN
                FirstReadLineContent = TempLine[0]
            ENDIF
            ; Read specified number of lines <added><20160730><dzliu>
            IF N_ELEMENTS(ReadValidLines) EQ 1 THEN BEGIN
                IF ReadValidLines LT N_ELEMENTS(TempLine) THEN BEGIN
                    TempLine = TempLine[0L:ReadValidLines-1L]
                ENDIF
            ENDIF
            ; Assumed there are no comment line in data block lines if set FastReadDataBlock
        ENDELSE
        ; 
        ;<20160716>;IF STRLEN(TempLine) LT HPosL THEN CONTINUE ; <Corrected><20131125><DzLiu>
        ;<20160716>IF i EQ 0 AND NOT KEYWORD_SET(IncludeHeader) THEN CONTINUE ; <Update><20131201><DzLiu>
        ; 
        ; Read data block line
        IF N_ELEMENTS(ColumnIndex) EQ 2 THEN BEGIN
            ; The input ColumnHeader is a string, 
            ; so the ColumnIndex is a two-element array containing the byte range of the target column, starting from 0. 
            HPosL = ColumnIndex[0]
            HPosR = ColumnIndex[1]
            NewHPosL = CrabStringSearch(TempLine,Splitter,HPosL,/REVERSE_SEARCH,/ONE_ELEMENT,/RETURN_FINAL_POS); & IF NewHPosL EQ -1 THEN NewHPosL=0
            NewHPosR = CrabStringSearch(TempLine,Splitter,HPosR,/ONE_ELEMENT,/RETURN_FINAL_POS); & IF NewHPosR EQ -1 THEN NewHPosR=STRLEN(TempLine)-1
            TempStr  = STRTRIM(STRMID(TempLine,NewHPosL,(NewHPosR-NewHPosL+1)),2)
            IF KEYWORD_SET(Verbose) THEN PRINT, 'Reading data block line ColumnIndex '+STRING(FORMAT='(I0," ",I0, " ")', ColumnIndex[0], ColumnIndex[1])+TempStr
            IF KEYWORD_SET(Verbose) THEN PRINT, 'Reading data block line ColumnIndex '+STRING(FORMAT='(I0," ",I0, " ")', NewHPosL, NewHPosR)+TempStr
            ColumnContents = [ColumnContents,TempStr]
            IF KEYWORD_SET(Verbose) AND TableLineId MOD 1000 EQ 0 THEN PRINT, "Read "+STRING(FORMAT='(I0)',TableLineId)+" lines."
        ENDIF ELSE BEGIN
            ; The input ColumnHeader is a number, 
            ; so the ColumnIndex is a single-element array containing the number index of target column, starting from 1. 
            IF NOT KEYWORD_SET(FastReadDataBlock) THEN BEGIN
                TempStrList = CrabStringSplit(TempLine,Splitter,RemoveBrackets=0)
                IF ColumnIndex-1 LE N_ELEMENTS(TempStrList)-1 THEN TempStr=TempStrList[ColumnIndex-1] ELSE TempStr=""
                ColumnContents = [ColumnContents,TempStr]
            ENDIF ELSE BEGIN
                TempStr = ((STRSPLIT(TempLine,/EXTRACT)).ToArray())[*,ColumnIndex-1] ; LIST 1st dimension is list item, so 2nd dimension is list item subscript. 
                ColumnContents = [ColumnContents,TempStr]
            ENDELSE
            IF KEYWORD_SET(Verbose) THEN PRINT, "CrabTableReadColumn: Read up to the number "+STRING(FORMAT='(I0)',TableLineId)+" line."
        ENDELSE
    ENDWHILE
    
;    Print,'TableHeadLine = ',TableHeadLine
;    Print,' '
;    Print,'ColumnHeader = ',ColumnHeader
;    Print,' '
;    Print,'ColumnContents = ',ColumnContents
    
    ;;; Close
    CLOSE, FileUnit
    CLOSE, FileUnit, /FORCE
    FREE_LUN, FileUnit
    FREE_LUN, FileUnit, /FORCE
    
    ;;; Return
    ReadValidLines = N_ELEMENTS(ColumnContents)
    IF KEYWORD_SET(Verbose) THEN PRINT, "CrabTableReadColumn: Read "+STRING(FORMAT='(I0)',ReadValidLines)+" valid lines."
    Return, ColumnContents
    
END