; 
; Print the content of a table to a file. 
; 
; Major modification: 
;     2014-07-14: changed almost the whole method
;     2014-10-30: TMPcode="I"; TMPwide=STRING(FORMAT='("0",I0)',TMPmax+2) changed to TMPwide=STRING(FORMAT='(I0)',TMPmax+2)
;                 add TMPsign
;     2014-12-11: consider the length of header
;                 do not print header when append
; 
; 
PRO CrabTablePrintC, FilePathOrFileUnit, Column1, Column2, Column3, Column4, Column5, Column6, Column7, Column8, $
                                                           Format=Format, Separator=Separator, APPEND=APPEND, Header=Header
    
    ; check input arguments (keywords are not included in N_PARAMS())
    IF N_PARAMS() LE 1 THEN BEGIN
        PRINT, 'Usage: CrabTablePrintC, "SomeTable.CSV", MyColumnData1, MyColumnData2, MyColumnData3'
    ENDIF
    
    IF N_ELEMENTS(Column1) EQ 0 THEN BEGIN
        PRINT, 'CrabTablePrintC: Error! Input column 1 is empty!'
        RETURN
    ENDIF
    
    ; Check input Format
    IF N_ELEMENTS(Format) EQ 0 THEN BEGIN
        ; Should decide column format by program!
        CCFcode = [] ; F
        CCFsign = [] ; +/-
        CCFwide = [] ; 10
        CCFprec = [] ; .3
        CCFhead = [] ; header
    ENDIF ELSE BEGIN
        ; Should use user-defined format!
        IF SIZE(Format,/TNAME) NE 'STRING' THEN BEGIN
            MESSAGE, 'CrabTablePrintC: Format should be a string!'
        ENDIF
        CCFcode = [] ; F
        CCFsign = [] ; +/-
        CCFwide = [] ; 10
        CCFprec = [] ; .3
        CCFhead = [] ; header
        ; Split Format for each column
        FOREACH FormatVar, Format DO BEGIN
            ; Split Format for each column
            IF STRMATCH(FormatVar,'*,*') THEN BEGIN 
                CCFtext = CrabStringSplit(FormatVar,Splitter=',',/RemoveBrackets,/Compress)
            ENDIF ELSE BEGIN
                CCFtext = CrabStringClean(FormatVar,TextsToRemove=['(',')'])
            ENDELSE
            FOR j=0,N_ELEMENTS(CCFtext)-1 DO BEGIN
                ; CCFtext could be  '5g10.3'  or  'F10.3'
                TMPtext = CCFtext[j]
                ; -- repeat format code
                TMPpeat = 1
                IF STREGEX(TMPtext,'^[0-9]+',/BOOLEAN) THEN BEGIN ; repeating format -- e.g. 3F10.3 = F10.3,F10.3,F10.3
                  TMPpeat = STREGEX(TMPtext,'^[0-9]+',/EXT)
                  TMPtext = STRMID(TMPtext,STRLEN(TMPpeat))
                ENDIF
                ; -- format code g/F/E/G/A
                TMPcode = ""
                IF STREGEX(TMPtext,'^[a-zA-Z]+',/BOOLEAN) THEN BEGIN
                  TMPcode = STREGEX(TMPtext,'^[a-zA-Z]+',/EXT)
                  TMPtext = STRMID(TMPtext,STRLEN(TMPcode))
                ENDIF
                ; -- format sign +/-
                TMPsign = ""
                IF STREGEX(TMPtext,'^[+-]+',/BOOLEAN) THEN BEGIN
                  TMPsign = STREGEX(TMPtext,'^[+-]+',/EXT)
                  TMPtext = STRMID(TMPtext,STRLEN(TMPsign))
                ENDIF
                ; -- format wide 10.
                TMPwide = ""
                IF STREGEX(TMPtext,'^[0-9]+',/BOOLEAN) THEN BEGIN
                  TMPwide = STREGEX(TMPtext,'^[0-9]+',/EXT)
                  TMPtext = STRMID(TMPtext,STRLEN(TMPwide))
                ENDIF
                ; -- format precise .3
                TMPprec = ""
                IF STRMATCH(TMPtext,'.*') THEN BEGIN
                  TMPtext = STRMID(TMPtext,1)
                  TMPprec = STREGEX(TMPtext,'^[0-9]+',/EXT) & TMPprec="."+TMPprec
                  TMPtext = STRMID(TMPtext,STRLEN(TMPtext))
                ENDIF
                ; -- save into array
                FOR k=0,TMPpeat-1 DO BEGIN
                  CCFcode = [CCFcode, TMPcode]
                  CCFsign = [CCFsign, TMPsign]
                  CCFwide = [CCFwide, TMPwide]
                  CCFprec = [CCFprec, TMPprec]
                ENDFOR
            ENDFOR
        ENDFOREACH
        ; PRINT, CCFcode, CCFsign, CCFwide, CCFprec
    ENDELSE
    
    
    ; Now we have a list of format as input
    
    
    
    ; 
    MaxRowCount = 0
    
    
    ; Then we analyse each column data!
    FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
        ; 
        ; Loop each column
        CCtitle = STRING(FORMAT='("Column",I0)',i)
        CColumn = scope_varfetch(CCtitle)
        CCcount = N_ELEMENTS(CColumn) & IF MaxRowCount LT CCcount THEN MaxRowCount=CCcount
        CHeader = scope_varname(scope_varfetch(CCtitle),LEVEL=-1) ; get the variable of this column in the caller scope
        ; 
        ; if format code is not provided, we need to decide column format, width, precision 
        ; IF i GT N_ELEMENTS(CCFcode) OR i GT N_ELEMENTS(CCFsign) OR i GT N_ELEMENTS(CCFwide) OR i GT N_ELEMENTS(CCFprec) THEN BEGIN
            TMPcode = ""
            TMPsign = "" ; & IF i EQ 1 THEN TMPsign="-"
            TMPwide = ""
            TMPprec = ""
            ; PRINT, SIZE(CColumn,/TNAME)
            IF SIZE(CColumn,/TNAME) EQ "STRING" THEN BEGIN
                TMPcode="A"
                TMPmax=0 & FOREACH TMPstr,CColumn DO IF TMPmax LT STRLEN(TMPstr) THEN TMPmax=STRLEN(TMPstr)
                IF TMPmax LT STRLEN(CHeader) THEN TMPmax=STRLEN(CHeader) ; consider the length of header
                TMPwide=STRING(FORMAT='(I0)',TMPmax+2)
            ENDIF
            IF SIZE(CColumn,/TNAME) EQ "LONG"   OR SIZE(CColumn,/TNAME) EQ "INT" THEN BEGIN
                TMPcode="I"
                TMPmax=MAX(ALOG10(ABS(CColumn)))+1+1 ; +1 to consider the - sign
                IF TMPmax LT STRLEN(CHeader) THEN TMPmax=STRLEN(CHeader) ; consider the length of header
                TMPwide=STRING(FORMAT='(I0)',TMPmax+2)
            ENDIF
            IF SIZE(CColumn,/TNAME) EQ "DOUBLE" OR SIZE(CColumn,/TNAME) EQ  "FLOAT"  THEN BEGIN
                TMPcode="F"
                TMPcoo=WHERE(CColumn NE 0.0,/NULL)
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPmax=MAX(ALOG10(ABS(CColumn[TMPcoo])))+1+1 ; +1 to consider the - sign
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPwoo=TMPmax ELSE TMPwoo=3 ; TMPwoo is the integer digits including sign mark
                IF SIZE(CColumn,/TNAME) EQ "DOUBLE" THEN TMPpoo=7 ELSE TMPpoo=4 ; TMPpoo is the fractional width, or the precision
                TMPcoo=WHERE(ABS(CColumn) LT 1E-4 OR ABS(CColumn) GE 1E4,/NULL)
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPpoo=7
                TMPcoo=WHERE(ABS(CColumn) LT 1E-6 OR ABS(CColumn) GE 1E6,/NULL)
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPpoo=6
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPcode="E"
                IF N_ELEMENTS(TMPcoo) GT 0 THEN TMPwoo=2+4
                TMPwoo=TMPwoo+1+TMPpoo ; +1 means including the dot
                IF TMPwoo LT STRLEN(CHeader) THEN TMPwoo=STRLEN(CHeader) ; consider the length of header
                TMPwide=STRING(FORMAT='(I0)',TMPwoo+2)
                TMPprec=STRING(FORMAT='(".",I0)',TMPpoo)
            ENDIF
        ;     PRINT, "Column "+STRING(FORMAT='(I0)',i)+": "+TMPcode+TMPsign+TMPwide+TMPprec
        ;     IF i GT N_ELEMENTS(CCFcode) THEN CCFcode = [CCFcode,TMPcode] ELSE CCFcode[i-1] = TMPcode
        ;     IF i GT N_ELEMENTS(CCFsign) THEN CCFsign = [CCFsign,TMPsign] ELSE CCFsign[i-1] = TMPsign
        ;     IF i GT N_ELEMENTS(CCFwide) THEN CCFwide = [CCFwide,TMPwide] ELSE CCFwide[i-1] = TMPwide
        ;     IF i GT N_ELEMENTS(CCFprec) THEN CCFprec = [CCFprec,TMPprec] ELSE CCFprec[i-1] = TMPprec
            IF i GT N_ELEMENTS(CCFcode) THEN CCFcode = [CCFcode,TMPcode] ELSE IF CCFcode[i-1] EQ "" THEN CCFcode[i-1] = TMPcode ELSE TMPcode = CCFcode[i-1]
            IF i GT N_ELEMENTS(CCFsign) THEN CCFsign = [CCFsign,TMPsign] ELSE IF CCFsign[i-1] EQ "" THEN CCFsign[i-1] = TMPsign ELSE TMPsign = CCFsign[i-1]
            IF i GT N_ELEMENTS(CCFwide) THEN CCFwide = [CCFwide,TMPwide] ELSE IF CCFwide[i-1] EQ "" THEN CCFwide[i-1] = TMPwide ELSE TMPwide = CCFwide[i-1]
            IF i GT N_ELEMENTS(CCFprec) THEN CCFprec = [CCFprec,TMPprec] ELSE IF CCFprec[i-1] EQ "" THEN CCFprec[i-1] = TMPprec ELSE TMPprec = CCFprec[i-1]
            PRINT, "Column "+STRING(FORMAT='(I0)',i)+": "+TMPcode+TMPsign+TMPwide+TMPprec 
        ; ENDIF ELSE BEGIN
        ;     TMPcode = CCFcode[i-1]
        ;     TMPsign = CCFsign[i-1]
        ;     TMPwide = CCFwide[i-1]
        ;     TMPprec = CCFprec[i-1]
        ;     PRINT, "Column "+STRING(FORMAT='(I0)',i)+": "+TMPcode+TMPsign+TMPwide+TMPprec
        ; ENDELSE
        
        
        ; get ColumnHeader
        IF i-1 LE N_ELEMENTS(Header)-1 THEN CHeader = Header[i-1]
        CCFhead = [CCFhead,CHeader]
        
        
    ENDFOR
    
    
    
    
    
    ; Now print header
    THeader = ""
    FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
        IF i EQ 1 THEN BEGIN
            THeader = "# " + STRING(FORMAT='(A'+CCFsign[i-1]+STRING(FORMAT='(I0)',FIX(CCFwide[i-1])-2)+')',CCFhead[i-1])
        ENDIF ELSE BEGIN
            THeader = THeader + STRING(FORMAT='(A'+CCFsign[i-1]+CCFwide[i-1]+')',CCFhead[i-1])
        ENDELSE
        ; add Separator
        IF N_ELEMENTS(Separator) GT 0 AND i NE N_PARAMS()-1 THEN BEGIN
            THeader = THeader + Separator
        ENDIF
    ENDFOR
    
    IF FilePathOrFileUnit NE '' THEN BEGIN
        ;<modified><20150702><dzliu>; create dir if dir non-exist
        IF NOT FILE_TEST(FILE_DIRNAME(FilePathOrFileUnit),/DIRECTORY) THEN FILE_MKDIR, FILE_DIRNAME(FilePathOrFileUnit)
        ;<modified><20150619><dzliu>; allow to input empty string as the FilePathOrFileUnit, so as to print directly on screen
        OPENW,  FileUnit, FilePathOrFileUnit, /GET_LUN, APPEND=APPEND
        IF NOT KEYWORD_SET(APPEND) THEN PRINTF, FileUnit, THeader ; do not print header when append
        IF NOT KEYWORD_SET(APPEND) THEN PRINTF, FileUnit, "#"
    ENDIF ELSE BEGIN
        PRINT, "#"
        PRINT, THeader
        PRINT, "#"
    ENDELSE
    
;   PRINT, THeader
;   PRINT, "#"
    
    
    
    TContents = MAKE_ARRAY(MaxRowCount,/STRING,VALUE='')
    FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
        ;
        ; Loop each column
        CCtitle = STRING(FORMAT='("Column",I0)',i)
        CColumn = scope_varfetch(CCtitle)
        ; 
        ; Fill up column
        CColumn = STRING(FORMAT='('+CCFcode[i-1]+CCFsign[i-1]+CCFwide[i-1]+CCFprec[i-1]+')',CColumn)
        IF N_ELEMENTS(CColumn) LT MaxRowCount THEN BEGIN
            CColumn = [CColumn,MAKE_ARRAY(MaxRowCount-N_ELEMENTS(CColumn),/STRING,VALUE=STRING(FORMAT='('+'A'+CCFwide[i-1]+')',' '))]
        ENDIF
        ; 
        ; Add Separator
        IF N_ELEMENTS(Separator) GT 0 AND i NE N_PARAMS()-1 THEN BEGIN
            FOR k=0,N_ELEMENTS(CColumn)-1 DO BEGIN
                CColumn[k] = CColumn[k] + Separator
            ENDFOR
        ENDIF
        ;
        ; Save column context
        TContents = TContents + CColumn
    ENDFOR
    
    
    ; Now print table content
    IF FilePathOrFileUnit NE '' THEN BEGIN
        PRINTF, FileUnit, FORMAT='(A)', TContents ; <New><20141217><DzLIU>
    ENDIF ELSE BEGIN
        PRINT, FORMAT='(A)', TContents            ; <added><20150619><dzliu>
        PRINT, " "                                ; <added><20150619><dzliu>
    ENDELSE
    
;    FOR j=0,MaxRowCount-1 DO BEGIN
;        TContext = ""
;        FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
;            ;
;            ; Loop each column
;            CCtitle = STRING(FORMAT='("Column",I0)',i)
;            CColumn = scope_varfetch(CCtitle)
;            ; 
;            ; Write column context
;            IF j LE N_ELEMENTS(CColumn)-1 THEN BEGIN
;                TContext = TContext + STRING(FORMAT='('+CCFcode[i-1]+CCFsign[i-1]+CCFwide[i-1]+CCFprec[i-1]+')',CColumn[j])
;            ENDIF ELSE BEGIN
;                ; if the column is empty -- <TODO> fill white spaces??? <TODO>
;                TContext = TContext + STRING(FORMAT='('+'A'+CCFwide[i-1]+')',' ')
;            ENDELSE
;        ENDFOR
;;       PRINT, TContext
;        PRINTF, FileUnit, TContext
;    ENDFOR ; Replaced by <New><20141217><DzLIU>
    
    
    ; 
    IF FilePathOrFileUnit NE '' THEN BEGIN
        CLOSE, FileUnit
        FREE_LUN, FileUnit
    ENDIF
    
    
    ; 
    RETURN
;    
;    
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    ; <TODO>
;    
;    ; check max row count and column format and column header
;    MaxRowCount = 0
;    ColumnFormat = []
;    ColumnDefStr = [] ; default value in string format
;    ColumnHeader = "#"
;    FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
;        ; get current column varname
;        CurrentColumnName = STRING(FORMAT='("Column",I0)',i)
;        ; get current column variable
;        CurrentColumn = scope_varfetch(CurrentColumnName)
;        ; get column array size
;        IF MaxRowCount LT N_ELEMENTS(CurrentColumn) THEN MaxRowCount=N_ELEMENTS(CurrentColumn)
;        ; decide column format
;        CurrentFormatForHeader = ""
;        IF N_ELEMENTS(Format) EQ 1 AND SIZE(Format,/TNAME) EQ "STRING" THEN BEGIN
;            ColumnFormat = Format[0]
;            CurrentDefStr = ""
;            CurrentWidth = 0
;            CurrentFormatForHeader = ""
;        ENDIF ELSE IF N_ELEMENTS(Format) GT 1 AND SIZE(Format,/TNAME) EQ "STRING" THEN BEGIN
;            IF i-1 LE N_ELEMENTS(Format)-1 THEN CurrentFormat=Format[i-1] ELSE CurrentFormat=Format[0]
;            CurrentDefStr = ""
;            CurrentWidth = 0
;            CurrentFormatForHeader = ""
;            IF NOT STRMATCH(CurrentFormat,'\(*\)') THEN CurrentFormat = '('+CurrentFormat+')'
;            ColumnFormat = [ ColumnFormat, CurrentFormat ]
;            ColumnDefStr = [ ColumnDefStr, CurrentDefStr ]
;        ENDIF ELSE BEGIN
;            ; automatically decide column format
;            IF SIZE(CurrentColumn,/TNAME) EQ "STRING" THEN BEGIN
;                CurrentWidth = MAX(STRLEN(CurrentColumn))
;                CurrentFtype = "A"
;                CurrentFormat = STRING(FORMAT='(A,I0)',CurrentFtype,CurrentWidth)
;                IF NOT STRMATCH(CurrentFormat,'\(*\)') THEN CurrentFormat = '('+CurrentFormat+')'
;                CurrentDefStr = STRING(FORMAT=CurrentFormat,' ')
;                IF i EQ 1 THEN CurrentWidth = CurrentWidth-1
;                CurrentFormatForHeader = STRING(FORMAT='("A",I0)',CurrentWidth)
;            ENDIF ELSE IF SIZE(CurrentColumn,/TNAME) EQ "FLOAT" OR SIZE(CurrentColumn,/TNAME) EQ "DOUBLE" THEN BEGIN
;;                CurrentFmean = MEAN(ABS(CurrentColumn),/DOUBLE,/NAN)
;;                IF CurrentFmean LT 1E-5 THEN BEGIN
;;                    CurrentFtype = "E"
;;                    CurrentWidth = 3+7+4+2
;;                    CurrentPreci = 7
;;                ENDIF ELSE IF CurrentFmean LT 1.0 THEN BEGIN
;;                    CurrentFtype = "F"
;;                    CurrentWidth = 3+7+2
;;                    CurrentPreci = 7
;;                ENDIF ELSE IF CurrentFmean LT 10.0 THEN BEGIN
;;                    CurrentFtype = "F"
;;                    CurrentWidth = 4+6+2
;;                    CurrentPreci = 6
;;                ENDIF ELSE IF CurrentFmean LT 100.0 THEN BEGIN
;;                    CurrentFtype = "F"
;;                    CurrentWidth = 4+6+2
;;                    CurrentPreci = 6
;;                ENDIF ELSE IF CurrentFmean LT 10.0 THEN BEGIN
;;                    CurrentFtype = "F"
;;                    CurrentWidth = 4+6+2
;;                    CurrentPreci = 6
;;                ENDIF ELSE IF CurrentFmean LT 10.0 THEN BEGIN
;;                    CurrentFtype = "F"
;;                    CurrentWidth = 4+6+2
;;                    CurrentPreci = 6
;;                ENDIF ELSE IF 
;                CurrentWidth = 15
;                CurrentPreci = 6
;                CurrentFtype = "G"
;                CurrentFormat = STRING(FORMAT='(A,I0,".",I0)',CurrentFtype,CurrentWidth,CurrentPreci)
;                CurrentDefStr = STRING(FORMAT='('+CurrentFormat+')',0.0)
;                IF i EQ 1 THEN CurrentWidth = CurrentWidth-1
;                CurrentFormatForHeader = STRING(FORMAT='("A",I0)',CurrentWidth)
;            ENDIF ELSE IF SIZE(CurrentColumn,/TNAME) EQ "INT" THEN BEGIN
;                CurrentWidth = 6
;                CurrentFtype = "I"
;                CurrentFormat = STRING(FORMAT='(A,I0)',CurrentFtype,CurrentWidth)
;                CurrentDefStr = STRING(FORMAT='('+CurrentFormat+')',0)
;                IF i EQ 1 THEN CurrentWidth = CurrentWidth-1
;                CurrentFormatForHeader = STRING(FORMAT='("A",I0)',CurrentWidth)
;            ENDIF ELSE IF SIZE(CurrentColumn,/TNAME) EQ "LONG" THEN BEGIN
;                CurrentWidth = 10
;                CurrentFtype = "I"
;                CurrentFormat = STRING(FORMAT='(A,I0)',CurrentFtype,CurrentWidth)
;                CurrentDefStr = STRING(FORMAT='('+CurrentFormat+')',0)
;                IF i EQ 1 THEN CurrentWidth = CurrentWidth-1
;                CurrentFormatForHeader = STRING(FORMAT='("A",I0)',CurrentWidth)
;            ENDIF
;            IF NOT STRMATCH(CurrentFormat,'\(*\)') THEN CurrentFormat = '('+CurrentFormat+')'
;            ColumnFormat = [ ColumnFormat, CurrentFormat ]
;            ColumnDefStr = [ ColumnDefStr, CurrentDefStr ]
;        ENDELSE
;        
;        ; get ColumnHeader
;        IF i-1 LE N_ELEMENTS(Header)-1 THEN BEGIN
;            ; get ColumnHeader from Input
;            CurrentHeader = Header[i-1]
;        ENDIF ELSE BEGIN
;            ; get ColumnHeader from varname of current column in the caller scope
;            CurrentHeader = scope_varname(scope_varfetch(CurrentColumnName),LEVEL=-1) ; get the variable of this column in the caller scope
;            ; <TODO> restore uppercase-lowercase
;            ; <TODO> 
;        ENDELSE
;        
;        ; decide format for ColumnHeader
;        IF CurrentFormatForHeader EQ "" THEN CurrentFormatForHeader = STRING(FORMAT='("(A",I0,")")',FIX(STRLEN(CurrentHeader)*1.5))
;        IF NOT STRMATCH(CurrentFormatForHeader,'\(*\)') THEN CurrentFormatForHeader = '('+CurrentFormatForHeader+')'
;        
;        ColumnHeader = ColumnHeader + STRING(FORMAT=CurrentFormatForHeader,CurrentHeader)
;        
;    ENDFOR
;    
;    ; loop each row and each column and printf
;    OPENW, FileUnit, FilePathOrFileUnit, /GET_LUN, APPEND=APPEND
;    PRINTF, FileUnit, ColumnHeader
;    PRINTF, FileUnit, "#"
;    FOR j=0,MaxRowCount-1 DO BEGIN
;        CurrentLine = ""
;        IF N_ELEMENTS(ColumnFormat) EQ 1 THEN BEGIN
;            CurrentFormat = ColumnFormat
;            IF NOT STRMATCH(CurrentFormat,'\(*\)') THEN CurrentFormat = '('+CurrentFormat+')'
;            IF N_PARAMS() EQ 1+1 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j])
;            IF N_PARAMS() EQ 1+2 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j])
;            IF N_PARAMS() EQ 1+3 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j])
;            IF N_PARAMS() EQ 1+4 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j], Column4[j])
;            IF N_PARAMS() EQ 1+5 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j], Column4[j], Column5[j])
;            IF N_PARAMS() EQ 1+6 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j], Column4[j], Column5[j], Column6[j])
;            IF N_PARAMS() EQ 1+7 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j], Column4[j], Column5[j], Column6[j], Column7[j])
;            IF N_PARAMS() EQ 1+8 THEN CurrentLine = STRING(FORMAT=CurrentFormat, Column1[j], Column2[j], Column3[j], Column4[j], Column5[j], Column6[j], Column7[j], Column8[j])
;        ENDIF ELSE BEGIN
;            FOR i=1,N_PARAMS()-1 DO BEGIN ; N_PARAMS() = 1 + NumCol
;                CurrentColumnName = STRING(FORMAT='("Column",I0)',i)
;                CurrentColumn = scope_varfetch(CurrentColumnName)
;                IF j LE N_ELEMENTS(CurrentColumn)-1 THEN BEGIN
;                    CurrentFormat = ColumnFormat[i-1]
;                    CurrentLine = CurrentLine + STRING(FORMAT=CurrentFormat,CurrentColumn[j])
;                ENDIF ELSE BEGIN
;                    CurrentDefStr = ColumnDefStr[i-1]
;                    CurrentLine = CurrentLine + CurrentDefStr
;                ENDELSE
;            ENDFOR
;        ENDELSE
;        PRINTF, FileUnit, CurrentLine
;    ENDFOR
;    CLOSE, FileUnit
;    FREE_LUN, FileUnit
;    
;    RETURN
;    
END
