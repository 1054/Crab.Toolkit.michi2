; ------------------------------------------------------------------------------------------------------
; CrabStringSplit --- 
;                 if input "[1,2,3,4,5.3]" then output ["1","2","3","4","5.3"]
;                 example, print, CrabStringSplit("STRPOS( Expression, Search_String [, Pos] [, /REVERSE_OFFSET] [, /REVERSE_SEARCH] )", " ")
; ------------------------------------------------------------------------------------------------------
FUNCTION CrabStringSplit, Str, SplitChars, Splitter=Splitter, RemoveBrackets=RemoveBrackets, Compress=Compress, $
                          ESCAPE=ESCAPE, REGEX=REGEX, FOLD_CASE=FOLD_CASE, PRESERVE_NULL=PRESERVE_NULL, Verbose=Verbose
    IF N_PARAMS() LT 1 THEN BEGIN
        MESSAGE,'CrabStringSplit requires at least 1 parameter: Str.'
        RETURN,[]
    ENDIF
    SizeOfStr  = SIZE(Str)
    IF SizeOfStr[N_ELEMENTS(SizeOfStr)-2] NE 7 THEN BEGIN
        MESSAGE,'Parameter Str must be of string type.'
        RETURN,[]
    ENDIF
    ;;; Default Setting
    ;<20160716>IF N_ELEMENTS(Splitter) LE 0       THEN Splitter=','
    ;<20160716>IF N_ELEMENTS(Compress) LE 0       THEN Compress=1
    IF N_ELEMENTS(SplitChars) EQ 0     THEN Splitter=[',',' '] ELSE Splitter = SplitChars
    IF N_ELEMENTS(RemoveBrackets) LE 0 THEN RemoveBrackets=1
    IF N_ELEMENTS(Splitter) EQ 1       THEN Splitter=[Splitter]
    ;<20160716>;;; Remove '='
    ;<20160716>StrCopied = STRCOMPRESS(STRMID(Str,strpos(Str,'=')+1),/REMOVE_ALL)
    StrCopied = Str
    ;;; Remove Brackets
    IF KEYWORD_SET(RemoveBrackets) THEN BEGIN
        StrCopied = CrabStringReplace(StrCopied,'=',' ')
        StrCopied = CrabStringReplace(StrCopied,'{',' ')
        StrCopied = CrabStringReplace(StrCopied,'}',' ')
        StrCopied = CrabStringReplace(StrCopied,'[',' ')
        StrCopied = CrabStringReplace(StrCopied,']',' ')
        StrCopied = CrabStringReplace(StrCopied,'(',' ')
        StrCopied = CrabStringReplace(StrCopied,')',' ')
    ENDIF
;    ;;; Remove Brackets
;    IF RemoveBrackets GE 1 THEN BEGIN
;        StrNew = StrCopied
;        ;;; Remove 3 Types of Brackets
;        BracketPairs = [ ['{','}'], ['[',']'], ['(',')'] ]
;        FOR j=0,2 DO BEGIN
;            Pos1 = STRPOS(StrCopied,BracketPairs[0,j])
;            Pos2 = STRPOS(StrCopied,BracketPairs[1,j])
;            IF Pos1 GE 0 AND Pos2 GE 0 AND Pos2 GT Pos1 THEN $
;                      StrNew = STRMID(StrCopied,Pos1+1,Pos2-Pos1-1)
;            IF STRPOS(StrNew,Splitter) GE 0 THEN StrCopied = StrNew
;        ENDFOR
;    ENDIF
    ;;; Compress
    ;<20160716>PRINT, StrCopied
    ;<20160716>IF KEYWORD_SET(Compress) THEN BEGIN
    ;<20160716>  ;IF Compress GE 1 AND STRPOS(Splitter,' ') LE -1 THEN BEGIN
    ;<20160716>    StrCopied = STRCOMPRESS(StrCopied)
    ;<20160716>  ;ENDIF
    ;<20160716>ENDIF
    ;<20160716>PRINT, StrCopied
    ;;; Extract
    ;<20160716>StrSplit   = STRSPLIT(StrCopied,Splitter,/EXTRACT,/REGEX)
    ;<20160716>RETURN,      StrSplit
    StrSplitted = []
    NumStrInput = N_ELEMENTS(StrCopied)
    IF NumStrInput EQ 1 THEN BEGIN
        StrSplitted = [StrCopied]
    ENDIF ELSE BEGIN
        StrSplitted = [StrCopied[0]]
        NumStrInput = 1
        PRINT, ""
        PRINT, "Warning! CrabStringSplit currently can only deal with one-element input Str!"
        PRINT, ""
        ;<20160716><TODO> ; Currently we can only deal with one-element input str
        ;<20160716><TODO> StrSplitted = MAKE_ARRAY
    ENDELSE
    ;;; Extract (speed-up)
    ;<20160716><TODO> FOR i=0,NumStrInput-1 DO BEGIN
        IF N_ELEMENTS(Splitter) EQ 1 THEN BEGIN
            StrSplitted = STRSPLIT(StrSplitted,Splitter,/EXTRACT, ESCAPE=ESCAPE, REGEX=REGEX, FOLD_CASE=FOLD_CASE, PRESERVE_NULL=PRESERVE_NULL)
        ENDIF ELSE BEGIN
            FOREACH TmpSplitter, Splitter DO BEGIN
                IF KEYWORD_SET(Compress) THEN StrCopied = STRCOMPRESS(StrCopied)
                iSplit = 0
                WHILE iSplit LT N_ELEMENTS(StrSplitted) DO BEGIN
                    TmpSplitted = STRSPLIT(StrSplitted[iSplit],TmpSplitter,/EXTRACT, ESCAPE=ESCAPE, REGEX=REGEX, FOLD_CASE=FOLD_CASE, PRESERVE_NULL=PRESERVE_NULL)
                    IF KEYWORD_SET(Verbose) THEN PRINT, 'Splitting "'+StrSplitted[iSplit]+'" by "'+TmpSplitter+'"'
                    IF KEYWORD_SET(Verbose) THEN PRINT, 'Splitted '+STRING(FORMAT='(I0)',N_ELEMENTS(TmpSplitted))+' "'+TmpSplitted+'"'
                    IF iSplit GT 0 THEN StrSplitted_1 = StrSplitted[0:iSplit-1] ELSE StrSplitted_1 = []
                    IF iSplit LT N_ELEMENTS(StrSplitted)-1 THEN StrSplitted_2 = StrSplitted[iSplit+1:N_ELEMENTS(StrSplitted)-1] ELSE StrSplitted_2 = []
                    StrSplitted = [StrSplitted_1, TmpSplitted, StrSplitted_2]
                    IF KEYWORD_SET(Verbose) THEN PRINT, 'Stored '+STRING(FORMAT='(I0)',N_ELEMENTS(StrSplitted))+' "'+StrSplitted+'"'
                    IF N_ELEMENTS(TmpSplitted) GT 1 THEN BEGIN
                        iSplit = iSplit
                    ENDIF ELSE BEGIN
                        iSplit = iSplit + 1
                    ENDELSE
                ENDWHILE
            ENDFOREACH
        ENDELSE
    ;<20160716><TODO> ENDFOR
    
    ;;; Return
    RETURN, StrSplitted
END