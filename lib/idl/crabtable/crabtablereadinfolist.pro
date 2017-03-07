FUNCTION CrabTableReadInfoList, FilePath, KeyNames, DOUBLE=DOUBLE, INT=INT, $
                                                    Equaller=Equaller, $
                                                    Splitter=Splitter, $
                                                    AllValues=AllValues, $
                                                    RemoveSuffix=RemoveSuffix, $
                                                    KeyNamePattern=KeyNamePattern, $
                                                    Verbose=Verbose, GiveWarning=GiveWarning
; 2012-12-15
; Read an ini file (information file), find the key values according to input key names. 
; e.g. "Source  = NGC6946"
;      "RA      = 123123123"
;      "Dec     = 123123123"
;      PRINT, CrabTableReadInfoList(ini,['Source','RA','Dec'])
;
    
    KeyValues = []
    
    FOREACH KeyName, KeyNames DO BEGIN
        KeyValue = CrabTableReadInfo(FilePath, KeyName, DOUBLE=DOUBLE, INT=INT, $
                                                        Equaller=Equaller, $
                                                        Splitter=Splitter, $
                                                        AllValues=AllValues, $
                                                        RemoveSuffix=RemoveSuffix, $
                                                        KeyNamePattern=KeyNamePattern, $
                                                        Verbose=Verbose, GiveWarning=GiveWarning)
         KeyValues = [KeyValues,KeyValue]
    ENDFOREACH
    
    RETURN, KeyValues
END