; 
; This code is for 
; computing chi^2 distribution 
; and plotting chi^2 distribution
; and outputing bestfit paarameter values and errors
; 
; The input files should have the same format for 
; michi2 code
; The first input is flux.dat file, 
; the second input is a list of lib files, 
; and the third input is the michi2 output chi2 file. 
; 
; The output of this code will be
; One EPS chi2 distribution figure for each parameter in each lib file, 
; and one histogram.txt and one bestfit.txt file for each of them. 
; 
; 
PRO PdChi2_v01, InputDAT, InputLIB, InputFIT, OutputName=OutputName
    
    ; resolve_all
    CD, C=CurrentPath
    IF NOT STRMATCH(!PATH,CurrentPath+"/idlcrab/crabarray:"+CurrentPath+"/idlcrab/crabstring:"+CurrentPath+"/idlcrab/crabtable:") THEN BEGIN
        !PATH = CurrentPath+"/idlcrab/crabarray:"+CurrentPath+"/idlcrab/crabstring:"+CurrentPath+"/idlcrab/crabtable:" + !PATH
    ENDIF
    resolve_all
    
    ; Welcome
    PRINT, "Welcome!"
    
    ; Check Command Line Argument
    InputCommandLineArguments = COMMAND_LINE_ARGS(Count=InputCommandLineArgCount)
    IF InputCommandLineArgCount GE 1 THEN BEGIN
        InputDAT = []
        InputLIB = []
        InputFIT = []
        PRINT, "Command line arguments: "
        FOR i = 0, InputCommandLineArgCount-1 DO BEGIN
            ;PRINT, InputCommandLineArguments[i]
            IF i EQ 0 THEN BEGIN
                InputDAT = InputCommandLineArguments[i]
                PRINT, "InputDAT = "+InputCommandLineArguments[i]
            ENDIF ELSE IF i NE InputCommandLineArgCount-1 THEN BEGIN
                InputLIB = [ InputLIB, InputCommandLineArguments[i] ]
                PRINT, "InputLIB += "+InputCommandLineArguments[i]
            ENDIF ELSE IF i EQ InputCommandLineArgCount-1 THEN BEGIN
                InputFIT = InputCommandLineArguments[i]
                PRINT, "InputFIT = "+InputCommandLineArguments[i]
            ENDIF
        ENDFOR
    ENDIF
    
    ; Check InputDAT
    IF N_ELEMENTS(InputDAT) EQ 0 THEN MESSAGE, "Error! InputDAT is not given!", /CONTINUE
    
    ; Check InputLIB
    IF N_ELEMENTS(InputLIB) EQ 0 THEN MESSAGE, "Error! InputLIB is not given!", /CONTINUE
    
    ; Check InputFIT
    IF N_ELEMENTS(InputFIT) EQ 0 THEN MESSAGE, "Error! InputFIT is not given!", /CONTINUE
    
    ; Check All
    IF N_ELEMENTS(InputDAT) EQ 0 OR $
       N_ELEMENTS(InputLIB) EQ 0 OR $
       N_ELEMENTS(InputFIT) EQ 0 THEN BEGIN
        PRINT, ""
        PRINT, "Usage: "
        PRINT, "    PdChi2_v01, InputDAT, InputLIB, InputFIT"
        PRINT, ""
        PRINT, "Example: "
        PRINT, "    PdChi2_v01, 'flux_obsframe.dat', ['FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED','MullaneyAGN.lib.SED','DL07.HiExCom.lib.SED','DL07.LoExCom.lib.SED'], 'fit_quadruple.out'"
        PRINT, ""
        PRINT, "Alternative:"
        PRINT, "    We can run from bash command line"
        PRINT, "    idl -e PdChi2_v01 -args flux_obsframe.dat FSPS.Padova.BaSeL.Z0.0190.EBV.lib.SED MullaneyAGN.lib.SED DL07.HiExCom.lib.SED DL07.LoExCom.lib.SED fit_quadruple.out"
        PRINT, ""
        RETURN
    ENDIF
    
    ; Reform InputDAT
    IF N_ELEMENTS(InputDAT) EQ 1 THEN PdChi2_DAT = InputDAT ELSE PdChi2_DAT = InputDAT[0]
    
    ; Reform InputLIB
    IF N_ELEMENTS(InputLIB) EQ 1 THEN PdChi2_LIB = [InputLIB] ELSE PdChi2_LIB = InputLIB
    
    ; Reform InputFIT
    IF N_ELEMENTS(InputFIT) EQ 1 THEN PdChi2_FIT = InputFIT ELSE PdChi2_FIT = InputFIT[0]
    
    ; Check file existennce
    Check_OK = 1
    IF FILE_TEST(PdChi2_DAT) EQ 0 THEN BEGIN
        MESSAGE, "Error! "+PdChi2_DAT+" is not found!", /CONTINUE
        Check_OK = 0
    ENDIF
    FOREACH Temp_LIB, PdChi2_LIB DO BEGIN 
        IF FILE_TEST(Temp_LIB) EQ 0 THEN BEGIN
            MESSAGE, "Error! "+Temp_LIB+" is not found!", /CONTINUE
            Check_OK = 0
        ENDIF
    ENDFOREACH
    IF FILE_TEST(PdChi2_FIT) EQ 0 THEN BEGIN
        MESSAGE, "Error! "+PdChi2_FIT+" is not found!", /CONTINUE
        Check_OK = 0
    ENDIF
    IF Check_OK EQ 0 THEN MESSAGE, "Error! Please check the above error information! Exit!"
    
    ; Read InputLIB, check how many parameters
    PdChi2_NPars = []
    FOREACH Temp_LIB, PdChi2_LIB DO BEGIN 
        PdChi2_NPars = [ PdChi2_NPars, FIX(CrabTableReadInfo(Temp_LIB,"# NPAR")) ]
    ENDFOREACH
    PRINT, "Number of LIB files = "+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_NPars))
    PRINT, "Number of LIB params = "+STRING(FORMAT='(I0)',TOTAL(PdChi2_NPars))
    
    ; Read FIT file (the chi2 output file from michi2)
    ; The FIT file starts with two columns i0 and chi2, 
    ; then for each LIB, there has two columns iLib and aLib
    ;;Verbose = 1 ;;<20160716> checked good no problem
    PdChi2_Col_i0 = LONG(CrabTableReadColumn(PdChi2_FIT,1,Verbose=Verbose,/FastReadDataBlock))
    PdChi2_Col_chi2 = DOUBLE(CrabTableReadColumn(PdChi2_FIT,2,Verbose=Verbose,/FastReadDataBlock))
    ;;CrabTablePrintC, "example/temp_i0.txt", PdChi2_Col_i0 ;;<20160716> checked good no problem
    ;;CrabTablePrintC, "example/temp_chi2.txt", PdChi2_Col_chi2 ;;<20160716> checked good no problem
    
    
    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ; CrabTableReadColumn TOO SLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ; Abort!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    
    ; Loop Parameters and Compute Chi2 Distribution
    FOR i=0, N_ELEMENTS(PdChi2_NPars)-1 DO BEGIN
        PdChi2_Col_iLib = LONG(CrabTableReadColumn(PdChi2_FIT,2+2*i+1))
        PdChi2_Col_aLib = DOUBLE(CrabTableReadColumn(PdChi2_FIT,2+2*i+1+1))
        FOR j=0, PdChi2_NPars[i]-1 DO BEGIN
            PRINT, "Analyzing LIB "+STRING(FORMAT='(I0)',i+1)+" PAR "+STRING(FORMAT='(I0)',j+1)
            IF i EQ 0 THEN BEGIN
                PdChi2_Col_tPar = CrabTableReadInfo(PdChi2_LIB[i],"# TPAR"+STRING(FORMAT='(I0)',j+1))
                PdChi2_Col_aPar = DOUBLE(CrabTableReadColumn(PdChi2_FIT,2+2*N_ELEMENTS(PdChi2_NPars)+j+1,Verbose=Verbose,/FastReadDataBlock))
            ENDIF ELSE BEGIN
                PdChi2_Col_tPar = CrabTableReadInfo(PdChi2_LIB[i],"# TPAR"+STRING(FORMAT='(I0)',j+1))
                PdChi2_Col_aPar = DOUBLE(CrabTableReadColumn(PdChi2_FIT,2+2*N_ELEMENTS(PdChi2_NPars)+TOTAL(PdChi2_NPars[0:i-1])+j+1,Verbose=Verbose,/FastReadDataBlock))
            ENDELSE
            PRINT, "    Read "+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Col_aPar))+" lines"
            PRINT, "    Read "+PdChi2_Col_tPar+" = [" + CrabStringPrintArray(PdChi2_Col_aPar[0:4],/NoBracket) + $
                   ",...," + CrabStringPrintArray(PdChi2_Col_aPar[N_ELEMENTS(PdChi2_Col_aPar)-5:N_ELEMENTS(PdChi2_Col_aPar)-1],/NoBracket) + "]"+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Col_aPar))+")"
            PRINT, "    Read "+      "chi2"   +" = [" + CrabStringPrintArray(PdChi2_Col_chi2[0:4],/NoBracket) + $
                   ",...," + CrabStringPrintArray(PdChi2_Col_chi2[N_ELEMENTS(PdChi2_Col_chi2)-5:N_ELEMENTS(PdChi2_Col_chi2)-1],/NoBracket) + "]"+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Col_chi2))+")"
            ; 
            ; Convert PAR value to log or linear
            PdChi2_Par_loga = 0
            IF STRMATCH(PdChi2_Col_tPar,"n_{H_2}") THEN BEGIN
                PdChi2_Col_aPar = ALOG10(PdChi2_Col_aPar)
                PdChi2_Col_aPar = LONG(PdChi2_Col_aPar*100) / 100.0D ; then make it only 2-digit precision
                PdChi2_Par_loga = 1
                PdChi2_Col_tPar = "log "+PdChi2_Col_tPar
            ENDIF
            ; 
            ; Prepare PAR grid
            ; (1) method one
            ;PdChi2_Par_Step = (MAX(PdChi2_Col_aPar)-MIN(PdChi2_Col_aPar))/10.0D
            ;PdChi2_Par_Hist = Histogram(PdChi2_Col_aPar, BINSIZE=PdChi2_Par_Step, LOCATIONS=PdChi2_Par_Grid)
            ;PdChi2_Par_Step = PdChi2_Par_Hist * 0.0D + PdChi2_Par_Step ; make it the same dimension as histogram grid
            ; (2) method two
            PdChi2_Par_Grid = PdChi2_Col_aPar[UNIQ(PdChi2_Col_aPar,SORT(PdChi2_Col_aPar))]
            IF N_ELEMENTS(PdChi2_Par_Grid) EQ 0 THEN MESSAGE, "Error! Could not make grid by UNIQ()!"
            IF N_ELEMENTS(PdChi2_Par_Grid) EQ 1 THEN BEGIN
                PdChi2_Par_Step = [0.0] ; compute the grid step, if only one grid point, then set step to 0
                PdChi2_Par_Grid = [PdChi2_Par_Grid]
                PdChi2_Par_Hist = [PdChi2_Par_Grid] * 0.0 ;<TODO> no need to compute histogram of PAR value
            ENDIF ELSE BEGIN
                PdChi2_Par_Step = PdChi2_Par_Grid[1:N_ELEMENTS(PdChi2_Par_Grid)-1] - PdChi2_Par_Grid[0:N_ELEMENTS(PdChi2_Par_Grid)-2] ; compute the grid step
                PdChi2_Par_Grid = PdChi2_Par_Grid[0:N_ELEMENTS(PdChi2_Par_Grid)-2]
                PdChi2_Par_Hist = PdChi2_Par_Grid * 0.0 ;<TODO> no need to compute histogram of PAR value
            ENDELSE
            ; 
            ; Print PAR grid
            PRINT, "    Grid "+PdChi2_Col_tPar+" = "+CrabStringPrintArray(PdChi2_Par_Grid)+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Par_Grid))+")"
            PRINT, "    Hist "+PdChi2_Col_tPar+" = "+CrabStringPrintArray(PdChi2_Par_Hist)+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Par_Hist))+")"
            ; 
            ; Check PAR grid
            IF N_ELEMENTS(PdChi2_Par_Grid) GT 30 THEN MESSAGE, "Error! Too much PAR grid points (>30)!"
            ; 
            ; Compute chi2 in each grid cell
            PdChi2_Chi_Grid = PdChi2_Par_Grid
            PdChi2_Chi_Hist = PdChi2_Par_Hist * 0.0D
            FOR k=0,N_ELEMENTS(PdChi2_Par_Grid)-1 DO BEGIN
                Temp_WER = WHERE(PdChi2_Col_aPar GE PdChi2_Par_Grid[k] AND PdChi2_Col_aPar LT PdChi2_Par_Grid[k]+PdChi2_Par_Step[k], /NULL) ; select chi2 data points in each grid cell
                IF k EQ N_ELEMENTS(PdChi2_Par_Grid)-1 THEN Temp_WER = [ Temp_WER, WHERE(PdChi2_Col_aPar EQ PdChi2_Par_Grid[k]+PdChi2_Par_Step[k], /NULL) ] ; the last grid cell, also consider the max value
                IF N_ELEMENTS(Temp_WER) GT 0 THEN BEGIN
                    PdChi2_Chi_Hist[k] = MIN(PdChi2_Col_chi2[Temp_WER])
                ENDIF
            ENDFOR
            PRINT, "    Grid "+"chi2"+" = "+CrabStringPrintArray(PdChi2_Chi_Grid)+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Chi_Grid))+")"
            PRINT, "    Hist "+"chi2"+" = "+CrabStringPrintArray(PdChi2_Chi_Hist)+" ("+STRING(FORMAT='(I0)',N_ELEMENTS(PdChi2_Chi_Hist))+")"
            FOR k=0,N_ELEMENTS(PdChi2_Par_Grid)-1 DO BEGIN
                PRINT, "    "+PdChi2_Col_tPar+" "+STRTRIM(STRING(PdChi2_Par_Grid[k]),2)+" "+STRTRIM(STRING(PdChi2_Par_Grid[k]+PdChi2_Par_Step[k]),2)+$
                       " chi2 "+STRTRIM(STRING(PdChi2_Chi_Grid[k]),2)+" histogram "+STRTRIM(STRING(PdChi2_Chi_Hist[k]),2)
            ENDFOR
            ; 
            ; Plot the chi2 PDF of PAR
            SET_PLOT, 'PS' & XSizeInCM=10.5 & YSizeInCM=7 ; & SET_FONT="NGC"
            IF N_ELEMENTS(OutputName) EQ 0 THEN OutputName="Output_PDCHI2"
            DEVICE, FILENAME=OutputName+"_"+"LIB"+"_"+STRING(FORMAT='(I0)',i+1)+"_PAR_"+STRING(FORMAT='(I0)',j+1)+".histogram.eps", $
                    /COLOR, BITS_PER_PIXEL=8, DECOMPOSED=1, /ENCAPSULATED, XSIZE=XSizeInCM, YSIZE=YSizeInCM
            IF N_ELEMENTS(SET_FONT) GT 0 THEN DEVICE, SET_FONT=SET_FONT, /TT_FONT
            IF N_ELEMENTS(PdChi2_Chi_Grid) EQ 1 THEN BEGIN ; if only one grid point
                PlotXRange = [ PdChi2_Chi_Grid[0]-1.0, PdChi2_Chi_Grid[0]+1.0 ]
                PlotYRange = [ PdChi2_Chi_Hist[0]-1.0, PdChi2_Chi_Hist[0]+1.0 ]
                PLOT, [PdChi2_Chi_Grid[0]-0.5,PdChi2_Chi_Grid[0]+0.5], [PdChi2_Chi_Hist[0],PdChi2_Chi_Hist[0]], /NODATA, XRANGE=PlotXRange, YRANGE=PlotYRange, $
                      THICK=3, XTHICK=3, YTHICK=3, XCHARSIZE=0.8, YCHARSIZE=0.8, CHARTHICK=2, $
                      XSTYLE=1, YSTYLE=1, FONT=1, POSITION=[0.13,0.18,0.95,0.95]
                OPLOT, [PdChi2_Chi_Grid[0]-0.5,PdChi2_Chi_Grid[0]+0.5], [PdChi2_Chi_Hist[0],PdChi2_Chi_Hist[0]], PSYM=10, THICK=4, Color=cgColor("blue")
            ENDIF ELSE BEGIN
                PlotXRange = [ 0.95*MIN(PdChi2_Chi_Grid), MAX(PdChi2_Chi_Grid)*1.05 ]
                PlotYRange = [ 0.75*MIN(PdChi2_Chi_Hist), MAX(PdChi2_Chi_Hist)*1.25 ]
                PLOT, PdChi2_Chi_Grid, PdChi2_Chi_Hist, /NODATA, XRANGE=PlotXRange, YRANGE=PlotYRange, $
                      THICK=3, XTHICK=3, YTHICK=3, XCHARSIZE=0.8, YCHARSIZE=0.8, CHARTHICK=2, $
                      XSTYLE=1, YSTYLE=1, FONT=1, POSITION=[0.13,0.18,0.95,0.95]
                OPLOT, PdChi2_Chi_Grid, PdChi2_Chi_Hist, PSYM=10, THICK=4, Color=cgColor("blue")
            ENDELSE
            XYOUTS, 0.55, 0.03, textoidl(PdChi2_Col_tPar), /NORMAL, ALIGNMENT=0.5, CHARTHICK=3, CHARSIZE=1.3
            XYOUTS, 0.05, 0.55, textoidl("\chi^{2}"), /NORMAL, ALIGNMENT=0.5, ORIENT=90, CHARTHICK=3, CHARSIZE=1.3
            ; 
            ; Overplot the minimum chi2 value
            PdChi2_Chi_Min = MIN(PdChi2_Col_chi2)
            PdChi2_Chi_Lim = PdChi2_Chi_Min + 4.5 ;<TODO> how to constructing confidence contours/regions?
                                                  ;       -- http://www.astro.sunysb.edu/metchev/PHY515/astrostatistics.html
                                                  ;       -- http://www.aip.de/groups/soe/local/numres/bookcpdf/c15-6.pdf (PDF page 9 table)
            PLOTS, PlotXRange, [PdChi2_Chi_Lim,PdChi2_Chi_Lim], LINESTYLE=2
            ; 
            ; Compute confidence region for the fitted N_ELEMENTS(PdChi2_NPars) parameters
            Temp_WER = WHERE(PdChi2_Chi_Hist GE PdChi2_Chi_Min AND PdChi2_Chi_Hist LE PdChi2_Chi_Lim, /NULL)
            ;<TODO> IF N_ELEMENTS(Temp_WER) EQ 1 THEN Temp_WER ;<TODO> what if only one bin selected? too narrow PAR grid
            PdChi2_Chi_Confidence_Region = PdChi2_Chi_Grid[Temp_WER]
            XYOUTS, PlotXRange[1]-0.01*(PlotXRange[1]-PlotXRange[0]), PdChi2_Chi_Lim+3.0*0.06*(PlotYRange[1]-PlotYRange[0]), $
                    '1'+'!Z(03C3)'+" confidence region", $ ; Unicode font chart http://www.unicode.org/charts/
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=1
            XYOUTS, PlotXRange[1]-0.01*(PlotXRange[1]-PlotXRange[0]), PdChi2_Chi_Lim+2.0*0.06*(PlotYRange[1]-PlotYRange[0]), $
                    "Mean = "+STRTRIM(STRING(MEAN(PdChi2_Chi_Confidence_Region)),2), $
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=1
            XYOUTS, PlotXRange[1]-0.01*(PlotXRange[1]-PlotXRange[0]), PdChi2_Chi_Lim+1.0*0.06*(PlotYRange[1]-PlotYRange[0]), $
                    "Error = "+STRTRIM(STRING((MAX(PdChi2_Chi_Confidence_Region)-MIN(PdChi2_Chi_Confidence_Region))/2.0),2), $
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=1
            ; 
            ; Write to histogram file
            OPENW, OutputLUN, OutputName+"_"+"LIB"+"_"+STRING(FORMAT='(I0)',i+1)+"_PAR_"+STRING(FORMAT='(I0)',j+1)+".histogram.txt", /GET_LUN
            PRINTF, OutputLUN, "# Current Time "+CrabStringCurrentTime()
            PRINTF, OutputLUN, "# "
            PRINTF, OutputLUN, "# "+PdChi2_Col_tPar+" "+"min_chi^2"
            PRINTF, OutputLUN, "# "
            FOR k=0,N_ELEMENTS(PdChi2_Chi_Hist)-1 DO BEGIN
                PRINTF, OutputLUN, FORMAT="(G20.10,G20.10)", PdChi2_Chi_Grid[k], PdChi2_Chi_Hist[k]
            ENDFOR
            CLOSE, OutputLUN
            FREE_LUN, OutputLUN
            ; 
            ; Write to bestfit file
            OPENW, OutputLUN, OutputName+"_"+"LIB"+"_"+STRING(FORMAT='(I0)',i+1)+"_PAR_"+STRING(FORMAT='(I0)',j+1)+".bestfit.txt", /GET_LUN
            PRINTF, OutputLUN, "# Current Time "+CrabStringCurrentTime()
            PRINTF, OutputLUN, "# "
            PRINTF, OutputLUN, "ILIB = "+STRING(FORMAT='(I0)',i+1)
            PRINTF, OutputLUN, 'TLIB = "'+PdChi2_LIB[i]+'"'
            PRINTF, OutputLUN, "IPAR = "+STRING(FORMAT='(I0)',j+1)
            PRINTF, OutputLUN, 'TPAR = "'+PdChi2_Col_tPar+'"'
            PRINTF, OutputLUN, "Min  = "+STRTRIM(STRING(MIN(PdChi2_Chi_Confidence_Region)),2)
            PRINTF, OutputLUN, "Max  = "+STRTRIM(STRING(MAX(PdChi2_Chi_Confidence_Region)),2)
            PRINTF, OutputLUN, "Mean = "+STRTRIM(STRING(MEAN(PdChi2_Chi_Confidence_Region)),2)
            PRINTF, OutputLUN, "Median = "+STRTRIM(STRING(MEDIAN(PdChi2_Chi_Confidence_Region)),2)
            PRINTF, OutputLUN, "Error = "+STRTRIM(STRING((MAX(PdChi2_Chi_Confidence_Region)-MIN(PdChi2_Chi_Confidence_Region))/2.0),2)
            CLOSE, OutputLUN
            FREE_LUN, OutputLUN
            ; 
            ; Close device
            DEVICE, /CLOSE
            SET_PLOT, "X"
            ;PlotDevice = Plot(PdChi2_Par_Grid, PdChi2_Chi_Hist, YRANGE=CrabMinMax(PdChi2_Chi_Hist), /STAIRSTEP, XTITLE=PdChi2_Col_tPar, YTITLE="chi2")
            
            ;<DEBUG>
            ;BREAK
        ENDFOR
        
        ;<DEBUG>
        ;BREAK
    ENDFOR
    ;; Prepare the format code for reading FIT file (the chi2 output file from michi2)
    ;PdChi2_FormatLibs = "L,D" ; the FIT file starts with two columns i0 and chi2, then for each LIB, there has two columns iN and aN
    ;PdChi2_FormatPars = "" ; then loop again for each LIB, there has N columns for N parameters
    ;FOR i=0, N_ELEMENTS(PdChi2_NPars)-1 DO BEGIN
    ;    PdChi2_FormatLibs += ",L,D"
    ;    FOR j=0, PdChi2_NPars[i]-1 DO BEGIN
    ;        PdChi2_FormatPars += ",D"
    ;    ENDFOR
    ;ENDFOR
    ;PdChi2_FormatFits = PdChi2_FormatLibs+PdChi2_FormatPars
    ;PRINT, "Format of FIT file = "+'('+PdChi2_FormatFits+')'
    
    ; Read FIT file
END