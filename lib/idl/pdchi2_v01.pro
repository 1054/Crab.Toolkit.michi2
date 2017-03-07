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
PRO PdChi2_v01, InputDAT, InputLIB, InputFIT, OutputName=OutputName, SET_FONT=SET_FONT, SET_XRANGE=SET_XRANGE, SET_YRANGE=SET_YRANGE, $
                                              Source=Source, Redshift=Redshift, Distance=Distance
    
    ; <TODO> 
    ;;Source = "ID12646"
    SET_FONT = 'NGC'
    
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
    
    ; Check input_redshift.sm
    IF N_ELEMENTS(Source) EQ 0 THEN BEGIN
        IF NOT FILE_TEST("input_redshift.sm") THEN MESSAGE, 'Error! "input_redshift.sm" was not found! Please prepare this file or give Redshift from command line!'
        Source = (CrabTableReadInfo("input_redshift.sm","set source"))
    ENDIF
    IF N_ELEMENTS(Redshift) EQ 0 THEN BEGIN
        IF NOT FILE_TEST("input_redshift.sm") THEN MESSAGE, 'Error! "input_redshift.sm" was not found! Please prepare this file or give Redshift from command line!'
        Redshift = Double(CrabTableReadInfo("input_redshift.sm","set z"))
    ENDIF
    IF N_ELEMENTS(Distance) EQ 0 THEN BEGIN
        IF NOT FILE_TEST("input_redshift.sm") THEN MESSAGE, 'Error! "input_redshift.sm" was not found! Please prepare this file or give Redshift from command line!'
        Distance = Double(CrabTableReadInfo("input_redshift.sm","set dL"))
    ENDIF
    
    ; Check OutputName
    IF N_ELEMENTS(OutputName) EQ 0 THEN OutputName="Output_PDCHI2"
    
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
    
    
    
    ; Prepare an array for best fitting i0 chi2 i1 a1 i2 a2 etc
    PdChi2_SED_Struct = { i0:0L, chi2:0D, $
                          w:PTR_NEW(/ALLOCATE), f:PTR_NEW(/ALLOCATE), $ ; total SED
                          iLib:REPLICATE(0L,N_ELEMENTS(PdChi2_NPars)), aLib:REPLICATE(0D,N_ELEMENTS(PdChi2_NPars)), $
                          tPar:PTRARR(/ALLOCATE,N_ELEMENTS(PdChi2_NPars)), aPar:PTRARR(/ALLOCATE,N_ELEMENTS(PdChi2_NPars)), $
                          wLib:PTRARR(/ALLOCATE,N_ELEMENTS(PdChi2_NPars)), fLib:PTRARR(/ALLOCATE,N_ELEMENTS(PdChi2_NPars)), $ ; We should not use REPLICATE to create PTR array! Use PTRARR!
                          Mstar:0D, Mdust:0D, Mdust_Cold:0D, Mdust_Warm:0D, LIR:0D, SFR:0D, U:0D, Tdust:0D, qIR:0D, LAGN:0D, fAGN:0D }
    PdChi2_SED_sorted_id = SORT(PdChi2_Col_chi2)
    PdChi2_SED_selected_id = PdChi2_SED_sorted_id[0:6] ;<TODO> select and plot 7 SEDs with best chi2
    PdChi2_SED_colors = [ '1c1ae4'xL,$
                          'b87e37'xL,$
                          '4aaf4d'xL,$
                          'a34e98'xL,$
                          '007fff'xL,$
                          '33ffff'xL,$
                          '2856a6'xL,$
                          'bf81f7'xL ] ; BGR
    PdChi2_SED_List = REPLICATE(PdChi2_SED_Struct,N_ELEMENTS(PdChi2_SED_selected_id))
    
    
    
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
            SET_PLOT, 'PS' & XSizeInCM=10.5 & YSizeInCM=7 &
            DEVICE, FILENAME=OutputName+"_"+"LIB"+"_"+STRING(FORMAT='(I0)',i+1)+"_PAR_"+STRING(FORMAT='(I0)',j+1)+".histogram.eps", $
                    /COLOR, BITS_PER_PIXEL=8, DECOMPOSED=1, /ENCAPSULATED, XSIZE=XSizeInCM, YSIZE=YSizeInCM
            IF N_ELEMENTS(SET_FONT) GT 0 THEN DEVICE, SET_FONT=SET_FONT, /TT_FONT
            IF N_ELEMENTS(SET_FONT) GT 0 THEN PlotFont=1
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
            XYOUTS, 0.05, 0.55, textoidl("\chi^{2}"), /NORMAL, ALIGNMENT=0.5, ORIENT=90, CHARTHICK=3, CHARSIZE=1.3 ; \chi ; '!Z(03C7)'+'!E2!N'
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
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=PlotFont
            XYOUTS, PlotXRange[1]-0.01*(PlotXRange[1]-PlotXRange[0]), PdChi2_Chi_Lim+2.0*0.06*(PlotYRange[1]-PlotYRange[0]), $
                    "Mean = "+STRTRIM(STRING(MEAN(PdChi2_Chi_Confidence_Region)),2), $
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=PlotFont
            XYOUTS, PlotXRange[1]-0.01*(PlotXRange[1]-PlotXRange[0]), PdChi2_Chi_Lim+1.0*0.06*(PlotYRange[1]-PlotYRange[0]), $
                    "Error = "+STRTRIM(STRING((MAX(PdChi2_Chi_Confidence_Region)-MIN(PdChi2_Chi_Confidence_Region))/2.0),2), $
                    ALIGNMENT=1.0, CHARTHICK=2, CHARSIZE=0.8, Color=cgColor("blue"), FONT=PlotFont
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
        
        ; Store into SED Struct
        ;Verbose = 1 ;<DEBUG><20160731>
        FOR sid=0, N_ELEMENTS(PdChi2_SED_selected_id)-1 DO BEGIN
            Temp_iw = FIX(CrabTableReadInfo(PdChi2_LIB[i],"# CVAR1"))
            Temp_if = FIX(CrabTableReadInfo(PdChi2_LIB[i],"# CVAR2"))
            Temp_Nw = FIX(CrabTableReadInfo(PdChi2_LIB[i],"# NVAR1")) ; Number of wavelength, i.e. NAXIS1
            Temp_iLib = PdChi2_Col_iLib[PdChi2_SED_selected_id[sid]]
            Temp_aLib = PdChi2_Col_aLib[PdChi2_SED_selected_id[sid]]
            Temp_chi2 = PdChi2_Col_chi2[PdChi2_SED_selected_id[sid]]
            ;PRINT, ""
            ;PRINT, "<TODO><DEBUG><20160731>"
            PRINT, "Analyzing LIB "+STRING(FORMAT='(I0)',i+1)+" SED "+STRING(FORMAT='(I0)',PdChi2_SED_selected_id[sid])+" which has a fitted chi2 of "+STRING(FORMAT='(G0)',Temp_chi2)
            PRINT, "    Reading SED Model from "+PdChi2_LIB[i]
            PRINT, "    Reading Column "+STRING(FORMAT='(I0)',Temp_iw)+" and "+STRING(FORMAT='(I0)',Temp_if)+" as wave (um) and flux (mJy)"
            PRINT, "    Reading "+STRING(FORMAT='(I0)',Temp_Nw)+" valid lines after skipping "+STRING(FORMAT='(I0)',Temp_iLib)+" valid lines"
            Temp_wLib = DOUBLE(CrabTableReadColumn(PdChi2_LIB[i],Temp_iw,SkipValidLines=Temp_iLib,ReadValidLines=Temp_Nw,FirstReadLineContent=Temp_sLib,Verbose=Verbose))
            Temp_fLib = DOUBLE(CrabTableReadColumn(PdChi2_LIB[i],Temp_if,SkipValidLines=Temp_iLib,ReadValidLines=Temp_Nw,FirstReadLineContent=Temp_sLib,Verbose=Verbose))
            Temp_aPar = MAKE_ARRAY(PdChi2_NPars[i],/DOUBLE)
            Temp_tPar = MAKE_ARRAY(PdChi2_NPars[i],/STRING)
            Temp_sPar = ''
            FOR j=0, PdChi2_NPars[i]-1 DO BEGIN
                Temp_iPar = FIX(CrabTableReadInfo(PdChi2_LIB[i],"# CPAR"+STRING(FORMAT='(I0)',j+1)))
                Temp_tPar[j] = CrabTableReadInfo(PdChi2_LIB[i],"# TPAR"+STRING(FORMAT='(I0)',j+1))
                Temp_aPar[j] = DOUBLE(CrabTableReadColumn(PdChi2_LIB[i],Temp_iPar,SkipValidLines=Temp_iLib,ReadValidLines=1,Verbose=Verbose))
                Temp_sPar = Temp_sPar + Temp_tPar[j] + "=" + STRING(FORMAT='(G0)',Temp_aPar[j])
                IF j LT PdChi2_NPars[i]-1 THEN Temp_sPar = Temp_sPar + ", "
            ENDFOR
            (PdChi2_SED_List[sid].i0) = PdChi2_SED_selected_id[sid]
            (PdChi2_SED_List[sid].chi2) = Temp_chi2
            (PdChi2_SED_List[sid].iLib[i]) = Temp_iLib
            (PdChi2_SED_List[sid].aLib[i]) = Temp_aLib
            *(PdChi2_SED_List[sid].tPar[i]) = Temp_tPar
            *(PdChi2_SED_List[sid].aPar[i]) = Temp_aPar
            *(PdChi2_SED_List[sid].wLib[i]) = Temp_wLib ; MAKE_ARRAY(N_ELEMENTS(Temp_wLib),/DOUBLE)
            *(PdChi2_SED_List[sid].fLib[i]) = Temp_fLib ; MAKE_ARRAY(N_ELEMENTS(Temp_fLib),/DOUBLE)
            
            FOR k=0,2 DO BEGIN
                PRINT, Temp_wLib[k], Temp_fLib[k], "     ", Temp_sPar
            ENDFOR
            FOR k=N_ELEMENTS(Temp_wLib)-4,N_ELEMENTS(Temp_wLib)-1 DO BEGIN
                PRINT, Temp_wLib[k], Temp_fLib[k], "     ", Temp_sPar
            ENDFOR
            PRINT, ""
            PRINT, ""
            
            ;BREAK
        ENDFOR
        
        ;<DEBUG>
        ;BREAK
    ENDFOR
    
    
    
    ; Save PdChi2_SED_List
    ;save, PdChi2_SED_List
    
    
    
    ; Plot the SEDs with best chi2
    SET_PLOT, 'PS' & XSizeInCM=10.5 & YSizeInCM=7 ; & SET_FONT="NGC"
    DEVICE, FILENAME=OutputName+"_"+"SED"+".eps", $
            /COLOR, BITS_PER_PIXEL=8, DECOMPOSED=1, /ENCAPSULATED, XSIZE=XSizeInCM, YSIZE=YSizeInCM
    IF N_ELEMENTS(SET_FONT) GT 0 THEN DEVICE, SET_FONT=SET_FONT, /TT_FONT
    IF N_ELEMENTS(SET_FONT) GT 0 THEN PlotFont=1
    IF N_ELEMENTS(SET_XRange) EQ 1 THEN PlotXRange=SET_XRange ELSE PlotXRange=([0.2,3e5])
    IF N_ELEMENTS(SET_YRange) EQ 1 THEN PlotYRange=SET_YRange ELSE PlotYRange=([1e-4,1e4])
    IF N_ELEMENTS(SET_XTitle) EQ 1 THEN PlotXTitle=SET_XTitle ELSE PlotXTitle='Observed Wavelength ['+'!Z(03BC)'+'m]' ; \mu
    IF N_ELEMENTS(SET_YTitle) EQ 1 THEN PlotYTitle=SET_YTitle ELSE PlotYTitle='Flux Density [mJy]'
    IF N_ELEMENTS(SET_Colors) EQ 1 THEN PlotColors=SET_Colors ELSE PlotColors=['green','yellow','red','blue','cyan'] ;<TODO> plot color
    PLOT, [0.0], [0.0], XRange=PlotXRange, YRange=PlotYRange, /NoData, XTitle=PlotXTitle, YTitle=PlotYTitle, FONT=PlotFont, $
                        /XLOG, /YLOG, XStyle=1, YStyle=1, XThick=2, YThick=2, Thick=2, XMinor=9, YMinor=9, XTickInterval=1, YTickInterval=1
    
    ; Plot each Library
    PRINT, ""
    PRINT, ""
    FOR sid=N_ELEMENTS(PdChi2_SED_selected_id)-1,0,-1 DO BEGIN
        ;<DEBUG><20160731>
        ;PRINT, PdChi2_SED_List[sid]
        ;HELP, PdChi2_SED_List[sid]
        ; 
        ; Print info
        PRINT, "Computing total flux for SED "+STRING(FORMAT='(I0)',PdChi2_SED_List[sid].i0)+" with chi2 "+STRING(FORMAT='(G0)',PdChi2_SED_List[sid].chi2)
        ; 
        ; Prepare to compute total flux
        Temp_wBin = 0.025D ;<TODO> interpolate interval
        Temp_wTot = 10^(CrabArrayIndGen(ALOG10(PlotXRange[0]),ALOG10(PlotXRange[1]),Temp_wBin)) 
        Temp_fTot = Temp_wTot * 0.0D
        ; 
        ; Loop each library and plot with different colors
        FOR i=0,N_ELEMENTS(PdChi2_NPars)-1 DO BEGIN
            
            ; prepare library coefficient, w, f, etc.
            Temp_aLib =  (PdChi2_SED_List[sid].aLib[i])
            Temp_wLib = *(PdChi2_SED_List[sid].wLib[i])
            Temp_fLib = *(PdChi2_SED_List[sid].fLib[i])
            
            ; <TODO> fix negative coefficient
            IF Temp_aLib LT 0.0 THEN Temp_aLib = 0.0D
            
            ; redshift rest-frame library w, normalize library f
            Temp_wLib = (Temp_wLib) * (1.0+Redshift)
            Temp_fLib = (Temp_fLib * Temp_aLib)
            
            ; interpolate to total flux wavelength grid
            Temp_iInterpol = WHERE((Temp_wTot GE MIN(Temp_wLib)) AND (Temp_wTot LE MAX(Temp_wLib)),/NULL)
            IF N_ELEMENTS(Temp_iInterpol) GT 0 AND TOTAL(Temp_fLib) GT 0 THEN BEGIN
                Temp_fInterpol = 10^(Interpol(ALOG10(Temp_fLib), ALOG10(Temp_wLib), ALOG10(Temp_wTot)))
                Temp_fTot[Temp_iInterpol] = Temp_fTot[Temp_iInterpol] + Temp_fInterpol[Temp_iInterpol]
                ;PRINT, "Interpolating tota flux ", CrabMinMax(Temp_wLib)
            ENDIF
            
            Temp_tPar = *(PdChi2_SED_List[sid].tPar[i])
            Temp_aPar = *(PdChi2_SED_List[sid].aPar[i])
            Temp_sPar = ''
            FOR j=0,PdChi2_NPars[i]-1 DO BEGIN
                Temp_sPar = Temp_sPar + Temp_tPar[j] + "=" + STRING(FORMAT='(G0)',Temp_aPar[j]) + ", "
                IF j LT PdChi2_NPars[i]-1 THEN Temp_sPar = Temp_sPar + ", "
            ENDFOR
            
            OPLOT, Temp_wLib, Temp_fLib, Color=cgColor(PlotColors[i]), LINESTYLE=1
            
            ;<DEBUG><20160731>
            ;PRINT, ""
            ;FOR k=0,2 DO BEGIN
            ;    PRINT, Temp_wLib[k], Temp_fLib[k], "     ", Temp_sPar
            ;ENDFOR
            ;FOR k=N_ELEMENTS(Temp_wLib)-4,N_ELEMENTS(Temp_wLib)-1 DO BEGIN
            ;    PRINT, Temp_wLib[k], Temp_fLib[k], "     ", Temp_sPar
            ;ENDFOR
            ;PRINT, ""
            
            ; ***********************************
            ; * Compute stellar mass 
            ; ***********************************
            IF STRMATCH(PdChi2_LIB[i],'FSPS.Padova.*.lib.SED') THEN BEGIN
                IF Temp_tPar[2] EQ 'Mass' THEN BEGIN
                    PdChi2_SED_List[sid].Mstar = Temp_aPar[2] * Temp_aLib / (1.0+Redshift) ; 1+z is because of rShift.sm
                    PRINT, "Computing stellar mass   "+STRING(FORMAT='(E0.6)', PdChi2_SED_List[sid].Mstar)+" [Msun] (z="+STRING(FORMAT='(F0.4)',Redshift)+")"
                ENDIF
            ENDIF
            
            ; ***********************************
            ; * Compute dust mass 
            ; ***********************************
            IF STRMATCH(PdChi2_LIB[i],'DL07.HiExCom.lib.SED') THEN BEGIN
                PdChi2_SED_List[sid].Mdust_Warm = Temp_aLib * Distance^2 / (1.0+Redshift) ; 1+z is because of rShift.sm
                PRINT, "Computing warm dust mass "+STRING(FORMAT='(E0.6)', Temp_aLib * Distance^2 / (1.0+Redshift))+" [Msun] (z="+STRING(FORMAT='(F0.4)',Redshift)+")"
            ENDIF
            IF STRMATCH(PdChi2_LIB[i],'DL07.LoExCom.lib.SED') THEN BEGIN
                PdChi2_SED_List[sid].Mdust_Cold = Temp_aLib * Distance^2 / (1.0+Redshift) ; 1+z is because of rShift.sm
                PRINT, "Computing cold dust mass "+STRING(FORMAT='(E0.6)', Temp_aLib * Distance^2 / (1.0+Redshift))+" [Msun] (z="+STRING(FORMAT='(F0.4)',Redshift)+")"
            ENDIF
            
            ; ***********************************
            ; * Compute AGN luminosity
            ; ***********************************
            IF STRMATCH(PdChi2_LIB[i],'MullaneyAGN.lib.SED') THEN BEGIN
                ;Temp_iTot = WHERE((Temp_wLib GE 3.0*(1.0+Redshift)) AND (Temp_wLib LE 1000.0*(1.0+Redshift)),/NULL)
                ;Temp_SAGN = TOTAL(Temp_fLib * (2.99792458e5/Temp_wLib) * Temp_wBin / ALOG10(exp(1)) )           ; integrated flux in mJy GHz -- daddi method
                ;Temp_SAGN = TOTAL(Temp_fLib[Temp_iTot] * (2.99792458e5/Temp_wLib[Temp_iTot]) * (10^Temp_wBin-(1.0/10^Temp_wBin))/2.0 ) ; integrated flux in mJy GHz -- dzliu method <TODO> Temp_wBin for wLib is not the same one as wTot!
                Temp_freq = CrabArrayIndGen(MIN(2.99792458e5/Temp_wLib),MAX(2.99792458e5/Temp_wLib),1.0)
                Temp_flog = Interpol(ALOG10(Temp_fLib),ALOG10(2.99792458e5/Temp_wLib),ALOG10(Temp_freq))
                Temp_SAGN = TOTAL(10^(Temp_flog)) ; mJy GHz -- dzliu stupid method
                ;PRINT, Temp_SAGN
                PdChi2_SED_List[sid].LAGN = Temp_SAGN / 1D20 * 4*!PI*Distance^2 * 9.52140D44 / 3.839D26 ; 1D20 converts mJy GHz to W m-2, 9.52140D44 converts Mpc^2 to m^2. 
                PRINT, "Computing AGN luminosity "+STRING(FORMAT='(E0.6)', PdChi2_SED_List[sid].LAGN)+" [Lsun] (z="+STRING(FORMAT='(F0.4)',Redshift)+")"
            ENDIF
            
        ENDFOR
        ; 
        ; Compute LIR and SFR
        IF 1 EQ 1 THEN BEGIN
            Temp_iTot = WHERE((Temp_wTot GE 8.0*(1.0+Redshift)) AND (Temp_wTot LE 1000.0*(1.0+Redshift)),/NULL)
            ;Temp_fTIR = TOTAL(Temp_fTot[Temp_iTot] * (2.99792458e5/Temp_wTot[Temp_iTot]) * Temp_wBin / ALOG10(exp(1)) )           ; integrated flux in mJy GHz -- daddi method
            ;Temp_fTIR = TOTAL(Temp_fTot[Temp_iTot] * (2.99792458e5/Temp_wTot[Temp_iTot]) * (10^Temp_wBin-(1.0/10^Temp_wBin))/2.0 ) ; integrated flux in mJy GHz
            Temp_freq = CrabArrayIndGen(MIN(2.99792458e5/Temp_wTot[Temp_iTot]),MAX(2.99792458e5/Temp_wTot[Temp_iTot]),1.0)
            Temp_flog = Interpol(ALOG10(Temp_fTot[Temp_iTot]),ALOG10(2.99792458e5/Temp_wTot[Temp_iTot]),ALOG10(Temp_freq))
            Temp_fTIR = TOTAL(10^(Temp_flog)) ; mJy GHz -- dzliu stupid method
            ;PRINT, Temp_fTIR
            PdChi2_SED_List[sid].LIR = Temp_fTIR / 1D20 * 4*!PI*Distance^2 * 9.52140D44 / 3.839D26 ; 1D20 converts mJy GHz to W m-2, 9.52140D44 converts Mpc^2 to m^2. 
            PdChi2_SED_List[sid].SFR = PdChi2_SED_List[sid].LIR / 1D10 ; Chabrier (2003) IMF
                PRINT, "Computing IR luminosity  "+STRING(FORMAT='(E0.6)', PdChi2_SED_List[sid].LIR)+" [Lsun] (z="+STRING(FORMAT='(F0.4)',Redshift)+")"
        ENDIF
        ; 
        ; Compute Mdust
        IF 1 EQ 1 THEN BEGIN
            PdChi2_SED_List[sid].Mdust = PdChi2_SED_List[sid].Mdust_Warm + PdChi2_SED_List[sid].Mdust_Cold
        ENDIF
        ; 
        ; Compute Umean
        IF 1 EQ 1 THEN BEGIN
            PdChi2_SED_List[sid].U = PdChi2_SED_List[sid].LIR / PdChi2_SED_List[sid].Mdust / 125.0D
        ENDIF
        ; 
        ; Compute Tdust
        IF 1 EQ 1 THEN BEGIN
            PdChi2_SED_List[sid].Tdust = 17.0D * (PdChi2_SED_List[sid].U)^(1.0D/6.0D)
        ENDIF
        ; 
        ; Add radio flux SED by computing L_IR_8_1000
        IF 1 EQ 1 THEN BEGIN
            Temp_iTot = WHERE((Temp_wTot GE 8) AND (Temp_wTot LE 1000),/NULL)
            Temp_fTIR = TOTAL(Temp_fTot[Temp_iTot] * (2.99792458e5/Temp_wTot[Temp_iTot]) * Temp_wBin / ALOG10(exp(1)) ) ; integrated flux in mJy GHz
            Temp_wRadio = Temp_wTot
            Temp_fRadio = Temp_fTIR*1e9 / 3.75e12 / 10^2.4 * ((2.99792458e5/Temp_wTot)/1.4)^(-0.80) ; # 3.75e12 is from arxiv.org/pdf/1005.1072, 10^2.4 is 10^qIR. 
            OPLOT, Temp_wRadio, Temp_fRadio, Color=cgColor(PlotColors[N_ELEMENTS(PlotColors)-1]), LINESTYLE=1
            Temp_fTot = Temp_fTot + Temp_fRadio
        ENDIF
        ; 
        ; Compute qIR
        IF 1 EQ 1 THEN BEGIN
            Temp_iTot = WHERE((Temp_wTot GE 8) AND (Temp_wTot LE 1000),/NULL)
            Temp_fTIR = TOTAL(Temp_fTot[Temp_iTot] * (2.99792458e5/Temp_wTot[Temp_iTot]) * Temp_wBin / ALOG10(exp(1)) ) ; integrated flux in mJy GHz
            Temp_wRadio = 2.99792458e5 / 1.4D ; 1.4GHz = 2e5um
            Temp_fRadio = 10^(Interpol(ALOG10(Temp_fTot),ALOG10(Temp_wTot),ALOG10(Temp_wRadio)))
            PdChi2_SED_List[sid].qIR = ALOG10(Temp_fTIR*1e9 / 3.75e12 / Temp_fRadio) ; # 3.75e12 is from arxiv.org/pdf/1005.1072
        ENDIF
        ; 
        ; Plot total flux SED
        IF sid EQ 0 THEN BEGIN
            OPLOT, Temp_wTot, Temp_fTot
        ENDIF ELSE BEGIN
            OPLOT, Temp_wTot, Temp_fTot, Color=PdChi2_SED_colors[sid]
        ENDELSE
        ; 
        ;<DEBUG><20160731> total flux
        ;PRINT, ""
        ;FOR k=0,2 DO BEGIN
        ;    PRINT, Temp_wTot[k], Temp_fTot[k]
        ;ENDFOR
        ;FOR k=N_ELEMENTS(Temp_wTot)-4,N_ELEMENTS(Temp_wTot)-1 DO BEGIN
        ;    PRINT, Temp_wTot[k], Temp_fTot[k]
        ;ENDFOR
        ;PRINT, ""
        ; 
        ; Overplot Legend if sid EQ 0
        IF sid EQ 0 THEN BEGIN
            PlotLegendX = 0.24
            PlotLegendY = 0.84
            PlotLegendDY = 0.045
            PlotLegendSize = 0.65
            PRINT, Source
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, Source                                                                        , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("z=",F0.4)',Redshift)                                         , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("dL=",I0," Mpc")',Distance)                                   , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("SFR=",F0.2," Msun/yr")',PdChi2_SED_List[sid].SFR)            , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("lgMstar=",F0.2," Msun")',ALOG10(PdChi2_SED_List[sid].Mstar)) , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("lgMdust=",F0.2," Msun")',ALOG10(PdChi2_SED_List[sid].Mdust)) , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("Umean=",F0.2)',PdChi2_SED_List[sid].U)                       , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("Tdust=",F0.2, " K")',PdChi2_SED_List[sid].Tdust)             , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("qIR=",F0.2)',PdChi2_SED_List[sid].qIR)                       , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
            XYOUTS, PlotLegendX, PlotLegendY, /NORMAL, STRING(FORMAT='("chi2=",F0.4)',PdChi2_SED_List[sid].chi2)                     , FONT=PlotFont, CharSize=PlotLegendSize & PlotLegendY=PlotLegendY-PlotLegendDY
        ENDIF
        
        ;;BREAK
    ENDFOR
    ; 
    ; Overplot the observed data points
    PlotSize = 0.45
    PlotThick = 1.5
    Obs_w  = DOUBLE(CrabTableReadColumn(PdChi2_DAT,1))
    Obs_f  = DOUBLE(CrabTableReadColumn(PdChi2_DAT,2))
    Obs_df = DOUBLE(CrabTableReadColumn(PdChi2_DAT,3))
    Obs_iDetected = WHERE(Obs_f GE 2.0D*Obs_df,/NULL)
    IF N_ELEMENTS(Obs_iDetected) GT 0 THEN BEGIN
        ;OPLOT, Obs_w[Obs_iDetected], Obs_f[Obs_iDetected]
        FOR k=0,N_ELEMENTS(Obs_iDetected)-1 DO BEGIN
            ; constraint S/N
            IF Obs_f[Obs_iDetected[k]]/Obs_df[Obs_iDetected[k]] GT 10.0D THEN BEGIN
                Obs_df[Obs_iDetected[k]] = Obs_f[Obs_iDetected[k]]/10.0D
            ENDIF
            Temp_f_high = ALOG10(Obs_f[Obs_iDetected[k]]) + Obs_df[Obs_iDetected[k]]/Obs_f[Obs_iDetected[k]]*1.08
            Temp_f_low = ALOG10(Obs_f[Obs_iDetected[k]]) - Obs_df[Obs_iDetected[k]]/Obs_f[Obs_iDetected[k]]*1.08
            Temp_f_high = 10^Temp_f_high
            Temp_f_low = 10^Temp_f_low
            ; plot a square symbol data point
            USERSYM, [-1.0,-1.0]*PlotSize, $
                     [+0.0,+1.0]*PlotSize, FILL=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]], PSYM=8
            USERSYM, [-1.0,+1.0]*PlotSize, $
                     [+1.0,+1.0]*PlotSize, FILL=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]], PSYM=8
            USERSYM, [+1.0,+1.0]*PlotSize, $
                     [+1.0,-1.0]*PlotSize, FILL=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]], PSYM=8
            USERSYM, [+1.0,-1.0]*PlotSize, $
                     [-1.0,-1.0]*PlotSize, FILL=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]], PSYM=8
            USERSYM, [-1.0,-1.0]*PlotSize, $
                     [-1.0,+0.0]*PlotSize, FILL=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]], PSYM=8
            ; and plot error bar
            USERSYM, [-1.2,+1.2]*PlotSize, $
                     [+0.0,+0.0]*PlotSize, FILL=0, Thick=PlotThick
            ;OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]+Obs_df[Obs_iDetected[k]]], PSYM=8
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Temp_f_high], PSYM=8
            USERSYM, [-1.2,+1.2]*PlotSize, $
                     [+0.0,+0.0]*PlotSize, FILL=0, Thick=PlotThick
            ;OPLOT, [Obs_w[Obs_iDetected[k]]], [Obs_f[Obs_iDetected[k]]-Obs_df[Obs_iDetected[k]]], PSYM=8
            OPLOT, [Obs_w[Obs_iDetected[k]]], [Temp_f_low], PSYM=8
            USERSYM, [+0.0,+0.0]*PlotSize, $
                     [+0.0,+0.0]*PlotSize, FILL=0, Thick=PlotThick
            ;OPLOT, [Obs_w[Obs_iDetected[k]],$
            ;        Obs_w[Obs_iDetected[k]]], $
            ;       [Obs_f[Obs_iDetected[k]]-Obs_df[Obs_iDetected[k]],$
            ;        Obs_f[Obs_iDetected[k]]+Obs_df[Obs_iDetected[k]]], LINESTYLE=0, Thick=PlotThick
            OPLOT, [Obs_w[Obs_iDetected[k]],$
                    Obs_w[Obs_iDetected[k]]], $
                   [Temp_f_high,$
                    Temp_f_low], LINESTYLE=0, Thick=PlotThick
        ENDFOR
        ;OPLOTERR, Obs_w[Obs_iDetected], Obs_f[Obs_iDetected], Obs_df[Obs_iDetected], HATLENGTH=0.02
    ENDIF
    Obs_iUndetected = WHERE(Obs_f LT 2.0D*Obs_df,/NULL)
    IF N_ELEMENTS(Obs_iUndetected) GT 0 THEN BEGIN
        USERSYM, [-1.2,+1.2]*PlotSize, $
                 [+0.0,+0.0]*PlotSize, FILL=0, Thick=PlotThick
        OPLOT, Obs_w[Obs_iUndetected], 3.0D*Obs_df[Obs_iUndetected], PSYM=8
        USERSYM, [+0.0,+0.0]*PlotSize, $
                 [+0.0,-2.6]*PlotSize, FILL=0, Thick=PlotThick
        OPLOT, Obs_w[Obs_iUndetected], 3.0D*Obs_df[Obs_iUndetected], PSYM=8
        USERSYM, [+1.0,+0.0,-1.0]*PlotSize, $
                 [-1.2,-2.6,-1.2]*PlotSize, FILL=0, Thick=PlotThick
        OPLOT, Obs_w[Obs_iUndetected], 3.0D*Obs_df[Obs_iUndetected], PSYM=8
    ENDIF
    ; 
    ; Close device
    DEVICE, /CLOSE
    SET_PLOT, "X"
    
    
    
END