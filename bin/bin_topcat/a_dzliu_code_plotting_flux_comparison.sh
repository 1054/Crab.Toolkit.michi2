#!/bin/bash
# 


if [[ $(type topcat 2>/dev/null | wc -l) -ne 1 ]]; then
    echo "Error! Topcat was not installed!"
    exit
fi


if [[ ! -f "Table_from_Shuowen.txt" ]]; then
    echo "Error! \"Table_from_Shuowen.txt\" was not found!"
    exit
fi

if [[ ! -f "../../SED_fitting_michi2_20180110/Results/best-fit_param_LIR_total.txt" ]]; then
    echo "Error! \"../../SED_fitting_michi2_20180110/Results/best-fit_param_LIR_total.txt\" was not found!"
    exit
fi


#DEBUG#
rm datatable_CrossMatched.txt
if [[ ! -f "datatable_CrossMatched.txt" ]]; then
topcat -stilts tmatchn nin=3 \
                        in1="Table_from_Shuowen.txt" ifmt1=ascii \
                        in2="../../SED_fitting_michi2_20180110/Results/best-fit_param_LIR_total.txt" ifmt2=ascii \
                        in3="Table_SNR.txt" ifmt3=ascii \
                        matcher=exact values1="Source" values2="Source" values3="Source" \
                        out="datatable_CrossMatched.txt" ofmt=ascii
fi


margin=(80 50 20 20) # left, bottom, right, top
topcat -stilts plot2plane \
                        xpix=500 ypix=400 \
                        insets="${margin[3]},${margin[0]},${margin[1]},${margin[2]}" \
                        xlabel="\Large \log \, L_{IR, \ Spdb \ (Shuowen)}" \
                        ylabel="\Large \log \, L_{IR, \ michi2 \ (full \ DL07)}" \
                        xlog=false \
                        ylog=false \
                        xmin=11.5 xmax=13.2 ymin=11.5 ymax=13.2 \
                        \
                        in="datatable_CrossMatched.txt" \
                        ifmt=ascii \
                        \
                        layer1=mark \
                        shape1=filled_circle \
                        size1=3 \
                        color1='#1564b2' \
                        shading1=flat \
                        x1="log10(SFR_IR*1e10)" \
                        y1="(LIR_total_best)" \
                        \
                        layer3=xyerror \
                        errorbar3=capped_lines \
                        color3='#1e90ff' \
                        shading3=flat \
                        x3="log10(SFR_IR*1e10)" \
                        y3="(LIR_total_median)" \
                        yerrlo3="(LIR_total_median)-(LIR_total_L68)" \
                        yerrhi3="(LIR_total_H68)-(LIR_total_median)" \
                        \
                        layer4=label \
                        color4='#0c3966' fontsize4=5 spacing4=1 \
                        x4="log10(SFR_IR*1e10)" \
                        y4="LIR_total_best+(RANDOM-0.5)/0.5*0.02" \
                        label4="substring(Source_1,5)" \
                        \
                        layer9=function \
                        fexpr9='(x)' \
                        color9=black \
                        antialias9=true \
                        thick9=1 \
                        leglabel9='1:1' \
                        \
                        legend=false \
                        seq="9,1,3,4" \
                        fontsize=16 \
                        texttype=latex \
                        omode=out \
                        out='Plot_comparison_L_IR.pdf'
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-usage.html
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-examples.html
                        # omode=swing


margin=(80 50 20 20) # left, bottom, right, top
topcat -stilts plot2plane \
                        xpix=500 ypix=400 \
                        insets="${margin[3]},${margin[0]},${margin[1]},${margin[2]}" \
                        xlabel="\Large S/N_{100-1000\,{\mu}m, \ Spdb}" \
                        ylabel="\large \log \, L_{IR \ (full \ DL07)} - \log \, L_{IR \ (Shuowen)}" \
                        xlog=true \
                        ylog=false \
                        xmin=1.0 xmax=100 ymin=-1.0 ymax=1.0 \
                        \
                        in="datatable_CrossMatched.txt" \
                        ifmt=ascii \
                        \
                        layer1=mark \
                        shape1=filled_circle \
                        size1=3 \
                        color1='#1564b2' \
                        shading1=flat \
                        x1="SN_IR" \
                        y1="(LIR_total_best) - log10(SFR_IR*1e10)" \
                        \
                        layer3=xyerror \
                        errorbar3=capped_lines \
                        color3='#1e90ff' \
                        shading3=flat \
                        x3="SN_IR" \
                        y3="(LIR_total_median) - log10(SFR_IR*1e10)" \
                        yerrlo3="(LIR_total_median)-(LIR_total_L68)" \
                        yerrhi3="(LIR_total_H68)-(LIR_total_median)" \
                        \
                        layer4=label \
                        color4='#0c3966' fontsize4=5 spacing4=1 \
                        x4="SN_IR" \
                        y4="(LIR_total_best) - log10(SFR_IR*1e10) + (RANDOM-0.5)/0.5*0.02" \
                        label4="substring(Source_1,5)" \
                        \
                        layer9=function \
                        fexpr9='(0)' \
                        color9=black \
                        antialias9=true \
                        thick9=1 \
                        leglabel9='Y=0' \
                        \
                        legend=false \
                        seq="9,1,3,4" \
                        fontsize=16 \
                        texttype=latex \
                        omode=out \
                        out='Plot_comparison_L_IR_diff.pdf'
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-usage.html
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-examples.html
                        # omode=swing


margin=(80 50 20 20) # left, bottom, right, top
topcat -stilts plot2plane \
                        xpix=500 ypix=400 \
                        insets="${margin[3]},${margin[0]},${margin[1]},${margin[2]}" \
                        xlabel="\Large S/N_{100-1000\,{\mu}m, \ Spdb}" \
                        ylabel="\large \log \, L_{IR \ (full \ DL07)} - \log \, L_{IR \ (Shuowen)}" \
                        xlog=true \
                        ylog=false \
                        xmin=1.0 xmax=100 ymin=-1.0 ymax=1.0 \
                        \
                        in="datatable_CrossMatched.txt" \
                        ifmt=ascii \
                        \
                        layer1=mark \
                        shape1=filled_circle \
                        size1=3 \
                        color1='#1564b2' \
                        shading1=flat \
                        x1="SN_IR" \
                        y1="(LIR_total_median) - log10(SFR_IR*1e10)" \
                        \
                        layer3=xyerror \
                        errorbar3=capped_lines \
                        color3='#1e90ff' \
                        shading3=flat \
                        x3="SN_IR" \
                        y3="(LIR_total_median) - log10(SFR_IR*1e10)" \
                        yerrlo3="(LIR_total_sigma)" \
                        yerrhi3="(LIR_total_sigma)" \
                        \
                        layer4=label \
                        color4='#0c3966' fontsize4=5 spacing4=1 \
                        x4="SN_IR" \
                        y4="(LIR_total_median) - log10(SFR_IR*1e10) + (RANDOM-0.5)/0.5*0.02" \
                        label4="substring(Source_1,5)" \
                        \
                        layer9=function \
                        fexpr9='(0)' \
                        color9=black \
                        antialias9=true \
                        thick9=1 \
                        leglabel9='Y=0' \
                        \
                        legend=false \
                        seq="9,1,3,4" \
                        fontsize=16 \
                        texttype=latex \
                        omode=out \
                        out='Plot_comparison_L_IR_total_median_diff.pdf'
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-usage.html
                        # http://www.star.bristol.ac.uk/~mbt/stilts/sun256/plot2plane-examples.html
                        # omode=swing

