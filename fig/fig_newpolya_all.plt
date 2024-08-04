# Gnuplot Template for multiplot 
# 5.2 patchlevel 7
# polyA profile at invididual position, based on "_histend3_C20.stat"
unset multiplot
reset
set term wxt font "arial, 5"
set term pdfcairo color font "Arial, 5" size 20cm,20cm solid linewidth 0.5

# 8 directories: 35S[1], 32-33S [2], 23S[3], 20-22S[4],	18S[5],	25S[6], 5.8S-7S[7], 27S[8]
# array S[14]: offset for starting position for each RNA. Need manually adjusted for each file to properly align the 0 position.
# If offset goes to negative, artifical empty blocks need to be added in "x_histend3_C20.stat"
# array S[14]: 23S[1], 22S[2], 21S[3], 20S[4], 35S[5], 33S[6], 32S[7], 27SA2[8], 27SA3[9], 27SB[10], 25S[11], 7S[12], 58S[13], 18S[14]
array DIR[11] =['../K31/', '../K32/', '../K33/', '../K34/', '../K35/', '../K36/', '../K37/', '../K39/', 'rex1 NOP7 rep2', '15'];array S[14] = [0, 0, 0, 129, 12, 5, 7, 3, 1, 5, 7, 33, 4, 1]



#unset key
set macros
# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set format x ''; unset xlabe5"
XTICS = "set format x; set xlabel ''"
NOYTICS = "set format y ''; unset ylabel"
#YTICS = "set format y '%g' ; set ylabel 'Fraction' offset 1,0"

## Histogram
#XTICS = "set xtics nomirror out scale 0 ('0' 0, '1-3' 1,'4-9' 2, '>10' 3); set xlabel 'Length'"
YTICS = "set format y '%g'; set ytics 0.2 offset 0,0 ; set mytic 2; set ylabel 'Fraction' offset 0,0"

POS = "at graph 0.0,1.1"
POS2 = "at graph 0.05,0.9"

XRANG = "set xrang [-1:50]; set xtic in 10 offset 0,0.4; set mxtic 5"

YRANGLOG = "set logscale y; set yrang [1e-5:1]; set ytic 10; set mytic 10"
YRANG = "unset logscale y; set yrang [0:0.2]; set ytic 0.2; set mytic 1"

YRANG_HIS = "unset logscale y; set yrang [0:1]"
# points 64 size = 1 column, data 62
XRANG_HIS = "set xrang [-1:63]; set xtic out 1 offset 0,0.4 rotate 90 nomirror font 'Arial, 4'; set mxtic 1"
# points 32 size = 0.5 column, data 30
XRANG_HIS2 = "set xrang [-1:31]; set xtic out 1 offset 0,0.4 rotate 90 nomirror font 'Arial, 4'; set mxtic 1"
# points 96 size = 1.5 column, 7S data = 89
XRANG_HIS3 = "set xrang [-1:95]; set xtic out 1 offset 0,0.4 rotate 90 nomirror font 'Arial, 4'; set mxtic 1"
# points 16 size = 0.25 column, data <= 14
XRANG_HIS4 = "set xrang [-1:15]; set xtic out 1 offset 0,0.4 rotate 90 nomirror font 'Arial, 4'; set mxtic 1"


TOT3D= "unset logscale xy;set xrange [-1e4:1e4]; set yrange [0:1e6]; stats dat1 using 2 nooutput; a = STATS_sum; stats dat2 using 2 nooutput; b = STATS_sum; stats dat3 using 2 nooutput; c = STATS_sum;  total = a + b + c"

TOT2D= "unset logscale xy;set xrange [-1e4:1e4]; set yrange [0:1e6]; stats dat1 using 2 nooutput; a = STATS_sum; stats dat2 using 2 nooutput; b = STATS_sum; total = a + b"

# global setting
set tics scale 0.4
set style data histograms
#set style histogram clustered gap 2
set style histogram columnstacked
set boxwidth 0.5 absolute
set style fill solid 1.0

print "fig_S".DIR[10]."_pA3_".DIR[9]."_b6.pdf"
set output "fig_S".DIR[10]."_pA3_".DIR[9]."_b6.pdf"
set multiplot title "Fig S".DIR[10].". ".DIR[9] font "Arial, 12"

# overall size
row_size = 0.12; col_size = 0.4; 

# adjust gap.  top_margin2 = top_margin (gap=0) 
top_margin = 0.95; left_margin = 0.1; 
top_margin2 = 1.01; left_margin2 = 0.1; 
# Although 8X4 panels are set here, you can use only 2x2. Just leave others
## Graph 1,1
set tmargin at screen top_margin - 0*row_size
set bmargin at screen top_margin2 - 1*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 2 "I" font "Arial, 12" at graph -0.1, 1.3
set label 1 "23S" @POS
set label 3 "A3" at 61,1.1 center
set arrow 3 from 61,1.0 to 61,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
dat1_1= DIR[3]."b6_23S_histend3_C20.stat"
plot for [IDX=S[1]:S[1]+61] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 3; unset arrow 3
unset label 2

## Graph 1,2
set tmargin at screen top_margin - 0*row_size
set bmargin at screen top_margin2 - 1*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "35S" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[1]."b6_35S_histend3_C20.stat"
plot for [IDX=S[5]:S[5]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4

## Graph 1,3
set tmargin at screen top_margin - 0*row_size
set bmargin at screen top_margin2 - 1*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 1,4
set tmargin at screen top_margin - 0*row_size
set bmargin at screen top_margin2 - 1*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 2,1
set tmargin at screen top_margin - 1*row_size
set bmargin at screen top_margin2 - 2*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 1 "22S" @POS
set label 3 "A3" at 61,1.1 center
set arrow 3 from 61,1.0 to 61,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
dat1_1= DIR[4]."b6_22S_histend3_C20.stat"
plot for [IDX=S[2]:S[2]+61] dat1_1 i IDX using 4 t columnhead(4) with histogram

## Graph 2,2
set tmargin at screen top_margin - 1*row_size
set bmargin at screen top_margin2 - 2*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "33S" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[2]."b6_33S_histend3_C20.stat"
plot for [IDX=S[6]:S[6]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4

## Graph 2,3
set tmargin at screen top_margin - 1*row_size
set bmargin at screen top_margin2 - 2*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 2,4
set tmargin at screen top_margin - 1*row_size
set bmargin at screen top_margin2 - 2*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 3,1
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 3*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 1 "21S" @POS
set label 3 "A3" at 61,1.1 center
set arrow 3 from 61,1.0 to 61,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
dat1_1= DIR[4]."b6_21S_histend3_C20.stat"
plot for [IDX=S[3]:S[3]+61] dat1_1 i IDX using 4 t columnhead(4) with histogram

## Graph 3,2
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 3*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "32S" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[2]."b6_32S_histend3_C20.stat"
plot for [IDX=S[7]:S[7]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4

## Graph 3,3
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 3*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 3,4
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 3*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 4,1
set tmargin at screen top_margin - 3*row_size
set bmargin at screen top_margin2 - 4*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 1 "20S" @POS
set label 3 "A2" at 37,1.1 center
set arrow 3 from 37,1.0 to 37,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
dat1_1= DIR[4]."b6_20S_histend3_C20.stat"
plot for [IDX=S[4]:S[4]+37] dat1_1 i IDX using 4 t columnhead(4) with histogram

## Graph 4,2
set tmargin at screen top_margin - 3*row_size
set bmargin at screen top_margin2 - 4*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "27SA2" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[8]."b6_27SA2_histend3_C20.stat"
plot for [IDX=S[8]:S[8]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4


## Graph 4,3
set tmargin at screen top_margin - 3*row_size
set bmargin at screen top_margin2 - 4*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 4,4
set tmargin at screen top_margin - 3*row_size
set bmargin at screen top_margin2 - 4*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 5,1
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 0.25*col_size
set label 1 "18S" @POS
set label 3 "D" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS4; @YRANG_HIS
dat1_1= DIR[5]."b6_18S_histend3_C20.stat"
plot for [IDX=S[14]:S[14]+10] dat1_1 i IDX using 4 t columnhead(4) with histogram

## Graph 5,2
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 0.50*col_size
set rmargin at screen left_margin2 + 1.0*col_size
set label 1 "25S" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[6]."b6_25S_histend3_C20.stat"
plot for [IDX=S[11]:S[11]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4


## Graph 5,3
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "27SA3" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[8]."b6_27SA3_histend3_C20.stat"
plot for [IDX=S[9]:S[9]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4

## Graph 5,4
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 6,1
set tmargin at screen top_margin - 5*row_size
set bmargin at screen top_margin2 - 6*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size

set label 1 "5.8S" @POS
set label 3 "E" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
dat1_1= DIR[7]."b6_58S_histend3_C20.stat"
plot for [IDX=S[13]:S[13]+60] dat1_1 i IDX using 4 t columnhead(4) with histogram


## Graph 6,2
set tmargin at screen top_margin - 5*row_size
set bmargin at screen top_margin2 - 6*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size
set label 1 "27SB" @POS
set label 3 "B2" at 10,1.1 center
set arrow 3 from 10,1.0 to 10,1.05 nohead
set label 4 "B0" at 25,1.1 center
set arrow 4 from 25,1.0 to 25,1.05 nohead
@XTICS; @NOYTICS
@XRANG_HIS2; @YRANG_HIS
dat1_1= DIR[8]."b6_27SB_histend3_C20.stat"
plot for [IDX=S[10]:S[10]+29] dat1_1 i IDX using 4 t columnhead(4) with histogram
unset label 4; unset arrow 4

## Graph 6,3
set tmargin at screen top_margin - 5*row_size
set bmargin at screen top_margin2 - 6*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 6,4
set tmargin at screen top_margin - 5*row_size
set bmargin at screen top_margin2 - 6*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 7,1
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1.5*col_size
set label 1 "7S" @POS
set label 3 "C2" at 27,1.1 center
set arrow 3 from 27,1.0 to 27,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS3; @YRANG_HIS
set label 3 "C2" at 88,1.1 center
set arrow 3 from 88,1.0 to 88,1.05 nohead
dat1_1= DIR[7]."b6_7S_histend3_C20.stat"  
plot for [IDX=S[12]-27:S[12]+61] dat1_1 i IDX using 4 t columnhead(4) with histogram

## Graph 7,2
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 1.05*col_size
set rmargin at screen left_margin2 + 1.55*col_size

## Graph 7,3
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 7,4
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 8,1
set tmargin at screen top_margin - 7*row_size
set bmargin at screen top_margin2 - 7.8*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1.5*col_size

set label 1 "" @POS
set label 3 "E" at 0,1.1 center
set arrow 3 from 0,1.0 to 0,1.05 nohead
set label 4 "C2" at 139,1.1 center
set arrow 4 from 139,1.0 to 139,1.05 nohead
@XTICS; @YTICS
@XRANG_HIS; @YRANG_HIS
set xrang [-1:140]
dat1_1= DIR[7]."b6_58S_histend3_C20.stat"
dat1_2= DIR[7]."b6_7S_histend3_C20.stat"  
#plot for [IDX=S[13]+10:S[13]+60] dat1_1 i IDX using 4 t columnhead(4) with histogram, \
#for [IDX=S[12]-27:S[12]+61] dat1_2 i IDX using 4 t columnhead(4) with histogram

## Graph 8,2
set tmargin at screen top_margin - 7*row_size
set bmargin at screen top_margin2 - 8*row_size
set lmargin at screen left_margin + 1*col_size
set rmargin at screen left_margin2 + 2*col_size

## Graph 8,3
set tmargin at screen top_margin - 7*row_size
set bmargin at screen top_margin2 - 8*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 8,4
set tmargin at screen top_margin - 7*row_size
set bmargin at screen top_margin2 - 8*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## END, return to terminal
unset multiplot
set output
set term pop
reset

