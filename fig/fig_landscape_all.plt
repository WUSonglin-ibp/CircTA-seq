# Gnuplot Template for multiplot 
# 5.2 patchlevel 7
unset multiplot
reset
set term wxt font "arial, 5"
set term pdfcairo color font "Arial, 6" size 20cm,30cm solid linewidth 0.5

# 8 directories: 35S[1], 32-33S [2], 23S[3], 20-22S[4],	18S[5],	25S[6], 5.8S-7S[7], 27S[8]
array DIR[10] =['../K31/', '../K32/', '../K33/', '../K34/', '../K35/', '../K36/', '../K37/', '../K39/', 'rex1 NOP7 rep2', '15']

set output "fig_S".DIR[10]."_landscape_".DIR[9]."_paper_b6.pdf"
print "fig_S".DIR[10]."_landscape_".DIR[9]."_paper_b6.pdf"
set multiplot title "Fig S".DIR[10].". ".DIR[9] font 'Arial, 12'

set colorsequence default
fn(v)=sprintf("%.3f",v)
# Points < shreshold will not disappear, but appear as ghost points, can be removed in AI
cf(v)=(v > 5 ? log10(v) : 1/0)

#unset key
set macros
# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set format x ''; unset xlabel"
XTICS = "set format x; set xlabel '' font 'arial, 5' offset 0,1"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%g' ; set ylabel '' font 'arial, 5' offset 1,0"

## Histogram
#XTICS = "set xtics nomirror out scale 0 ('0' 0, '1-3' 1,'4-9' 2, '>10' 3); set xlabel 'Length'"
#YTICS = "set format y '%.1f'; set ytics scale 0.5 0, 0.2, 0.8; set ylabel 'Fraction'"

UNK = "b6_UNK_coord.dat.sort"
HIS3END = '_histend3.stat'

POS = "font 'arial, 12' at graph 0.6,0.9"
POS2 = "font 'arial, 12' at graph -0.05,1.1"

XRANG = "set xrang [-20:20]; set xtic 10 offset 0,0.2; set mxtic 5"
XRANG_ITS1 = "set xrang [-160:20]; set xtic 10 offset 0,0.2; set mxtic 5"
XRANG_ITS1A = "set xrang [-15:300]; set xtic 10 offset 0,0.2; set mxtic 5"
XRANG_ITS1S = "set xrang [-15:5]; set xtic 5 offset 0,0.2; set mxtic 5"
XRANG_ITS2L = "set xrang [-20:150]; set xtic 10 offset 0,0.2; set mxtic 5"
XRANG_ITS2R = "set xrang [-100:10]; set xtic 10 offset 0,0.2; set mxtic 5"
XRANG_3ETS = "set xrang [-5:20]; set xtic 5 offset 0,0.2; set mxtic 5"
XRANG_A1 = "set xrang [-10:10]; set xtic 5 offset 0,0.2; set mxtic 5"
XRANG_PA = "set xrang [-1:30]; set xtic 5 offset 0,0.2; set mxtic 5"

YRANG = "unset logscale y; set yrang [-0.05:1]; set ytic 0.2; set mytic 1"
YRANGLOG = "set logscale y; set yrang [1e-6:1]; set ytic 10; set mytic 10"
YRANG_PA = "unset logscale y; set yrang [-0.01:0.15]; set ytic 0.05; set mytic 5"
YRANG_HIS = "unset logscale y; set yrang [0:1]; set ytic 0.2; set mytic 2"
XRANG_HIS = "set xrang [-1:14]; set xtic out 1 offset 0,0.4 rotate 90 nomirror; set mxtic 1"

#YRANG = "unset logscale y; set yrang [0:1]; set ytic 0.2; set mytic 1"


TOT3D= "unset logscale xy; set xrange [-1e4:1e4]; set yrange [0:1e9]; stats dat1 using 2 nooutput; a = STATS_sum; stats dat2 using 2 nooutput; b = STATS_sum; stats dat3 using 2 nooutput; c = STATS_sum;  total = a + b + c"

TOT2D= "unset logscale xy; set xrange [-1e4:1e4]; set yrange [0:1e9];  stats dat1 using 2 nooutput; a = STATS_sum; stats dat2 using 2 nooutput; b = STATS_sum; total = a + b;"

# global setting
set tics scale 1.0
set style histogram columnstacked
set boxwidth 0.5 absolute

# headmap
set palette positive
set cbrange [1:4]
#set palette defined (0 'white', 4 'black')
set palette rgb 7,5,15
#set xtics rotate by 90 offset 0, -1
set style fill transparent solid 0.5 noborder

set grid x2tics 
set grid y2tics 
set y2tics out 
set x2tics out 

# overall size
row_size = 0.2; col_size = 0.19; 

# adjust gap.  top_margin2 = top_margin (gap=0) 
top_margin = 0.95; left_margin = 0.1; 
#top_margin2 = 1.0; left_margin2 = 0.06; 
top_margin2 = 0.95; left_margin2 = 0.1; 
# Although 8X4 panels are set here, you can use only 2x2. Just leave others
## Graph 1,1
set tmargin at screen top_margin - 0.1*row_size
set bmargin at screen top_margin2 - 0.4*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 0.2*col_size
set label 10 "A" font 'arial, 12' at graph -0.1,1.3
@XTICS; @YTICS
unset colorbox
set x2tics ("5'" 1, 'A0' 611, 'A1' 701, 'A2' 2713,'A4' 2749, 'A3' 2785, 'B1' 2862,'C2' 3159, 'C1' 3254) offset 0,-0.5 font 'arial, 5'
set y2tics ('D' 2500, 'A2' 2712, 'A3' 2784, 'E' 3019, 'C2' 3158, 'C1' 3253, 'B2' 6647) offset -0.7,0 font 'arial, 5'
set y2tics ('D' 2500, 'A2' 2712, 'A3' 2784, 'E' 3019, 'C2' 3158, 'C1' 3253, '' 6647) offset -0.7,0 font 'arial, 5'
set xrang [-10:10]; set yrang [6640:6670]; 
set xtics 50 nomirror font 'arial, 5'; set mxtics 10; set ytics 6640,10,6670 nomirror font 'arial, 5'; set mytics 5;
set arrow 1 from 1, 6647 to -13, 6657 front nohead
set label 1 "35S" font 'arial, 5' front at -13, 6657 right
plot DIR[1].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle
unset label 10

## Graph 1,2
set tmargin at screen top_margin - 0.1*row_size
set bmargin at screen top_margin2 - 0.4*row_size
set lmargin at screen left_margin + 0.25*col_size
set rmargin at screen left_margin2 + 1.35*col_size
set label 10 "" @POS2
@XTICS; @NOYTICS
set xrang [600:710]; set yrang [6640:6670]; 
set xtics 50 nomirror font 'arial, 5'; set mxtics 10; set ytics 6640,10,6670 nomirror font 'arial, 5'; set mytics 5;
set arrow 1 from 611, 6647 to 621, 6657 front nohead
set label 1 "33S" font 'arial, 5' front at 621, 6657 left
set arrow 2 from 701, 6647 to 691, 6657 front nohead
set label 2 "32S" font 'arial, 5' front at 691, 6657 right
plot DIR[2].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle

## Graph 1,3
set tmargin at screen top_margin - 0.1*row_size
set bmargin at screen top_margin2 - 0.4*row_size
set lmargin at screen left_margin + 1.55*col_size
set rmargin at screen left_margin2 + 3.25*col_size
set label 10 "" @POS2
set xrang [2705:2875]; set yrang [6640:6670]; 
@XTICS; @NOYTICS
set arrow 1 from 2713, 6647 to 2723, 6657 front nohead
set label 1 "27SA2" font 'arial, 5' front at 2723, 6657 left
set arrow 2 from 2785, 6647 to 2795, 6657 front nohead
set label 2 "27SA3" font 'arial, 5' front at 2795, 6657 left
set arrow 3 from 2862, 6647 to 2847, 6657 front nohead
set label 3 "27SB" font 'arial, 5' front at 2847, 6657 right
plot DIR[8].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle

## Graph 1,4
set tmargin at screen top_margin - 0.1*row_size
set bmargin at screen top_margin2 - 0.4*row_size
set lmargin at screen left_margin + 3.3*col_size
set rmargin at screen left_margin2 + 4.4*col_size
set label 10 "" @POS2
@XTICS; @NOYTICS
set xrang [3150:3260]; set yrang [6640:6670]; 
set xtics 50 nomirror; set mxtics 10; set ytics 6640,10,6670 nomirror; set mytics 5;
set y2tics out ('D' 2500, '' 2712, 'A3' 2784, 'E' 3019, '' 3158, 'C1' 3253, 'B2' 6647) offset -0.7,0 font 'arial, 5'
set arrow 1 from 3246, 6647 to 3236, 6657 front nohead
set label 1 "25S'" font 'arial, 5' front at 3236, 6657 right
set arrow 2 from 3254, 6647 to 3264, 6657 front nohead
set label 2 "25S" font 'arial, 5' front at 3264, 6657 left
set arrow 3 from 3159, 6647 to 3169, 6657 front nohead
set label 3 "26S" font 'arial, 5' front at 3169, 6657 left

plot DIR[6].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle

## Graph 2,1
set tmargin at screen top_margin - 0.5*row_size
set bmargin at screen top_margin2 - 1.1*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 0.2*col_size
set label 10 "" @POS2
@XTICS; @YTICS
set xrang [-10:10]; set yrang [2480:2800]; 
set xtics 50 nomirror; set mxtics 10; set ytics 100 nomirror; set mytics 10;
set y2tics out ('' 2500, '' 2712, '' 2784, 'E' 3019, '' 3158, 'C1' 3253, 'B2' 6647) offset -0.7,0 font 'arial, 5'
set ylabel "3' end position" font 'Arial, 7' offset 0,0
set arrow 1 from 1, 2784 to -13, 2754 front nohead
set label 1 "23S" font 'arial, 5' front at -13, 2754 right
plot DIR[3].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle
unset ylabel

## Graph 2,2
set tmargin at screen top_margin - 0.5*row_size
set bmargin at screen top_margin2 - 1.1*row_size
set lmargin at screen left_margin + 0.25*col_size
set rmargin at screen left_margin2 + 1.35*col_size
set label 10 "" @POS2
@XTICS; @NOYTICS
set xrang [600:710]; set yrang [2480:2800]; 
set y2tics ('D' 2500, 'A2' 2712, 'A3' 2784, 'E' 3019, 'C2' 3158, 'C1' 3253, 'B2' 6647) offset -0.7,0 font 'arial, 5'
set xlabel "5' end position" font 'Arial, 7' offset 0,0.2
set arrow 1 from 611, 2784 to 621, 2814 front nohead
set label 1 "22S" font 'arial, 5' front at 621, 2814 left
set arrow 3 from 701, 2784 to 695, 2814 front nohead
set label 3 "21S" font 'arial, 5' front at 695, 2814 right
set arrow 4 from 701, 2712 to 691, 2742 front nohead
set label 4 "20S" font 'arial, 5' front at 691, 2742 right
set arrow 5 from 701, 2500 to 691, 2530 front nohead
set label 5 "18S" font 'arial, 5' front at 691, 2530 right
set arrow 6 from 696, 2784 to 686, 2814 front nohead
set label 6 "21S'" font 'arial, 5' front at 686, 2814 right
plot DIR[5].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle
unset xlabel

## Graph 2,3
set tmargin at screen top_margin - 0.5*row_size
set bmargin at screen top_margin2 - 0.8*row_size
set lmargin at screen left_margin + 1.55*col_size
set rmargin at screen left_margin2 + 3.25*col_size
set label 10 "" @POS2
@XTICS; @YTICS
set xrang [2705:2875]; set yrang [3000:3170]; 

set colorbox user origin graph 1.1, 0.0 size 0.01, 0.06
set cbtics ('1' 0, '10' 1, '100' 2, '1000' 3, '10000' 4, '100000' 5)

set arrow 5 from 2862, 3158 to 2872, 3188 front nohead
set label 5 "7S" font 'arial, 5' front at 2872, 3188 left
set arrow 6 from 2862, 3019 to 2872, 3049 front nohead
set label 6 "5.8S" font 'arial, 5' front at 2872, 3049 left

plot DIR[7].UNK using ($1):($2):(2):(cf($3)) notitle with circles palette z, x notitle
unset colorbox

## Graph 2,4
set tmargin at screen top_margin - 1.2*row_size
set bmargin at screen top_margin2 - 1.7*row_size
set lmargin at screen left_margin + 3.3*col_size
set rmargin at screen left_margin2 + 4.4*col_size


## Graph 3,1
set tmargin at screen top_margin - 1.3*row_size
set bmargin at screen top_margin2 - 1.8*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 1 "B" @POS2
#set label 2 "5'ETS" at graph 0.02, 0.9
@XTICS; @YTICS
@XRANG_A1
@YRANG
set style fill solid 1.0
set x2tics ('A1' 0, "A1'" -5) offset 0,-0.5 font 'arial, 6'
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
plot \
DIR[2]."b6_32S_5EXT.dat" using (-$1):3 t "32S" with lines lw 1, \
DIR[4]."b6_21S_5EXT.dat" using (-$1):3 t "21S" with lines lw 1, \
DIR[4]."b6_20S_5EXT.dat" using (-$1):3 t "20S" with lines lw 1, \
DIR[5]."b6_18S_5EXT.dat" using (-$1):3 t "18S" with lines lw 1

## Graph 3,2
set tmargin at screen top_margin - 1.3*row_size
set bmargin at screen top_margin2 - 1.8*row_size
set lmargin at screen left_margin + 1.2*col_size
set rmargin at screen left_margin2 + 4.2*col_size

## Graph 3,3
set tmargin at screen top_margin - 1.5*row_size
set bmargin at screen top_margin2 - 1.8*row_size
set lmargin at screen left_margin + 1.2*col_size
set rmargin at screen left_margin2 + 4.2*col_size


## Graph 3,4
set tmargin at screen top_margin - 1.4*row_size
set bmargin at screen top_margin2 - 2.4*row_size
set lmargin at screen left_margin + 3.5*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 4,1
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 2.5*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 2*col_size
set label 1 "C" @POS2
#set label 2 "ITS1" at graph 0.1, 0.9
@XTICS; @YTICS
dat1 = DIR[8]."b6_27SA2_5EXT.dat"; dat2 = DIR[8]."b6_27SA3_5EXT.dat";  dat3 = DIR[8]."b6_27SB_5EXT.dat"; @TOT3D; total1 = total
@XRANG_ITS1
@YRANGLOG
#set key at graph 0.8,0.9
set key center top
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
set x2tics ('B1' 0, 'A3' -77,  'A4' -113, 'A2' -149) offset 0,-0.5 font 'arial, 6'
#set format y '10^%T'
plot \
DIR[8]."b6_27SA2_5EXT.dat" using (-$1-77-72):($2/total1) t "27S"  w lines  lc 1, \
DIR[8]."b6_27SA3_5EXT.dat" every ::4 using (-$1-77):($2/total1) notitle  w lines lc 1, \
DIR[8]."b6_27SB_5EXT.dat"  using (-$1):($2/total1) notitle w lines lc 1
set key default

## Graph 4,2
set tmargin at screen top_margin - 2*row_size
set bmargin at screen top_margin2 - 2.5*row_size
set lmargin at screen left_margin + 2.3*col_size
set rmargin at screen left_margin2 + 3.3*col_size
set label 1 "D" @POS2
#set label 2 "ITS1" at graph 0.02,0.9
@XTICS; @YTICS
@XRANG_ITS1S
@YRANG
set x2tics ('B1' 0, 'B1L' -6) offset 0,-0.5 font 'arial, 6'
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
set key center top
plot \
DIR[8]."b6_27SB_5EXT.dat" using (-$1):($3) t "27SB" w lines lc 1,  \
DIR[7]."b6_7S_5EXT.dat" using (-$1):($3) t "7S" w lines lc 3, \
DIR[7]."b6_58S_5EXT.dat" using (-$1):($3) t "5.8S" w lines lc 4

set key default

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
set tmargin at screen top_margin - 2.7*row_size
set bmargin at screen top_margin2 - 3.2*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 2*col_size
set label 1 "E" @POS2
#set label 2 "ITS2" at graph 0.02,0.9
dat1 = DIR[7]."b6_7S_3EXT.dat"; dat2 = DIR[7]."b6_58S_3EXT.dat";  @TOT2D; total1 = total
@XTICS; @YTICS
@XRANG_ITS2L
@YRANGLOG
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
set x2tics ('E' 0, '+7' 7, '+31' 31, 'C2' 139) offset 0,-0.5 font 'arial, 6'
set key center top
plot \
DIR[7]."b6_7S_3EXT.dat" every ::5 using (139+$1):($2/total1) t "7S" w lines lt 1 lc 1, \
DIR[7]."b6_58S_3EXT.dat" using ($1):($2/total1) notitle w lines lt 1 lc 1

## Graph 5,2
set tmargin at screen top_margin - 2.7*row_size
set bmargin at screen top_margin2 - 3.2*row_size
set lmargin at screen left_margin + 2.3*col_size
set rmargin at screen left_margin2 + 4.3*col_size
set label 1 "F" @POS2
@XTICS; @YTICS
dat1 = DIR[6]."b6_26S_5EXT.dat"; dat2 = DIR[6]."b6_25S_5EXT.dat"; @TOT2D; total1 = total
@XRANG_ITS2R
@YRANGLOG
set x2tics ('' 0, '-2' -2, '-7' -7, 'C2' -95) offset 0,-0.5 font 'arial, 6'
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
set label 13 "C1" at 0,2.5 left font 'arial, 6'
set key center top
plot \
DIR[6]."b6_26S_5EXT.dat" using (-95-$1):($2/total1) t "26S" w lines lt 1 lc 1, \
DIR[6]."b6_25S_5EXT.dat" using (-$1):($2/total1) notitle w lines lt 1 lc 1
set key default
unset label 13

## Graph 5,3
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 2*col_size
set rmargin at screen left_margin2 + 3*col_size

## Graph 5,4
set tmargin at screen top_margin - 4*row_size
set bmargin at screen top_margin2 - 5*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 6,1
set tmargin at screen top_margin - 3.4*row_size
set bmargin at screen top_margin2 - 3.9*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size
set label 1 "G" @POS2
@XTICS; @YTICS
@XRANG_3ETS
@YRANG
set ylabel 'Fraction' font 'Arial, 7' offset 0,0; set xlabel 'Relative position' font 'Arial, 7' offset 0,0.2
set x2tics ('B2' 0, 'B0' 15) offset 0,-0.5 font 'arial, 6'

plot \
DIR[8]."b6_27SA2_3EXT.dat"  using ($1):($3) t "27SA2" w lines dt 1, \
DIR[2]."b6_33S_3EXT.dat"  using ($1):($3) t "33S" w lines dt 1, \
DIR[1]."b6_35S_3EXT.dat" using ($1):($3) t "35S" w lines dt 1, \
DIR[2]."b6_32S_3EXT.dat" using ($1):($3) t "32S" w lines dt 1,\
DIR[8]."b6_27SA3_3EXT.dat"  using ($1):($3) t "27SA3" w lines dt 1, \
DIR[8]."b6_27SB_3EXT.dat" using ($1):($3) t "27SB" w lines dt 1, \
DIR[6]."b6_26S_3EXT.dat" using ($1):($3) t "26S" w lines, \
DIR[6]."b6_25S_3EXT.dat" using ($1):($3) t "25S" w lines
set key default
unset label 2

## Graph 6,2
set tmargin at screen top_margin - 3.4*row_size
set bmargin at screen top_margin2 - 3.9*row_size
set lmargin at screen left_margin + 1.2*col_size
set rmargin at screen left_margin2 + 2.2*col_size


## Graph 6,3
set tmargin at screen top_margin - 3.4*row_size
set bmargin at screen top_margin2 - 3.9*row_size
set lmargin at screen left_margin + 1.4*col_size
set rmargin at screen left_margin2 + 2.4*col_size
set label 1 "H" @POS2
@NOXTICS; @YTICS
@XRANG_HIS
@YRANG_HIS
unset x2tics
set xtics out ("35S" 0, "33S" 1,"32S" 2,"23S" 3,"22S" 4,"21S" 5,"20S" 6,"27SA2" 7,"27SA3" 8,"27SB" 9,"7S" 10, "5.8S" 11, "18S" 12, "25S" 13) font 'arial, 6'
set ylabel 'Fraction' font 'Arial, 7' offset 0,0

plot \
DIR[1]."b6_35S_PAB.dat" using 3  with histograms, \
DIR[2]."b6_33S_PAB.dat" using 3  with histograms, \
DIR[2]."b6_32S_PAB.dat" using 3  with histograms, \
DIR[3]."b6_23S_PAB.dat" using 3  with histograms, \
DIR[4]."b6_22S_PAB.dat" using 3  with histograms, \
DIR[4]."b6_21S_PAB.dat" using 3  with histograms, \
DIR[4]."b6_20S_PAB.dat" using 3  with histograms, \
DIR[8]."b6_27SA2_PAB.dat" using 3  with histograms, \
DIR[8]."b6_27SA3_PAB.dat" using 3  with histograms, \
DIR[8]."b6_27SB_PAB.dat" using 3  with histograms, \
DIR[7]."b6_7S_PAB.dat" using 3  with histograms, \
DIR[7]."b6_58S_PAB.dat" using 3  with histograms, \
DIR[5]."b6_18S_PAB.dat" using 3  with histograms, \
DIR[6]."b6_25S_PAB.dat" using 3  with histograms, \

## Graph 6,4
set tmargin at screen top_margin - 5*row_size
set bmargin at screen top_margin2 - 6*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 7,1
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size

## Graph 7,2
set tmargin at screen top_margin - 4.1*row_size
set bmargin at screen top_margin2 - 4.6*row_size
set lmargin at screen left_margin + 2.4*col_size
set rmargin at screen left_margin2 + 3.4*col_size

## Graph 7,3
set tmargin at screen top_margin - 4.1*row_size
set bmargin at screen top_margin2 - 4.6*row_size
set lmargin at screen left_margin + 2.4*col_size
set rmargin at screen left_margin2 + 3.4*col_size

## Graph 7,4
set tmargin at screen top_margin - 6*row_size
set bmargin at screen top_margin2 - 7*row_size
set lmargin at screen left_margin + 3*col_size
set rmargin at screen left_margin2 + 4*col_size

## Graph 8,1
set tmargin at screen top_margin - 7*row_size
set bmargin at screen top_margin2 - 8*row_size
set lmargin at screen left_margin + 0*col_size
set rmargin at screen left_margin2 + 1*col_size

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


