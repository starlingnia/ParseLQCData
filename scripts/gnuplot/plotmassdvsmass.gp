# 来源: /Users/junxiongnie/code/ana/dat/ccor/plotmassdvsmass.gp (原样复制, 仅把硬编码输入文件名参数化)
set datafile separator ","
set terminal pdfcairo enhanced font "Helvetica,12"
set xlabel "mass"
set ylabel "delta mass(MeV) at 220 MeV"
set title "Symmetry Mass Differences at T=220MeV"
set xrange [0:0.008]
set grid
set output "deltamassvsm.pdf"

set zeroaxis linetype 1 linecolor "black"
plot \
input_fname using (($1==1)?$2:1/0):($3*2640):($4*2640) with yerrorbars pt 6 lc rgb "#69b3a2" title "V-A", \
input_fname using (($1==2)?$2:1/0):($3*2640):($4*2640) with yerrorbars pt 8 lc rgb "#64b5f6" title "T-X", \
input_fname using (($1==3)?$2:1/0):($3*2640):($4*2640) with yerrorbars pt 4 lc rgb "#ba68c8" title "P-S", \
input_fname using (($1==4)?$2:1/0):($3*2640):($4*2640) with yerrorbars pt 10 lc rgb "#ffb74d" title "A-X"



