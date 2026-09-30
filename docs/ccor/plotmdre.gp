set datafile separator ","
set terminal pdfcairo enhanced font "Helvetica,12"
set xlabel "temperature(MeV)"
set ylabel "mass(MeV)"
set xrange [125:225]
set yrange [0:1350]
set grid
# 1. 设置一个矩形对象
#    'from 155, graph 0'  : 起始点 x=155, y=图表底部
#    'to 158, graph 1'    : 结束点 x=158, y=图表顶部
set object 1 rectangle from 155, graph 0 to 158, graph 1 fillstyle solid noborder fillcolor "gray"



# 3. (重要) 把这个矩形放到 "后面"，防止它遮住你的数据线
set object 1 behind
set pointsize 0.6
set key right bottom
do for [mass in "0.002 0.0035 0.007 0.0120"] {
    set output sprintf("massvtem_%.4f.pdf", real(mass))
    set title sprintf("mass vs. temperature (mass=%.4f)", real(mass))
    plot \
        input_fname using (($2==1 && abs($3-mass)<1e-6)?$1-1.75:1/0):(($2==1 && abs($3-mass)<1e-6)?$4:1/0):(($2==1 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 11 lc rgb "#9467bd" title "V", \
        input_fname using (($2==2 && abs($3-mass)<1e-6)?$1-1.25:1/0):(($2==2 && abs($3-mass)<1e-6)?$4:1/0):(($2==2 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 7 lc rgb "#1f77b4" title "A", \
        input_fname using (($2==3 && abs($3-mass)<1e-6)?$1+0.25:1/0):(($2==3 && abs($3-mass)<1e-6)?$4:1/0):(($2==3 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 13 lc rgb "#d62728" title "Tt", \
        input_fname using (($2==4 && abs($3-mass)<1e-6)?$1+1.25:1/0):(($2==4 && abs($3-mass)<1e-6)?$4:1/0):(($2==4 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 2 lc rgb "#8c564b" title "Xt", \
        input_fname using (($2==5 && abs($3-mass)<1e-6)?$1-0.25:1/0):(($2==5 && abs($3-mass)<1e-6)?$4:1/0):(($2==5 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 9 lc rgb "#2ca02c" title "S", \
        input_fname using (($2==6 && abs($3-mass)<1e-6)?$1-0.75:1/0):(($2==6 && abs($3-mass)<1e-6)?$4:1/0):(($2==6 && abs($3-mass)<1e-6)?$5:1/0) with yerrorbars pt 5 lc rgb "#ff7f0e" title "Ps",\
        2*pi*x
}
