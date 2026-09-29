# 来源: /Users/junxiongnie/code/ana/dat/ccor/plotmd.gp (原样复制, 仅把硬编码输入文件名参数化)
set terminal pdfcairo enhanced font "Helvetica,12"
set output output_fname
set title symt
set xlabel "temperature(MeV) "
set ylabel ytle
set xrange [140:225]
set grid
# 1. 设置一个矩形对象
#    'from 155, graph 0'  : 起始点 x=155, y=图表底部
#    'to 158, graph 1'    : 结束点 x=158, y=图表顶部
set object 1 rectangle from 155, graph 0 to 158, graph 1  fillstyle solid noborder fillcolor "gray" 



# 3. (重要) 把这个矩形放到 "后面"，防止它遮住你的数据线
set object 1 behind
set zeroaxis linetype 1 linecolor "black"
plot \
input_fname using ($2==7?$1:1/0):($2==7?$3:1/0):($2==7?$4:1/0) with yerrorbars pt 10 lc rgb "#ffb74d" title "mass 0.0120", \
input_fname using ($2==4?$1:1/0):($2==4?$3:1/0):($2==4?$4:1/0) with yerrorbars pt 8 lc rgb "#64b5f6" title "mass 0.0070", \
input_fname using ($2==2?$1:1/0):($2==2?$3:1/0):($2==2?$4:1/0) with yerrorbars pt 4 lc rgb "#ba68c8" title "mass 0.0035", \
input_fname using ($2==1?$1:1/0):($2==1?$3:1/0):($2==1?$4:1/0) with yerrorbars pt 6 lc rgb "#69b3a2" title "mass 0.0020"