set datafile separator ","
set terminal pdfcairo enhanced font "Helvetica,12"
set xlabel "mass"
set ylabel "delta mass(MeV) at 204 MeV"
set title "Symmetry Mass Differences at T=204 MeV"
set xrange [0:0.014]
set grid
set output "deltamassvsm.pdf"

set zeroaxis linetype 1 linecolor "black"
set key right top
plot \
"deltamass_nt12.csv" using (($1==1)?$2:1/0):($3*2453):($4*2453) with yerrorbars pt 6 lc rgb "#69b3a2" title "V-A (SU(2)_L × SU(2)_R)", \
"deltamass_nt12.csv" using (($1==2)?$2:1/0):($3*2453):($4*2453) with yerrorbars pt 8 lc rgb "#64b5f6" title "T-X (U(1)_A)", \
"deltamass_nt12.csv" using (($1==3)?$2:1/0):($3*2453):($4*2453) with yerrorbars pt 4 lc rgb "#ba68c8" title "P-S (U(1)_A)", \
"deltamass_nt12.csv" using (($1==4)?$2:1/0):($3*2453):($4*2453) with yerrorbars pt 10 lc rgb "#ffb74d" title "A-X (SU(2)_{CS})"
