set xlabel "Iterações"
set xlabel font ",15" offset 0,-1.5
set ylabel "Ocorrências do menor caminho"
set ylabel font ",15"
set ylabel offset -2,0

set xtics font ",12"
set xtics offset -1,0
set ytics font ",12"
#unset xtics

set key box left
set key font ",11"

set grid #linewidth 1
set arrow 1 from 6675,1526 to 10000,1526 nohead lc rgb '#696969' dashtype 2
set arrow 2 from 6675,-150 to 6675,1526 nohead lc rgb '#696969' dashtype 2
set arrow 3 from 7460,-240 to 7460,1526 nohead lc rgb '#696969' dashtype 2
set arrow 4 from 8320,-310 to 8320,1526 nohead lc rgb '#696969' dashtype 2
set label "1526" at 10050,1526 font ",10"
set label "6675" at 6325,-190 font ",10"
set label "7460" at 7110,-280 font ",10"
set label "8320" at 7980,-350 font ",10"

plot "100 formigas" with lines, "90 formigas" with lines, "80 formigas" with lines, "70 formigas" with lines

