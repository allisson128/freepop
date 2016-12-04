set xlabel "Iterações"
set xlabel font ",15" offset 0,0
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
set arrow from 0,500 to 4182,500 nohead lc rgb '#696969' dashtype 2
set arrow from 2935,0 to 2935,500 nohead lc rgb '#696969' dashtype 2
set arrow from 3211,0 to 3211,500 nohead lc rgb '#696969' dashtype 2
set arrow from 3617,0 to 3617,500 nohead lc rgb '#696969' dashtype 2
set arrow from 4182,0 to 4182,500 nohead lc rgb '#696969' dashtype 2
#set label "2935" at 2550,-50 font ",10"
#set label "3211" at 2951,-50 font ",10"
#set label "3617" at 3317,-50 font ",10"
#set label "4182" at 3942,-50 font ",10"



plot "100 formigas" with lines, "90 formigas" with lines, "80 formigas" with lines, "70 formigas" with lines

