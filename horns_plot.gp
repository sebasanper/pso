unset key
set yrange [-10:3600]
set xrange [-10:400]
do for [i=0:1499]{
plot '3horns_pso.dat' u ($2):($3):1 every :::i::i pt 7 w labels textcolor rgb 'blue', 3600
pause 0.5
}
