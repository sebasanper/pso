unset key
set xrange [-1000:7000]
set yrange [-1000:5000]
do for [i=0:30]{
#plot 'r5_layout_jensen.dat' u 2:3:1 every :::i::i w labels, 3907, 0, -3907/412*x+3907, 3907/417*(5457-x)
plot 'ref1_random_layout_ainslie.dat' u 2:3:1 every :::i::i w labels, 3907, 0, -3907/412*x+3907, 3907/417*(5457-x)
pause -1
}

# pt 7 lc rgb 'red'
