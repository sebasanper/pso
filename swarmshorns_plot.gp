unset key
set yrange [-10:4000]
set xrange [-10:5500]
do for [i=0:99]{
plot 'ref1_random_layout_ainslie.dat' u ($2):($3):1 every :::i::i pt 7 w labels textcolor rgb 'blue',-3907./412. * x + 3907., 3907./417.*(5457. - x), 3907., 0
pause -1
}
