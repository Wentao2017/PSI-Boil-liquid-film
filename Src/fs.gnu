set xlabel "Time (s)"
set ylabel "Wave elevation (m)"
set gri
plot "fs.out" u 2:4 t "Elevation #1" w l,"" u 2:6 t "Elevation #2" w l,"" u 2:8 t "Elevation #3" w l, "" u 2:10 t "Elevation #4" w l
