set xlabel "Time (s)" 
set ylabel "Superficial velocity of liquid (m/s)" 
set y2label "Superficial velocity of gas (m/s)" 
set y2tic auto
plot "vel.out" u 2:3 t "liquid" w l,"vel.out"  u 2:4 axis x1y2 w l t "gas"
