set xlabel "Time (s)"
set ylabel "Void fraction in entire domain"
plot "vol.out" u 2:($3/($3+$4)) w l
