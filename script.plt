a = 0.003703
b = -1.288406
c = 176.578090
f(x) = a*x*x + b*x + c
h = 185
set term png
set output "output.png"
set xrange[0:200]
set label at h, f(h) "" point pointtype 7 pointsize 2
plot "data.dat", f(x)
