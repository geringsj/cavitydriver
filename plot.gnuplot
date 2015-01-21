set term wxt 0
set autoscale
set xlabel "time"
set ylabel "velocity"
set style fill transparent solid 0.2 noborder
plot \
	"plot.data" using 1:2 title 'E[U@(120,5)]' with lines, \
		"plot.data" using 1:($2+$8):($2-$8) title 'V[U@(120,5)]' with filledcurves, \
	"plot.data" using 1:3 title 'E[V@(120,5)]' with lines, \
		"plot.data" using 1:($3+$9):($3-$9) title 'V[V@(120,5)]' with filledcurves, \
	"plot.data" using 1:4 title 'E[U@(64,64)]' with lines, \
		"plot.data" using 1:($4+$10):($4-$10) title 'V[U@(64,64)]' with filledcurves, \
	"plot.data" using 1:5 title 'E[V@(64,64)]' with lines, \
		"plot.data" using 1:($5+$11):($5-$11) title 'V[V@(64,64)]' with filledcurves, \
	"plot.data" using 1:6 title 'E[U@(5,120)]' with lines, \
		"plot.data" using 1:($6+$12):($6-$12) title 'V[U@(5,120)]' with filledcurves, \
	"plot.data" using 1:7 title 'E[V@(5,120)]' with lines, \
		"plot.data" using 1:($7+$13):($7-$13) title 'V[V@(5,120)]' with filledcurves
set term png truecolor medium size 1920,1080
set output "all_in_one.png"
replot

set term wxt 1
set autoscale
set xlabel "time"
set ylabel "velocity"
set style fill transparent solid 0.2 noborder
plot \
	"plot.data" using 1:2 title 'E[U@(120,5)]' with lines, \
		"plot.data" using 1:($2+$8):($2-$8) title 'V[U@(120,5)]' with filledcurves, \
	"plot.data" using 1:3 title 'E[V@(120,5)]' with lines, \
		"plot.data" using 1:($3+$9):($3-$9) title 'V[V@(120,5)]' with filledcurves
set term png truecolor medium size 1920,1080
set output "UV_120_5.png"
replot

set term wxt 2
set autoscale
set xlabel "time"
set ylabel "velocity"
set style fill transparent solid 0.2 noborder
plot \
	"plot.data" using 1:4 title 'E[U@(64,64)]' with lines, \
		"plot.data" using 1:($4+$10):($4-$10) title 'V[U@(64,64)]' with filledcurves, \
	"plot.data" using 1:5 title 'E[V@(64,64)]' with lines, \
		"plot.data" using 1:($5+$11):($5-$11) title 'V[V@(64,64)]' with filledcurves
set term png truecolor medium size 1920,1080
set output "UV_64_64.png"
replot

set term wxt 3
set autoscale
set xlabel "time"
set ylabel "velocity"
set style fill transparent solid 0.2 noborder
plot \
	"plot.data" using 1:6 title 'E[U@(5,120)]' with lines, \
		"plot.data" using 1:($6+$12):($6-$12) title 'V[U@(5,120)]' with filledcurves, \
	"plot.data" using 1:7 title 'E[V@(5,120)]' with lines, \
		"plot.data" using 1:($7+$13):($7-$13) title 'V[V@(5,120)]' with filledcurves
set term png truecolor medium size 1920,1080
set output "UV_5_120.png"
replot


set term wxt 4
set autoscale
set xlabel "Re"
set ylabel "P[Re]"
Norm(x,m,s) = 1./(sqrt(2*pi)*s) * exp( -(x-m)**2 / (2*s*s) )
plot [500:2500] Norm(x,1500,(2000-1000)/6), \
	"reynolds_numbers.txt" using 1:(Norm($1,1500,(2000-1000)/6))
set term png truecolor medium size 1920,1080
set output "reynolds_numbers_samples.png"
replot

