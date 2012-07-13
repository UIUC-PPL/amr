set terminal post color font "Helvetica,24"
set border 1+2+4+8 lw 1.3

set style line 1 pt 1 lt 1 lc rgb "red" lw 4
set style line 2 pt 3 lt 1 lc rgb "#ff8800" lw 4
set style line 3 pt 4 lt 1 lc rgb "blue" lw 4
set style line 4 pt 6 lt 3 lc rgb "#006400" lw 4

set logscale x
set yrange [-100:5000]
set xrange [10:3000]

set xtics autofreq nomirror (16,32,64,128,256,512,1024,2048)

set bars 4.0
set style fill empty
f(x)=x+4
plot 'QDcandle.320' using 9:5:4:8:7 with candlesticks ls 1 lw 2 whiskerbars,\
     'QDcandle.320' using 9:6:6:6:6 with candlesticks ls 1 lw 2 lt -1 notitle,\
     'QDcandle.340' using 9:5:4:8:7 with candlesticks ls 2 lw 2 whiskerbars,\
     'QDcandle.340' using 9:6:6:6:6 with candlesticks ls 2 lw 2 lt -1 notitle,\
     'QDcandle.300' using 9:5:4:8:7 with candlesticks ls 3 lw 2 whiskerbars,\
     'QDcandle.300' using 9:6:6:6:6 with candlesticks ls 3 lw 2 lt -1 notitle