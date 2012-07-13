#!/opt/local/bin/perl

use strict;
use warnings;

my ($output, $d1,$d2) = @ARGV;

my $out = <<OUT;
set encoding utf8
set terminal post color font "Helvetica,16"
set border 1+2+4+8 lw 1.3

set pointsize 1

set style line 1 pt 1 lt 2 lc rgb "red" lw 4
set style line 2 pt 3 lt 5 lc rgb "#ff8800" lw 4
set style line 3 pt 4 lt 8 lc rgb "blue" lw 4
set style line 4 pt 6 lt 3 lc rgb "#006400" lw 4
# set style line 4 pt 6 lt 4 lc rgb "#ff8800" lw 4
# set style line 5 pt 8 lt 6 lc rgb "#b03060" lw 4
# set style line 6 pt 2 lt 5 lc rgb "#8299AA" lw 4
# set style line 7 pt 5 lt 7 lc rgb "#000000" lw 4
# set style line 8 pt 9 lt 8 lc rgb "brown" lw 4
set style line 10 pt 0 lw 2 lc rgb "black"

set key top left
#set grid mytics ytics

set logscale xy

set xlabel 'Number of Cores'
set ylabel 'Timesteps / sec'
set xtics autofreq nomirror (16,32,64,128,256,512,1024,2048,4096,8192,16384,32768)
set ytics autofreq 10
set mytics 5
set yrange [1:10000]
set xrange [8:40000]

plot '$output/$d1' using 1:2 with linespoints ls 1 title "Cray XK6",\\
     '$output/$d1' using 1:3 with linespoints ls 10 title "",\\
     '$output/$d2' using 1:2 with linespoints ls 2 title "IBM BG/Q",\\
     '$output/$d2' using 1:3 with linespoints ls 10 title ""
OUT

open FILE, ">", "temp.gnuplot";
print FILE $out;
close FILE;
