#!/bin/bash

set -e
OUTPUT=$1

mkdir  -p $OUTPUT

touch $OUTPUT/status

for i in {18,20,22,26,28,30,34,36,38}
do
    ./scalingPlot.pl titan${i}0 data/AMRAdvection$i* >> $OUTPUT/status
done

for i in {102,103,104}
do
    echo data/AMRAdvection${i}{0..7}*
    ./scalingPlot.pl titan${i}0 data/AMRAdvection${i}{0..7}* >> $OUTPUT/status
done

for i in {18,20,22,26,28,30,34,36,102,103,104} #38,
do
    mv *titan${i}0 $OUTPUT
done

./scalingPlot.pl bgq9 data/strong.9/*output   >> $OUTPUT/status
./scalingPlot.pl bgq10 data/strong.10/*output >> $OUTPUT/status
./scalingPlot.pl bgq11 data/strong.11/*output >> $OUTPUT/status

./scalingPlot.pl bgq12 data/mid.512.9/*output >> $OUTPUT/status
./scalingPlot.pl bgq13 data/mid.512.10/*output >> $OUTPUT/status
./scalingPlot.pl bgq14 data/mid.512.11/*output >> $OUTPUT/status

./scalingPlot.pl bgq15 data/mid.1024.9/*output >> $OUTPUT/status
./scalingPlot.pl bgq16 data/mid.1024.10/*output >> $OUTPUT/status
./scalingPlot.pl bgq17 data/mid.1024.11/*output >> $OUTPUT/status

./scalingPlot.pl bgq18 data/remesh.9/*output >> $OUTPUT/status
./scalingPlot.pl bgq19 data/remesh.10/*output >> $OUTPUT/status
./scalingPlot.pl bgq20 data/remesh.11/*output >> $OUTPUT/status

mv *bgq{9..20} $OUTPUT

# box and whisker plots for TD
./qdCandlestick.pl $OUTPUT QDcandle.titan{180,200,220} titan
gnuplot temp.gnuplot > $OUTPUT/qdCandlestickTitan.ps
#./qdCandlestick.pl $OUTPUT QDcandle.bgq{9,10,11}
./qdCandlestick.pl $OUTPUT QDcandle.bgq{18,19,20} bgq
gnuplot temp.gnuplot > $OUTPUT/qdCandlestickBGQ.ps

# strong scaling (ts/sec) paired by machine
./tsPerSec.pl $OUTPUT timestepPerSec.titan180 timestepPerSec.bgq9 timestepPerSec.titan260 timestepPerSec.bgq12 timestepPerSec.titan340 timestepPerSec.bgq15
gnuplot temp.gnuplot > $OUTPUT/timestepPerSec9.ps
./tsPerSec.pl $OUTPUT timestepPerSec.titan200 timestepPerSec.bgq10 timestepPerSec.titan280 timestepPerSec.bgq13 timestepPerSec.titan360 timestepPerSec.bgq16
gnuplot temp.gnuplot > $OUTPUT/timestepPerSec10.ps
./tsPerSec.pl $OUTPUT timestepPerSec.titan220 timestepPerSec.bgq11 timestepPerSec.titan300 timestepPerSec.bgq14 timestepPerSec.titan380 timestepPerSec.bgq17
gnuplot temp.gnuplot > $OUTPUT/timestepPerSec11.ps

# histograms for remeshing
# ./remeshHistogram.pl $OUTPUT histoRemesh.32.titan1020
# gnuplot temp.gnuplot > $OUTPUT/histoRemesh32-1020.ps
# ./remeshHistogram.pl $OUTPUT histoRemesh.256.titan1020
# gnuplot temp.gnuplot > $OUTPUT/histoRemesh256-1020.ps

# box and whisker plots for TD, (titan is ordered a bit unintuitively: depth 9,11,10)
./remeshLatency.pl $OUTPUT RMcandle.titan{1020,1040,1030} titan
gnuplot temp.gnuplot > $OUTPUT/remeshLatencyTitan.ps
./remeshLatency.pl $OUTPUT RMcandle.bgq{18,19,20} bgq
gnuplot temp.gnuplot > $OUTPUT/remeshLatencyBGQ.ps

./remeshMedian.pl $OUTPUT RMcandle.titan{1020,1040,1030} titan 1100
gnuplot temp.gnuplot > $OUTPUT/remeshMedianTitan.ps
./remeshMedian.pl $OUTPUT RMcandle.bgq{18,19,20} bgq 2400
gnuplot temp.gnuplot > $OUTPUT/remeshMedianBGQ.ps


# step time plots
./stepTime.pl $OUTPUT stepTime titan180 left 150
gnuplot temp.gnuplot > $OUTPUT/stepTimeTitanDepth9.ps
./stepTime.pl $OUTPUT stepTime titan200 left 150
gnuplot temp.gnuplot > $OUTPUT/stepTimeTitanDepth10.ps
./stepTime.pl $OUTPUT stepTime titan220 right 150
gnuplot temp.gnuplot > $OUTPUT/stepTimeTitanDepth11.ps

./stepTime.pl $OUTPUT stepTime bgq9 left 800
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth9.ps
./stepTime.pl $OUTPUT stepTime bgq10 left 800
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth10.ps
./stepTime.pl $OUTPUT stepTime bgq11 right 800
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth11.ps

./stepTime.pl $OUTPUT stepTime bgq12 left 1600
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth12.ps
./stepTime.pl $OUTPUT stepTime bgq13 left 1600
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth13.ps
./stepTime.pl $OUTPUT stepTime bgq14 right 1600
gnuplot temp.gnuplot > $OUTPUT/stepTimeBGQDepth14.ps


rm temp.gnuplot

for i in $OUTPUT/*ps
do
    echo $i
    ps2pdf $i $i.pdf
    BNAME=`basename $i .ps`
    pdfcrop $i.pdf $OUTPUT/$BNAME.pdf
    rm $i.pdf
done