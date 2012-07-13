#!/bin/bash

set -e
OUTPUT=$1

mkdir  -p $OUTPUT

touch $OUTPUT/status
./scalingPlot.pl titan180 data/AMRAdvection18* >> $OUTPUT/status
./scalingPlot.pl titan200 data/AMRAdvection20* >> $OUTPUT/status
./scalingPlot.pl titan220 data/AMRAdvection22* >> $OUTPUT/status
mv *titan{180,200,220} $OUTPUT

./scalingPlot.pl bgq9 data/strong.9/*output   >> $OUTPUT/status
./scalingPlot.pl bgq10 data/strong.10/*output >> $OUTPUT/status
./scalingPlot.pl bgq11 data/strong.11/*output >> $OUTPUT/status
mv *bgq{9,10,11} $OUTPUT

./qdCandlestick.pl $OUTPUT QDcandle.titan{180,200,220}
gnuplot temp.gnuplot > $OUTPUT/qdCandlestickTitan.ps
./qdCandlestick.pl $OUTPUT QDcandle.bgq{9,10,11}
gnuplot temp.gnuplot > $OUTPUT/qdCandlestickBGQ.ps

rm temp.gnuplot

for i in $OUTPUT/*ps
do
    echo $i
    ps2pdf $i $i.pdf
    BNAME=`basename $i .ps`
    pdfcrop $i.pdf $OUTPUT/$BNAME.pdf
    rm $i.pdf
done