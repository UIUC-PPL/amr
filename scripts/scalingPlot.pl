#!/usr/bin/perl

use strict;
use warnings;

my @files = @ARGV;
my %scaling;
my $bestThroughput = 0.0;

for my $file (@files) {
    open FILE, "<", $file;

    my $proc = -1;
    my $time = 0.0;

    for (<FILE>) {
        if (/(\d+) processors/) {
            $proc = $1;
        }
        if (/simulation time: ([0-9.]+)/) {
            $time = $1;
        }
    }

    if ($bestThroughput == 0.0 or
        $bestThroughput > $time * $proc) {
        $bestThroughput = $time * $proc;
    }

    $scaling{$proc} = $time;

    print STDERR "$file: $proc cores, duration $time\n";
    #print "$proc $time\n";

    close FILE;
}

my $perfect = 0.0;

for my $key (sort {$a <=> $b} (keys %scaling)) {
    next if ($key == -1);
    my $val = $scaling{$key};
    if ($perfect == 0.0) {$perfect = $val * $key;}
    my $ideal = $perfect / $key;
    print "$key $val $ideal\n";
}
