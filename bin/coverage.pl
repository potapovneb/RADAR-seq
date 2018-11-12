#!/usr/bin/perl -w

################################################################################
# RADAR-seq: detecting and locating DNA damage on a genome-wide scale
# Copyright (C) 2018 New England Biolabs, Inc.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
################################################################################

use strict;
use Getopt::Long;
use Text::ParseWords;

### command-line options
my $o_edist  = 30;
my $o_ndist  = 30;
my $o_nsites =  5;
my @o_bases  = ();
my %o_bases  = ();
my $o_nsfile = "";
my $o_msfile = "";

GetOptions(
    "bases=s"  => \@o_bases,
    "edist=f"  => \$o_edist,
    "ndist=f"  => \$o_ndist,
    "nsites=f" => \$o_nsites,
    "nsfile=s" => \$o_nsfile,
    "msfile=s" => \$o_msfile,
    );

for my $x ( split(/;/, join(";", @o_bases)) ) {$o_bases{$x} = 1;}

### command-line usage
if ( @ARGV == 0 )
{
    print "usage: $0 [options] reference.fasta methylation_sites.csv modifications.csv.gz\n\n";
    print "options:\n";
    print "  --bases\tA or C or A,C\n";
    print "  --edist\textension distance\n";
    print "  --ndist\tnicking site distance\n";
    print "  --nsites\tmodifiable sites tolerated in extension\n";
    print "  --nsfile\tnicking site file\n";
    print "  --msfile\tmethylation site file\n";
    exit;
}

my @modfiles = @ARGV;

my %cov = ();

for my $modfile ( @modfiles )
{
    print STDERR "  Reading '$modfile'...\n";

    open(IPD, "gzip -cd $modfile |") || die "Can't open '$modfile'";

    while ( my $line = <IPD> )
    {
	### skip CSV header line
	if ( substr($line, 0, 7) eq "refName" )
	{
	    next;
	}

	chomp($line);

	my ($rname,
	    $pos,
	    $strand,
	    $base,
	    $score,
	    $tMean,
	    $tErr,
	    $modelPrediction,
	    $ipdRatio,
	    $coverage,
	    $frac,
	    $fracLow,
	    $fracUpp) = split(/,/, $line);

	$rname =~ s/^\"|\"$//g;

	$cov{$rname}{$pos}{$strand}++;
    }
}

print join(",", qw/Reference Position N0 N1/), "\n";

for my $rname ( keys %cov )
{
    for my $pos ( sort { $a <=> $b } keys %{$cov{$rname}} )
    {
	my @values = ( $rname, $pos );

	for my $strand ( 0, 1 )
	{
	    if ( exists $cov{$rname}{$pos}{$strand} )
	    {
		push @values, $cov{$rname}{$pos}{$strand};
	    }
	    else
	    {
		push @values, 0;
	    }
	}

	print join(",", @values), "\n";
    }
}
