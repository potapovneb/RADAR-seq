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

my $o_prefix = "test";

GetOptions( "prefix=s" => \$o_prefix );

### command-line arguments
if ( @ARGV == 0 )
{
    print "usage: $0 [options] input.bam\n\n";
    exit;
}

my $bamfile = shift @ARGV;

# data structures for parsing CIGAR string
my %digit = (
    "0" => 1,
    "1" => 1,
    "2" => 1,
    "3" => 1,
    "4" => 1,
    "5" => 1,
    "6" => 1,
    "7" => 1,
    "8" => 1,
    "9" => 1 );

my %operator = (
    "M" => 1,
    "I" => 1,
    "D" => 1,
    "=" => 1,
    "X" => 1 );

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### extract mapping information
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

open(T1, ">", "${o_prefix}.mapping.csv") || die "Can't open '${o_prefix}.mapping.csv'";

print T1 "Movie,ZMW,Coord,Reference,Flag,MAPQ,AlnStart,AlnEnd\n";

my %data = ();

open(BAM, "samtools view $bamfile |") || die "Can't open '$bamfile'";

while ( my $line = <BAM> )
{
    my $s1 = ""; # aligned reference
    my $s2 = ""; # aligned read
    my $bq = ""; # base quality

    ### parse fields
    my ($qname,
	$flag,
	$rname,
	$pos,
	$mapq,
	$cigar,
	$rnext,
	$pnext,
	$tlen,
	$read,
	$qual,
	@tag) = split(/\t/,$line);

    ### fwd/rev strand
    my $strand = ( $flag & 0x10 ) ? 1 : 0;

    ### skip unmapped reads
    next if ( $flag & 0x4 );

    ### skip secondary alignment
    next if ( $flag & 0x100 );

    ### skip supplementary alignment
    next if ( $flag & 0x800 );

    # ----- process CIGAR ------------------------------------------------------

    my $p1  = $pos - 1; # start of aligned reference
    my $p2  = 0;        # start of aligned read
    my $num = 0;        # length of the block to operate
    my $len = 0;        # read length
    
    my @cigar = split(//,$cigar);
    
    while ( @cigar )
    {
	my $c = shift @cigar;
	
	if ( exists $digit{$c} )
	{
	    $num .= $c;
	    next;
	}
	
	if ( $c eq "M" || $c eq "=" || $c eq "X" )
	{
	    $len += $num;
	}
	elsif ( $c eq "D" )
	{
	    $len += $num;
	}

	$num = "";
    }

    my ($movie, $zmw, @oth) = split(/\//, $qname);

    print T1 join(",",
		  $movie,
		  $zmw,
		  ( @oth ? $oth[0] : "NA" ),
		  "\"$rname\"",
		  $flag,
		  $mapq,
		  $pos,
		  $pos - 1 + $len,
	), "\n";
    
    $qname = "$movie/$zmw";
    
    if( ! exists $data{$rname}{$qname}{"AlnStart"} || $data{$rname}{$qname}{"AlnStart"} > $pos )
    {
	$data{$rname}{$qname}{"AlnStart"} = $pos;
    }
    
    if( ! exists $data{$rname}{$qname}{"AlnEnd"} || $data{$rname}{$qname}{"AlnEnd"} < ($pos - 1 + $len) )
    {
	$data{$rname}{$qname}{"AlnEnd"} = ($pos - 1 + $len);
    }
}

close(T1);


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### arrange non-overlapping ZMW reads into layers
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

open(T2, ">", "${o_prefix}.strata.csv") || die "Can't open '${o_prefix}.strata.csv'";

print T2 "Movie,ZMW,AlnStart,AlnEnd,Stratum\n";

my %stats = ();

for my $rname ( sort keys %data )
{
    ### sort ZMW reads by AlnStart
    my @sorted = sort { $data{$rname}{$a}{"AlnStart"} <=> $data{$rname}{$b}{"AlnStart"} } keys %{$data{$rname}};

    ### arrange reads into layers
    my @strata = ();

    for ( my $i = 0; $i < @sorted; $i++ )
    {
	my $qname = $sorted[$i];
	
	my @found = ();
	
	### find existing layers that can accomodate read
	for ( my $k = 0; $k < @strata; $k++ )
	{
	    if ( $data{$rname}{$qname}{"AlnStart"} > $strata[$k] + 500 )
	    {
		my $gap = $data{$rname}{$qname}{"AlnStart"} - $strata[$k];
		
		push @found, {"strata" => $k, "gap" => $gap};
	    }
	}
	
	### otherwise, add new layer
	if ( ! @found )
	{
	    my $n = scalar(@strata);
	    push @found, {"strata" => $n, "gap" => 0};
	}

	### find layer that will introduce minimum gap
	my @sorted_found = sort { $$a{"gap"} <=> $$b{"gap"} } @found;
	my $resulting_strata = $sorted_found[0]{"strata"};

	### update right-most position in the layer
	$strata[$resulting_strata] = $data{$rname}{$qname}{"AlnEnd"};
	
	### print assigned layer for each read
	my ($movie, $zmw) = split(/\//,$qname);
	
    	print T2 join(",",
		      $movie,
		      $zmw,
		      $data{$rname}{$qname}{"AlnStart"},
		      $data{$rname}{$qname}{"AlnEnd"},
		      $resulting_strata + 1,
    	    ), "\n";
	
	### update number of read per layer
	$stats{$resulting_strata + 1}++;
    }
}

close(T2);


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### number of read per layer
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

open(T3, ">", "${o_prefix}.stats.csv") || die "Can't open '${o_prefix}.stats.csv'";

print T3 "Stratum,Count\n";

for my $stratum ( sort { $a <=> $b } keys %stats )
{
    print T3 join(",", $stratum, $stats{$stratum}), "\n";
}

close(T3);
