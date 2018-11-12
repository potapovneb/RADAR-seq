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

if( @ARGV == 0 )
{
    print "usage: $0 [options] coverage.csv patches.csv\n\n";
    exit;
}

my $covfile = shift @ARGV;
my $patfile = shift @ARGV;

my %data = ();

### load coverage data
my @head = ();

open(COVERAGE,$covfile) || die "Can't open '$covfile'";

while( my $line = <COVERAGE> )
{
    chomp($line);
    
    my @tokens = split(/,/,$line,-1);
    
    if( @head == 0 )
    {
	@head = @tokens;
    }
    else
    {
	my %entry = ();

	for( my $i = 0; $i < @tokens; $i++ )
	{
	    $entry{$head[$i]} = $tokens[$i];
	}

	my $rname = $entry{"Reference"};

 	$data{$rname}{0}{"bases"} += $entry{"N0"};
 	$data{$rname}{1}{"bases"} += $entry{"N1"};
    }
}

close(COVERAGE);

### init patch counts
for my $rname ( keys %data )
{
    for my $strand ( keys %{$data{$rname}} )
    {
	$data{$rname}{$strand}{"patches"} = 0;
    }
}

### load patch data
@head = ();

open(CSV,$patfile) || die "Can't open '$patfile'";

while( my $line = <CSV> )
{
    chomp($line);
    
    # my @tokens = quotewords(',',0,$line);
    my @tokens = split(/,/,$line,-1);
    
    if( @head == 0 )
    {
	@head = @tokens;
    }
    else
    {
	my %entry = ();

	for( my $i = 0; $i < @tokens; $i++ )
	{
	    $entry{$head[$i]} = $tokens[$i];
	}

	my $rname = $entry{"Reference"};
	my $strand = $entry{"Strand"};

	if( $entry{"SVM_Prediction"} == 1 )
	{
	    $data{$rname}{$strand}{"patches"}++;
	}
    }
}

close(CSV);

my $bases_total = 0;
my $patches_total = 0;

print join(",","Reference", "Strand", "Bases", "Patches", "Frequency"), "\n";

for my $rname ( sort keys %data )
{
    my $bases = 0;
    my $patches = 0;

    for my $strand ( sort { $a <=> $b } keys %{$data{$rname}} )
    {
	printf( "%s,%i,%i,%i,%f\n",
		$rname,
		$strand,
		$data{$rname}{$strand}{"bases"},
		$data{$rname}{$strand}{"patches"},
		$data{$rname}{$strand}{"patches"} / $data{$rname}{$strand}{"bases"},
	    );

	$bases += $data{$rname}{$strand}{"bases"};
	$patches += $data{$rname}{$strand}{"patches"};
    }

    printf( "%s,%s,%i,%i,%f\n",
	    $rname,
	    "*",
	    $bases,
	    $patches,
	    $patches / $bases,
	);

    print ",\n";

    $bases_total += $bases;
    $patches_total += $patches;
}

printf( "%s,%s,%i,%i,%f\n",
	"*",
	"*",
	$bases_total,
	$patches_total,
	$patches_total / $bases_total,
    );
