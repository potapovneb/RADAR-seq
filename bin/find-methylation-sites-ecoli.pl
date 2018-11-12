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

if( @ARGV == 0 )
{
    print "usage: $0 genome.fasta\n\n";
    exit;
}

my $file = shift @ARGV;

### load reference sequence
my $refs = load_references($file);

print join(",",
	   "Reference",
	   "Position",
	   "Strand",
	   "Motif",
	   "Offset"), "\n";

for my $name ( keys %$refs )
{
    my $seq = uc($$refs{$name});
    my $n = length($seq);

    for( my $i = 0; $i < length($seq); $i++ )
    {
	if( $i + 4 < $n )
	{
	    my $motif = substr($seq,$i,4);

	    #    0123
	    #     *
	    # 5' GATC 3'
	    # 5' CTAG 3'
	    #      *

	    if( $motif eq "GATC" )
	    {
		print join(",", $name, ($i+1) + 1, 0, $motif, 1), "\n";
		print join(",", $name, ($i+1) + 2, 1, $motif, 2), "\n";
	    }
	}

	if( $i + 5 < $n )
	{
	    my $motif = substr($seq,$i,5);

	    #    01234
	    #     *
	    # 5' CCAGG 3'
	    # 5' GGTCC 3'
	    #       *

	    if( $motif eq "CCAGG" )
	    {
		print join(",", $name, ($i+1) + 1, 0, $motif, 1), "\n";
		print join(",", $name, ($i+1) + 3, 1, $motif, 3), "\n";
	    }

	    #    01234
	    #     *
	    # 5' CCTGG 3'
	    # 5' GGACC 3'
	    #       *

	    if( $motif eq "CCTGG" )
	    {
		print join(",", $name, ($i+1) + 1, 0, $motif, 1), "\n";
		print join(",", $name, ($i+1) + 3, 1, $motif, 3), "\n";
	    }
	}

	if( $i + 13 < $n )
	{
	    my $motif = substr($seq,$i,13);

	    # 0123456789012
	    #  *
	    # AACNNNNNNGTGC
	    # TTGNNNNNNCACG
	    #           *

	    if( substr($motif,0,3) eq "AAC" && substr($motif,9,4) eq "GTGC" )
	    {
		print join(",", $name, ($i+1) +  1, 0, $motif,  1), "\n";
		print join(",", $name, ($i+1) + 10, 1, $motif, 10), "\n";
	    }

	    # 0123456789012
	    #   *
	    # GCACNNNNNNGTT
	    # CGTGNNNNNNCAA
	    #            *

	    if( substr($motif,0,4) eq "GCAC" && substr($motif,9,3) eq "GTT" )
 	    {
		print join(",", $name, ($i+1) +  2, 0, $motif,  2), "\n";
		print join(",", $name, ($i+1) + 11, 1, $motif, 11), "\n";
	    }
	}
    }
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    my $name = "";
    
    open(FA, $file) || die "Can't open '$file'";
    
    while ( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line, 0, 1) eq ">" )
	{
	    $name = substr($line, 1);
	}
	else
	{
	    $refs{$name} .= $line;
	}
    }
    
    close(FA);
    
    return \%refs;
}
