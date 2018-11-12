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
    # Methylation sites:
    #           *
    # 5'... GTGAAG ...3'     5'... CTTCAC ...3'
    # 3'... CACTTC ...5'     3'... GAAGTG ...5'
    #                               *

    # Methylation sites:
    #           *
    # 5'... TTCAAG ...3'     5'... CTTGAA ...3'
    # 3'... AAGTTC ...5'     3'... GAACTT ...5'
    #                               *

    my $seq = uc($$refs{$name});
    
    for( my $i = 0; $i < length($seq) - 6; $i++ )
    {
	my $motif = substr($seq,$i,6);

	if( $motif eq "GTGAAG" )
	{
	    print join(",", $name, ($i+4) + 1, 0, $motif, 4), "\n";
	}
	elsif( $motif eq "CTTCAC" )
	{
	    print join(",", $name, ($i+1) + 1, 1, $motif, 1), "\n";
	}
	elsif( $motif eq "TTCAAG" )
	{
	    print join(",", $name, ($i+4) + 1, 0, $motif, 4), "\n";
	}
	elsif( $motif eq "CTTGAA" )
	{
	    print join(",", $name, ($i+1) + 1, 1, $motif, 1), "\n";
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
