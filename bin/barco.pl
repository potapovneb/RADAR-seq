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

my $o_bc = "";

GetOptions( "bc=s" => \$o_bc );

my $bamfile = shift @ARGV;

my %data = ();

open(BAM,"samtools view -h $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
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

    if( $o_bc ne "" && index($line,"@") == 0 )
    {
	print $line;
    }

    my ($movie,$zmw,$range) = split(/\//,$qname);
    my $bc = "";

    ### parse PacBio-specific tags
    for my $tag ( @tag )
    {
	my ($name,$type,$value) = split(/:/,$tag);
	
	if( $name eq "bc" )
	{
	    $bc = $value;

	    $data{$bc}{"$movie/$zmw"} = 1;

	    if( $o_bc ne "" && $bc eq $o_bc )
	    {
		print $line;
	    }

	    last;
	}
    }
}    

close(BAM);

print STDERR join("\t", "Barcode", "N", "Type", "B1", "B2"), "\n";

for my $bc ( keys %data )
{
    my ($type,$b1,$b2) = split(/,/,$bc);

    my $n = scalar(keys %{$data{$bc}});

    print STDERR join("\t", $bc, $n, $type, $b1, $b2), "\n";
}
