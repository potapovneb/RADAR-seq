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

my $bamfile = shift @ARGV;
my $wl_file = shift @ARGV;

### load whitelist
my %wl = ();

open(WL, $wl_file) || die "Can't open '$wl_file'";

while( my $line = <WL> )
{
    chomp($line);

    my ($movie, $zmw) = split(/\//, $line);

    $wl{$movie}{$zmw} = 1;
}

close(WL);

### filter aligned reads
open(BAM,"samtools view -h $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    if( index($line, "@") == 0 )
    {
	print $line;
    }
    else
    {
	my @tokens = split(/\t/, $line);
	
	my ($movie,$zmw,$coord) = split(/\//, $tokens[0]);
	
	if ( exists $wl{$movie}{$zmw} )
	{
	    print $line;
	}
    }
}

close(BAM);
