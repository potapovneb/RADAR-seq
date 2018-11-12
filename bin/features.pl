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

my $o_class    = "";
my @o_features = ();
my $o_mtotal   = -1;
my $o_logfile  = "";

if ( @ARGV == 0 )
{
    print "usage: $0 features.csv\n";
    print "\n";
    print "options:\n";
    print "  --class\n";
    print "  --features\n";
    print "  --mtotal\n";
    print "  --logfile\n";
    print "\n";
    exit;
}

GetOptions(
    "class=f"    => \$o_class,
    "features=s" => \@o_features,
    "mtotal=f"   => \$o_mtotal,
    "logfile=s"  => \$o_logfile);

my $file = shift @ARGV;

my @default_features = qw/Length mA mC mG mT mTotal A C G T Total iA iC iG iT iTtotal Modifiable Fraction/;

my @features = @o_features ? split(/;/, join(";", @o_features)) : @default_features;

my @head = ();

open(CSV, $file) || die "Can't open '$file'";

while ( my $line = <CSV> )
{
    chomp($line);
    
    my @tokens = split(/,/, $line);
    
    if ( @head == 0 )
    {
	@head = @tokens;

	if ( $o_logfile ne "" )
	{
	    open(LOG, ">", $o_logfile) || die "Can't write '$o_logfile'";
	    
	    print LOG join(",", @head), "\n";
	}
    }
    else
    {
	my %entry = ();
	
	for ( my $i = 0; $i < @tokens; $i++ )
	{
	    $entry{$head[$i]} = $tokens[$i];
	}

	### ignore patches with less than N modified bases
	if ( $o_mtotal != -1 && $entry{"mTotal"} < $o_mtotal )
	{
	    next;
	}

	### assign class
	my @values = $o_class ne "" ? $o_class : $entry{"Class"};

	for ( my $i = 0; $i < @features; $i++ )
	{
	    push @values, sprintf("%i:%s", $i + 1, $entry{$features[$i]});
	}

	### print formatted features
	print join(" ", @values), "\n";
	
	if ( $o_logfile ne "" )
	{
	    print LOG join(",", map { $entry{$_} } @head), "\n";
	}
    }
}

close(CSV);

if ( $o_logfile ne "" )
{
    close(LOG);
}
