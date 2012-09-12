#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: test.pl
#
#        USAGE: ./test.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Christoph Kaempf (CK), kaempf@bioinf.uni-leipzig.de
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 19.07.2012 11:03:53
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my %hash = (1 => "eins");
my @array = ("1", "2");

$hash{"ARRAY_REF"} = \@array;

my @rarray = @$hash{"ARRAY_REF"};

print scalar( @rarray )."\n"; 
# @array = ();
$hash->{"ARRAY_REF"} = [];
print ref %hash;
print scalar( @$hash{"ARRAY_REF"} )." ".scalar( @array ) ."\n"; 
