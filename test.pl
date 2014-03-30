#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: blat.pl
#
#        USAGE: ./blat.pl
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
#      CREATED: 05.08.2012 15:28:44
#     REVISION: ---
#===============================================================================

## Loading modules and initializing variables ##
use strict;
use warnings;
use utf8;

use feature "switch";
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Image::Magick;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Math::BigFloat;
use Pod::Usage;

my $dbn = "(((.[.{.<..)))]}>";
my @dbnArr = split("", $dbn);

my %opening_brackets = ( 0 => '(', 1 => '[', 2 => '{', 3 => '<' );
my %closing_brackets = ( 0 => ')', 1 => ']', 2 => '}', 3 => '>' );

my $i = 3;
my $incorrect_brackets = "";
foreach my $key (keys(%opening_brackets)){
    $incorrect_brackets .= '\\'.$opening_brackets{$key} unless ( $key == $i );
}
my $re_open_bracket = '\\'.$opening_brackets{$i}.'['
    .$incorrect_brackets.'\.\)\]\}\>]*$';
print("$re_open_bracket\n");
$dbn =~ /(?=$re_open_bracket)/g;
print("Last opening bracket: $opening_brackets{$i} at position ". ($-[0]+1)."\n");
my $re_closing_bracket = '\\'.$closing_brackets{$i}.'.*$';
print("$re_closing_bracket\n");
$dbn =~ /(?=$re_closing_bracket)/;
print("Last closing bracket: $closing_brackets{$i} at position ". ($-[0]+1) ."\n");
