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
#use RDF::Trine::Parser;
#use RDF::Helper;
#use Path::Class;
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir); 
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;
#require RNAprobing::BLASTresult;
require RNAprobing::RNAupFile;


my %STANDARD = (
    POWER => Math::BigFloat->new(10),  # Watt used for gel run
    TIME => Math::BigFloat->new(60),   # duration of gel run
    LDNA => Math::BigFloat->new(10000), # longest separable DNA fargment in bp
    LDIST => Math::BigFloat->new(10),  # traveled distance of LDNA fragment in mm
    SDNA => Math::BigFloat->new(1000),  # smallest separable DNA fragment in bp
    SDIST => Math::BigFloat->new(50) );   # traveled distance of SDNA fragment in mm

my $wanderstrecke = Math::BigFloat->bzero();
my $length = Math::BigFloat->new(4000);
my $power = Math::BigFloat->new('10');
my $time = Math::BigFloat->new('40');
my $SLOPE = '100'; # 

&calculate_wanderlust(\%STANDARD, 4000);

my $q = Image::Magick->new;
$q->Set(size=>'100x100', 'pixel[50,50]'=>'255,0,0');
$q->ReadImage('xc:white');

my $colour = 155;
my $gel = Image::Magick->new;
$gel->Set(size => '1000x1000' );
$gel->ReadImage('canvas:black');
$gel->Draw(primitive => "rectangle", points => "100,500 200,600", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
$gel->Write('png:out.png');
my $out = $gel->Append(stack => 'true');



sub calculate_wanderlust{
    my ($standard, $frag_length) = @_;
    my $fl = Math::BigFloat->new($frag_length);
    my $test = Math::BigFloat->new($frag_length);
    print $test."\n";
    $test->blog(10);
    print $test."\n";

    my $y = $fl->blog();
    my $y1 = $standard->{LDNA}->blog();
    my $x1 = $standard->{LDIST};
    my $y2 = $standard->{SDNA}->blog();
    my $x2 = $standard->{SDIST};
    my $slope =  ($y2 - $y1) / ( $x2 - $x1 );
    print '$slope: '.$slope.'=('.$y2.' - '.$y1.') / ('.$x2.' - '.$x1.")\n";
    my $x = (($y - $y1) + ($slope * $x1)) / $slope;
    $x->precision(-2);
    print '$x: '.$x.'=(('.$y.' - '.$y1.') + ('.$slope.' * '.$x1.')) / '.$slope."\n";
    print $slope.' = '. ($y - $y1) / ($x - $x1)."\n";
}


# 
#`blat database.fa query.fa -q=rna -out=blast8 out.psl`;
#system("blat", "-q=rna", "-minIdentity=98", "-out=blast8", "database.fa", "query.fa", "out.blast_out");
