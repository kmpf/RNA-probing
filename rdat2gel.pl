#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: rdat2gel.pl
#
#        USAGE: ./rdat2gel.pl  
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
#      CREATED: 05.08.2012 15:26:53
#     REVISION: ---
#===============================================================================

## Loading modules and initializing variables ##
use strict;
use warnings;
use utf8;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Image::Magick;
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Path::Class;
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir); 
require RNAprobing::RDATFile;

## Definition of constants describing properties of the gel

# Types of gel chamber defines size of the gel
my %AGAROSE_GEL_CHAMBER = (
    width => 1000,      # width of gel
    height => 1000,     # height of gel
    band_height => 10,  # 
    ltr_space => 33,  # space to the left, at the top and to the right
    lane_space => 26);
my %PAA_GEL_CHAMBER = (

);

# Types of combs define  
my %TEN_WELL_COMB = (
    nr_wells => 10,
    well_size_percent => 7);




## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $geltype = "Agarose";
my $verbose = 0;
my $output_format = "png";
my $labeled_end = 5;
# Gel parameter


GetOptions(
	"file|f=s" => \@files,
    "geltype|g=s" => \$geltype,
	"verbose|v+" => \$verbose,
    "output-format|o=s" => \$output_format,
    "labeled-end|l=s" => \$labeled_end);

###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
my $this_file = __FILE__;
$this_file =~ s/scripts/RNAprobing/g;
my $log4perl_conf = file(dirname($this_file), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");


my @rdat_files = &checkFiles(@files);

####### Application section #######
my @rdat_objects = ();
foreach my $rdat_file ( @rdat_files ) {
    my $rdat_object = RNAprobing::RDATFile->new();
    $rdat_object->read_file($rdat_file);
    push (@rdat_objects, $rdat_object);

    print $rdat_object->serialize_rdat_version()."\n";   # funktioniert
    print $rdat_object->serialize_name()."\n";           # funktioniert
    print $rdat_object->serialize_sequence()."\n";       # funktioniert
    print $rdat_object->serialize_structure()."\n";      # funktioniert
    print $rdat_object->serialize_offset()."\n";         # funktioniert
    print $rdat_object->serialize_seqpos()."\n";  # funktioniert
    print $rdat_object->serialize_mutpos()."\n";  # funktioniert
    print $rdat_object->serialize_annotation()."\n";
    print $rdat_object->serialize_comment()."\n";
    print $rdat_object->serialize_annotation_data()."\n";
    print $rdat_object->serialize_reactivity()."\n";
    print $rdat_object->serialize_reactivity_error()."\n";

    my $gel = ($geltype eq "Agarose") ? \%AGAROSE_GEL : \%PAA_GEL ;

    &make_gel($rdat_object, $gel);   
}

###############################################################
##              
##              Subroutine section
##
###############################################################

###############################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $WARN
## -- 1 => $INFO
## -- 2 => $DEBUG
## 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    print $verbose;
    SELECT:{
	    if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
	    else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
    }
    return $logger;
}

###############################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################

sub checkFiles {
    my @testfiles = shift;
    my @checkedfiles = ();
    my $logger = get_logger("rdat2gel");
    foreach (@testfiles) {
        ## Check if files are readable
	    if ( -r $_){
		    push(@checkedfiles, $_);
            $logger->info("$_ is readable.")
	    } else {
            $logger->error("$_ is not readable.");
        }
    }
    return @checkedfiles;
}

###############################################################
##
## &make_gel($sequence, $offset, @seq_pos, @reactivity);
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################


sub make_gel {
    # get the RDAT object that was used to call us
    my ($rdat_object, $gel) = @_;

    my $logger = get_logger("rdat2gel");
    my $sequence = $rdat_object->sequence();
    my $offset = $rdat_object->offset();
    my @seq_pos = @{ $rdat_object->seqpos() };
    my $reactivity = $rdat_object->scaled_reactivity() ;
    
    # log what we know
    $logger->info("Sequence length: ".length($sequence));
    $logger->info("Offset: ".$offset);
    $logger->info("# of sequence positions: ".scalar(@seq_pos));
    $logger->info('Reactivities: '.Dumper($reactivity) );

    # $gel is the image object everything should end up in
    my $gel = Image::Magick->new(width => ${$gel}{width}, height => ${$gel}{height} );
    my $ 
    my $lane = Image::Magick->new(width=>${$gel}{lane_width});
    
    my $band = Image::Magick->new(size => ${$gel}{lane_width}."x".${$gel}{band_height});
    for ( my $j = 0; $j < scalar(@{ $reactivity}); $j++ ) {
        $logger->info("Lane number: ".$j);
        $logger->info("# of reactivties for lane $j: ". scalar( @{ $reactivity}[$j] ));
        my $lane_reactivities = $reactivity->[$j];
        print '@{$lane_reactivities}: '.Dumper(@{$lane_reactivities})."\n";

        for ( my $i = 0; $i < scalar(@{$lane_reactivities}); $i++ ){
            my ($nucleotide, $colour, $nuc_reactivity, $perlmagick_error);

            if ( $labeled_end == 5 ) {
                $nucleotide = $i;
            } elsif ( $labeled_end == 3 ){
                $nucleotide = scalar(@{$lane_reactivities}) - $i;
            } else {
                $logger->error("Value $labeled_end not allowed for option \"-labeled-end\".") && exit;
            }

            $nuc_reactivity = $lane_reactivities->[$nucleotide];
            if ($nuc_reactivity < 0 ) {
                $nuc_reactivity = 0;
                $logger->error("There is a nucleotide with reactivity < 0");
            }
            if ($nuc_reactivity > 1 ) {
                $nuc_reactivity = 1;
                $logger->error("There is a nucleotide with reactivity > 1");
            }
            $colour = sprintf("%d", 255 * (1 - $nuc_reactivity) );
            $logger->info("RGB of Nucleotide $nucleotide is rgb($colour, $colour, $colour)");

            $perlmagick_error = $band->Read('xc:white');
            warn $perlmagick_error if $perlmagick_error;
            $perlmagick_error = $band->Draw(primitive => "line", points => "0,15 $BAND_WIDTH,15", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
#            $perlmagick_error = $band->Draw(primitive => "line", points => "0,15 $BAND_WIDTH,15", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
            warn $perlmagick_error if $perlmagick_error;
            push (@{$lane}, @{$band});
            @{$band} = ();
            
#            $lane->Blur(sigma =>'3', radius => '3');
        }
        push ( @$gel, $lane->Append(stack=>'true') );
        @$lane = ();
    }
    my $out = $gel->Append(stack => 'false');
    $out->Write('png:out.png');
}

###############################################################
##
## &calculate_wanderlust($tandard, $frag_length);
## - Performs file checks and returns an array with all succesfully checked files
## - $standard = Hash reference
## - $frag_length = length of fragment for which migration distance is to be calculated
##
###############################################################


sub calculate_wanderlust{
    my ($standard, $frag_length) = @_;
    my $fl = Math::BigFloat->new($frag_length);
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
