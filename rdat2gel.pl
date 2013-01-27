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
use Math::BigFloat;
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
    top_space => 100    # free space at the top of the gel

#    ltr_space => 33     # space to the left, at the top and to the right
);
my %PAA_GEL_CHAMBER = (
    width => 1000,
    height => 5000,
    top_space => 100
);

# Type of comb defines lane number and size
my %TEN_WELL_COMB = (
    nr_wells => 10,
    lane_width => 70,
    inter_lane_space => 26
);


# Type of detection method describes the colors of background and bands
my %EtBr = (
    background => "canvas:black",
    bands => "255"
);

my %RADIOGRAPHY = (
    background => "canvas:white",
    bands => "0"
);

# Type of gel describes the migration properties of a DNA fragment
my %TWO_PERCENT_AGAROSE = (
    POWER => Math::BigFloat->new(10),  # Watt used for gel run
    TIME => Math::BigFloat->new(60),   # duration of gel run
    LDNA => Math::BigFloat->new(150), # longest separable DNA fargment in bp
    LDIST => Math::BigFloat->new(100),  # traveled distance of LDNA fragment in pixels
    SDNA => Math::BigFloat->new(1),  # smallest separable DNA fragment in bp
    SDIST => Math::BigFloat->new(900) );   # traveled distance of SDNA fragment in pixels

my $BAND_HEIGHT = 10;


## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $gel_chamber_type = "Agarose";
my $gel_type = "2%";
my $comb_type = "10";
my $detection_type = "EtBr";
my $verbose = 0;
my $output_format = "png";
my $labeled_end = "5'";
# Gel parameter


GetOptions(
	"file|f=s" => \@files,
    "geltype|g=s" => \$gel_type,
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

    print $rdat_object->serialize_rdat_version()."\n";  # funktioniert
    print $rdat_object->serialize_name()."\n";          # funktioniert
    print $rdat_object->serialize_sequence()."\n";      # funktioniert
    print $rdat_object->serialize_structure()."\n";     # funktioniert
    print $rdat_object->serialize_offset()."\n";        # funktioniert
    print $rdat_object->serialize_seqpos()."\n";        # funktioniert
    print $rdat_object->serialize_mutpos()."\n";        # funktioniert
    print $rdat_object->serialize_annotation()."\n";
    print $rdat_object->serialize_comment()."\n";
    print $rdat_object->serialize_annotation_data()."\n";
    print $rdat_object->serialize_reactivity()."\n";
    print $rdat_object->serialize_reactivity_error()."\n";
}

my $gel = ($gel_chamber_type eq "Agarose") ? \%AGAROSE_GEL_CHAMBER : \%PAA_GEL_CHAMBER ;
my $comb;
if ($comb_type eq "10") {
    $comb = \%TEN_WELL_COMB;
} elsif ($comb_type ) {

}
my $detection = ($detection_type eq "EtBr") ? \%EtBr : \%RADIOGRAPHY;
my $standard;
if ($gel_type eq "2%") {
    $standard = \%TWO_PERCENT_AGAROSE;
} #elsif () {}
&make_gel(\@rdat_objects, $gel, $comb, $detection, $standard);   


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
## &make_gel($rdat_object, $gel, $comb, $detection);   
## - prints the gel picture
## - $rdat_object = array of RNAprobing::RDATFile object ref
## - $gel = hash ref
## - $comb = hash ref
## - $detection = hash ref
##
###############################################################


sub make_gel {
    # get all referemces which were used to call us
    my ($rdat_objects, $gel, $comb, $detection, $standard) = @_;

    my $logger = get_logger("rdat2gel");
    my $left_space = ($gel->{width} - ($comb->{nr_wells} * $comb->{lane_width} + 
        ($comb->{nr_wells} - 1) * $comb->{inter_lane_space}) ) / 2;
    $logger->info('$left_space: '.$left_space);


    # $gel is the image object everything should end up in
    my $gel_image = Image::Magick->new;
    $gel_image->Set(size => "$gel->{width}x$gel->{height}" );
#    $gel_image->Set(width => $gel->{width}, height => $gel->{height} );
    $gel_image->ReadImage($detection->{background});

#    my $lane = Image::Magick->new(width => ${$comb}{lane_width});
#    my $band = Image::Magick->new(size => ${$gel}{lane_width}."x".${$gel}{band_height});

    my $well_nr = 0;
    for (my  $j = 0; $j < scalar(@{$rdat_objects}); $j++) {
        my $rdat_object = ${$rdat_objects}[$j];
        my $sequence = $rdat_object->sequence();
        $logger->info("Sequence length: ".length($sequence));

#        $logger->info("Offset: ".$offset);
#        my $offset = $rdat_object->offset();

        my @seq_pos = @{ $rdat_object->seqpos() };
        $logger->info("# of sequence positions: ".scalar(@seq_pos));

        foreach my $reactivity (@{ $rdat_object->scaled_reactivity() }) {
            exit  if ( $well_nr == $comb->{nr_wells}); # Abbruch wenn zu viele wells benutzt werden
            $logger->info("Reactivity[$well_nr]: ".Dumper($reactivity) );
            $logger->info("# of reactivties for lane $well_nr: ". scalar( @{$reactivity} ));
            my $x = $left_space + $well_nr * $comb->{lane_width} + 
                    $well_nr * $comb->{inter_lane_space};
            for (my $i = 0; $i < scalar(@seq_pos); $i++) {
                my $frag_length;
                if ( $labeled_end eq "5'" ) {
                    $frag_length = $i + 1;
                } elsif ( $labeled_end eq "3'" ){
                    $frag_length = scalar(@seq_pos) - $i;
                } else {
                    $logger->error("Value $labeled_end not allowed for option \"-labeled-end\".") && exit;
                }                
#                my $mig_dist = &calculate_wanderlust($standard, $frag_length);
                my $mig_dist = &calculate_wanderlust($frag_length);
                my $y_start = $mig_dist + $gel->{top_space};
                #next if ($y > $gel->{height});
                $logger->info("Migration distance of $frag_length bp is $y_start px.");
                my $colour = sprintf("%d", 255 * (1 - $reactivity->[$frag_length-1]) );
                $logger->info("RGB of Nucleotide $frag_length is rgb($colour, $colour, $colour)");
                my $y_end = $y_start + $comb->{lane_width};
                $logger->info();
                $gel_image->Draw(primitive => "line", points => "$x,$y_start $x,$y_end",
                           stroke => "rgb($colour, $colour, $colour)", strokewidth => '5');
            }
            $well_nr++;
        }
    }
    $logger->info('$left_space: '.$left_space);
    $gel_image->Write('png:out.png');

#           $colour = sprintf("%d", 255 * (1 - $nuc_reactivity) );

#            $perlmagick_error = $band->Read('xc:white');
#            warn $perlmagick_error if $perlmagick_error;
#            $perlmagick_error = $band->Draw(primitive => "line", points => "0,15 $BAND_WIDTH,15", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
#            $perlmagick_error = $band->Draw(primitive => "line", points => "0,15 $BAND_WIDTH,15", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
#            warn $perlmagick_error if $perlmagick_error;
#            push (@{$lane}, @{$band});
#            @{$band} = ();           
#            $lane->Blur(sigma =>'3', radius => '3');
#        }
#        push ( @$gel, $lane->Append(stack=>'true') );
#        @$lane = ();
#    }
#    my $out = $gel->Append(stack => 'false');
#    $out->Write('png:out.png');
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
#    my ($standard, $frag_length) = @_;
    my ($frag_length) = @_;
    my $fl = Math::BigFloat->new($frag_length);
    my $y = $fl->blog();
    my $logger = get_logger("rdat2gel");
#    $logger->info($standard->{LDNA});
#    $logger->info($y);
#    $logger->info($standard->{LDNA}->copy()->blog());
    my $y1 = $standard->{LDNA}->copy()->blog();
    my $x1 = $standard->{LDIST};
    my $y2 = $standard->{SDNA}->copy()->blog();
    my $x2 = $standard->{SDIST};
    my $slope =  ($y2 - $y1) / ( $x2 - $x1 );
#    print '$slope: '.$slope.'=('.$y2.' - '.$y1.') / ('.$x2.' - '.$x1.")\n";
    my $x = (($y - $y1) + ($slope * $x1)) / $slope;
#    print '$x: '.$x.'=(('.$y.' - '.$y1.') + ('.$slope.' * '.$x1.')) / '.$slope."\n";
#    print $slope.' = '. ($y - $y1) / ($x - $x1)."\n";
    $x->precision(0);
    return $x;
}
