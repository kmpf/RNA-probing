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
    type => "Agarose",
    width => 1000,      # width of gel
    height => 1000,     # height of gel
    top_space => 100    # free space at the top of the gel
);

my %PAA_GEL_CHAMBER = (
    type => "PAA",
    width => 1800,
    height => 3050,
    top_space => 50
);

my %PAA_PROBING_GEL_CHAMBER = (
    type => "PAA",
    width => 1800,
    height => 6050,
    top_space => 50
);

# Type of comb defines lane number and size
my %TEN_WELL_COMB = (
    nr_wells => 10,
    lane_width => 70,
    inter_lane_space => 26
);


# Type of detection method describes the colors of background and bands
my %EtBr = (
    type => "EtBr",
    background => "canvas:black",
    bands => "255"
);

my %RADIOGRAPHY = (
    type => "Radiography",
    background => "canvas:white",
    bands => "0"
);

# Type of standard describes the migration properties of a DNA fragment
my %ONE_PERCENT_AGAROSE = (
    POWER => Math::BigFloat->new(10),   # Watt used for gel run
    TIME => Math::BigFloat->new(60),    # duration of gel run
    LDNA => Math::BigFloat->new(2000),  # longest separable DNA fargment in bp
    LDIST => Math::BigFloat->new(10),   # traveled distance of LDNA fragment in pixels
    SDNA => Math::BigFloat->new(200),    # smallest separable DNA fragment in bp
    SDIST => Math::BigFloat->new(900)   # traveled distance of SDNA fragment in pixels
);

my %TWO_PERCENT_AGAROSE = (
    POWER => Math::BigFloat->new(10),   # Watt used for gel run
    TIME => Math::BigFloat->new(60),    # duration of gel run
    LDNA => Math::BigFloat->new(1000),  # longest separable DNA fargment in bp
    LDIST => Math::BigFloat->new(10),   # traveled distance of LDNA fragment in pixels
    SDNA => Math::BigFloat->new(50),    # smallest separable DNA fragment in bp
    SDIST => Math::BigFloat->new(900)   # traveled distance of SDNA fragment in pixels
);

my %TEN_PERCENT_PAA = (
    POWER => Math::BigFloat->new(10),   # Watt used for gel run
    TIME => Math::BigFloat->new(60),    # duration of gel run
    LDNA => Math::BigFloat->new(130),  # longest separable DNA fargment in bp
    LDIST => Math::BigFloat->new(10),   # traveled distance of LDNA fragment in pixels
    SDNA => Math::BigFloat->new(1),    # smallest separable DNA fragment in bp
    SDIST => Math::BigFloat->new(2950)   # traveled distance of SDNA fragment in pixels
);

my $BAND_HEIGHT = 10;


## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $gel_chamber_type = "PAA";
my $gel_type = "10%";
my $comb_type = "10";
my $detection_type = "Radiography";
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
}

my ($gel, $comb, $detection, $standard);

if ($gel_chamber_type eq "PAA" && $gel_type eq "10%") {
    $gel = \%PAA_GEL_CHAMBER;
    $standard = \%TEN_PERCENT_PAA;
} elsif ($gel_chamber_type eq "PAA_probing" && $gel_type eq "10%") {
    $gel = \%PAA_PROBING_GEL_CHAMBER;
    $standard = \%TEN_PERCENT_PAA;
} elsif ($gel_chamber_type eq "Agarose" && $gel_type eq "2%") {
    $gel = \%AGAROSE_GEL_CHAMBER;
    $standard = \%TWO_PERCENT_AGAROSE;
} else {
    die $logger->error("Gel type and gel chamber are incompatible!");
} 
 
if ($comb_type eq "10") {
    $comb = \%TEN_WELL_COMB;
} elsif ($comb_type ) {

}

if ($detection_type eq "EtBr") {
    $detection = \%EtBr;
} elsif ($detection_type eq "Radiography") {
    $detection = \%RADIOGRAPHY;
} else {
    die $logger->error("Detection method not provided!");
}

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
    my $logger = get_logger();
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

    my $logger = get_logger();
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

    for (my $i = 0; $i < $comb->{nr_wells}; $i++) {
        my $x1 = $left_space + $i * $comb->{lane_width} + $i * $comb->{inter_lane_space};
        my $x2 = $x1 + $comb->{lane_width};
        my $y1 = $gel->{top_space} - 20;
        my $y2 = $gel->{top_space};
        my $colour = $detection->{bands};
        $gel_image->Draw(primitive => "rectangle", points => "$x1,$y1 $x2,$y2",
                           stroke => "rgb($colour, $colour, $colour)", strokewidth => '3');   
    }

    my $well_nr = 0;
    my $filename = "";
    for (my  $j = 0; $j < scalar(@{$rdat_objects}); $j++) {
        next  if ( $well_nr == $comb->{nr_wells}); # Abbruch wenn zu viele wells benutzt werden
        my $rdat_object = ${$rdat_objects}[$j];
        $filename .= $rdat_object->filename();
        my $sequence = $rdat_object->sequence();
        $logger->info("Sequence length: ".length($sequence));

#        $logger->info("Offset: ".$offset);
#        my $offset = $rdat_object->offset();

        my @seq_pos = @{ $rdat_object->seqpos() };
        $logger->info("# of sequence positions: ".scalar(@seq_pos));

        my $slope = &calculate_slope($standard);
        $logger->info("Slope of mig_dist~log(length): $slope");

        foreach my $reactivity (@{ $rdat_object->scaled_reactivity() }) {
            next if ( $well_nr == $comb->{nr_wells}); # Abbruch wenn zu viele wells benutzt werden
            $logger->info("Reactivity[$well_nr]: ".Dumper($reactivity) );
            $logger->info("# of reactivties for lane $well_nr: ". scalar( @{$reactivity} ));
            my $x_start = $left_space + $well_nr * $comb->{lane_width} + 
                    $well_nr * $comb->{inter_lane_space};
            for (my $i = 0; $i < scalar(@seq_pos); $i++) {
                my $frag_length;
                if ( $labeled_end eq "5'" ) {
                    $frag_length = $i + 1;
                } elsif ( $labeled_end eq "3'" ){
                    # 3' part of seq until base $i
                    $frag_length = scalar(@seq_pos) - $i;
                } else {
                    $logger->error("Value $labeled_end not allowed for option \"-labeled-end\".") && exit;
                }                
                my $mig_dist = &calculate_wanderlust($standard, $slope, $frag_length);
                my $y = $mig_dist + $gel->{top_space};
                next if ($y > $gel->{height});
#                $logger->info("Migration distance of $frag_length bp is $y_start px.");
                my $colour;
                if ($detection->{type} eq "EtBr") {
                    $colour = sprintf("%d", 255 * $reactivity->[$i] ); # wie berechne ich die Farbe gegeben die gel information
                } elsif ($detection->{type} eq "Radiography") {
                    $colour = sprintf("%d", 255 * (1 - $reactivity->[$i]) ); # wie berechne ich die Farbe gegeben die gel information
                }
#                $logger->info("RGB of Nucleotide $frag_length is rgb($colour, $colour, $colour)");
                my $x_end = $x_start + $comb->{lane_width};
                $logger->info("Band Position[$frag_length]: $x_start,$y $x_end,$y");
                $gel_image->Draw(primitive => "line", points => "$x_start,$y $x_end,$y",
                           stroke => "rgb($colour, $colour, $colour)", strokewidth => '5');
            }
            $well_nr++;
        }
    }
    $logger->info('$left_space: '.$left_space);
    $gel_image->Blur(sigma =>'5', radius => '5');
    $filename =~ s/\.rdat//g;
    $gel_image->Write("png:$filename.png");

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
    my ($standard, $slope, $frag_length) = @_;
    my $fl = Math::BigFloat->new($frag_length);
    my $y = $fl->blog();
    my $y1 = $standard->{LDNA}->copy()->blog();
    my $x1 = $standard->{LDIST};
    my $x = (($y - $y1) + ($slope * $x1)) / $slope;
    $x->precision(0);
    $logger->info("Fragment der Länge $frag_length wandert $x pixel");
    return $x;
}

sub calculate_slope{
    my ($standard) = @_;
    my $y1 = $standard->{LDNA}->copy()->blog();
    my $x1 = $standard->{LDIST};
    my $y2 = $standard->{SDNA}->copy()->blog();
    my $x2 = $standard->{SDIST};
    my $slope =  ($y2 - $y1) / ( $x2 - $x1 );
    return $slope;
}

