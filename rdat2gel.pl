#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: rdat2gel.pl
#
#        USAGE: ./rdat2gel.pl  
#
#  DESCRIPTION: This script takes any number of provided RDAT files and creates
#               a gel picture from them. It is still subject to further develop-
#               ment. Better support for different combs and chamber types should
#               be implemented.
#
#      OPTIONS: -f, --file           Can be used multiple times naming a file
#                                    everytime it is called
#               -g, --geltype        [PAA, Agarose] 
#               -v, --verbose        Increases the verbosity level
#               -o, --output-format  Not used at the moment
#               -l, --labelled-end   Not used at the moment
#
# REQUIREMENTS: ImageMagick Perl bindings
#         BUGS: 
#        NOTES: 
#       AUTHOR: Christoph Kaempf (CK), kaempf@bioinf.uni-leipzig.de
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 05.08.2012 15:26:53
#     REVISION: 
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
push(@INC, $module_dir);

## Definition of constants describing properties of the gel

# Types of gel chamber defines size of the gel
my %CHAMBERS = (
    small => {
        height => 720,     # height of gel
        top_space => 120,    # free space at the top of the gel
        left_space => 300,
        right_space => 10,
        lane_width => 200,
        inter_lane_space => 70
    },
    large => {
        height => 6050,
        top_space => 120,
        left_space => 300,
        right_space => 10,
        lane_width => 200,
        inter_lane_space => 70
    }
    );

# Type of detection method describes the colors of background and bands
my %DETECTION = (
    EtBr => {
	type => "EtBr",
	background => "canvas:black",
	bands => "1" # high intensity bands should be 255
    },
    Radiography => {
	type => "Radiography",
	background => "canvas:white",
	bands => "0" #  high intensity bands should be 0 
    }
    );

# Type of standard describes the migration properties of a DNA fragment
my %GEL_TYPES = (
    onePercentAgarose => {
	type => "Agarose",
	percentage => Math::BigFloat->new(1),
	POWER => Math::BigFloat->new(10),   # Watt used for gel run
	TIME => Math::BigFloat->new(60),    # duration of gel run
	LDNA => Math::BigFloat->new(2000),  # longest separable DNA fargment in bp
	LDIST => Math::BigFloat->new(10),   # traveled distance of LDNA fragment in pixels
	SDNA => Math::BigFloat->new(200),    # smallest separable DNA fragment in bp
	SDIST => Math::BigFloat->new(900)   # traveled distance of SDNA fragment in pixels
    },
    twoPercentAgarose => {
	type => "Agarose",
	percentage => Math::BigFloat->new(2),
	POWER => Math::BigFloat->new(10),   # Watt used for gel run
	TIME => Math::BigFloat->new(60),    # duration of gel run
	LDNA => Math::BigFloat->new(1000),  # longest separable DNA fargment in bp
	LDIST => Math::BigFloat->new(10),   # traveled distance of LDNA fragment in pixels
	SDNA => Math::BigFloat->new(50),    # smallest separable DNA fragment in bp
	SDIST => Math::BigFloat->new(900)   # traveled distance of SDNA fragment in pixels
    },
    tenPercentPaa => {
	type => "PAA",
	percentage => Math::BigFloat->new(10),
	POWER => Math::BigFloat->new(10),   # Watt used for gel run
	TIME => Math::BigFloat->new(60),    # duration of gel run
	LDNA => Math::BigFloat->new(300),  # longest separable DNA fargment in bp
	LDIST => Math::BigFloat->new(500),   # traveled distance of LDNA fragment in pixels
	SDNA => Math::BigFloat->new(10),    # smallest separable DNA fragment in bp
	SDIST => Math::BigFloat->new(6000)   # traveled distance of SDNA fragment in pixels
    }
    );

my $BAND_HEIGHT = 20;

my @ENDS = ("five-prime", "three-prime");

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
# set defaults
my $gel_size = "large";
my $gel_type = "tenPercentPaa";
my $detection_type = "Radiography";
my $verbose = 0;
my $output_format = "png";
my $labelled_end = "five-prime";
# Gel parameter
my $help = 0;
my $man = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "file|f=s" => \@files,
    "gel-type=s" => \$gel_type, # onePercentAgarose || twoPercentAgarose || tenPercentPaa
    "gel-size=s" => \$gel_size, # small || large
    "detection=s" => \$detection_type, # EtBr || Radiography
    "verbose|v+" => \$verbose,
    "output-format|o=s" => \$output_format,
    "labelled-end|l=s" => \$labelled_end ); # five-prime || three-prime

if ( $help || scalar(@files) == 0 ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
} elsif ($man) {
    pod2usage( { -verbose => 2});
}

###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
my $log4perl_conf = file(dirname(__FILE__), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");
# Get the logger
my $logger_name = "RNAprobing";
my $logger = get_logger($logger_name);

&configureLogger($verbose);

$logger->info("++++ ".__FILE__." has been started. ++++");

# require RNAprobing classes just after logger initialization
require RNAprobing::RDATFile;

####### Application section #######
my $rdat_files = &checkFiles(\@files);

my @rdat_objects = ();
foreach my $rdat_file ( @{$rdat_files} ) {
    $logger->info($rdat_file);
    my $rdat_object = RNAprobing::RDATFile->new();
    $rdat_object->read_file($rdat_file);
    push (@rdat_objects, $rdat_object);
}

my ($gel, $gel_sizes, $detection);

if ( grep($_ eq $gel_type, keys(%GEL_TYPES)) == 1 ){
    $gel = $GEL_TYPES{$gel_type};
} else {
    die $logger->error("Unknown gel value \"". $gel_type 
		       ."\", should be one of: ". join(", ", keys(%GEL_TYPES))."!");
}

if ( grep($_ eq $gel_size, keys(%CHAMBERS)) == 1 ) {
    $gel_sizes = $CHAMBERS{$gel_size};
} else {
    die $logger->error("Unknown gel size value \"". $gel_size 
		       ."\", should be one of: ". join(", ", keys(%CHAMBERS))."!");
}

if ( grep($_ eq $detection_type, keys(%DETECTION)) == 1 ) {
    $detection = $DETECTION{$detection_type};
} else {
    die $logger->error("Unknown detection value \"". $detection_type
		       ."\", should be one of: ". join(", ", keys(%DETECTION))."!");
}

if ( !(grep($_ eq $labelled_end, values(@ENDS)) == 1) ) {
    die $logger->error("Unknown labelled end value \"". $labelled_end
		       ."\", should be one of: ". join(", ", values(@ENDS))."!");
}

&make_gel(\@rdat_objects, $gel, $gel_sizes, $detection);

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
## -- 0 => $ERROR
## -- 1 => $WARN
## -- 2 => $INFO
## -- >=3 => $DEBUG
## 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    $logger->info("Verbosity level: $verbose");
    SELECT:{
	    if ($verbose == 0){$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
	    if ($verbose == 1){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    else { $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
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
    my ($testfiles) = shift;
    my @checkedfiles = ();
    my $logger = get_logger();
    foreach  (@{$testfiles}) {
        ## Check if files are readable
	    if ( -r $_){
		push(@checkedfiles, $_);
		$logger->info("$_ is readable.")
	    } else {
		$logger->warn("$_ is not readable.");
	    }
    }
    return \@checkedfiles;
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
    # get all references which were used to call us
    my ($rdat_objects, $gel, $gel_sizes, $detection) = @_;
    my $logger = get_logger();

    # calculate number of lanes
    my $nr_of_lanes = 0;
    foreach my $rdat_object (@{$rdat_objects}) {
        my $data = $rdat_object->data();
        $nr_of_lanes += scalar( @{$data->indices()} );
    }

    # define all necessary local variables for spaces
    my $top_space = $gel_sizes->{top_space};
    my $left_space = $gel_sizes->{left_space};
    my $right_space = $gel_sizes->{right_space};
    my $lane_width = $gel_sizes->{lane_width};
    my $inter_lane_space = $gel_sizes->{inter_lane_space};
    my $image_width = $left_space + ($nr_of_lanes * $lane_width) + 
        ( ($nr_of_lanes - 1) * $inter_lane_space ) + $right_space;
    my $image_height = $gel_sizes->{height};

    $logger->info('$left_space: '.$left_space);
    $logger->info('$right_space: '.$right_space);
    $logger->info('$lane_width: '.$lane_width);
    $logger->info('$inter_lane_space: '.$inter_lane_space);
    $logger->info('$image_width: '.$image_width);
    $logger->info('$image_height: '.$image_height);

    # define $gel_image is the image object everything should end up in
    my $gel_image = Image::Magick->new;
    $gel_image->Set(size => $image_width."x".$image_height);
    $gel_image->ReadImage($detection->{background});

    my $text_colour = "";
    if ($detection->{type} eq "EtBr") {
        $text_colour = 255;
    } elsif ($detection->{type} eq "Radiography") {
        $text_colour = 0;
    }

    # Draw the lane numbers for each lane we've got
    for (my $i = 0; $i < $nr_of_lanes; $i++) {
        # Calculate position
        my $x1 = $left_space + $i * $lane_width + $i * $inter_lane_space +
            int($lane_width / 2);
        my $y1 = $top_space - 20;
        # Draw text
        my $lane_nr = $i + 1;
        $gel_image->Draw(
            primitive => "text",
            points => "$x1,$y1",
            pointsize => 100,
            fill => "rgb($text_colour, $text_colour, $text_colour)",
            text => "$lane_nr");
    }

    # Draw the distinct bands ...
    my $well_nr = 0;
    my $imagename = "";
    my $longest_rna_on_gel = &find_longest_rna_on_gel($rdat_objects);
    
    # ... for each RDAT file ...
    foreach my $rdat_object  (@{$rdat_objects}) {
        my ($filename, $directories) = fileparse($rdat_object->filename());
        $filename =~ s/\.rdat$//g;
        $imagename .=  $filename . "-";
        my $data = $rdat_object->data();
        $logger->info("Band: $well_nr");
        # ... which might contain several probing tracks
        foreach my $index ( @{$data->indices()} ) {
            my $pos_reac =  $data->seqpos_scaled_reactivity_map($index);
            $logger->info(Dumper($pos_reac));
            # $rna_length is the number of nucleotides we have reactivities for
            my $rna_length = scalar(@{$rdat_object->seqpos()});
            for( my $i = 0; $i < $rna_length; $i++ ) {
                # $pos is essentially a key for all the hashmaps created ...
                my $pos = $rdat_object->seqpos()->[$i];
                # ... but it is not the length of a fragment that is just $i
              
                next if ($pos <= 0);
                $logger->info($pos);
                my $x_start = $left_space + $well_nr * $lane_width + 
                    $well_nr * $inter_lane_space;
                # Calculate fragment length of migrating fragment
                my $frag_length = 
                    # $i is used for the fragment length calculation as I said
                    &calculate_fragment_length_by_label($i, $rna_length, $labelled_end);
                my $mig_dist = 
                    &calculate_fragment_migration($gel, $gel_sizes, $frag_length, $longest_rna_on_gel);
                my $y = $mig_dist + $top_space;
                next if ($y > $image_height);

                # Calcualte the colour of the bands given the detection type
                my $colour = "";
                if ($detection->{type} eq "EtBr") {
                    $colour = sprintf("%d", 255 * $pos_reac->{$pos} );
                } elsif ($detection->{type} eq "Radiography") {
                    $colour = sprintf("%d", 255 * (1 - $pos_reac->{$pos}) );
                }

                $logger->info("Nuc: $pos / Colour: $colour");
                my $x_end = $x_start + $lane_width;
                $logger->info("Band Position[$frag_length]: ".
                              "$x_start,$y $x_end,$y");
                $gel_image->Draw(primitive => "line",
                                 points => "$x_start,$y $x_end,$y",
                                 stroke => "rgb($colour, $colour, $colour)",
                                 strokewidth => "$BAND_HEIGHT");
                $logger->info("Nucleotide[$pos]: ".$rdat_object->offset_sequence_map()->{$pos} );
                # Mark every 10th nucleotide on the gel
                if ( $pos % 10 == 0 ) {
                    my $nuc_at_pos = uc($rdat_object->offset_sequence_map()->{$pos});
                    my $x_text = 10;
                    my $y_text = $y - 10;
                    $logger->info("Text Position: $x_text,$y_text");
                    $gel_image->Draw(
                        primitive => "text",
                        points => "$x_text,$y_text",
                        pointsize => 100,
                        fill => "rgb($text_colour, $text_colour, $text_colour)",
                        text => "$nuc_at_pos$pos");
                    $gel_image->Draw(
                        primitive => "line",
                        points => ($x_text).",$y ".($x_text + ($left_space - 50) ).",$y",
                        stroke => "rgb($text_colour, $text_colour, $text_colour)",
                        strokewidth => "10");
                }
            }
            $well_nr++;
        }
    }
    $imagename =~ s/-$//g;
    $gel_image->Write("png:$imagename.png");
}

#########################################################################################
##
## &calculate_fragment_migration($gel, $chamber, $frag_length, $logest_rna);
## - Performs file checks and returns an array with all succesfully checked files
## - $standard = Hash reference
## - $frag_length = length of fragment for which migration distance is to be
##                  calculated
##
################################################################################


sub calculate_fragment_migration{
    my ($gel, $gel_sizes, $frag_length, $longest_rna) = @_;
    my $logger = get_logger();
    # adjust your RNA lengths
    my ($b_frag_length, $b_longest_rna, $b_shortest_rna) = undef;
    $b_frag_length = &adjuste_frag_length($frag_length);
    $b_longest_rna = &adjuste_frag_length($longest_rna);
    $b_shortest_rna = &adjuste_frag_length(1);
    $logger->info("Adjusted frag_length: ". $b_frag_length);
    $logger->info("Adjusted longest_rna: ". $b_longest_rna);
    $logger->info("Adjusted shortest_rna: ". $b_shortest_rna);
    # determine and scale gel concentration
    my $b_gel_conc = Math::BigFloat->bnan();
    if ($gel->{type} eq "Agarose") {
        $b_gel_conc = $gel->{percentage}->copy()->bdiv(20);
    } elsif ($gel->{type} eq "PAA") {
        $b_gel_conc = $gel->{percentage}->copy()->bdiv(20);
    } else {
        $logger->error("Gel type unknown");
        exit(1);
    }
    $logger->info("Gel concentration: " . $b_gel_conc); 

    # calculate exponential function for:
    my $b_exp_frag_length = &calculate_e_function($b_gel_conc->copy(),
                                                  $b_frag_length->copy());
#    $logger->info("Exponential frag_length: ".$b_exp_frag_length);
    my $b_exp_longest_rna = &calculate_e_function($b_gel_conc->copy(),
                                                  $b_longest_rna->copy());
#    $logger->info("Exponential longest_rna: ".$b_exp_longest_rna);
    my $b_exp_shortest_rna = &calculate_e_function($b_gel_conc->copy(),
                                                   $b_shortest_rna->copy());
#    $logger->info("Exponential shortest_rna: ".$b_exp_shortest_rna);
    
    # calculate scale factor and scaled exp. function for frag_length
    my $scale_factor = $b_exp_shortest_rna->bsub($b_exp_longest_rna)->copy();
    my $scaled_exp_frag_length = $b_exp_frag_length->bsub($b_exp_longest_rna)->copy();
    
    # now we can calculate the migration (0 <= mig <= 1)
    my $migration = $scaled_exp_frag_length->bdiv($scale_factor)->copy;
#    $logger->info($migration ."=". $scaled_exp_frag_length ."/". $scale_factor);
#    $migration->bsub( Math::BigFloat->bone() );

    $logger->info("Scaled migration: ".$migration);
#    exit(1);
    my $lane_length = $gel_sizes->{height} - (2 * $gel_sizes->{top_space});
    my $max_migration = Math::BigFloat->new($lane_length);
    $logger->info("Max migration: ". $max_migration);
    $migration->bmul($max_migration);
    $logger->info("Migration: ".$migration);
    $migration->precision(0);
    return $migration;
}

sub adjuste_frag_length {
    my ($frag_length) = @_;
    my $b_frag_length = Math::BigFloat->new($frag_length);
    $b_frag_length->badd("0");
    $b_frag_length->bmul("0.01");
    return $b_frag_length->copy();
}

sub calculate_e_function {
    my ($b_gel_conc, $b_frag_length) = @_;
    my $lambda = Math::BigFloat->bone('-')
        ->bmul($b_gel_conc->bmul($b_frag_length))->copy();
    my $exponential = $lambda->bexp()->copy();
    return $exponential->copy();
}

sub calculate_fragment_length_by_label {
    my ($pos_in_rna, $total_rna_length, $labelled_end) = @_;
    # Calculate fragment length of migrating fragment
    my $frag_length;
    if ( $labelled_end eq "five-prime" ) {
        $frag_length = $pos_in_rna + 1;
    } elsif ( $labelled_end eq "three-prime" ){
        # 3' part of seq until base $pos_in_rna
        $frag_length = $total_rna_length - $pos_in_rna;
    } else {
        $logger->error("Value $labelled_end not allowed for ".
                       "option \"--labelled-end\". Should be ".
                       "'five-prime' or 'three-prime'.");
        exit(1);
    }
    return $frag_length;
}

sub find_longest_rna_on_gel{
    my ($rdat_objects) = @_;
    my $longest_rna_on_gel = 0;
    foreach my $rdat_object (@{$rdat_objects}) {
	my $data = $rdat_object->data();
	foreach my $index ( @{$data->indices()} ){
	    my $rna_length = scalar( @{$rdat_object->seqpos()} );
	    if ($rna_length > $longest_rna_on_gel){
		$longest_rna_on_gel = $rna_length;
	    }
	}
    }
    return $longest_rna_on_gel;
}

#sub calculate_slope{
#    my ($standard) = @_;
#    my $y1 = $standard->{LDNA}->copy()->blog();
#    my $x1 = $standard->{LDIST};
#    my $y2 = $standard->{SDNA}->copy()->blog();
#    my $x2 = $standard->{SDIST};
#    my $slope =  ($y2 - $y1) / ( $x2 - $x1 );
#    return $slope;
#}

__END__

=head1 NAME

rdat2gel.pl - creates a gel picture from RNA probing information in a RDAT file

=head1 SYNOPSIS

rdat2gel.pl [options] -f,--file F<rdat-file>

=head1 DESCRIPTION

This script creates a RDAT file, containing the results of I<in silico> probing experiment.
To perform a probing reaction it needs to be provided with a RNA sequence stored in a FASTA file and the reactivity rules in a special file format that looks like this:


=head1 OPTIONS

=over 8

=item B<-h, --help>       

Display help message

=item B<-m, --man>

Display whole man page

=item B<-f, --file>

RDAT file containing RNA probing information (mandatory). Can be given multiple times.

=item b<--gel-type [onePercentAgarose, twoPercentAgarose, tenPercentPaa]>

Indicates the type of the gel, this sets some variables concerning the migration of fragments.

=item b<--gel-size [small, large]>

Sets the size of the gel.

=item b<--detection [EtBr, Radiography]>

Can be either EtBr or Radiography. EtBr leads to bright bands on dark background. 
Radiography leads to dark bands on bright background.

=item b<--labelled-end [five-prime, three-prime]>

The resulting fragment after modification/cleavage is labelled on either the 5' or 3' end.

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (most verbose if used 3 or more times) 

=back

=cut
