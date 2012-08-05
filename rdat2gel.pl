#!/usr/bin/env perl

## Loading modules and initializing variables ##
use strict;
use warnings;
use feature "switch";
use lib "/homes/olymp/christoph.kaempf/bin/scripts";
require ProbingRNA::RDATFile;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Image::Magick;
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Path::Class;

## Definition of constants used for Image::Magick

my $BAND_HEIGHT = 15;
my $BAND_WIDTH = 30;


## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $verbose = 0;
my $output_format = "png";
my $labeled_end = 5;
GetOptions(
	"file|f=s" => \@files,
	"verbose|v+" => \$verbose,
    "output-format|o=s" => \$output_format,
    "labeled-end|l=s" => \$labeled_end);

# Configuration file for the logger
my $log4perl_conf = file(dirname(__FILE__), "rdat2gel.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger = &configureLogger($verbose, "rdat2gel");

my @rdat_files = &checkFiles(@files);
# my @gel_files = ();
# foreach my $rdat_file ( @rdat_files ) {
#     if ($rdat_file =~ /rdat$/) {
#         $rdat_file =~ /(*)rdat$/$1/;
#         push(@gel_files, $rdat_file);
#     }
# 
# }


####### Application section #######
my @rdat_objects = ();
foreach my $rdat_file ( @rdat_files ) {
    my $rdat_object = ProbingRNA::RDATFile->new($rdat_file);
    push (@rdat_objects, $rdat_object);

    print $rdat_object->rdat_version()."\n";   # funktioniert
    print $rdat_object->name()."\n";           # funktioniert
    print $rdat_object->sequence()."\n";       # funktioniert
    print $rdat_object->structure()."\n";      # funktioniert
    print $rdat_object->offset()."\n";         # funktioniert
    print join(" ", @{ $rdat_object->seqpos()})."\n";  # funktioniert
    print join(" ", @{ $rdat_object->mutpos()})."\n";  # funktioniert
    
    my $annotation = $rdat_object->annotation();
    foreach my $key (keys %{ $annotation } ) {
        print $key.":".$annotation->{$key}."\n";
    }
    print $rdat_object->comment()."\n";
    
    my $annotation_data = $rdat_object->annotation_data();
    for (my $i = 0; $i < scalar( @{ $annotation_data } ); $i++ ) {
        foreach my $key (keys %{ $annotation_data->[$i]}){
            print $key.":".$annotation_data->[$i]->{$key}."|".ref(\$annotation_data->[$i]->{$key})."\n";
        }
    }
    foreach my $aref ($rdat_object->reactivity()) {
        print join( " ", @$aref )."\n";
    }
    my $reactivity_error = $rdat_object->reactivity_error();
    print "reactivity_error:";
    print join( " ",  @{ $reactivity_error->[0]} )."\n";

    &make_gel($rdat_object);

    
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
    my $rdat_object = shift;

    my $logger = get_logger("rdat2gel");
    my $sequence = $rdat_object->sequence();
    my $offset = $rdat_object->offset();
    my @seq_pos = @{ $rdat_object->seqpos() };
    my @reactivity = @{ $rdat_object->reactivity() };

    # log what we know
    $logger->info("Sequence length: ".length($sequence));
    $logger->info("Offset: ".$offset);
    $logger->info("# of sequence positions: ".scalar(@seq_pos));

    # the image object everything should end up in
    my $gel = Image::Magick->new;
    my $lane = Image::Magick->new;
    my $band = Image::Magick->new(size => $BAND_WIDTH."x".$BAND_HEIGHT);
    for ( my $j = 0; $j < scalar(@reactivity); $j++ ) {
        $logger->info("Lane number: ".$j);
        $logger->info("# of reactivties for lane $j: ". scalar( @{ $reactivity[$j]} ));
        my @lane_reactivities = @{ $reactivity[$j]};


        for ( my $i = 0; $i < scalar(@lane_reactivities); $i++ ){
            my ($nucleotide, $colour, $nuc_reactivity, $perlmagick_error);

            if ( $labeled_end == 5 ) {
                $nucleotide = $i;
            } elsif ( $labeled_end == 3 ){
                $nucleotide = scalar(@lane_reactivities) - $i;
            } else {
                $logger->error("Value $labeled_end not allowed for option \"-labeled-end\".") && exit;
            }

            $nuc_reactivity = $lane_reactivities[$nucleotide];
            $nuc_reactivity = 0 if ($nuc_reactivity < 0 );
            $nuc_reactivity = 1 if ($nuc_reactivity > 1 );
            $colour = sprintf("%d", 255 * (1 - $nuc_reactivity) );
            $logger->info("RGB of Nucleotide $nucleotide is rgb($colour, $colour, $colour)");

            $perlmagick_error = $band->Read('xc:white');
            warn $perlmagick_error if $perlmagick_error;
            $perlmagick_error = $band->Draw(primitive => "line", points => "0,15 $BAND_WIDTH,15", stroke => "rgb($colour, $colour, $colour)", strokewidth => '10');
            warn $perlmagick_error if $perlmagick_error;
            push (@$lane, @$band);
            @$band = ();
            
#            $lane->Blur(sigma =>'3', radius => '3');
        }
        push ( @$gel, $lane->Append(stack=>'true') );
        @$lane = ();
    }
    my $out = $gel->Append(stack => 'false');
    $out->Write('png:out.png');
}