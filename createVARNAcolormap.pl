#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: createVARNAcolormap.pl
#
#        USAGE: ./createVARNAcolormap.pl
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
#      CREATED: 05.08.2012 15:29:50
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

use feature "switch";
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Math::BigFloat;
use Pod::Usage;
use Path::Class;
my $module_dir = dirname(__FILE__);
push(@INC, $module_dir); 

###############################################################################
#
# Options section
#
################################################################################
my $help = 0;
my $man = 0;
my $off_file = "";
my $rdat_file = "";
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "off=s" => \$off_file, # mandatory
    "rdat=s" => \$rdat_file, # mandatory
    "verbose|v+" => \$verbose # optional
    );

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( $off_file eq "" || $rdat_file eq "" ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
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
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");

# require RNAprobing classes just after logger initialization
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;


################################################################################
#                 
# Files to objects
#                 
################################################################################

# check if input files exist
$off_file = &checkFiles($off_file) if ( $off_file ne "" );
$rdat_file = &checkFiles($rdat_file) if ( $rdat_file ne "" );

# creation of RDAT file object and information extraction
my $rdat_object     = RNAprobing::RDATFile->new($rdat_file);
my $rdat_filename   = fileparse($rdat_object->filename());
my $rdat_seq_startpos = $rdat_object->seq_startpos();
#my $rdat_reactivity = $rdat_object->data()->reactivity();
#my $rdat_seqpos_reactivity = $rdat_object->seqpos_reactivity_map();
my @rdat_seq        = split(//, $rdat_object->sequence());
my @rdat_seqpos     = @{$rdat_object->seqpos()};

# creation of OFF file object and information extraction
my $off_object      = RNAprobing::OFFFile->new($off_file);

###############################################################################
#
# Map the sequences according to the BLAT/BLAST output
#                 
###############################################################################

my $off_sequence = $off_object->sequence();
my $rdat_sequence = $rdat_object->sequence();
my $rdat_offset = $rdat_object->offset();
my ($start_position, $end_position); # position where both sequences start to be equal


if ($rdat_sequence =~ /$off_sequence/) {
    $start_position = $-[0];
    $end_position = $start_position + length($off_sequence);
} else {
    $logger->error("Can't match sequence in ".$off_object->filename().
        " to sequence in ".$rdat_object->filename());
    exit(1);
}

foreach my $rdat_index ( @{$rdat_object->data()->indices()} ) {
    my $color_map_filename = $rdat_filename."_".$rdat_index.".color_map";
    $color_map_filename =~ s/\.rdat//g;
    open(my $color_map, ">", $color_map_filename) or die "Couldn't open file $color_map_filename. Error: $!";
    for (my $i = $start_position; $i <= $end_position; $i++){
        my @sequence = split("", $rdat_sequence);
        print($color_map $rdat_object->data()->seqpos_reactivity_map($rdat_index)->{1+$rdat_offset+$i}."\n");
    }
}
exit(0);

###############################################################################
##              
##              Subroutine section
##
###############################################################################

###############################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $ERROR
## -- 1 => $WARN
## -- 2 => $INFO
## -- >2 => $DEBUG
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

###############################################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################################

sub checkFiles {
    my $test_file = shift;
    my $checked_file = "";
    my $logger = get_logger();
    # Check if files are readable
    if ( -r $test_file){
        $checked_file = $test_file;
        $logger->info("$test_file is readable.");
    } else {
        $logger->error("$test_file is not readable.");
        exit;
    }
    return $checked_file;
}

__END__

=head1 NAME

createVARNAcolormap.pl - Creates a file which can be used as a color map for the VARNA RNA visulalization tool.

=head1 SYNOPSIS

createVARNAcolormap.pl --off </path/to/off-file> --rdat </path/to/rdat-file>

=head1 DESCRIPTION

Creates a file which can be used as a color map for the VARNA RNA visulalization tool.

=head1 OPTIONS

=over 8

=item B<-h, --help>       

Display help message

=item B<-m, --man>

Display whole man page

=item B<--off>

DBN or OFF file containing the sequence and structure information of a RNA (normally extracted from NDB RNAML files); mandatory

=item B<--rdat>

RDAT file containing the probing information for a RNA (normally downloaded from the RMDB); mandatory

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (verbosest if used 3 or more times) 

=back

=cut
