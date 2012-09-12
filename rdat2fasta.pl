#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: rdat2fasta.pl
#
#        USAGE: ./rdat2fasta.pl  
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
use File::Basename;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $help = 0;
my $verbose = 0;
my $to_dna = 0;
GetOptions(
	"file|f=s" => \@files,
    "toDNA|t" => \$to_dna,
	"verbose|v+" => \$verbose);

if ( $help || scalar(@files) == 0 ){
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



my @rdat_files = &checkFiles(@files);
foreach my $rdat_file ( @rdat_files ) {
    my $fasta_file = $rdat_file;
    $fasta_file =~ s/rdat$/fa/g;
    open(FASTAFILE, ">", $fasta_file) or die("Can't open $fasta_file.");
    my $rdat_object = RNAprobing::RDATFile->new($rdat_file);
    print FASTAFILE ">".$rdat_file."\n";
    my $sequence = $rdat_object->sequence();
    if ($to_dna) {
        $sequence =~ s/U/T/g;
        $sequence =~ s/u/t/g;
    }
    print FASTAFILE $sequence."\n";
    close(FASTAFILE);
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
    my $logger = get_logger("rdat2fasta");
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
    my $logger = get_logger("rdat2fasta");
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

__END__


=head1 NAME

rdat2fasta.pl - Creating .fasta files fom .rdat files

=head1 SYNOPSIS

rdat2fasta.pl -f=</path/to/file> -v -v -v -t


=head1 OPTIONS

=over 4

=item -f, --file=</path/to/file>

RDAT file(s) to be converted to FASTA files

=item -t, --toDNA

if set convert RNA sequences are converted into DNA sequences

=item -v, --verbose

verbosity level increases by multiple times option given

=item -h, --help

prints this help page

=back


=back
=head1 DESCRIPTION
B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut

