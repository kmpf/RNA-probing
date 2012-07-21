#!/usr/bin/env perl

## Loading modules and initializing variables ##
use strict;
use warnings;
use lib "/home/hubert/bin/scripts";
require ProbingRNA::RDATFile;
use File::Basename;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;


## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my $verbose = 0;
my $to_dna = 0;
GetOptions(
	"file|f=s" => \@files,
    "toDNA" => \$to_dna,
	"verbose|v+" => \$verbose);

print "TO DNA: $to_dna\n";

# Configuration file for the logger
my $log4perl_conf = file(dirname(__FILE__), "rdat2fasta.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger = &configureLogger($verbose);



my @rdat_files = &checkFiles(@files);
foreach my $rdat_file ( @rdat_files ) {
    my $fasta_file = $rdat_file;
    $fasta_file =~ s/rdat$/fa/g;
    open(FASTAFILE, ">", $fasta_file) or die("Can't open $fasta_file.");
    my $rdat_object = ProbingRNA::RDATFile->new($rdat_file);
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
