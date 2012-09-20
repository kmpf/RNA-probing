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

## Loading modules and initializing variables ##
use strict;
use warnings;
use utf8;
use File::Basename;
use Getopt::Long;
use LWP::Simple;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my $database = "";
my $query = "";
my $help = 0;
my $verbose = 0;

GetOptions(
    "help|h"    => \$help,
# this has to be a DNA sequence so in case it's not we have to translate it
	"database|d=s"  => \$database, # mostly just one file
# in our case a text file that lists all .fa files to compare against the database
# this is due to the fact that its more easy for the query to be a RNA sequence the
    "query|q=s" => \$query, 
	"verbose|v+" => \$verbose);

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

###############################################################
##
##              Application section
## 
###############################################################

my %replace = (
    U => "U",
    u => "t"
);

my $regex = join( "|", keys %replace);
$regex = qr/$regex/;

my ($directories, $filename, $suffix) = fileparse($database);
my $dna_database = $directories.$filename."dna.fa";
open(my $database_file, "<" , $database);

my ($content, $is_rna, $db) = (undef, undef, $database);

while ( <$database_file> ) {
    if ( $_ =~ /[Uu]/) {
        $is_rna = 1;
        $db = $dna_database;
        $logger->info("$database contains a RNA sequence.");
        $_ =~ s/($regex)/$replace{$1}/g;
    }
    $content .= $_;
}

close( $database_file );

if ( $is_rna ) {
    open(my $dna_database_file, ">", $dna_database);
    print $dna_database $content;
    close( $dna_database_file );
}

system("bash", "blat", "-q=rna", "-minIdentity=98", "-out=blast8", $db, $query, $file.".blast_out") == 0 
        or die ("Error during execution of BLAT search.\n $!");

###############################################################
##
##              Subroutine section
##
###############################################################

###############################################################
#
# Invocation - &configureLogger($verbosityLevel)
# - Configures and initialzes the Logger
# - $verbosityLevel = scalar value that sets log level
# -- 0 => $WARN
# -- 1 => $INFO
# -- 2 => $DEBUG
# 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger = get_logger("RNAprobing");
    $logger->info("Verbosity level: $verbose");
    SELECT:{
        if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
        if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
        if ($verbose >= 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
        else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
    }
    return $logger;
}

