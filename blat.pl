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
use File::Find;
use Getopt::Long;
use LWP::Simple;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my $database = "";
my $query_list = "";
my $help = 0;
my @query_directories = ();
my @query_files = ();
my $verbose = 0;

GetOptions(
    "help|h"    => \$help,
# this has to be a DNA sequence so in case it's not we have to translate it
	"database|d=s"  => \$database,  # mostly just one file containing a single sequence
    "query|q=s" => \$query_list,         # in our case a text file that lists all .fa files to compare against the database
                                    # this is due to the fact that its more easy for the query to be a RNA sequence the
    "query-directory|u=s" => \@query_directories,
    "query-file|f=s" => \@query_files,
	"verbose|v+" => \$verbose);

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
$logger->info("++++ ".$this_file." has been started. ++++");

###############################################################
##
##              Application section
## 
###############################################################

# check @query_directories for .fa files and add them to @query_files
find( \&wanted, @query_directories);

# add all .fa files from $query_list to @query_files
if ( $query_list ) {
    open(my $query_list_fh, "<", $query_list ) or die ("Can't open $query_list.");
    while ( <$query_list_fh> ) {
        push(@query_files, $_) if ( $_ =~ /\.fa$/ && (-f $_) && ( -r $_) );
    }
    close($query_list_fh);
}

my %replace = (
    U => "T",
    u => "t"
);

my $regex = join( "|", keys %replace);
$regex = qr/$regex/;

my($filename, $directories, $suffix) = fileparse($database);
my $dna_database = $directories.$filename."dna.fa";
my ($content, $is_rna, $db) = (undef, undef, $database);

open(my $database_fh, "<", $database) or die ("Can't open $database.");
while ( <$database_fh> ) {
    if ( $_ =~ /[Uu]/) {
        $is_rna = 1;
        $db = $dna_database;
        $logger->info("$database contains a RNA sequence.");
        $_ =~ s/($regex)/$replace{$1}/g;
    }
    $content .= $_;
}
close( $database_fh );

if ( $is_rna ) {
    open(my $dna_database_fh, ">", $dna_database);
    print $dna_database $content;
    close( $dna_database_fh );
}

my $query = $db.".query.list";
open( my $query_fh,  ">", $query) or die ("Can't open $query.");
print $query_fh join("\n", @query_files);
close($query_fh);

my $result_file = $db.".blast.out";
system( "blat", "-q=rna", "-minIdentity=98", "-out=blast8", $db, $query, $result_file.".blast_out" ) == 0 
        or die ("Error during execution of BLAT search.\n $!");

###############################################################
##
##              Subroutine section
##
###############################################################

###############################################################################
##              
##  Invocation - \&wanted
##      - Subroutine used by File::Find
##
###############################################################################

sub wanted {
    my $logger = get_logger("RNAprobing");
    $logger->info("$_ is a directory") if ( -d $_ );
    $logger->info("$_ isn't a directory") unless ( -d $_ );
    if ( $_ =~ /\.fa$/ && ( -f $_ ) && ( -r $_ ) ) {
        $logger->info($File::Find::name." is a .fa file");
        # @files is a global variable
        push(@query_files, $File::Find::name);
    } 
}

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

