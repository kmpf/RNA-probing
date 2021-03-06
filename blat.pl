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
	"database|d=s"  => \$database,
# mostly just one file containing a single sequence
    "query|q=s" => \$query_list,
# a text file that lists all .fa files to compare against the database, this is
# due to the fact that its more easy for the query to be a RNA sequence then
    "query-directory|u=s" => \@query_directories,
    "query-file|f=s" => \@query_files,
	"verbose|v+" => \$verbose);

if ( $help || $database eq "" || scalar(@query_files) == 0 && scalar(@query_directories) == 0){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
}

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

# what to replace to get from RNA to DNA
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
    next if ( $_ =~ /^>/);
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

my $result_file = $db;
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

__END__

=head1 NAME

blat.pl - Using BLAT to find similar sequences between a single, so called, database file and different query sequences

=head1 SYNOPSIS

 blat.pl -d </path/to/database.fa> -u </path/to/fa-files/>  -v -v -v
 blat.pl --database </path/to/database.fa> --query-directory </path/to/fa-files/> 

=head1 OPTIONS

=over 4

=item -d, --database=</path/to/database-file>

The database is a single fasta file containing a DNA sequence.

=item -q, --query=</path/to/file.list>

The file given lists the path to a fasta file per line. The sequence contained in every file is used as query.

=item -u, --query-directory=</path/to/fasta-directory/>

Every fasta file (ending on .fa) in this directory and all its subdirectories is used as a query sequence. The option can be given multiple times.

=item -f, --query-file=</path/to/file.fa>

A single fasta file can be given to this option. The option can be given multiple times.

=item -v, --verbose

This option increases the verbosity level everytime it is given.


=head1 DESCRIPTION

This script starts a BLAT sequence similarirty search. It compares the database sequence with the sequences given in the "query" files. All query files found are written as a list into a file named 
"<database>.query.list". And the result is written into a file named "<database>.blast_out".

 This perl script relies on the Log::Log4perl, Path::Class and XML::Simple modules. Please make sure the according packages are around on your system.

=head1 ARGUMENTS

=cut
