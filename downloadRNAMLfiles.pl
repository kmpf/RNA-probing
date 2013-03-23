#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: downloadRNAMLfiles.pl
#
#        USAGE: ./downloadRNAMLfiles.pl
#
#  DESCRIPTION: This script downloads RNAML files from the Rutgers Nucleic
#               Acid Database given a file where every line begins with the
#               4-character long NDB ID. Such a file can be downloaded from
#               the NDB. 
#
#      OPTIONS: -f, --file      File containing a NDB ID (4 characters) per
#                               line 
#               -d, --dir       Directory in which to store all downloaded
#                               RNAML files
#               -v, --verbose   More verbose output
#               -h, --help      Display help message
#               -m, --man       Display man page
#
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
use LWP::Simple;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my $file = "";
my $store_dir = "";
my $help = 0;
my $man = 0;
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "file|f=s" => \$file,
    "dir|d=s" => \$store_dir,
    "verbose|v+" => \$verbose) or pod2usage( -verbose => 1 && exit;

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ($file eq "") {
    pod2usage( -verbose => 1 ) && exit;
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

my $base_url = "http://ndbserver.rutgers.edu/atlas/";
my $nmr = "nmr/structures/";
my $xray = "xray/structures/";
my $rnaml_ext = "-rna1.xml";
if ( -e $store_dir) {
    logger.info("Directory $store_dir already exists.");
elsif (! (-e $store_dir)){
    logger.debug("Directory $store_dir does not exist.")
    mkdir($store_directory);
    logger.debug("Directory $store_dir has been created.")
}
my $rnaml_files_list = "rnaml_files.list";

open(my $ndb_entries_file, "<", $file) or die("Can't open $file.");
while ( <$ndb_entries_file> ) {
    next if ( $_ =~ /^ID/ );
    $_ =~ /^(\w{4,})/; 
    my $ndb_id = $1;
    my $nmr_url = $base_url . $nmr . "id/". lc($ndb_id) . "/" . uc($ndb_id) . $rnaml_ext;
    my @ndb_id_split = split('', uc($ndb_id) );
    my $xray_url = $base_url . $xray . $ndb_id_split[0] ."/". lc($ndb_id) . "/" . uc($ndb_id) . $rnaml_ext;
    $file = $store_directory . uc($ndb_id) . "-rna1.xml";
    # try downloading as if it is a NMR structure
    if ( getstore($nmr_url, $file) == 200 ) {
        open (my $list_file, ">>", $rnaml_files_list);
        print $list_file uc($ndb_id) ."-rna1.xml";
        close $list_file;
        next;    
    }
    # try downloading as if it is a Xray structure
    if ( getstore($xray_url, $file) == 200 ) {
        open (my $list_file, ">>", $rnaml_files_list);
        print $list_file uc($ndb_id) ."-rna1.xml";
        close $list_file;
        next;
    } 
}

close( $ndb_entries_file );


__END__


=head1 NAME

downloadRNAMLfiles.pl - downloads RNAML files from the Rutgers NDB

=head1 DESCRIPTION

This script downloads RNAML files from the Rutgers Nucleic Acid Database given a file where every line begins with the 4-character long NDB ID. Such a file can be downloaded from the NDB. It creates a file listing all downloaded RNAML files, called 'rnaml_files.list'.

=head1 SYNOPSIS

downloadRNAMLfiles.pl -f=</path/to/file> -d </path/to/dir> -v -v -v


=head1 OPTIONS

=over 8

=item B<-f, --file>

File that holds a NDB ID (4 characters long) at the beginning of each line

=item B<-d, --dir>

Directory in which the RNAML files are stored.

=item B<-h, --help>

Print help message

=item B<-m, --man>

Display man page

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (highest level if used 3 or more times) 

=back

=cut

