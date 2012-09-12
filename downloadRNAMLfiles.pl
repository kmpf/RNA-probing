#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: downloadRNAMLfiles.pl
#
#        USAGE: ./downloadRNAMLfiles.pl
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
my $help = 0;
my $verbose = 0;

GetOptions(
    "help|h"    => \$help,
	"file|f=s"  => \$file,
	"verbose|v+" => \$verbose);

my $base_url = "http://ndbserver.rutgers.edu/atlas/";
my $nmr = "nmr/structures/";
my $xray = "xray/structures/";

mkdir("./RNAML_files_from_NDB");

open(my $ndb_entries_file, ">", $file) or die("Can't open $file.");
while ( <$ndb_entries_file> ) {
    next if ( $_ =~ /ID/ );
    $_ =~ /^(\w{4,})/; 
    my $ndb_id = $1;
    my $nmr_url = $base_url . $nmr . "id/". lc($ndb_id) . "/" . uc($ndb_id) . "-rna1.xml";
    my @ndb_id_split = split('', uc($ndb_id) );
    my $xray_url = $base_url . $xray . $ndb_id_split[0] ."/". lc($ndb_id) . "/" . uc($ndb_id) . "-rna1.xml";
    $file = "./RNAML_files_from_NDB/".uc($ndb_id)."-rna1.xml";
    getstore($nmr_url, $file); # try downloading as if it is a NMR structure
    getstore($xray_url, $file); # try downloading as if it is a Xray structure
}

