#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: createRNAprobingModel.pl
#
#        USAGE: ./createRNAprobingModel.pl
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
use Image::Magick;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Math::BigFloat;
use Pod::Usage;
use RDF::Trine::Parser;
use RDF::Helper;
use Path::Class;
my $module_dir = dirname(__FILE__);
# $module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir); 
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;
require RNAprobing::BLASTresult;
require RNAprobing::RNAupFile;

###############################################################################
#
# Options section
#
################################################################################
my $help = 0;
my $rdf_file ="";
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "rdf=s"=> \$rdf_file,
    "verbose|v+" => \$verbose);

if ( $help ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
} elsif ( $rdf_file eq "" ){
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
$logger->info("++++ ".__FILE__." has been started. ++++");

my $rdf = RDF::Helper->new(
    BaseInterface => 'RDF::Trine',
    namespaces => {
        bioinf => "http://www.bioinf.uni-leipzig.de/~kaempf/RNAprobing.owl#",
        dcterms => 'http://purl.org/dc/terms/',
        rdfs => "http://www.w3.org/2000/01/rdf-schema#",
        rdf => "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        xsd => "http://www.w3.org/2001/XMLSchema#",
        '#default' => "http://purl.org/rss/1.0/",
        },
    ExpandQNames => 1);

# HowTo parse a file into a model

my $parser = RDF::Trine::Parser->new( 'rdfxml' );
my $base_uri = 'http://www.bioinf.uni-leipzig.de/~kaempf/RNAprobing.owl#';
$parser->parse_file_into_model( $base_uri, $rdf_file, $rdf->model() );

# Write the reduced model to file

my $reduced_rdf = $rdf_file;
$reduced_rdf =~ s/\.rdf$/\.reduced\.rdf/g;
open (my $rdf_fh, ">", $reduced_rdf)  or die("Can't open $reduced_rdf.");
my $string = $rdf->serialize( format => 'rdfxml');
# my $string = $rdf->serialize( format => 'turtle');
# $logger->info("$string");
print $rdf_fh $string;
close $rdf_fh;

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
    print Dumper($logger);
    SELECT:{
	    if ($verbose == 0){$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
	    if ($verbose == 1){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    else { $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
    }
    return $logger;
}

__END__

=h1

clearRNAprobingModel.pl - Removes http://www.bioinf.uni-leipzig.de/~kaempf/RNAprobing.owl#hasExperimentalUnpairedProbability from RDF model

=head1 SYNOPSIS

clearRNAprobingModel.pl --rdf=</path/to/rdf-model> --pos=</path/to/pos-sparql-query> --neg--neg=</path/to/neg-sparql-query> -v -v -v

=head1 OPTIONS

=over 4

=item --rdf=</path/to/rdf-model>

RDF file containing the RDF model

=item -v, --verbose

Verbosity level increases by multiple times option given

=item -h, --help

Prints this help page

=back

=cut

