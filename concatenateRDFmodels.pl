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

################################################################################
#
# Options section
#
################################################################################
my $help = 0;
my $man = 0;
my @rdf_files = ();
my $rdf_out = "";
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "rdf=s"=> \@rdf_files,
    "out=s" => \$rdf_out,
    "verbose|v+" => \$verbose);

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( $rdf_out eq "" || scalar(@rdf_files) <= 2 ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
}

###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
# my $this_file = __FILE__;
# $this_file =~ s/scripts/RNAprobing/g;
my $log4perl_conf = file(dirname(__FILE__), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");

# create the model

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
foreach my $rdf_file (@rdf_files) {
    $parser->parse_file_into_model( $base_uri, $rdf_file, $rdf->model() );
}
$logger->info($rdf_out);
open (my $rdf_out_file, ">", $rdf_out)  or die("Can't open $rdf_out.");
my $string = $rdf->serialize( format => 'rdfxml');
#$logger->info($string);
print $rdf_out_file $string;
close $rdf_out_file;

###############################################################################
##              
##              Subroutine section
##
###############################################################################



###############################################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $WARN
## -- 1 => $INFO
## -- 2 => $DEBUG
## 
###############################################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    SELECT:{
	    if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
	    else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
    }
    return $logger;
}


__END__

=head1 NAME

concatenateRDFmodels.pl - Takes multiple RDF/XML models as input and writes one new model containing them all.

=head1 SYNOPSIS

concatenateRDFmodels.pl --rdf F<rdat-file> --out F<out.rdf> [-v -v -v] 

=head1 DESCRIPTION

This script concatenates all RDF/XML models given. The newly created model is then written to an user-specified file.

=head1 OPTIONS

=over 8

=item B<-h, --help>   

Display help message

=item B<-m, --man>

Display whole man page

=item B<--rdf /path/to/in.rdf>

All RDF/XML models that will be concatenated are given to the script via this option. Must be given at least two times.

=item B<--out /path/to/out.rdf>

The file where the new RDF/XML model is written to.

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (verbosest if used 3 or more times) 

=back

=cut
