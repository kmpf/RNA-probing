#!/usr/bin/env perl
#===============================================================================
#
#         FILE: SPARQLqueryRDF.pl
#
#        USAGE: ./SPARQLqueryRDF.pl  
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

## Loading modules and initializing variables ##
use strict;
use warnings;
use feature "switch";
use Data::Dumper;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Image::Magick;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Math::BigFloat;
use Pod::Usage;
use RDF::Query;
use RDF::Trine::Parser;
use RDF::Helper;
use Path::Class;
use Pod::Usage;
use Scalar::Util qw(blessed);
# load my modules
my $module_dir = dirname(__FILE__);
# $module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;
require RNAprobing::BLASTresult;
require RNAprobing::RNAupFile;


my $rdf_file = "";
my $pos = "";
my $neg = "";
my $owl_file =file(dirname(__FILE__), "RNAprobing.owl");
my $verbose = 0;
my $help = "";
my $man = 0;

GetOptions(
    "rdf=s" => \$rdf_file,
    "pos=s" => \$pos,
    "neg=s" => \$neg,
    "owl=s" => \$owl_file,
    "verbose|v+" => \$verbose,
    "help|h" => \$help,
    "man|m" => \$man) or pod2usage(-verbose => 1) && exit;

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( $rdf_file eq "" ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
}

###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
my $log4perl_conf = file(dirname(__FILE__), "RNAprobing.log.conf");

pod2usage(-verbose => 1) && exit if ( $help );
pod2usage(-verbose => 1) && exit if ( ($rdf_file eq "") || ($pos eq "") );
pod2usage(-verbose => 2) && exit if ( $man );


# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");



# Configure RDF::Helper
my $rdf = RDF::Helper->new (
    BaseInterface => 'RDF::Trine',
    namespaces => {
        bioinf => "http://www.bioinf.uni-leipzig.de/~kaempf/RNAprobing.owl#",
        rdfs => "http://www.w3.org/2000/01/rdf-schema#",
        rdf => "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        xsd => "http://www.w3.org/2001/XMLSchema#",
        '#default' => "http://purl.org/rss/1.0/",
        },
    ExpandQNames => 1);

my $model = $rdf->model();
my $base_uri = 'http://purl.org/rss/1.0/';
my $parser = RDF::Trine::Parser->new( 'rdfxml' );
$parser->parse_file_into_model( $base_uri, $rdf_file, $model );

my $pos_sparql_query = "";
my $neg_sparql_query = "";

open( my $pos_sparql_fh, "<", $pos) or die "Couldn't open file $pos. Error: $!";
while (<$pos_sparql_fh>) {
    $pos_sparql_query .= $_;
}
close $pos_sparql_fh;

# &test($rdf);
# my $positive_list = &query_model( $pos_sparql_query, $rdf );
# my $negative_list= &query_model( $neg_sparql_query, $rdf );
my @positive_list = ();
my @negative_list = ();


if ( $pos_sparql_query =~ /[Ss][Ee][Ll][Ee][Cc][Tt]\s/ ) {
    # SPARQL SELECT Query
    my $query = RDF::Query->new( $pos_sparql_query );
    my $iterator = $query->execute( $model );
    while (my $row = $iterator->next) {
      # $row is a HASHref containing variable name -> RDF Term bindings
      my @vars = keys %$row;
      foreach my $key ( @vars) {
        my $pos_element = $row->{ $key }->as_string;
        $pos_element =~ s/[<>]/"/g;
        $logger->info($pos_element);
        push(@positive_list, $pos_element);
      }
    }
}

open( my $neg_sparql_fh, "<", $neg) or die "Couldn't open file $neg. Error: $!";
while (<$neg_sparql_fh>) {
    $neg_sparql_query .= $_;
}
close $neg_sparql_fh;

if ( $neg_sparql_query =~ /[Ss][Ee][Ll][Ee][Cc][Tt]\s/ ) {
    # SPARQL SELECT Query
    my $query = RDF::Query->new( $neg_sparql_query );
    my $iterator = $query->execute( $model );
    while (my $row = $iterator->next) {
      # $row is a HASHref containing variable name -> RDF Term bindings
      my @vars = keys %$row;
      foreach my $key ( @vars) {
        my $neg_element = $row->{ $key }->as_string;
        $neg_element =~ s/[<>]/"/g;
        $logger->info($neg_element);
        push(@negative_list, $neg_element);
      }
    }
}

my $pos_query_name = fileparse( $pos );
$pos_query_name =~ s/\.sparql$//g;
my $conf_file = $rdf_file;
$conf_file =~ s/\.rdf$/\.$pos_query_name\.conf/g;

my $conf_content = "";
$conf_content .= "// knowledge sources\n".
                 "ks1.type = \"OWL File\"\n".
                 "ks1.fileName = \"$rdf_file\"\n".
#                 "ks2.type = \"OWL File\"\n".
#                 "ks2.fileName = \"$owl_file\"\n\n".
                 "// reasoner\n".
                 "reasoner.type = \"fast instance checker\"\n".
                 "reasoner.sources = { ks1 }\n\n".
                 "// learning problem\n".
                 "lp.type = \"posNegStandard\"\n".
                 "lp.accuracyMethod = \"fmeasure\"\n\n".
                 "// learning algorithm\n".
                 "h.type =\"celoe_heuristic\"\n".
                 "// h.expansionPenaltyFactor = 0.2\n".
                 "h.expansionPenaltyFactor = 0.01\n\n".
                 "op.type = \"rho\"\n".
                 "op.useCardinalityRestrictions = true\n".
                 "op.useNegation = true\n\n".
                 "alg.type = \"celoe\"\n".
                 "// alg.nrOfThreads = 4\n".
                 "alg.maxExecutionTimeInSeconds = 60\n".
                 "alg.noisePercentage = 30\n\n".
                 "lp.positiveExamples = {".join(",", @positive_list)."}\n".
                 "lp.negativeExamples = {".join(",", @negative_list)."}\n";

$logger->debug($conf_content);
open(my $conf_fh, ">", $conf_file) or die "Couldn't open file $conf_file. Error: $!";
print $conf_fh $conf_content;
close $conf_fh;



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
## -- 0 => $ERROR
## -- 1 => $WARN
## -- 2 => $INFO
## -- >2 => $DEBUG
## 
###############################################################################

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
## &query_model($sparql_query, $model)
## - queries $model with $sparql_query
## 
###############################################################################

sub query_model {
    my ($sparql_query, $model) = @_;
#    print Dumper($model);
    my $logger = get_logger();
#    &test($model);
    my $result_list = "";

    if ( $sparql_query =~ /[Ss][Ee][Ll][Ee][Cc][Tt]\s/ ) {
        # SPARQL SELECT Query
        my $query = RDF::Query->new( $sparql_query );
        my $iterator = $query->execute( $model );
        while (my $row = $iterator->next) {
          # $row is a HASHref containing variable name -> RDF Term bindings
          my @vars = keys %$row;
          foreach my $key ( @vars) {
            $result_list .= $row->{ $key }->as_string.",";
          }
        }
    }

    if ( $sparql_query =~ /[Cc][Oo][Nn][Ss][Tt][Rr][Uu][Cc][Tt]\s/ ) {
        # SPARQL CONSTRUCT/DESCRIBE Query
        my $query = RDF::Query->new( $sparql_query );
        my $iterator = $query->execute( $model );
        while (my $st = $iterator->next) {
          # $st is a RDF::Trine::Statement object representing an RDF triple
          $result_list .= $st->as_string.","; 
          print $st->as_string;
        }
    }
    return $result_list;
}

__END__


=head1 NAME

makeDLLconf.pl - Querys a RDF model

=head1 SYNOPSIS

makeDLLconf.pl --rdf=</path/to/rdf-model> --pos=</path/to/pos-sparql-query> --neg--neg=</path/to/neg-sparql-query> -v -v -v

=head1 OPTIONS

=over 4

=item --rdf=</path/to/rdf-model>

RDF file containing the RDF model

=item --pos=</path/to/pos-sparql-query>

File containing a SPARQL query that returns the positive set of nodes from the RDF model

=item --neg=</path/to/neg-sparql-query>

File containing a SPARQL query that returns the negative set of nodes from the RDF model

=item -v, --verbose

Verbosity level increases by multiple times option given

=item -h, --help

Prints this help page

=item -m, --man

Prints the complete man page

=back


=back
=head1 DESCRIPTION
B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut

