#!/usr/bin/env perl
#===============================================================================
#
#         FILE: createDLLconf.pl
#
#        USAGE: ./createDLLconf.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: --rdf           File containing the RDF model
#               --pos           File containing a SPQRQL query defining the
#                               positive set
#               --neg           File containing a SPARQL query defining the
#                               negative set
#               -v, --verbose   More verbose output
#               -h, --help      Display help message
#               -m, --man       Display man page
#
# REQUIREMENTS: RDF::Helper, RDF::Query, RDF::Trine::Parser
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
use Getopt::Long;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Math::BigFloat;
use Pod::Usage;
use RDF::Query;
use RDF::Trine::Parser;
use RDF::Helper;
use Path::Class;
use Scalar::Util qw(blessed);
# load my modules
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;
require RNAprobing::BLASTresult;
require RNAprobing::RNAupFile;

my $rdf_file = "";
my $pos = "";
my $neg = "";
my $verbose = 0;
my $help = "";
my $man = 0;

GetOptions(
    "rdf=s" => \$rdf_file,
    "pos=s" => \$pos,
    "neg=s" => \$neg,
    "verbose|v+" => \$verbose,
    "help|h" => \$help,
    "man|m" => \$man) or pod2usage(-verbose => 1) && exit;

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( ($rdf_file eq "") && ($pos eq "") ) {
    pod2usage( -verbose => 1 ) && exit ;
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

open( my $neg_sparql_fh, "<", $neg) or die "Couldn't open file $neg. Error: $!";
while (<$neg_sparql_fh>) {
    $neg_sparql_query .= $_;
}
close $neg_sparql_fh;


&test($rdf);
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
        push(@positive_list, $pos_element);
      }
    }
} else {
    $logger->error("Could not use given SPARQL query.");
    exit;
}


my $pos_query_name = fileparse( $pos );
$pos_query_name =~ s/\.sparql$//g;
my $conf_file = $rdf_file;
$conf_file =~ s/\.rdf$/\.$pos_query_name\.conf/g;


open(my $conf_fh, ">", $conf_file) or die "Couldn't open file $conf_file. Error: $!";
print $conf_fh "// knowledge source
ks.type = \"OWL File\"
ks.fileName = \"$rdf_file\"
lp.positiveExamples = {".join(",", @positive_list)."}
lp.negativeExamples = {".join(",", @negative_list)."}";

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
    &test($model);
    my $result_list = ();


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

sub test{
    my $model = @_;
    my $logger = get_logger();
    $logger->info( blessed($model) );
    return $model;
}

__END__


=head1 NAME

SPARQLqueryRDF.pl - Querys a RDF model

=head1 DESCRIPTION

This script takes a RDF model and one or two files containing SPARQL queries. Mandatory describing positive and/or negative  will read the given input file(s) and do something
useful with the contents thereof.

=head1 SYNOPSIS

SPARQLqueryRDF.pl -f=</path/to/file> -v -v -v -t


=head1 OPTIONS

=over 8

=item B<--rdf >

RDF file containing the RDF model

=item B<--pos> 

File containing a SPARQL query defining the positive set

=item B<--neg>

File containing a SPARQL query defining the negative set

=item B<-h, --help>

Print help message

=item B<-m, --man>

Display man page

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (highest level if used 3 or more times) 

=back

=cut

