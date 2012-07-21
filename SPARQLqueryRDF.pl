#!/usr/bin/env perl

## Loading modules and initializing variables ##
use strict;
use warnings;
use feature "switch";
use lib "/home/hubert/bin/scripts";
require ProbingRNA::RDATFile;
require ProbingRNA::OFFFile;
require ProbingRNA::BLASTresult;
require ProbingRNA::RNAupFile;
use Data::Dumper;
use File::Basename;
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

my $rdf_file = "";
my $sparql_file = "";

GetOptions(
    "rdf=s" => \$rdf_file,
    "sparql=s" => \$sparql_file);


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

open(SPARQL_QUERY, "<$sparql_file") or die "Couldn't open file $sparql_file. Error: $!";
my $sparql_query = "";
while (<SPARQL_QUERY>) {
    $sparql_query .= $_;
}
print $sparql_query;

if ( $sparql_query =~ /[Ss][Ee][Ll][Ee][Cc][Tt]\s/ ) {
    # SPARQL SELECT Query
    my $query = RDF::Query->new( $sparql_query );
    my $iterator = $query->execute( $model );
    while (my $row = $iterator->next) {
      # $row is a HASHref containing variable name -> RDF Term bindings
      my @vars = keys %$row;
      foreach my $key ( @vars){
          print $row->{ $key }->as_string."\n";
      }
    }
}

if ( $sparql_query =~ /[Cc][Oo][Nn][Ss][Tt][Rr][Uu][Cc][Tt]\s/ ) {
    # SPARQL CONSTRUCT/DESCRIBE Query
    my $query = RDF::Query->new( $sparql_query );
    my $iterator = $query->execute( $model );
    while (my $st = $iterator->next) {
      # $st is a RDF::Trine::Statement object representing an RDF triple
      print $st->as_string;
    } 
}

