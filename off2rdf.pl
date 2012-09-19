#!/usr/bin/env perl
#===============================================================================
#
#         FILE: off2rdf.pl
#
#        USAGE: ./off2rdf.pl  
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
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;
require RNAprobing::OFFFile;
require RNAprobing::BLASTresult;
require RNAprobing::RNAupFile;
 
###############################################################################
#
# Options section
#
###############################################################################
my $off_file = "";
my $verbose = 0;

GetOptions(
    "off=s" => \$off_file,
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

###############################################################################
#                 
# Files to objects
#                 
###############################################################################

$off_file = &checkFiles($off_file);

# creation of OFF file object and information extraction
my $off_object      = RNAprobing::OFFFile->new($off_file);
my @off_seq         = split(//, $off_object->sequence());

###############################################################################
#
# Creating hashes for all necessary information
#   - Keys are almost always the position indices of the query sequence!
#                 
###############################################################################

&generate_rdf_model( \@off_seq, $off_object );

###############################################################################
##              
##              Subroutine section
##
###############################################################################

###############################################################################
##
## &generate_rdf_model( \@off_seq, $off_object )
## - 
## - 
##
###############################################################################

sub generate_rdf_model {
    my $off_seq = shift;
    my $off_object = shift;
    my $off_query_start = 0;
    my $off_query_end = scalar( @{ $off_seq } );
    my $logger = get_logger();

    my $rdf_filename = fileparse($off_object->filename());
    $rdf_filename =~ s/off$/rdf/g;

    my $off_id = $off_object->fasta_id();

    # Configure RDF::Helper
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
    my $model = $rdf->model();
    $parser->parse_file_into_model( $base_uri, '/home/hubert/Data/ontologies/RNAprobing.owl', $model );

    # define classes

    my $ribonucleotide_uri = 'bioinf:Ribonucleotide';
    my $adenosine_uri = 'bioinf:Adenosine';
    my $cytidine_uri = 'bioinf:Cytidine';
    my $guanosine_uri = 'bioinf:Guanosine';
    my $uridine_uri = 'bioinf:Uridine';
    my $hoogsteen_edge_uri = 'bioinf:HoogsteenEdge';
    my $sugar_edge_uri = 'bioinf:SugarEdge';
    my $watson_crick_edge_uri = 'bioinf:WatsonCrickEdge';

    # define object properties
    
    my $has_edge_uri = 'bioinf:hasEdge';
    my $has_tertiary_interaction_with_uri = 'bioinf:hasTertiaryInteractionWith';
    my $is_three_prime_of_uri = 'bioinf:isThreePrimeOf';
    my $is_five_prime_of_uri = 'bioinf:isFivePrimeOf';
    my $is_paired_with_uri = 'bioinf:isPairedWith';

    # define data properties

    my $has_glycosidic_bond_orientation_uri = 'bioinf:hasGlycosidicBondOrientation';

    # populate the RDF graph

    foreach my $querypos ( $off_query_start .. $off_query_end ) {
        my $nuc_uri = 'bioinf:'.$off_id.'/'.$off_seq[$querypos-1].$querypos;
        $logger->debug("Nucleotide: $nuc_uri");

        # isThreePrimeOf

        if ( $querypos < $off_query_end ) {
            my $three_prime_nuc_uri = 'bioinf:'.$off_id.'/'.$off_seq[$querypos].($querypos+1);
            $rdf->assert_resource($three_prime_nuc_uri, $is_three_prime_of_uri, $nuc_uri);
        }

        # isFivePrimeOf

        if ( $querypos > $off_query_start ) {
            my $five_prime_nuc_uri =  'bioinf:'.$off_id.'/'.$off_seq[$querypos-2].($querypos-1);
            $rdf->assert_resource($five_prime_nuc_uri , $is_five_prime_of_uri, $nuc_uri);
        }

        # type of nucleotide at position $querypos

        $rdf->assert_resource($nuc_uri, 'rdf:type', $adenosine_uri) if ( $off_seq[$querypos-1] eq 'A');
        $rdf->assert_resource($nuc_uri, 'rdf:type', $cytidine_uri) if ( $off_seq[$querypos-1] eq 'C');
        $rdf->assert_resource($nuc_uri, 'rdf:type', $guanosine_uri) if ( $off_seq[$querypos-1] eq 'G');
        $rdf->assert_resource($nuc_uri, 'rdf:type', $uridine_uri) if ( $off_seq[$querypos-1] eq 'U');
        $rdf->assert_resource($nuc_uri, 'rdf:type', $ribonucleotide_uri) if ( $off_seq[$querypos-1] eq 'N');
    }

    # insert all Leontis-Westhof base pairs and tertiary interactions
    for ( my $i = 0; $i < scalar( @{ $off_object->edges() } ); $i++) {
        my @edge = @{$off_object->edges()->[$i]};
        $logger->debug($edge[0]);
        $logger->debug(scalar( @{ $off_object->edges() } ));
        if ( $edge[0] ne "LW" ){
            $logger->info("Unknown notation $edge[0] trying next line.");
            next;
        }
        my $five_prime_nuc_uri = 'bioinf:'.$off_id.'/'.$off_seq[$edge[1]-1].$edge[1];
        my $three_prime_nuc_uri = 'bioinf:'.$off_id.'/'.$off_seq[$edge[2]-1].$edge[2];
        if ( $edge[3] =~ /^[ct!][WEHST][WEHST]$/ ) {
            $edge[3] =~ /^([ct!])([WEHST])([WEHST])$/;
            my $orientation = $1;
            my $five_prime_edge = $2;
            my $three_prime_edge = $3;

            # insert TertiaryInteraction

            if ( $orientation eq '!' && $five_prime_edge eq 'T' && $three_prime_edge eq 'T') {
                $rdf->assert_resource($five_prime_nuc_uri, $has_tertiary_interaction_with_uri, $three_prime_nuc_uri );
            } elsif ($orientation eq 'c' || $orientation eq 't' ) {

                # insert Watson-Crick, Hoogsteen or Sugar edges

                &insert_edge( $rdf, $five_prime_nuc_uri, $five_prime_edge );
                &insert_edge( $rdf, $three_prime_nuc_uri, $three_prime_edge );
                my $glycosidic_bond_orientation;
                if ( $orientation eq "c" ){
                    $glycosidic_bond_orientation = 'bioinf:Cis';
                } elsif ($orientation eq "t" ){
                    $glycosidic_bond_orientation = 'bioinf:Trans';
                }
                $rdf->assert_resource($five_prime_nuc_uri.$five_prime_edge, $is_paired_with_uri, $three_prime_nuc_uri.$three_prime_edge );
                $rdf->assert_resource($three_prime_nuc_uri.$three_prime_edge, $is_paired_with_uri, $five_prime_nuc_uri.$five_prime_edge );
                $rdf->assert_resource($five_prime_nuc_uri.$five_prime_edge, $has_glycosidic_bond_orientation_uri, $glycosidic_bond_orientation );
                $rdf->assert_resource($three_prime_nuc_uri.$three_prime_edge, $has_glycosidic_bond_orientation_uri, $glycosidic_bond_orientation );
            }
        }
        else {
            $logger->error("The edge description $edge[3] can not be parsed.");
        }
    }

    
    open (my $fh, ">", $rdf_filename);
    my $string = $rdf->serialize( format => 'rdfxml');
#    my $string = $rdf->serialize( format => 'turtle');
    print ($fh $string);
    close($fh);

    return;
}

###############################################################################
##
## &insert_edge( $rdf, $nuc_uri, $edge_type );
## * 
## * 
## * 
##
###############################################################################

sub insert_edge {
    my $rdf_model = shift;
    my $nuc_uri = shift;
    my $edge_type = shift;
    my $edge_uri = $nuc_uri.$edge_type;
    my $logger = get_logger();

    # edge types
    my $hoogsteen_edge_uri = 'bioinf:HoogsteenEdge';
    my $sugar_edge_uri = 'bioinf:SugarEdge';
    my $std_watson_crick_edge_uri = 'bioinf:StandardWatsonCrickEdge';
    my $non_std_watson_crick_edge_uri = 'bioinf:NonStandardWatsonCrickEdge';

    # object properties
    my $has_edge_uri = 'bioinf:hasEdge';

    if ( $edge_type eq 'W') {
        $rdf_model->assert_resource($edge_uri, 'rdf:type', $std_watson_crick_edge_uri);
        $rdf_model->assert_resource($nuc_uri, $has_edge_uri, $edge_uri);
    } elsif ($edge_type eq 'E') {
        $rdf_model->assert_resource($edge_uri, 'rdf:type', $non_std_watson_crick_edge_uri);
        $rdf_model->assert_resource($nuc_uri, $has_edge_uri, $edge_uri);
    } elsif ($edge_type eq 'H') {
        $rdf_model->assert_resource($edge_uri, 'rdf:type', $hoogsteen_edge_uri);
        $rdf_model->assert_resource($nuc_uri, $has_edge_uri, $edge_uri);
    } elsif ($edge_type eq 'S') {
        $rdf_model->assert_resource($edge_uri, 'rdf:type', $sugar_edge_uri);
        $rdf_model->assert_resource($nuc_uri, $has_edge_uri, $edge_uri);
    } else {
        $logger->error("Edge type $edge_type for $nuc_uri unknown.");
        exit 0;
    }
    return $rdf_model;
}



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

###############################################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################################

sub checkFiles {
    my $test_file = shift;
    my $checked_file = "";
    my $logger = get_logger();
    # Check if files are readable
    if ( -r $test_file){
        $checked_file = $test_file;
        $logger->info("$test_file is readable.");
    } else {
        $logger->error("$test_file is not readable.");
        exit;
    }
    return $checked_file;
}


# Output should be a RDAT file

__END__

=h1

This script requires the modules RDF::Helper, Moose, Namespace::Autoclean (???)

