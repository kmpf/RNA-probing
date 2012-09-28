#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: RNAprobing.pl
#
#        USAGE: ./RNAprobing.pl  
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
################################################################################
my $blast_result = "";
my $csv = 0;
my $help = 0;
my $off_file = "";
my $rdat_file = "";
my $rnaup_file = "";
my $rdf_out = 0;
my $owl_file ="";
my $verbose = 0;

GetOptions(
    "blastResults=s" => \$blast_result,
    "csv" => \$csv,
    "off=s" => \$off_file,
    "rdat=s" => \$rdat_file,
    "rnaup=s" => \$rnaup_file,
    "rdf" => \$rdf_out,
    "owl=s"=> \$owl_file,
    "verbose|v+" => \$verbose);

if ( $help || $off_file eq "" || $rdat_file eq "" || $rnaup_file eq "" || $blast_result eq "" ){
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

###############################################################################
#                 
# Files to objects
#                 
###############################################################################

$blast_result = &checkFiles($blast_result);
$off_file = &checkFiles($off_file);
$rdat_file = &checkFiles($rdat_file);
$owl_file = &checkFiles($owl_file);

# creation of RDAT file object and information extraction
my $rdat_object     = RNAprobing::RDATFile->new($rdat_file);
my $rdat_filename   = fileparse($rdat_object->filename());
my $rdat_offset     = $rdat_object->offset();
my $rdat_reactivity = $rdat_object->reactivity();
my $rdat_seqpos_reactivity = $rdat_object->seqpos_reactivity_map();
my @rdat_seq        = split(//, $rdat_object->sequence());
my @rdat_seqpos     = @{$rdat_object->seqpos()};

# creation of OFF file object and information extraction
my $off_object      = RNAprobing::OFFFile->new($off_file);
my @off_seq         = split(//, $off_object->sequence());

# creation of BLASTresult object and information extraction
my $blast_result_object = RNAprobing::BLASTresult->new($blast_result);

# creation of RNAup file object and information extraction
my $rnaup_object    = RNAprobing::RNAupFile->new($rnaup_file);

###############################################################################
#                 
# Test if the files map correctly
#                 
###############################################################################

# this is the array index to identify the correct line in the BLAST output
my $blast_array_index = -1;

for ( my $i = 0; $i < scalar( @{$blast_result_object->query_id()} ); $i++ ) {
    $blast_array_index = $i if ( ${$blast_result_object->query_id()}[$i] eq $off_object->fasta_id() );
}

if ( $blast_array_index == -1 ) {
    die "No line in ".$blast_result_object->filename()." found with query ID: ".$off_object->fasta_id()."\n";
}

if ( ${$blast_result_object->subject_id()}[$blast_array_index] eq $rdat_filename ) {
    $logger->debug("The subject ID ".${$blast_result_object->subject_id()}[$blast_array_index]." from BLAST output 
        is identical with the .radt filename ".$rdat_filename.".");
} else {
    $logger->error("The subject ID ".${$blast_result_object->subject_id()}[$blast_array_index]." in ".$blast_result_object->filename()
        ." is not identical with the .rdat filename ".$rdat_filename." in ".$rdat_object->filename() );
    exit 1;
}

                 
###############################################################################
#
# Map the sequences according to the BLAT/BLAST output
#                 
###############################################################################

# get the query and subject start and stop values from the BLAT BLAST output
my $off_query_start     = ${$blast_result_object->query_start()}[$blast_array_index];
my $off_query_end       = ${$blast_result_object->query_end()}[$blast_array_index];
my $rdat_subject_start  = ${$blast_result_object->subject_start()}[$blast_array_index];
my $rdat_subject_end    = ${$blast_result_object->subject_end()}[$blast_array_index];

$logger->debug(".off query start: $off_query_start");
$logger->debug(".off query end: $off_query_end");
$logger->debug(".rdat subject start: $rdat_subject_start");
$logger->debug(".rdat subject end: $rdat_subject_end");

# calculate the start and end positions in the rdat_offset number system
my $offset_subject_start = $rdat_offset + $rdat_subject_start;
my $offset_subject_end = $rdat_offset + $rdat_subject_end;

# If the match is not fully covered with probing data
# fit the *_start and *_end to the probing data.
# And adjust:
# - $off_query_start
# - $off_query_end
# - $rdat_subject_start
# - $rdat_subject_end

if ( $rdat_seqpos[0] > $offset_subject_start ) {
    $offset_subject_start = $rdat_seqpos[0];
    my $diff = ($offset_subject_start - $rdat_offset) - $rdat_subject_start ;
    $rdat_subject_start += $diff;
    $off_query_start += $diff;
}

if ( $rdat_seqpos[$#rdat_seqpos] < $offset_subject_end ) {
    $offset_subject_end = $rdat_seqpos[$#rdat_seqpos];
    my $diff = ($offset_subject_end - $rdat_offset) - $rdat_subject_end ;
    $rdat_subject_end += $diff;
    $off_query_end += $diff;
}

$logger->debug("New .rdat subject start: $rdat_subject_start");
$logger->debug("New .rdat subject end: $rdat_subject_end");
$logger->debug("New .off query start: $off_query_start");
$logger->debug("New .off query end: $off_query_end");
$logger->debug("Start position of reactivity: $offset_subject_start");
$logger->debug("End position of reactivity: $offset_subject_end");


if ( $rdat_subject_end - $rdat_subject_start != $off_query_end - $off_query_start ||
     $rdat_subject_end - $rdat_subject_start != $offset_subject_end - $offset_subject_start ||
     $off_query_end - $off_query_start != $offset_subject_end - $offset_subject_start ) {
    die $logger->error("Subject and Query sequences are of different length.");
}



# get the with BLAT matched sequences
my @match_in_off = @off_seq[($off_query_start - 1) .. ($off_query_end - 1)];
my @rnaup_values = @{$rnaup_object->rnaup_values()}[($off_query_start - 1) .. ($off_query_end - 1)];
my @match_in_rdat = @rdat_seq[($rdat_subject_start - 1) .. ($rdat_subject_end - 1)];

if ( uc(join("",@match_in_off)) ne uc(join("",@match_in_rdat)) ){
    $logger->error("Sequence in the RDAT file ".$rdat_object->filename());
    $logger->error(join("",@match_in_rdat));
    $logger->error("is unequal with sequence in OFF file ".$off_object->filename());
    $logger->error(join("",@match_in_off));
    die;
} 

###############################################################################
#
# Creating hashes for all necessary information
#   - Keys are almost always the position indices of the query sequence!
#                 
###############################################################################

#  %querypos_rnaup_energies - map position to calculated rnaup energies
my %querypos_rnaup_energies = ();
my %querypos_rnaup_prob = ();
for ( my ($i, $j) = ($off_query_start, 0) ; $i <= $off_query_end; ($i++, $j++) ) {
    $querypos_rnaup_energies{$i} = $rnaup_values[$j];
    $querypos_rnaup_prob{$i} = &calculate_probability_from_rnaup_energies($querypos_rnaup_energies{$i});
}


# %querypos_subjectpos - maps query position to subject position
my %rdat_seqpos_subjectpos = ();
# %subjectpos_querypos - maps subject position to query position
my %rdat_seqpos_querypos = ();
for (my ($i, $j) = ($offset_subject_start, 0); $i <= $offset_subject_end; ($i++, $j++) ) {
    $rdat_seqpos_querypos{ $i } = $off_query_start + $j;
    $rdat_seqpos_subjectpos{ $i } = $rdat_subject_start + $j;
}

# $querypos_reactivity - reference of an array with each element a hash mapping query positions to reactivity
my $querypos_reactivity = &create_querypos_reactivity( $rdat_seqpos_reactivity, \%rdat_seqpos_querypos );
my $querypos_reactivity_prob = &calculate_probability_from_experimental_reactivtity($querypos_reactivity);

if ( $csv ){
    $logger->info("CSV file containing RNAup and reactivity is going to be created.");
    my $csv_filename = fileparse($off_object->filename());
    $csv_filename =~ s/off$/csv/g;
    &create_csv_file( $csv_filename, \%querypos_rnaup_energies, \%querypos_rnaup_prob,
        $querypos_reactivity, $querypos_reactivity_prob, $off_query_start, $off_query_end);
}


if ( $rdf_out ){
    $logger->info("RDF/XML model is going to be created.");
    &generate_rdf_model( \@off_seq, $off_query_start, $off_query_end, \%querypos_rnaup_prob, $querypos_reactivity_prob, $rdat_object, $off_object, $owl_file );
}



###############################################################################
##              
##              Subroutine section
##
###############################################################################

###############################################################################
#
#     &create_csv_file( $csv_filename, \%querypos_rnaup_energies, \%querypos_rnaup_prob,
#        $querypos_reactivity, $querypos_reactivity_prob, $off_query_start, $off_query_end);
# Create .csv file containing
#   1. column: nucleotide position (numbering as in .rdat)
#   2. column: calculated melting energy by RNAup
#   3. column: probability calculated from RNAup energy
#   4. column: reactivity
#                 
###############################################################################

sub create_csv_file {
    my $csv_filename = shift;
    my $querypos_rnaup_energies = shift;
    my $querypos_rnaup_prob = shift;
    my $querypos_reactivity = shift;
    my $querypos_reactivity_prob = shift;
    my $off_query_start = shift;
    my $off_query_end = shift;


    open(CSV, ">$csv_filename") or die "Couldn't open file $csv_filename. Error: $!";

    # generate header line
    print(CSV "pos,rnaup_energy,rnaup_prob");
    for (my $index = 0; $index < scalar(@{$querypos_reactivity} ); $index++ ) {
        print(CSV ",reactivity_$index,reactivity_prob_$index");
    }
    print(CSV "\n");

    # generate the data lines
    for ( my $i = $off_query_start ; $i <= $off_query_end; $i++ ) {
        print( CSV $i."," );
        print( CSV $querypos_rnaup_energies->{$i}."," );
        print( CSV $querypos_rnaup_prob->{$i}."," );
        for (my $j = 0; $j < scalar( @{$querypos_reactivity} ); $j++ ) {
            print( CSV $querypos_reactivity->[$j]->{$i}.",");
            print( CSV $querypos_reactivity_prob->[$j]->{$i});
            print( CSV ",") if ( $j + 1 != scalar(@{$querypos_reactivity}) );
        }
        print( CSV "\n");
    }
    close(CSV);
}

###############################################################################
##
## &calculate_probability_from_rnaup_energies( $rnaup_energy, $temperature )
## - the output of RNAup is a list of energies in kcal/mol needed
##   to unpair a certain portion of the RNA structure
## - this function converts those energies into probabilities
##
###############################################################################

sub calculate_probability_from_rnaup_energies {
    # Boltzmann constant in kcal/mol/K
    my $BOLTZMANN_CONST = 0.0019872041;
    # temperature (default: 37 °C = 310.15 K)
    my $temperature = 310.15;
    # By RNAup predicted energy needed to unpair nucleotide(s)
    my $rnaup_energy = shift;
    # Temperature if available
    $temperature = shift if ( scalar(@_) >= 1 && $_[0] =~ /[\d\.]+/ );
    my $exponent = -( $rnaup_energy/($BOLTZMANN_CONST * $temperature) );
    Math::BigFloat->precision(-5);
    my $up_probability = Math::BigFloat->bexp($exponent);
    return $up_probability->bstr();

}

###############################################################################
##
## &calculate_probability_from_experimental_reactivtity( $querypos_reactivity)
## - the output of RNAup is a list of energies in kcal/mol needed
##   to unpair a certain portion of the RNA structure
## - this function converts those energies into probabilities
##
###############################################################################

sub calculate_probability_from_experimental_reactivtity {
    my $querypos_reactivity = shift;
    my @querypos_reactivity_prob = ();
    my $logger = get_logger();

    for (my $i = 0; $i < @{$querypos_reactivity}; $i++) {
        my $max_reactivity = max(values(%{$querypos_reactivity->[$i]}));
        $logger->debug($max_reactivity);
        my $min_reactivity = min(values(%{$querypos_reactivity->[$i]}));
        $logger->debug($min_reactivity);
        my $reactivity_span = $max_reactivity - $min_reactivity;
        foreach my $key (keys(%{$querypos_reactivity->[$i]})) {
            my $str = ( $querypos_reactivity->[$i]->{$key} - $min_reactivity ) / $reactivity_span;
            Math::BigFloat->precision(-5);
            my $prob = Math::BigFloat->new( $str );
            $querypos_reactivity_prob[$i]{$key} = $prob->bstr();
        }
    }
#    $logger->debug( Dumper(@querypos_reactivity_prob) );
    return \@querypos_reactivity_prob;
}

###############################################################################
##
## &generate_rdf_model( \@off_seq, $off_query_start, $off_query_end, \%querypos_rnaup_prob, $querypos_reactivity_prob, $rdat_object, $off_object )
## - 
## - 
##
###############################################################################

sub generate_rdf_model {
    my $off_seq = shift;
    my $off_query_start = shift;
    my $off_query_end = shift;
    my $querypos_rnaup_prob = shift;
    my $querypos_reactivity_prob = shift;
    my $rdat_object = shift;
    my $off_object = shift;
    my $owl_file = shift;
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
    $parser->parse_file_into_model( $base_uri, $owl_file, $model );

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

    my $has_exp_unpaired_prob_uri = 'bioinf:hasExperimentalUnpairedProbability';
    my $has_rnaup_unpaired_prob_uri = 'bioinf:hasRNAupUnpairedProbability';
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

        # RNAupProbability for nucleotide at $querypos

        my $rna_up = $rdf->new_literal($querypos_rnaup_prob->{$querypos}, '', 'xsd:decimal');
        $rdf->assert_literal($nuc_uri, $has_rnaup_unpaired_prob_uri, $rna_up);

        # ExperimentalUnpairedProbability for nucleotide at $querypos

        for (my $i = 0; $i < scalar( @{$rdat_object->mutpos()} ); $i++ ) {
            next unless ( $rdat_object->mutpos()->[$i] eq "WT" );
            my $react_prob = $rdf->new_literal($querypos_reactivity_prob->[$i]->{$querypos}, '', 'xsd:decimal');
            $rdf->assert_literal($nuc_uri, $has_exp_unpaired_prob_uri, $react_prob);
        }
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
        if ( $edge[3] =~ /^[ct!][WHST][WHST]$/ ) {
            $edge[3] =~ /^([ct!])([WHST])([WHST])$/;
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

    
    open (my $rdf_fh, ">", $rdf_filename)  or die("Can't open $rdf_filename.");
    my $string = $rdf->serialize( format => 'rdfxml');
#    my $string = $rdf->serialize( format => 'turtle');
    $logger->info("$string");
    print $rdf_fh $string;

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
    my $watson_crick_edge_uri = 'bioinf:WatsonCrickEdge';

    # object properties
    my $has_edge_uri = 'bioinf:hasEdge';

    if ( $edge_type eq 'W') {
        $rdf_model->assert_resource($edge_uri, 'rdf:type', $watson_crick_edge_uri);
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
## &create_querypos_reactivity( $rdat_seqpos_reactivity, \%rdat_seqpos_querypos )
## * returns an array with one entry per experiment
## * each entry contains a hash
## ** the hash contains a reactivity value for every position in the query
##
###############################################################################

sub create_querypos_reactivity {
    my $rdat_seqpos_reactivity = shift;
    my $rdat_seqpos_querypos = shift;
    my @querypos_reactivity = ();
    my $logger = get_logger();

    for (my $j = 0; $j < scalar( @{$rdat_seqpos_reactivity} ); $j++ ) {
        my %pos_reac = ();
        foreach my $key (keys %{$rdat_seqpos_querypos}) {
            $pos_reac{ $rdat_seqpos_querypos{$key} } = $rdat_seqpos_reactivity->[$j]->{$key};
        }
        push( @querypos_reactivity, \%pos_reac );        
    }
    return \@querypos_reactivity;
}

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

