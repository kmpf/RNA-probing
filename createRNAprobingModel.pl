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
#$module_dir =~ s/scripts$/RNAprobing/g;
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
my $man = 0;
my $off_file = "";
my $rdat_file = "";
my $rnaup_file = "";
my $rdf_out = 0;
my $owl_file =file(dirname(__FILE__), "RNAprobing.owl");
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "blastResults=s" => \$blast_result, # mandatory
    "csv" => \$csv, # optional
    "off=s" => \$off_file, # mandatory
    "rdat=s" => \$rdat_file, # mandatory
    "rnaup=s" => \$rnaup_file, # optional
    "rdf" => \$rdf_out, # optional
    "owl=s"=> \$owl_file, 
    "verbose|v+" => \$verbose, # optional
    );

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( $off_file eq "" || $rdat_file eq "" || $blast_result eq "" ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
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

################################################################################
#                 
# Files to objects
#                 
################################################################################

# check if input files exist
$blast_result = &checkFiles($blast_result) if ( $blast_result ne "" );
$off_file = &checkFiles($off_file) if ( $off_file ne "" );
$rdat_file = &checkFiles($rdat_file) if ( $rdat_file ne "" );
$owl_file = &checkFiles($owl_file) if ( $owl_file ne "" );

# creation of RDAT file object and information extraction
my $rdat_object     = RNAprobing::RDATFile->new($rdat_file);
my $rdat_filename   = fileparse($rdat_object->filename());
my $rdat_seq_startpos = $rdat_object->seq_startpos();
#my $rdat_reactivity = $rdat_object->data()->reactivity();
#my $rdat_seqpos_reactivity = $rdat_object->seqpos_reactivity_map();
my @rdat_seq        = split(//, $rdat_object->sequence());
my @rdat_seqpos     = @{$rdat_object->seqpos()};

# creation of OFF file object and information extraction
my $off_object      = RNAprobing::OFFFile->new($off_file);

# creation of BLASTresult object and information extraction
my $blast_result_object = RNAprobing::BLASTresult->new($blast_result);

# creation of RNAup file object and information extraction
my $rnaup_object    = RNAprobing::RNAupFile->new($rnaup_file);

################################################################################
#                 
# Test if the BLAST mapping is correct
#                 
################################################################################

# this is the array index to identify the correct line in the BLAST output
my $blast_index = -1;

for ( my $i = 0; $i < scalar( @{$blast_result_object->subject_id()} ); $i++ ) {
    $blast_index = $i 
    if ( $blast_result_object->subject_id()->[$i] eq $off_object->fasta_id() );
}

if ( $blast_index == -1 ) {
    die("No line in ".$blast_result_object->filename().
    " found with query ID: ".$off_object->fasta_id() );
}

my $query_id = $blast_result_object->query_id()->[$blast_index];

if ( $query_id eq $rdat_filename ) {
    $logger->debug( "The subject ID ".$query_id." from BLAST output 
        is identical with the .radt filename ".$rdat_filename.".");
} else {
    die("The subject ID ".$query_id." in "
           .$blast_result_object->filename()
           ." is not identical with the .rdat filename "
           .$rdat_filename." in ".$rdat_object->filename() );
}

if ( $blast_result_object->mismatches()->[$blast_index] != 0 ||
     $blast_result_object->gap_openings()->[$blast_index] != 0 ) {
    die("The mapping of ".$blast_result_object->subject_id()->[$blast_index]
    ."onto ".$blast_result_object->query_id()->[$blast_index]
    ." contains gaps or mismatches.");
}
                 
###############################################################################
#
# Map the sequences according to the BLAT/BLAST output
#                 
###############################################################################

# get the query and subject start and stop values from the BLAT BLAST output
my $query_start_rdat     = $blast_result_object->query_start()->[$blast_index];
my $query_end_rdat       = $blast_result_object->query_end()->[$blast_index];
my $subject_start_off    = $blast_result_object->subject_start()->[$blast_index];
my $subject_end_off      = $blast_result_object->subject_end()->[$blast_index];

$logger->debug(".rdat query start: $query_start_rdat");
$logger->debug(".rdat query end: $query_end_rdat");
my $l = $query_end_rdat - $query_start_rdat + 1;
$logger->debug("Length of match in RDAT: ". $l );
$logger->debug(".off subject start: $subject_start_off");
$logger->debug(".off subject end: $subject_end_off");
$l = $subject_end_off - $subject_start_off + 1;
$logger->debug("Length of match in OFF: ". $l);

## calculate the start and end positions in the rdat_offset number system
# "- 1" is needed cause BLAST output has 1-indexed sequence positions
# such that match start at position 1 gives a addition of 0
my $seq_startpos = $rdat_object->seq_startpos() + $query_start_rdat - 1;
my $seq_endpos = $rdat_object->seq_startpos() + $query_end_rdat - 1;

# If the match is not fully covered with probing data
# fit the *_start and *_end to the probing data.
# And adjust:
# - $query_start_rdat
# - $query_end_rdat
# - $subject_start_off
# - $subject_end_off

if ( $rdat_seqpos[0] > $seq_startpos ) {
    my $diff = $rdat_seqpos[0] - $seq_startpos;
    $seq_startpos += $diff;
    $subject_start_off += $diff;
    $query_start_rdat += $diff;
}

if ( $rdat_seqpos[$#rdat_seqpos] < $seq_endpos ) {
    my $diff = $seq_endpos - $rdat_seqpos[$#rdat_seqpos];
    $seq_endpos -= $diff;
    $subject_end_off -= $diff;
    $query_end_rdat -= $diff;
}

$logger->debug("New .rdat query start: $query_start_rdat");
$logger->debug("New .rdat query end: $query_end_rdat");
$l = $query_end_rdat - $query_start_rdat + 1;
$logger->debug("New length of match in RDAT: ". $l);
$logger->debug("New .off subject start: $subject_start_off");
$logger->debug("New .off subject end: $subject_end_off");
$l = $subject_end_off - $subject_start_off + 1;
$logger->debug("New length of match in OFF: ". $l);
$logger->debug("RDAT start position of reactivity: $seq_startpos");
$logger->debug("RDAT end position of reactivity: $seq_endpos");


###############################################################################
#
# Creating hashes for all necessary information
#   - Keys are the 1-indexed positions relating to the sequence contained in 
#     OFF file
#   - Values are scalars
#     1. %pos_reac -> reactivity at position $key
#                  -> is a hash of hashes
#     2. %pos_seq -> nucleotide at position $key
#                 
###############################################################################
my %pos_reac = ();
my %pos_seq = ();

# create one hash per $index
foreach my $index ( @{$rdat_object->data()->indices()} ) {
    $pos_reac{$index} = {}
}

for (my $i = 0; $i < ($seq_endpos - $seq_startpos);  $i++ ) {
    # fill hash with keys eqaul sequence position in OFF file and
    # related reactivity from RDAT file
    $logger->info("$i. OFF: ".$off_object->sequence_one_indexed_map()->{$subject_start_off + $i}." ".($subject_start_off + $i) );
    $logger->info("$i. RDAT: ".$rdat_object->offset_sequence_map()->{$seq_startpos + $i}." ".($seq_startpos + $i) );
    if ( $off_object->sequence_one_indexed_map()->{$subject_start_off + $i} eq
        $rdat_object->offset_sequence_map()->{$seq_startpos + $i}) {
        $pos_seq{$subject_start_off + $i} = 
            $rdat_object->offset_sequence_map()->{$seq_startpos + $i};
    } else {
    die("Sequence mapping went wrong.");
    }
    foreach my $index ( @{$rdat_object->data()->indices()} ) {
        $pos_reac{$index}{$subject_start_off + $i} = 
            $rdat_object->seqpos_scaled_reactivity_map($index)->{$seq_startpos + $i};
    }
}
###############################################################################
#
# Calling the subs with the appropriate information
#                 
###############################################################################


#if ( $csv ){
#    $logger->info("CSV file containing RNAup and reactivity is going to be created.");
#    my $csv_filename = fileparse($off_object->filename());
#    $csv_filename =~ s/off$/csv/g;
#    &create_csv_file( $csv_filename, \%querypos_rnaup_energies, \%querypos_rnaup_prob,
#        $querypos_reactivity, $querypos_reactivity_prob, $query_start_rdat, $query_end_rdat,
#        \@match_in_off, \@match_in_rdat );
#}


if ( $rdf_out ){
    $logger->info("RDF/XML model is going to be created.");
    &generate_rdf_model(\%pos_seq, \%pos_reac, $subject_start_off,
            $subject_end_off, $off_object, $rdat_object, $owl_file); 
}



###############################################################################
##              
##              Subroutine section
##
###############################################################################

###############################################################################
#
#     &create_csv_file( $csv_filename, \%querypos_rnaup_energies, \%querypos_rnaup_prob,
#        $querypos_reactivity, $querypos_reactivity_prob, $query_start_rdat, $query_end_rdat);
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
    my $query_start_rdat = shift;
    my $query_end_rdat = shift;
    my $match_in_off = shift; 
    my $match_in_rdat = shift;

    my $output = "";


    # generate header line
    $output = "pos,rnaup_energy,rnaup_prob,off_seq,rdat_seq";
    for (my $index = 0; $index < scalar(@{$querypos_reactivity} ); $index++ ) {
        $output .= ",reactivity_$index,reactivity_prob_$index";
    }
    $output .= "\n";

    # generate the data lines
    for ( my $i = $query_start_rdat ; $i <= $query_end_rdat; $i++ ) {
        $output .= $i.",";
        $output .= $querypos_rnaup_energies->{$i}.",";
        $output .= $querypos_rnaup_prob->{$i}.",";
        $output .= $match_in_off->[$i - 1].",";
        $output .= $match_in_rdat->[$i - 1].",";
        for (my $j = 0; $j < scalar( @{$querypos_reactivity} ); $j++ ) {
            $output .= $querypos_reactivity->[$j]->{$i}.",";
            $output .= $querypos_reactivity_prob->[$j]->{$i};
            $output .= "," if ( $j + 1 != scalar(@{$querypos_reactivity}) );
        }
        $output .= "\n";
    }
    open(my $csv, ">", $csv_filename) or die "Couldn't open file $csv_filename. Error: $!";
    print($csv $output);
    close($csv);
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
## &generate_rdf_model(\%pos_seq, \%pos_reac, $query_start_rdat, 
##                     $query_end_rdat, $off_object, $rdat_object)
##
###############################################################################

sub generate_rdf_model {
    my ($pos_seq, $pos_reac, $seq_start, $seq_end, $off_object,
    $rdat_object, $owl_file) = @_;
    my $logger = get_logger();
#    $logger->debug(Dumper($pos_seq));
#    $logger->debug(Dumper($pos_reac));
#    $logger->debug(Dumper($seq_start));
#    $logger->debug(Dumper($seq_end));
#    $logger->debug(Dumper($off_object));
#    $logger->debug(Dumper($rdat_object));
#    $logger->debug(Dumper($owl_file));

    my $rdat_name = fileparse($rdat_object->filename());
    $rdat_name =~ s/\.rdat$//g;
    $logger->debug($rdat_name);
    my $off_name = fileparse($off_object->filename());
    $off_name =~ s/\.off$//g;
    $logger->debug($off_name);

    my $rdf_id = $rdat_name."-".$off_name;
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


    # For each elemnet in @{$rdat_object->mutpos()} create a new RDF::Helper
    # object, to get only one experimentalUnpairedProbability per nucleotide

    # Save a RDF::Helper objects in @rdf for each "WT" element in
    # @{$rdat_object->mutpos()} 
    my @rdf_models = [];

    # Configure RDF::Helper

    for (my $i = 0; $i < scalar( @{$rdat_object->mutpos()} ); $i++ ) {
        next if ( $rdat_object->mutpos()->[$i] ne "WT" );

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
        $parser->parse_file_into_model( $base_uri, $owl_file->stringify, $model );
        # populate the RDF graph

        foreach my $querypos ( keys(%{$pos_seq}) ) {
        # each nucleotide gets its node and ... 
            my $nuc_uri = 'bioinf:'.$rdf_id.'/'.$pos_seq{$querypos}.$querypos;
            $logger->debug("Nucleotide: $nuc_uri");

        # ... gets connected to its neighbours ...
            # isThreePrimeOf

            if ( $querypos < $seq_end-1 ) {
                $logger->info($querypos." < ".$seq_end);
                my $three_prime_nuc_uri = 'bioinf:'.$rdf_id.'/'
                    .$pos_seq{$querypos+1}.($querypos+1);
                $rdf->assert_resource($three_prime_nuc_uri,
                      $is_three_prime_of_uri, $nuc_uri);
            }

            # isFivePrimeOf

            if ( $querypos > $seq_start ) {
                my $five_prime_nuc_uri =  'bioinf:'.$rdf_id.'/'
                    .$pos_seq{$querypos-1}.($querypos-1);
                $rdf->assert_resource($five_prime_nuc_uri,
                      $is_five_prime_of_uri, $nuc_uri);
            }

        # ... and gets a type ...
            # type of nucleotide at position $querypos

            $rdf->assert_resource($nuc_uri, 'rdf:type', $adenosine_uri) 
                if ( $pos_seq{$querypos} eq 'A');
            $rdf->assert_resource($nuc_uri, 'rdf:type', $cytidine_uri)
                if ( $pos_seq{$querypos} eq 'C');
            $rdf->assert_resource($nuc_uri, 'rdf:type', $guanosine_uri)
                if ( $pos_seq{$querypos} eq 'G');
            $rdf->assert_resource($nuc_uri, 'rdf:type', $uridine_uri)
                if ( $pos_seq{$querypos} eq 'U');
            $rdf->assert_resource($nuc_uri, 'rdf:type', $ribonucleotide_uri)
                if ( $pos_seq{$querypos} eq 'N');

            # ... and gets ALL its edges.
            # insert Watson-Crick, Hoogsteen or Sugar edges
            &insert_edge( $rdf, $nuc_uri, 'W' );
            &insert_edge( $rdf, $nuc_uri, 'H' );
            &insert_edge( $rdf, $nuc_uri, 'S' );

            # RNAupProbability for nucleotide at $querypos

#            my $rna_up = $rdf->new_literal($querypos_rnaup_prob->{$querypos},
#                       '', 'xsd:decimal');
#            $rdf->assert_literal($nuc_uri, $has_rnaup_unpaired_prob_uri, $rna_up);

            # ExperimentalUnpairedProbability for nucleotide at $querypos

            # add 1 because $i is 0-indexed but the keys are 1-indexed
            my $reactivity = $pos_reac{$i + 1}->{$querypos};
            my $react_prob = $rdf->new_literal($reactivity, '', 'xsd:decimal');
            $rdf->assert_literal($nuc_uri, $has_exp_unpaired_prob_uri
                 , $react_prob);
        }
        
    # insert all Leontis-Westhof base pairs and tertiary interactions

        for ( my $j = 0; $j < scalar( @{ $off_object->edges() } ); $j++) {

            my @edge = @{$off_object->edges()->[$i]};
            next if ($edge[1] < $seq_start || $seq_end < $edge[2] );

            if ( $edge[0] ne "LW" ) {
                $logger->info("Unknown notation $edge[0] trying next line.");
                next;
            }
            my $five_prime_nuc_uri = 'bioinf:'.$rdf_id.'/'
                .$pos_seq{$edge[1]}.$edge[1];
            my $three_prime_nuc_uri = 'bioinf:'.$rdf_id.'/'
                .$pos_seq{$edge[2]}.$edge[2];
            if ( $edge[3] =~ /^[ct!][WEHST][WEHST]$/ ) {
                $edge[3] =~ /^([ct!])([WEHST])([WEHST])$/;
                my $orientation = $1;
                my $five_prime_edge = $2;
                my $three_prime_edge = $3;

                # insert TertiaryInteraction
                if ( $orientation eq '!' && $five_prime_edge eq 'T' 
                     && $three_prime_edge eq 'T') {
                    $rdf->assert_resource($five_prime_nuc_uri
                      , $has_tertiary_interaction_with_uri
                      , $three_prime_nuc_uri );
                } elsif ($orientation eq 'c' || $orientation eq 't' ) {
                    my $glycosidic_bond_orientation = ();
                    if ( $orientation eq "c" ){
                        $glycosidic_bond_orientation = 'bioinf:Cis';
                    } elsif ($orientation eq "t" ){
                        $glycosidic_bond_orientation = 'bioinf:Trans';
                    }
                    next unless ( $five_prime_edge =~ /WEHS/ && 
                        $three_prime_edge =~ /WEHS/ );
                    $rdf->assert_resource($five_prime_nuc_uri.$five_prime_edge,
                      $is_paired_with_uri,
                      $three_prime_nuc_uri.$three_prime_edge);
                    $rdf->assert_resource($three_prime_nuc_uri.$three_prime_edge,
                      $is_paired_with_uri,
                      $five_prime_nuc_uri.$five_prime_edge );
                    $rdf->assert_resource($five_prime_nuc_uri.$five_prime_edge,
                      $has_glycosidic_bond_orientation_uri,
                      $glycosidic_bond_orientation );
                    $rdf->assert_resource($three_prime_nuc_uri.$three_prime_edge,
                      $has_glycosidic_bond_orientation_uri,
                      $glycosidic_bond_orientation );
                }
            }
            else {
                $logger->error("The edge description $edge[3] can not be parsed.");
            }
        }

#        my $anno = $rdat_object->annotation_data()->[$i];
#        print Dumper($anno );
#        print Dumper($anno->chemicals() );

       my $rdf_filename = $rdf_id.".$i.rdf";
       $logger->info("RDF file: $rdf_filename");
       open (my $rdf_fh, ">", $rdf_filename) or die("Can't open $rdf_filename.");
       my $string = $rdf->serialize( format => 'rdfxml');
       # my $string = $rdf->serialize( format => 'turtle');
       # $logger->info("$string");
       print $rdf_fh $string;
       close $rdf_fh;
       $logger->info("Wrote $rdf_filename");

    }

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
    my ($rdf_model, $nuc_uri, $edge_type) = @_;
    my $edge_uri = $nuc_uri.$edge_type;
    my $logger = get_logger();

    # edge types
    my $hoogsteen_edge_uri = 'bioinf:HoogsteenEdge';
    my $sugar_edge_uri = 'bioinf:SugarEdge';
    my $watson_crick_edge_uri = 'bioinf:WatsonCrickEdge';

    # object properties
    my $has_edge_uri = 'bioinf:hasEdge';

    if ( $edge_type eq 'W' || $edge_type eq 'E') {
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
            $pos_reac{ $rdat_seqpos_querypos->{$key} } = $rdat_seqpos_reactivity->[$j]->{$key};
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

=head1 NAME

createRNAprobingModel.pl - combines the structure information of an OFF file with the probing information of a RDAT file, given a mapping in BLAST format 

=head1 SYNOPSIS

rdat2gel.pl [options] -f,--file F<rdat-file>

=head1 DESCRIPTION

This script creates a RDAT file, containing the results of I<in silico> probing experiment.
To perform a probing reaction it needs to be provided with a RNA sequence stored in a FASTA file and the reactivity rules in a special file format that looks like this:


=head1 OPTIONS

=over 8
    "rnaup=s" => \$rnaup_file, #optional

=item B<-h, --help>       

Display help message

=item B<-m, --man>

Display whole man page

=item B<--blastResults>

File containing the results of an alignment of the sequences from the RDAT file (see --rdat) and the OFF file (see --off).

=item B<--off>

OFF file containing the sequence and structure information of a RNA (normally extracted from NDB RNAML files)

=item B<--rdat>

RDAT file containing the probing information for a RNA (normally downloaded from the RMDB)

=item B<--rdf>

Flag that if set leads to the generation of the RDF model

=item B<--owl>

OWL file which will be included into the RDF model. By default uses "RNAprobing.owl" which is contained in the RNAprobing repository.

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (verbosest if used 3 or more times) 

=back

This script requires the modules RDF::Helper, Moose, Namespace::Autoclean (???)


=cut
