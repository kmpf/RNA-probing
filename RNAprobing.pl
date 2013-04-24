#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: RNAprobing.pl
#
#        USAGE: ./RNAprobing.pl  
#
#  DESCRIPTION: this script generates a RDAT file given a RNA sequence in FASTA
#               format and a file describing the reactivities of a chemical 
#
#      OPTIONS: -h, --help     Display help message
#               --fasta        Fasta file containing RNA sequence to be probed
#               --chemical     Chemical file describing the reactivities of the
#                              probing reagent
#               -v, --verbose  
#               
# REQUIREMENTS: RNAlib Perl bindings
#         BUGS: 
#        NOTES: 
#       AUTHOR: Christoph Kaempf (CK), kaempf@bioinf.uni-leipzig.de
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 12.09.2012 12:53:10
#     REVISION: 
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;
use RNA;
my $module_dir = dirname(__FILE__);
push(@INC, $module_dir);

################################################################################
#
# Options section
#
################################################################################
my $help = 0;
my $man = 0;
my $rdat_file;
my $fasta_file;
my $chemical_file;
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "fasta=s" => \$fasta_file,
    "chemical=s" => \$chemical_file,
    "verbose|v+" => \$verbose);

if ( $help || !( defined $fasta_file && defined $chemical_file) ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
} elsif ($man) {
    pod2usage( { -verbose => 2});
}


###############################################################################
#                 
# Logger initiation
#                 
###############################################################################
my $this_file = __FILE__;
my $log4perl_conf = file(dirname($this_file), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the loggerperl file
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");

# require RNAprobing classes just after logger initialization
require RNAprobing::Chemical;
require RNAprobing::RDATFile;
require RNAprobing::RDATFile::Annotation;
require RNAprobing::RDATFile::Name;
require RNAprobing::OFFFile;

###############################################################################
#                 
# Program logic
#                 
###############################################################################

# Input needed is a FASTA file and a reactivity file
my $fasta = RNAprobing::OFFFile->new($fasta_file);
$logger->debug("Loaded Fasta file ".$fasta_file);
my $chemical = RNAprobing::Chemical->new($chemical_file);

my $sample_size = 2; # number of stochastic sampled RNA structures to probe 
                     # could be an option
my $seq = $fasta->sequence();
my @probing_profile = (0) x length($seq);
my $length = length(join("", @probing_profile));
my @structures = &stochastic_sampling($seq, $sample_size);
my @structure_description = &db2sd(@structures);
foreach (@structures){
    print( "Structure: ".$_."\n");
}
foreach (@structure_description){
    print( "Structure description: ".$_."\n");
}

@probing_profile = &simulate_probing(\@structure_description, \@probing_profile,
                                     $seq, $chemical);

# === Result output ===
# Should be logged instead of printed

print("=== Results ===\n");
for (my $i = 0; $i < scalar(@structures); $i++) {
    print("$i. Structure:\n$structures[$i]\n$structure_description[$i]\n");
}
print(join("", @probing_profile)."\n");

# === Assemble a RDAT file ===

my $rdat_file_name = $fasta->fasta_id()."_".$chemical->probe_name().".rdat";
$logger->info("++++ ".$rdat_file_name." ++++");
my $rdat_out = RNAprobing::RDATFile->new($rdat_file_name);

my $construct_name = $fasta->fasta_id()."_in_silico_probed_using_".$chemical->probe_name();

my $construct_object = $rdat_out->name($construct_name);

$construct_object->sequence($seq);
my ($struct, $mfe) = RNA::fold($seq);  # predict mfe structure of $seq
$construct_object->structure($struct);
$construct_object->offset(0);
my @seqpos = (1 .. length($seq));
$construct_object->seqpos( \@seqpos );
$construct_object->reactivity(\@probing_profile, 1); # 1 <- index
print Dumper($construct_object);
$rdat_out->write_file();

################################################################################
################################################################################
##
##                              Subroutines
##
################################################################################
################################################################################

################################################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $ERROR
## -- 1 => $WARN
## -- 2 => $INFO
## -- >2 => $DEBUG
## 
################################################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    $logger->info("Verbosity level: $verbose");
    # print Dumper($logger);
    SELECT:{
	    if ($verbose == 0){$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
	    if ($verbose == 1){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    else { $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
    }
    return $logger;
}

################################################################################
#
#             Stochastic sampling of RNA secondary structures
#
################################################################################

sub stochastic_sampling {
    my ($seq, $sample_size) = @_;
    # compute partition function and pair pobabilities
    my $structure;
    $RNA::st_back = 1;
    my $gfe = RNA::pf_fold($seq, $structure);
    my @structures;

    for (my $i = 0; $i < $sample_size; $i++){
        push( @structures, RNA::pbacktrack($seq) );
    }

    return @structures;
}

################################################################################
#
#           Convert Dot-Bracket into structure description 
#
################################################################################

sub db2sd {

my (@structures) = @_;
my @structure_description = ();
foreach (@structures) {

    my @db = split('', $_);

    my @dots = (-1);
    my @opening_br = (-1);

    # at first we expect every nucleotide to be unpaired "U"
    my @struc_dec = ("U") x length($_);

    # fill arrays with positions
    for (my $i = 0; $i < scalar(@db); $i++) {
        push(@dots, $i) if ( $db[$i] eq "." );
        push(@opening_br, $i) if ( $db[$i] eq "(" );

        # closing bracket found, so lets classify enclosed unpaired nucleotides
        # if there are any
        # if clause leads to errors, if @db contains substrings like "()"
        if ( $db[$i] eq ")" && $dots[$#dots] > $opening_br[$#opening_br] ) {

            # declare enclosing bracket positions
            my $op_br_pos = pop(@opening_br);
            my $cl_br_pos = $i;

            # Declare paired nucleotides
            $struc_dec[$op_br_pos] = "P";
            $struc_dec[$cl_br_pos] = "P";


            my @enclosed_dots = ();
            # count the number of stems enclosed by the unpaired nucleotides
            my $enclosed_stems = 0;

            # collect all enclosed nucleotides in @enclosed_dots
            while ( $dots[$#dots] > $op_br_pos ) {
                push(@enclosed_dots, pop(@dots));
                # stop if you reached the last enclosed nucleotide
                last if ($dots[$#dots] < $op_br_pos);
                # if found non-continous numbers we found an enclosed stem
                if ( $dots[$#dots]+1 < $enclosed_dots[$#enclosed_dots]) {
                    $enclosed_stems++;
                }
            }
#           print("Enclosed unpaired nucleotides: "                    .join(",", @enclosed_dots)."\n");

            # Hairpin or bulge detected
            if ( $enclosed_stems == 0 ) {
                # Hairpin detected if the enclosed dots reach from opening to
                # closing bracket
                if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {            
                    foreach (@enclosed_dots) { $struc_dec[$_] = "H" }
#                    print("Hairpin found.\n");
                }
                # bulge detected if the enclosed dots just touch one bracket
                elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                        $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "B" }
#                    print("Bulge found.\n");
                }
            }
            # Multi loop with two included stems or an interior loop detected
            elsif ( $enclosed_stems == 1 ) {
                # Interior loop detected if the enclosed dots reach from opening
                # to closing bracket
                if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "I" }
#                    print("Interior loop found.\n");
                }
                # Multi loop with two stems detected if the enclosed dots just
                # touch one bracket
                elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                        $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
#                    print("Multi loop found.\n");
                }
            }
            # Multi loop detected if more than one stem is enclosed
            elsif( $enclosed_stems > 1 ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
#                    print("Multi loop found.\n");
            }
        }
        # No enclosed unpaired nucleotides so lets declare the base pair
        elsif ( $db[$i] eq ")" ) {
            $struc_dec[$i] = "P";
            my $rel_open_br = pop(@opening_br);
            $struc_dec[$rel_open_br] = "P";
        }
    }
    my $str_desc = join("",@struc_dec);
#    print("$_\n$str_desc\n");
    push(@structure_description, join("",@struc_dec));
}
return @structure_description;

}

################################################################################
#
#           Probe sequence and secondary structure
#
################################################################################

sub simulate_probing {
    my ($structure_description, $probing_profile, $seq, $chemical) = @_;

    for ( my $i = 0; $i < scalar(@{$chemical->probe_reac()}); $i++ ) {
        my $prob_reac = ${$chemical->probe_reac()}[$i];
        my $prob_seq = ${$chemical->probe_seq()}[$i];
        my $prob_str = ${$chemical->probe_str()}[$i];
        my $prob_cut = ${$chemical->probe_cut()}[$i];

        # create regex from probing sequence
        $prob_seq =~ s/N/[ACGTU]/g;

        # create regex from probing structure
        $prob_str =~ s/U/[HBIMn]/g;
        $prob_str =~ s/n/U/g;

        foreach my $str_desc (@{$structure_description}) {
            # Find the modification/cut point
            $prob_cut =~ /\|/;
            my $cut_pos = $-[0];


            my @seq_matches = &match_all_positions($prob_seq, $seq, $cut_pos);
            my @str_matches = &match_all_positions($prob_str, $str_desc, $cut_pos);


            if (scalar(@seq_matches) > scalar(@str_matches) ) {
                my $more_matches = join(",",@seq_matches);
                foreach (@str_matches) {
                    ${probing_profile}[$_] += $prob_reac if ($more_matches =~ /,$_,/);
                }
            }
            else {
                 my $more_matches = join(",",@str_matches);
                foreach (@seq_matches) {
                    ${probing_profile}[$_] += $prob_reac  if ($more_matches =~ /,$_,/);
                }
            }

    #        print(join("",@{$probing_profile})."\n");
    #        print("$seq\n");
        }
    }
    return @{$probing_profile};
}

sub match_all_positions {
    my ($regex, $string, $cut_pos) = @_;
    my @ret;
    while ($string =~ /(?=$regex)/g) {
        push(@ret,  $-[0] + $cut_pos);
    }
    return @ret
}

__END__

=head1 NAME

sample - Using GetOpt::Long and Pod::Usage

=head1 SYNOPSIS

RNAprobing.pl [options] --fasta F<fasta-file> --chemical F<chemical-file>

=head1 DESCRIPTION

This script creates a RDAT file, containing the results of I<in silico> probing experiment.
To perform a probing reaction it needs to be provided with a RNA sequence stored in a FASTA file and the reactivity rules in a special file format that looks like this:

=head2 EXAMPLE REACTIVITY FILE 

    > U2 nuclease # Name of reagent

    1.0 # modification strength
    A # specific RNA sequence
    U # specific sec. structure
    | # modification point
    
    0.9
    G # specific RNA sequence
    U # specific sec. structure
    | # modification point
    
    0.2
    C # specific RNA sequence
    U # specific sec. structure
    | # modification point
    
    0.1
    U # specific RNA sequence
    U # specific sec. structure
    | # modification point


=head1 OPTIONS

=over 8

=item B<-h, --help>       

Display help message

=item B<--fasta>

Fasta file containing RNA sequence to be probed

=item B<--chemical>

Chemical file describing the reactivities of the probing reagent

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (highest level if used 3 or more times) 

=back

=cut
