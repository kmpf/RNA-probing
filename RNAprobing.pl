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
use Scalar::Util;
use Scalar::Util::Numeric;
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
my $dbn_file;
my $samples = 1000;
my $offset = 0;
my ($seqpos_begin, $seqpos_end);
my $verbose = 0;

GetOptions(
    "help|h" => \$help,
    "man|m" => \$man,
    "fasta|f=s" => \$fasta_file,
    "chemical|c=s" => \$chemical_file,
    "dbn|d=s" => \$dbn_file,
    "samples=i" => \$samples,
    "offset|o=i" => \$offset,
    "begin|b=i" => \$seqpos_begin,
    "end|e=i" => \$seqpos_end,
    "verbose|v+" => \$verbose);

if ( $help ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
} elsif ( !(defined $fasta_file || defined $dbn_file) 
          || (defined $fasta_file && defined $dbn_file) ) {
    pod2usage( { -verbose => 1,
                 -message => "Either --fasta or --dbn must be set. Not both:\n"});
} elsif ( !(defined $chemical_file) ) {
    pod2usage( { -verbose => 1,
                 -message => "Please set --chemical:\n"});
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
require RNAprobing::RDATFile::Data;
require RNAprobing::OFFFile;

###############################################################################
#                 
# Program logic
#                 
###############################################################################

# Input needed is a FASTA file and a reactivity file
my ($fasta, $chemical);
if (defined $fasta_file) {
    $fasta = RNAprobing::OFFFile->new($fasta_file);
    $fasta->read_file();
    $logger->debug("Loaded fasta file ".$fasta_file);
}
if (defined $dbn_file) {
    $fasta = RNAprobing::OFFFile->new($dbn_file);
    $fasta->read_file();
    $logger->debug("Loaded dot-bracket notation file ".$dbn_file);
}
if (defined $chemical_file) {
    $chemical = RNAprobing::Chemical->new($chemical_file);
    $chemical->read_file();
    $logger->debug("Loaded chemical file ".$chemical_file);
}
## Sanity check of input parameter: --samples
my $sample_size; # number of stochastically sampled RNA structures to probe 
if (! defined($samples)) {
    $logger->error("--samples not used. Please provide positive integer ".
                   "value via '--samples' option.");
    exit 1;
} elsif (Scalar::Util::looks_like_number($samples) && 
        Scalar::Util::Numeric::isint($samples) ) {
    if ($samples <= 0) {
        $logger->error("Sample size is set to ".$samples.". Please provide ".
                       "positive integer value via '--samples' option.");
        exit 1;
    } else {
        $sample_size = $samples;
        $logger->debug("--samples has value: ".$sample_size);
    }
} else {
    $logger->error("Sample size must be a positive integer. ".
                   "Set via '--samples' option.");
    exit 1;
}

## Sanity check of input parameters:

# --offset
if ( ! (Scalar::Util::Numeric::isint($offset)) ) {
    $logger->error("Please provide valid integer value via '--offset' option'");
    exit 1;
}
# --begin
if (defined $seqpos_begin && 
    ! (Scalar::Util::Numeric::isint($seqpos_begin)) ) {
    $logger->error("Please provide valid integer value via '--begin' option'");
    exit 1;
}
# --end
if (defined $seqpos_end && 
    ! (Scalar::Util::Numeric::isint($seqpos_end)) ) {
    $logger->error("Please provide integer value via '--end' option'");
    exit 1;
}

# === Perform the probing ===

my ($seq, @structures);
if ( defined $fasta_file ) {
    $seq = $fasta->sequence();
    @structures = &stochastic_sampling($seq, $sample_size);
} elsif ( defined $dbn_file ){
    $seq = $fasta->sequence();
    @structures = ($fasta->structure) x $sample_size;
}
$seq = uc($seq);
$logger->debug($seq);
my @probing_profile = (0) x length($seq);
#my $length = scalar(@probing_profile);
# my @structure_description = &dot_bracket_to_structure_description(@structures);
@probing_profile = &simulate_probing(\@structures, \@probing_profile, 
                                     $seq, $chemical);
#print(Dumper($bp_per_structure[1]));

# === Log Result if needed ===
# Should be logged instead of printed

$logger->info("=== Results ===");
#for (my $i = 0; $i < scalar(@structures); $i++) {
#    $logger->info("$i. Structure:\n$structures[$i]\n$structure_description[$i]");
#}
$logger->info(join(",", @probing_profile));

# === Fiddle around with offset, seqpos_begin and seqpos_start ===
my ($reactivity_begin, $reactivity_end);
if (defined($seqpos_begin) && defined($seqpos_end) &&
    # EITHER $seqpos_begin and $seqpos_end are set via options ...
    $offset < $seqpos_begin && $seqpos_begin < $seqpos_end && 
    $seqpos_end < length($seq)+$offset
    # ... and they obey the rules ...
    ) {
    # ... new $reactivity_begin and $reactivity_end of REACTIVITY are calculated, ...
    ## $offset minus $seqpos_* is fine, but for instance first nucleotide would be 1
    ## so we need to substract 1 to get a fine array index
    $reactivity_begin = abs($offset - $seqpos_begin ) - 1;
    $reactivity_end = abs($offset - $seqpos_end) - 1;
} else {
    # ... OR $seqpos_begin and $seqpos_end are set here.
    ## f**k you RDAT indices
    ## sequence enumeration starts at OFFSET plus 1
    ## and SEQPOS is 1-indexed  
    $seqpos_begin = 1 + $offset;
    ## length() starts counting at 1 so no need for extra addition
    $seqpos_end = length($seq) + $offset;
    ## $reactivity_begin and reactivity_end are array indices
    ## so substract one from the $seqpos_* values
    $reactivity_begin = $seqpos_begin - 1;
    $reactivity_end = $seqpos_end - 1;
}

# Fill @reactivity with the correct values depending on $offset and $seqpos_begin
# and $seqpos_end 
my @reactivity = @probing_profile[$reactivity_begin..$reactivity_end];

# === Assemble the RDAT file ===

my $fasta_id = $fasta->fasta_id();
$fasta_id =~ s/\.\w*$//g;
my $rdat_file_name = $fasta_id."_".$chemical->probe_name()."_".$sample_size.".rdat";
$rdat_file_name =~ s/\s//g;
$logger->info("++++ ".$rdat_file_name." ++++");
my $rdat_out = RNAprobing::RDATFile->new($rdat_file_name);

my $construct_name = $fasta->fasta_id()."_in_silico_probed_using_".$chemical->probe_name();

$rdat_out->name($construct_name);
$rdat_out->sequence($seq);

my ($struct, $mfe) = RNA::fold($seq);  # predict mfe structure of $seq
$rdat_out->structure($struct);
$rdat_out->offset($offset);
$rdat_out->seqpos([$seqpos_begin..$seqpos_end]);
$rdat_out->data()->reactivity(1, \@reactivity); # 1 is the index of the DATA line
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
    SELECT:{
	    if ($verbose == 0){
	           $logger->level($ERROR);
	           $logger->debug("Log level is ERROR");
	           last SELECT;
        } elsif ($verbose == 1){ 
               $logger->level($WARN);
               $logger->debug("Log level is WARN");
               last SELECT;
	    } elsif ($verbose == 2){
	           $logger->level($INFO);
	           $logger->debug("Log level is INFO");
	           last SELECT;}
	    else {
	           $logger->level($DEBUG);
	           $logger->debug("Log level is DEBUG");
	           last SELECT;
	    }
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

sub dot_bracket_to_structure_description {
    my ($structure) = @_;
    my $structure_description = "";
       
    my @db = split('', $structure);
    
    my @dots = (-1);
    my @opening_br = (-1);
    
    # at first we expect every nucleotide to be unpaired "U"
    my @struc_dec = ("U") x length($structure);
    
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
                # stop if the last enclosed nucleotide is reached
                last if ($dots[$#dots] < $op_br_pos);
                # if found non-continous numbers we found an enclosed stem
                if ( $dots[$#dots]+1 < $enclosed_dots[$#enclosed_dots]) {
                    $enclosed_stems++;
                }
            }
            
            # Hairpin or bulge detected
            if ( $enclosed_stems == 0 ) {
                # Hairpin detected if the enclosed dots reach from opening to
                # closing bracket
                if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {            
                    foreach (@enclosed_dots) { $struc_dec[$_] = "H" }
                }
                # bulge detected if the enclosed dots just touch one bracket
                elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                        $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "B" }
                }
            }
            # Multi loop with two included stems or an interior loop detected
            elsif ( $enclosed_stems == 1 ) {
                # Interior loop detected if the enclosed dots reach from opening
                # to closing bracket
                if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "I" }
                }
                # Multi loop with two stems detected if the enclosed dots just
                # touch one bracket
                elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                        $cl_br_pos - 1 == $enclosed_dots[0] ) {
                    foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
                }
            }
            # Multi loop detected if more than one stem is enclosed
            elsif( $enclosed_stems > 1 ) {
                foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
            }
        }
        # No enclosed unpaired nucleotides so lets declare the base pair
        elsif ( $db[$i] eq ")" ) {
            $struc_dec[$i] = "P";
            my $rel_open_br = pop(@opening_br);
            $struc_dec[$rel_open_br] = "P";
        }
    }
    $structure_description = join("",@struc_dec);
    return $structure_description;
}

sub create_bp_hash {
    my ($structure) = @_;

    my %base_base_pairings = ();
    my @db = split('', $structure);
    my @dots = (-1);
    my @opening_br = (-1);
    my @closing_br = (-1);
    # fill arrays with positions
    for (my $i = 0; $i < scalar(@db); $i++) {
        push(@dots, $i) if ( $db[$i] eq "." );
        push(@opening_br, $i) if ( $db[$i] eq "(" );
        if ( $db[$i] eq ")" ){
            my $five_prime_pos = pop(@opening_br);
            my $three_prime_pos = $i;
            $base_base_pairings{$three_prime_pos} = $five_prime_pos;
            $base_base_pairings{$five_prime_pos} = $three_prime_pos;
        }
    }
    return %base_base_pairings;
}


################################################################################
#
#           Probe sequence and secondary structure
#
################################################################################

sub simulate_probing {
    my ($structures, $probing_profile, $seq, $chemical) = @_;

    for ( my $i = 0; $i < scalar(@{$chemical->probe_reac()}); $i++ ) {
        my $prob_reac = ${$chemical->probe_reac()}[$i];
        my $prob_seq = @{$chemical->probe_seq()}[$i];
        my $prob_str = @{$chemical->probe_str()}[$i];
        my $prob_mod = @{$chemical->probe_mod()}[$i];
        my $prob_mod_pos = @{$chemical->probe_mod_pos()}[$i];
        
        my %block_connect = &find_str_block_connect($prob_str);

        # because @prob_seq, @prob_str, @prob_mod, and @prob_mod_pos have equal
        # length I can iterate over one of them

        my (@seq_regexp, @str_regexp);
        for ( my $i=0; $i < scalar(@{$prob_seq}); $i++) {
            # create regexp for sequence given in chemical file
            my $seq_re = '';
            my %letter2nucleotides = &single_letter_codes_for_nucleotides();
            foreach my $nuc ( split('', $prob_seq->[$i]) ) {
                $seq_re .= '['.$letter2nucleotides{$nuc}.']';
            }
            push(@seq_regexp, $seq_re);
        }

        for ( my $i=0; $i < scalar(@{$prob_str}); $i++) {
            # create regexp for structure given in chemical file
            my $str_re = '';
            my %letter2structure = &letter_codes_for_structure_elements();
            foreach my $str ( split('', $prob_str->[$i]) ) {
                $str_re .= '['.$letter2structure{$str}.']';
            }
            push(@str_regexp, $str_re);
        }            
        
        foreach my $structure (@{$structures}) {
            my $str_desc = 
                &dot_bracket_to_structure_description($structure);
            my (@seq_matches, @str_matches);
            # find all matching positions for sequence regexp
            foreach my $seq_re (@seq_regexp) {
                my @matches = &match_all_positions($seq_re, $seq);
                push(@seq_matches, \@matches);
            }
            # find all matching positions for the structure regexp
            foreach my $str_re (@str_regexp) {
                my @matches = &match_all_positions($str_re, $str_desc);
                push(@str_matches, \@matches);
            }

            my @block_match;
            # find identical matches in @seq_matches and @str_matches
            for (my $j = 0; $j < scalar(@seq_matches); $j++ ) {
                my @matching_positions =();
                my ($str_el, $seq_el) = (0, 0);
                my $seq_block_matches = $seq_matches[$j];
                my $str_block_matches = $str_matches[$j];

                while ( $seq_el < scalar(@{$seq_block_matches}) && 
                        $str_el < scalar(@{$str_block_matches}) ) {
                    if ( $seq_block_matches->[$seq_el] < 
                         $str_block_matches->[$str_el] ) {
                        $seq_el++;
                    } elsif ( $seq_block_matches->[$seq_el] > 
                              $str_block_matches->[$str_el] ) {
                        $str_el++;
                    } elsif ($seq_block_matches->[$seq_el] == 
                             $str_block_matches->[$str_el] ) {
                        push(@matching_positions, $seq_block_matches->[$seq_el]);
                        $str_el++;
                        $seq_el++;
                    }
                }
                push(@block_match, \@matching_positions);
#                print("Matches: ".join(",", @matching_positions)."\n");
            }

            # check if found matches for each block are connected by base pairs

            my %bp_per_structure = &create_bp_hash($structure);

            ## At first check if all requirements are fulfilled ...

            # variables for later use
            my $matching_blocks = 0;
            my @matching_positions;

            # ... each block must have matched ...             
            for ( my $j = 0; $j < scalar(@block_match); $j++) {
                # check if there has something been found for block $j
                $matching_blocks++ if (@{$block_match[$j]});
            }
                # ... and only matches which obey the rules are passed on ...

            if ( %block_connect ) {
                foreach my $key ( keys(%block_connect) ){
                    my @key_pos = split(",", $key);
                    my @val_pos = split(",", $block_connect{$key});
                    foreach my $key_offset ( @{$block_match[$key_pos[0]]} ){
                        my $key_base_pos = $key_offset + $key_pos[1];
                        foreach my $val_offset (@{$block_match[$val_pos[0]]}){
                            my $val_base_pos = $val_offset + $val_pos[1];
                            if ( $bp_per_structure{$key_base_pos} == 
                                 $val_base_pos ) {
                                # need to calculate the modification point
                                foreach my $mod_pos (@{$prob_mod_pos->[$key_pos[0]]} ) {
 #                                   print(Dumper($prob_mod_pos));
#                                    print("Modification point: $mod_pos\n");
                                    my $prob_pos = $key_offset + $mod_pos;
#                                    print("Found a valid bp: $key_base_pos"
#                                          ."/$val_base_pos\n");
#                                    print("Probed position: $prob_pos\n");
                                    push(@matching_positions, $prob_pos);
                                }
                            }
                        }
                    }
                }
            } else {
                for ( my $j = 0; $j < scalar(@block_match); $j++) {
                    foreach my $offset ( @{$block_match[$j]} ) {
                        foreach my $mod_pos (@{$prob_mod_pos->[$j]}) {
                            my $prob_pos = $offset + $mod_pos;
                            push(@matching_positions, $prob_pos);
                        }
                    }
                }
            }
#            print(Dumper(\%block_connect));
#            print("Found matches for "
#                  .$matching_blocks."/".scalar(@{$prob_seq})
#                  ." blocks\n");
            
            ## Finally modify the probing profile if all blocks matched
            if ( $matching_blocks == scalar(@{$prob_seq}) ) {
                foreach my $prob_pos (@matching_positions) {
                    $probing_profile->[$prob_pos] += $prob_reac;
                    $logger->info( "Probing profile:\n".
                           join(",",@{$probing_profile})."\n" );
                }
            }
        }
    }
    return @{$probing_profile};
}

sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /(?=$regex)/g) {
        push(@ret,  $-[0]);
    }
    return @ret;
}

sub find_str_block_connect {
    my ($prob_str) = @_;

    my (@op_br, @cl_br, %block_connect);
    for (my $i=0; $i < scalar(@{$prob_str}); $i++ ) {
        my @str = split('', $prob_str->[$i]);
        for (my $j=0; $j < scalar(@str); $j++ ) {
            if ( $str[$j] eq '(' ) {
                # create unique key for opening bracket:
                # * key = block,position
                my $op_br_pos = join(",", $i, $j);
                # save position on stack
                push( @op_br,  $op_br_pos);
            }
            if ( $str[$j] eq ')' ) {
                # create unique key for closing bracket
                # * key = block,position
                my $cl_br_pos = join( ",", $i, $j);
                # get associated opening bracket
                my $op_br_pos = pop( @op_br );
                # fill hash with coordinates
                $block_connect{$op_br_pos} = $cl_br_pos;
                $block_connect{$cl_br_pos} = $op_br_pos;
            }
        }
    }
    return %block_connect;
}

sub single_letter_codes_for_nucleotides {
    my %letter2nucleotides = ( 
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        U => 'U',
        M => 'AC',
        R => 'AG',
        W => 'ATU',
        S => 'CG',
        Y => 'CTU',
        K => 'GTU',
        V => 'ACG',
        H => 'ACTU',
        D => 'AGTU',
        B => 'CGTU',
        N => 'ACGTU' );
    return %letter2nucleotides;
}


sub letter_codes_for_structure_elements {
    my %letter2structure = (
        'B' => 'B',
        'H' => 'H',
        'I' => 'I',
        'M' => 'M',
        'P' => 'P',
        'U' => 'HBIMU',
        '(' => 'P',
        ')' => 'P' );
    return %letter2structure;
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

=item B<-f, --fasta>

Fasta file containing RNA sequence to be probed

=item B<-c, --chemical>

Chemical file describing the reactivities of the probing reagent

=item B<--samples>

Sets value for Boltzmann ensemble sampling. Expects positive integer value as parameter. Default=1000

=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (highest level if used 3 or more times) 

=back

=cut
