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
my $length = scalar(@probing_profile);
my @structure_description = &dot_bracket_to_structure_description(@structures);
my @bp_per_structure = &create_bp_hash(\@structures);
@probing_profile = &simulate_probing(\@structures, \@structure_description, 
                                     \@probing_profile, $seq, $chemical);
print(Dumper($bp_per_structure[1]));

# === Log Result if needed ===
# Should be logged instead of printed

$logger->info("=== Results ===");
for (my $i = 0; $i < scalar(@structures); $i++) {
    $logger->info("$i. Structure:\n$structures[$i]\n$structure_description[$i]");
}
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
        my $str_desc = join("",@struc_dec);
        push(@structure_description, join("",@struc_dec));
    }
    return @structure_description;
}

sub create_bp_hash {
    my ($structures) = @_;
    my @base_base_per_structure = ();
    foreach my $structure (@{$structures}) {
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
        push(@base_base_per_structure, \%base_base_pairings);
    }
    return @base_base_per_structure;
}


################################################################################
#
#           Probe sequence and secondary structure
#
################################################################################

sub simulate_probing {
    my ($bp_per_structure, $structure_description, 
        $probing_profile, $seq, $chemical) = 
        @_;

    for ( my $i = 0; $i < scalar(@{$chemical->probe_reac()}); $i++ ) {
        my $prob_reac = ${$chemical->probe_reac()}[$i];
        my $prob_seq = ${$chemical->probe_seq()}[$i];
        my $prob_str = ${$chemical->probe_str()}[$i];
        my $prob_cut = ${$chemical->probe_cut()}[$i];

        # create regex from probing sequence
        
        my $seq_regex = '';
        my %letter2nucleotides = &single_letter_codes_for_nucleotides();
        foreach my $nuc ( split('', $prob_seq) ) {
            $seq_regex .= '['.$letter2nucleotides{$nuc}.']';
        }

        # create regex from probing structure
        my $str_regex = $prob_str;
        $str_regex =~ s/U/[HBIMn]/g;
        $str_regex =~ s/n/U/g;

        foreach my $str_desc (@{$structure_description}) {
            # Find the modification/cut point
            $prob_cut =~ /\|/;
            my $cut_pos = $-[0];

            my @seq_matches = &match_all_positions($seq_regex, $seq, $cut_pos);
            my @str_matches = &match_all_positions($str_regex, $str_desc, $cut_pos);

            if (scalar(@seq_matches) > scalar(@str_matches) ) {
                my $more_matches = join(",",@seq_matches);
                foreach (@str_matches) {
                    ${probing_profile}[$_] += $prob_reac if ($more_matches =~ /^$_,|,$_,|,$_$/);
                }
            }
            else {
                 my $more_matches = join(",",@str_matches);
                foreach (@seq_matches) {
                    ${probing_profile}[$_] += $prob_reac  if ($more_matches =~ /^$_,|,$_,|,$_$/);
                }
            }
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
    return @ret;
}


sub single_letter_codes_for_nucleotides {
    my %letter2nucleotides = ( A => 'A',
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
                               N => 'ACGTU'
        );
    return %letter2nucleotides;
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
