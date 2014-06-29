package RNAprobing::RDATFile;

use strict;
use warnings;
use Carp;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;
require RNAprobing::RDATFile::Annotation;
require RNAprobing::RDATFile::Data;

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

sub new {
    my ($classname, $filename) = @_;
    my $self = {};

    &filename($self, $filename);
    &rdat_version($self);
    &name($self);
    &sequence($self);
    &structure($self);
    &offset($self);
    &seqpos($self);
    &annotation($self);
    &comment($self);
    &mutpos($self);
    &data($self);

    bless $self, $classname;
    if ( defined $filename ) {
        if ( -e $filename ) {
            $self->read_file($filename);
        } else {
            die "File ".$filename." does not exist.";
        }
    }
    return $self;
}
################################################################################
################################################################################
##
##  Subroutine section
##
################################################################################
################################################################################

################################################################################
##
##  Subroutine that reads in a RDAT file
##
################################################################################

sub read_file {
    my ($self, $filename) = @_;
    $self->filename( $filename ) if ( defined $filename );
    
    open ( my $rdat_file , "<", $self->filename() ) or 
	croak "Couldn't open file $self->filename(). Error: $!";

    my $lines = "";

    while (my $line = <$rdat_file>) {
        chomp( $line );
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        my @words = split(/\s+/, $line, 2);
        my @split = split( /\s+/, $words[1]);

        # Construct section lines
        $self->rdat_version( $words[1] ) if ($words[0] =~ /^RDAT_VERSION$/);
        $self->name( $words[1] ) if ($words[0] =~ /^NAME$/);
        $self->sequence( $words[1] ) if ($words[0] =~ /^SEQUENCE$/);
        $self->structure( $words[1] ) if ($words[0] =~ /^STRUCTURE$/);
        $self->offset( $words[1] ) if ($words[0] =~ /^OFFSET$/);
        $self->seqpos( [@split] ) if ($words[0] =~ /^SEQPOS$/);
        $self->annotation( $words[1] ) if ($words[0] =~ /^ANNOTATION$/);
        $self->comment( $words[1]." " ) if ($words[0] =~ /^COMMENT$/);
        $self->xsel( [@split] ) if ($words[0] =~ /^XSEL$/);
        $self->mutpos( [@split] ) if ($words[0] =~ /^MUTPOS$/);

        # Collect data section lines
        if ($words[0] =~ /^[\w]*:(\d+)$/) {
            $lines .= $line."\n";
        }
    }
    close($rdat_file);
    $self->data($lines);
    foreach my $index ( @{$self->data()->indices()} ){
        $self->seqpos_reactivity_map( $index, $self->seqpos() );
        $self->seqpos_scaled_reactivity_map( $index, $self->seqpos() );
        $self->seqpos_non_negative_reactivity_map( $index, $self->seqpos() );
    }

    if ( $self->rdat_version() == 0.24 ) {
        $logger->info("Correct version.");
    } else {
        $logger->error("Incorrect Version. There may occur errors while parsing ".$filename);
    }

    return $self;
}

################################################################################
##
##  Write file subroutine
##
################################################################################

sub write_file {
    my ($self, $filename) = @_;
    $self->filename($filename) if (defined $filename);
    my $file_content = $self->serialize_rdat_version().
        $self->serialize_name().
	$self->serialize_sequence().
        $self->serialize_structure().
        $self->serialize_offset().
        $self->serialize_seqpos().
        $self->serialize_annotation().
        $self->serialize_comment().
        $self->serialize_xsel().
        $self->serialize_mutpos().
	$self->serialize_data();

    open( my $rdat_file, ">", $self->filename() ) or 
	croak( "Couldn't open file". $self->filename() );
    print($rdat_file  $file_content);
    close($rdat_file);
}


# Construct section

sub seqpos_reactivity_map {
    my ($self, $index, $seqpos) = @_;
    if ( defined $index ){
	if ( ref $seqpos eq "ARRAY" ){
	    return $self->data()->seqpos_reactivity_map($index, $seqpos);
	}
	return $self->data()->seqpos_reactivity_map($index);
    }
    return $self->data()->seqpos_reactivity_map();
}

sub seqpos_scaled_reactivity_map {
    my ($self, $index, $seqpos) = @_;
    if ( defined $index ){
	    if ( ref $seqpos eq "ARRAY" ){
	        return $self->data()->seqpos_scaled_reactivity_map($index, $seqpos);
	    }
	    return $self->data()->seqpos_scaled_reactivity_map($index);
    }
    return $self->data()->seqpos_scaled_reactivity_map();
}

sub seqpos_non_negative_reactivity_map {
    my ($self, $index, $seqpos) = @_;
    if ( defined $index ){
	    if ( ref $seqpos eq "ARRAY" ){
	        return $self->data()->seqpos_non_negative_reactivity_map($index, $seqpos);
	    }
	    return $self->data()->seqpos_non_negative_reactivity_map($index);
    }
    return $self->data()->seqpos_non_negative_reactivity_map();
}

################################################################################
##
##  Getter/Setter subroutines 
##
################################################################################

sub filename {
    my ($self, $filename) = @_;
    my $method_key = "FILENAME";
    if ( defined $filename ){
        $self->{$method_key} = $filename;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

#
# ==== RDAT General section ====
#

sub rdat_version {
    my ($self, $rdat_version) = @_;
    my $method_key = "RDAT_VERSION";
    
    if ( defined $rdat_version){
        $self->{$method_key} = $rdat_version;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "0.24";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_rdat_version {
    my $self = shift;
    my $line = "RDAT_VERSION ".$self->rdat_version()."\n";
    return $line;
}

#
# ==== RDAT Construct section ==== 
#

# NAME
# The name of the construct in the following format: "RNA system name, Source"
# where source can be the organism or gene where the system is usually studied
# (in some cases, one can be more relevant than the other).. If the RNA is of
# synthetic origin, a source must still be specified.

sub name {
    my ($self, $name) = @_;
    my $method_key = "NAME";

    if ( defined $name ){
        $self->{$method_key} = $name;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_name {
    my $self = shift;
    if ( $self->name() ){
        my $line = "NAME ".$self->name()."\n";
        return $line;
    } else {
        $logger->error("Missing mandatory NAME entry.");
        $logger->error("RDAT output is erroneous.");
    }
}

# SEQUENCE
# The sequence of the RNA molecule studied. The sequence of interest should be
# in upper case, while auxiliary sequences (e.g. primer binding sites or barcodes)
# are in lower case. Degenerate sequences describing populations of sequences
# probed can be described with IUPAC codes (e.g. AGGAAAGGNNRNNAAGAA). This is
# particularly useful when an experiment compares different sequences that are
# expected to fold in the same structure, such as those found in the EteRNA
# project.

sub sequence {
    my ($self, $sequence) = @_;
    my $method_key = "SEQUENCE";

    if ( defined $sequence){
        $self->{$method_key} = $sequence;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_sequence {
    my $self = shift;
    if ( $self->sequence() ){
        my $line = "SEQUENCE ".$self->sequence()."\n";
        return $line;
    } else {
        $logger->error("Missing mandatory SEQUENCE entry.");
        $logger->error("RDAT output is erroneous.");
    }
}

# STRUCTURE
# A string describing the best known secondary structure model of the RNA of
# interest, in dot bracket notation (e.g. "((....))"). Pseudoknots can be
# included using "{}" curly braces (e.g "((..{.{))}.}").

sub structure {
    my ($self, $structure) = @_;
    my $method_key = "STRUCTURE";
    # might check for correct notation here
    if ( defined $structure){
        $self->{$method_key} = $structure;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar (hopefully in dot-bracket notation)
}

sub serialize_structure {
    my $self = shift;
    if ( $self->structure() ){
        my $line = "STRUCTURE ".$self->structure()."\n";
        return $line;
    } else {
        $logger->error("Missing mandatory STRUCTURE entry.");
        $logger->error("RDAT output is erroneous.");
    }
}

# OFFSET
# Sometimes, the sequence studied is a sub-sequence from a bigger sequence (such
# as the P4-P6 domain of the Tetrahymena ribozyme group I intron) and the bigger
# sequence numbering scheme will want to be preserved. The offset then specifies
# the number such that Pos = OrPos - OFFSET, where Pos is the position in the
# subsequence and OrPos is the same position but in the numbering scheme of the
# bigger sequence. Note that positions are always considered to be 1-indexed.

sub offset {
    my ($self, $offset) = @_;
    my $method_key = "OFFSET";

    if ( defined $offset){
        $self->{$method_key} = $offset;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_offset {
    my $self = shift;
    my $line = "OFFSET ".$self->offset()."\n";
    return $line;
}

# SEQPOS
# A list of sequence positions and their respective residues that were probed in
# the experiment. If OFFSET is defined, each position in SEQPOS must have OFFSET
# already added to it (thus, a position OrPos defined in SEQPOS refers to position
# Pos = OrPos - OFFSET in the sequence defined in SEQUENCE.). Positions in SEQPOS
# are 1-indexed.

sub seqpos {
    my ($self, $seqpos) = @_;
    my $method_key = "SEQPOS";

    if ( defined $seqpos){
        $self->{$method_key} = $seqpos;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_seqpos{
    my $self = shift;
    my $line = "SEQPOS ".join( " ", @{ $self->seqpos() } )."\n";
    return $line;
}

# ANNOTATION
# Annotations for this construct section. See the Annotations section below for
# more details.

sub annotation {
    my ($self, $annotation) = @_;
    my $method_key = "ANNOTATION";

    if ( defined $annotation){
        $self->{$method_key} = RNAprobing::RDATFile::Annotation->new($annotation);
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = ();
    }
    return $self->{$method_key} ; # returns a RNAprobing::RDATFile::Annotation
                                  # object reference or an undefined value
}

sub serialize_annotation{
    my $self = shift;
    my $annotation = $self->annotation();
    my $line = "";
    if (defined $annotation ) {
        $line .= "ANNOTATION". $annotation->serialize_annotation() ."\n";
    }
    return $line;
}

# COMMENT
# Free text comments for the construct.

sub comment {
    my ($self, $comment) = @_;
    my $method_key = "COMMENT";

    if ( defined $comment){
        $self->{$method_key} = $comment;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_comment{
    my $self = shift;
    my $line = "";
    if ( $self->comment() ) {
        $line = "COMMENT ".$self->comment()."\n";
    }
    return $line;
}

# XSEL
# List of tab-separated integers that specify the place in the traces used for
# peak calling in the sequence assignment step of the analysis, used globally.
# This and XSEL_REFINE (see Data Section description below) are mutually
# exclusive and when provided with both, XSEL will be ignored.

sub xsel {
    my ($self, $xsel) = @_;
    my $method_key = "XSEL";

    if ( defined $xsel){
        $self->{$method_key} = $xsel;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_xsel{
    my $self = shift;
    my $line = "";
    if ( @{$self->xsel()} ) {
        $line = "XSEL\t".join( "\t", @{ $self->xsel() } )."\n";
    }
    return $line;
}

# MUTPOS
# List of tab-separated values with size of the number of data lanes that
# indicate which samples have been mutated and where (only single point mutations
# are accepted in MUTPOS). WT indicates that there is no mutation on the sample.

sub mutpos {
    my ($self, $mutpos) = @_;
    my $method_key = "MUTPOS";

    if ( defined $mutpos){
        $self->{$method_key} = $mutpos;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_mutpos{
    my $self = shift;
    my $line = "";
    if ( @{$self->mutpos()} ) {
        $line = "MUTPOS\t".join( "\t", @{ $self->mutpos() } )."\n";
    }
    return $line;
}

#
# ==== RDAT Data section ====
#


# Everything from the data section in an RDAT file is going into a new 
# RNAprobing::RDATFile::Data object which can hold all information given there

sub data {
    my ($self, $lines) = @_;
    my $method_key = "DATA";

    if ( defined $lines){
        $self->{$method_key} = RNAprobing::RDATFile::Data->new($lines);
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = RNAprobing::RDATFile::Data->new();
    }
    return $self->{$method_key}; # returns a RNAprobing::RDATFile::Data object
    # or a undefined value
}

sub serialize_data{
    my ($self) = @_;
    my $line = "";
    if ( $self->data() ) {
        $line = $self->data()->serialize_data()."\n";
    }
    return $line;
}


################################################################################
##
##  Mappings for RDAT contained information 
##
################################################################################

sub seq_startpos {
    my ($self) = @_;
    my $startpos = "";
    # 1 is used because sequence is 1-indexed
    if ( $self->offset() =~ /[+-]?\d+/ ) {
        $startpos = 1 + $self->offset();
        $logger->debug("Reactivity start position given OFFSET: $startpos");
    } else {
        undef $startpos;
    }
    return $startpos;
}

sub seq_endpos {
    my ($self) = @_;
    my $endpos = "";
    if ( defined $self->sequence() && $self->offset() =~ /[+-]?\d+/ ) {
	    $endpos = length($self->sequence()) + $self->offset();
	    $logger->debug("Reactivity end position given OFFSET: $endpos");
    } else {
    	undef $endpos;
    }
    return $endpos;
}

# this mapping creates a hash whose keys start at OFFSET and cover the
# whole sequence, such that the keys of this hash and those of SEQPOS are
# mapping towards the same base

sub offset_sequence_map {
    my ($self) = @_;
    my $method_key = "OFFSET_SEQUENCE_MAP";
    return $self->{$method_key} if ( defined $self->{$method_key} );
    if ( defined $self->offset() && defined $self->sequence() ){
	    my $offset = $self->offset();
	    my @sequence = split(//, $self->sequence());
	    # 1 is used because sequence is 1-indexed
	    my $startpos = $self->seq_startpos();
	    my $endpos = $self->seq_endpos();
	    $logger->info("OFFSET = ".$offset);
	    $logger->info("Sequence length = ".scalar(@sequence));
	    $logger->info("Start position given OFFSET = ".$startpos);
	    $logger->info("End position given OFFSET = ".$endpos);
	    my $off_seq_map = {};
	    my $seq_index = 0;
	    for (my $i = $startpos; $i <= $endpos; $i++ ) {
	        $off_seq_map->{$i} = $sequence[$seq_index];
	        $seq_index++;
	    }
        $self->{$method_key} = $off_seq_map;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = {};
    }
    return $self->{$method_key}; # returns a hash reference
}


# this mapping creates a hash whose keys are the positions of the nucleotides in
# the RNA sequence and map these onto the values starting with OFFSET.
sub ntpos_offset_map {
    my ($self) = @_;
    my $method_key = "NTPOS_OFFSET_MAP";
    return $self->{$method_key} if ( defined $self->{$method_key} );
    if ( defined $self->offset() && defined $self->sequence() ){
        my $offset = $self->offset();
        my $sequence = $self->sequence();
        my $ntpos_offset_map = {};
        for (my $i = 0; $i < length($sequence); $i++) {
            $ntpos_offset_map->{$i} = $self->offset() + $i;
        }
        $self->{$method_key} = $ntpos_offset_map;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = {};
    }
    return $self->{method_key}; # returns a hash reference
}

1;
