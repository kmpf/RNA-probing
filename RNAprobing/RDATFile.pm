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

    &annotation_data($self);
    &reactivity($self);
    &reactivity_error($self);
    &xsel_refine($self);
    &seqpos_reactivity_map($self);
    &scaled_reactivity($self);

    bless $self, $classname;

    $self->read_file($filename) if ($filename ne "");

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
    $self->filename( $filename );
    
    open ( my $rdat_file , "<", $self->filename() ) or croak "Couldn't open file $self->filename(). Error: $!";

    my ($name, $lines) = "";

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
    
    my $data_object = RNAprobing::RDATFile::Data->new($lines);

    $self->seqpos_reactivity_map( $self->seqpos(), $self->reactivity() );
    $self->scaled_reactivity( $self->reactivity() );

    if ( $self->rdat_version() == 0.24 ) {
        $logger->info("Correct version.");
    } else {
        $logger->error("Incorrect Version. There may occur errors while parsing this file.");
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
        $self->serialize_annotation_data().
        $self->serialize_reactivity().
        $self->serialize_reactivity_error().
        $self->serialize_xsel_refine().
#        $self->serialize_data_seqpos().
        $self->serialize_trace().
        $self->serialize_reads();

    open( my $rdat_file, ">", $self->filename() ) or 
	croak( "Couldn't open file". $self->filename() );
    print($rdat_file  $file_content);
    close($rdat_file);
    $logger->info($file_content);
}


# Construct section

sub seqpos_reactivity_map {
    my ($self, $seqpos, $reactivity) = @_;
    my $method_key = "SEQPOS_REACTIVITY_MAP";

    if ( defined $reactivity && defined $seqpos ){
        foreach my $reac_entry ( @{$reactivity}){
            push( @{$self->{$method_key}}, 
                $self->_create_seqpos_reactivity_hash($seqpos, $reac_entry) );
        }
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference to an array of
                                 # hashes or just an empty array
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


################################################################################
##
##  Subroutine that creates SEQPOS => REACTIVITY hash
##
################################################################################

sub _create_seqpos_reactivity_hash {
    my ($self, $seqpos, $reac_entry) = (@_);
#    my @seqpos_reactivity = ();
    my %pos_reac = ();

    if (defined $reac_entry && defined $seqpos ){
        if ( scalar( @{$seqpos} ) == scalar( @{$reac_entry} ) ) {
            for ( my $i = 0; $i < scalar( @{$seqpos} ); $i++ ) {
                $pos_reac{ $seqpos->[$i] } = $reac_entry->[$i];
            }
        } else {
            $logger->error("Line SEQPOS in ".$self->filename()." has ".scalar(@{ $seqpos })." entries.");
            $logger->error("Line REACTIVITY in ".$self->filename()." has ".scalar(@{ $reac_entry })." entries.");
            $logger->error("Both lines should have an identical number of entries.");
            $logger->error("Check your .rdat file for consistency!");
        }
    }
    return \%pos_reac;
}

sub _href_to_string{
    my ($self, $href) = @_;
    my $line = "";
    foreach my $key ( keys %{ $href } ) {
        $line = " ".$key.":".$href->{$key};
    }
    return $line;
}

1;
