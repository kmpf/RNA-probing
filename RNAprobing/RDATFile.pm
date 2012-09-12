package ProbingRNA::RDATFile;

use strict;
use warnings;
use Carp;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

sub new {
    my ($classname, $filename) = @_;
    my $self = {};
    print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
    my $logger = get_logger();

    $self->{"FILENAME"} = $filename;
    $self->{"RDAT_VERSION"} = "";
    $self->{"NAME"} = "";
    $self->{"SEQUENCE"} = "";
    $self->{"STRUCTURE"} = "";
    $self->{"OFFSET"} = "";
    $self->{"SEQPOS"} = [];
    $self->{"MUTPOS"} = [];
    $self->{"ANNOTATION"} = [];
    $self->{"COMMENT"} = "";
    $self->{"ANNOTATION_DATA"} = [];
    $self->{"ANNOTATION_DATA"} = [];
    $self->{"REACTIVITY"} = [];
    $self->{"REACTIVITY_ERROR"} = [];
    $self->{"SEQPOS_REACTIVITY_MAP"} = [];
    
    open ( my $rdat_file , "<", $self->{FILENAME} ) or croak "Couldn't open file $self->{FILENAME}. Error: $!";
    while (my $line = <$rdat_file>) {
        chomp( $line );
        my $string = '^#|^\s*$';
        next if ($line =~ /$string/ );
        my @words = split(/\s+/, $line, 2);
        my @split = split( /\s+/, $words[1]);
        $self->{"RDAT_VERSION"} = $words[1]                     if ($words[0] =~ /^RDAT_VERSION$/);
        $self->{"NAME"} = $words[1]                             if ($words[0] =~ /^NAME$/);
        $self->{"SEQUENCE"} = $words[1]                         if ($words[0] =~ /^SEQUENCE$/);
        $self->{"STRUCTURE"} = $words[1]                        if ($words[0] =~ /^STRUCTURE$/);
        $self->{"OFFSET"} = $words[1]                           if ($words[0] =~ /^OFFSET$/);
        $self->{"SEQPOS"} = [@split]                            if ($words[0] =~ /^SEQPOS$/);
        $self->{"MUTPOS"} = [@split]                            if ($words[0] =~ /^MUTPOS$/);
        $self->{"ANNOTATION"} = &_parse_annotation($words[1])   if ($words[0] =~ /^ANNOTATION$/);
        $self->{"COMMENT"} .= $words[1]." "                     if ($words[0] =~ /^COMMENT$/);
        push(@{ $self->{"ANNOTATION_DATA"} }, 
            &_parse_annotation($words[1]) )                     if ($words[0] =~ /^ANNOTATION_DATA:(\d+)$/);
        push(@{ $self->{"REACTIVITY"} }, [@split])              if ($words[0] =~ /^REACTIVITY:(\d+)$/);
        push(@{ $self->{"REACTIVITY_ERROR"} }, [@split])        if ($words[0] =~ /^REACTIVITY_ERROR:(\d+)$/);
        
    }
    close($rdat_file);

    
    $self->{"SEQPOS_REACTIVITY_MAP"} = &_create_seqpos_reactivity_hash($self);

    if ($self->{"RDAT_VERSION"} == 0.24 ) {
        $logger->info("Correct version.");
    } else {
        $logger->error("Incorrect Version. There may occur errors while parsing this file.");
    }

    bless( $self, $classname );
    return $self;
}

###############################################################
##
##  Subroutine section
##
###############################################################

###############################################################
##
##  Subroutine that parses ANNOTATION and ANNOTATION_DATA lines
##
###############################################################

sub _parse_annotation {
    my ($line) = @_;
    my %annotation = ();
    my @chemicals = (); # array for all chemicals that have been used
    my @processings = (); # array for all processings that have been done

    my @words = split(/\s+/, $line);
    for (my $i = 0; $i < scalar(@words); $i++) {
        my $word = $words[$i];
        $annotation{"EXPERIMENT_TYPE"} = $1     if ($word =~ /^experimentType:(\S+)/ && splice(@words, $i, 1));
        $annotation{"MODIFIER"} = $1            if ($word =~ /^modifier:(\S+)/ && splice(@words, $i, 1) );
        $annotation{"MUTATION"} = $1            if ($word =~ /^mutation:(\S+)/ && splice(@words, $i, 1) );
        push (@processings, $1)                 if ($word =~ /^processing:(\S+)/ && splice(@words, $i, 1) );
        $annotation{"TEMPERATURE"} = $1         if ($word =~ /^temperature:(\S+)/ && splice(@words, $i, 1) );
        push (@chemicals, my %chemical = ("CHEMICAL" => $1, "CHEMICAL_CONCETRATION" => $2) )
                                                if ($word =~ /^chemical:(\S+):(\S+)/ && splice(@words, $i, 1) );

    }
    print "Undigestible rest: ". join(" ",@words) ."\n" if ( scalar(@words) != 0 );
    $annotation{"CHEMICALS"} = [@chemicals]       if ( scalar(@chemicals) > 0 );
    $annotation{"PROCESSING"} = [@processings]    if ( scalar(@processings) > 0 );
    return \%annotation;
}

###############################################################
##
##  Subroutine that creates SEQPOS => REACTIVITY hash
##
###############################################################

sub _create_seqpos_reactivity_hash {
    my $self = shift;
    my @seqpos_reactivity = ();
    my $logger = get_logger();
    for (my $j = 0; $j < scalar( @{$self->{"REACTIVITY"}} ); $j++ ) {
        if ( scalar( @{$self->{"SEQPOS"}} ) == scalar( @{$self->{"REACTIVITY"}[$j]} ) ) {
            my %pos_reac = ();
            for ( my $i = 0; $i < scalar(@{ $self->{"SEQPOS"} }); $i++ ) {
                $pos_reac{ $self->{"SEQPOS"}[$i] } = $self->{"REACTIVITY"}[$j][$i];
            }
            push( @seqpos_reactivity, \%pos_reac );
        } else {
            $logger->error("Line SEQPOS in ".$self->{"FILENAME"}." has ".scalar(@{ $self->{"SEQPOS"} })." entries.");
            $logger->error("Line REACTIVITY in ".$self->{"FILENAME"}." has ".scalar(@{ $self->{"REACTIVITY"} })." entries.");
            $logger->error("Both lines should have an identical number of entries.");
            $logger->error("Check your .rdat file for consistency!");
        }
    }
    return \@seqpos_reactivity;
}


###############################################################
##
##  Write subroutine
##
###############################################################

sub write_rdat_file {
    my $self = shift;
    my $logger = get_logger();
    open( my $rdat_file, ">", $self->filename() ) or croak( "Couldn't open file". $self->filename() );
    print($rdat_file, ( $self->serialize_rdat_version(),
    $self->serialize_name(),
    $self->serialize_sequence(),
    $self->serialize_structure(),
    $self->serialize_offset(),
    $self->serialize_seqpos(),
    $self->serialize_mutpos(),
    $self->serialize_annotation(),
    $self->serialize_comment(),
    $self->serialize_annotation_data(),
    $self->serialize_reactivity(),
    $self->serialize_reactivity_error() ) );
    close($rdat_file);
}

###############################################################
##
##  Getter/Setter subroutines 
##
###############################################################

sub filename {
    my ($self, $filename) = @_;
    if ( defined $filename){
        $self->{"FILENAME"} = $filename;
    }
    return $self->{"FILENAME"};
}

sub rdat_version {
    my ($self, $rdat_version) = @_;
    if ( defined $rdat_version){
        $self->{"RDAT_VERSION"} = $rdat_version;
    }
    return $self->{"RDAT_VERSION"}; # returns a scalar
}

sub serialize_rdat_version {
    my $self = shift;
    my $line = "RDAT_VERSION ".$self->rdat_version()."\n";
    return $line;
}

sub name {
    my ($self, $name) = @_;
    if ( defined $name){
        $self->{"NAME"} = $name;
    }
    return $self->{"NAME"}; # returns a scalar
}

sub serialize_name {
    my $self = shift;
    my $line = "NAME ".$self->name()."\n";
    return $line;
}

sub sequence {
    my ($self, $sequence) = @_;
    if ( defined $sequence){
        $self->{"SEQUENCE"} = $sequence;
    }
    return $self->{"SEQUENCE"}; # returns a scalar
}

sub serialize_sequence {
    my $self = shift;
    my $line = "SEQUENCE ".join( "", @{ $self->sequence() } )."\n";
    return $line;
}

sub structure {
    my ($self, $structure) = @_;
    if ( defined $structure){
        $self->{"STRUCTURE"} = $structure;
    }
    return $self->{"STRUCTURE"}; # returns a scalar (hopefully in dot-bracket notation)
}

sub serialize_structure {
    my $self = shift;
    my $line = "STRUCTURE ".join( "", @{ $self->structure() } )."\n";
    return $line;
}

sub offset {
    my ($self, $offset) = @_;
    if ( defined $offset){
        $self->{"OFFSET"} = $offset;
    }
    return $self->{"OFFSET"}; # returns a scalar
}

sub serialize_offset {
    my $self = shift;
    my $line = "OFFSET ".$self->offset()."\n";
    return $line;
}

sub seqpos {
    my ($self, $seqpos) = @_;
    if ( defined $seqpos){
        $self->{"SEQPOS"} = $seqpos;
    }
    return $self->{"SEQPOS"}; # returns an array reference
}

sub serialize_seqpos{
    my $self = shift;
    my $line = "SEQPOS ".join( " ", @{ $self->seqpos() } )."\n";
    return $line;
}

sub mutpos {
    my ($self, $mutpos) = @_;
    if ( defined $mutpos){
        $self->{"MUTPOS"} = $mutpos;
    }
    return $self->{"MUTPOS"}; # returns an array reference
}

sub serialize_mutpos{
    my $self = shift;
    my $line = "MUTPOS ".join( " ", @{ $self->mutpos() } )."\n";
    return $line;
}

sub annotation {
    my ($self, $annotation) = @_;
    if ( defined $annotation){
        $self->{"ANNOTATION"} = $annotation;
    }
    return $self->{"ANNOTATION"} ; # returns a hash reference
}

sub serialize_annotation{
    my $self = shift;
    my $line = "ANNOTATION".$self->_href_to_string( $self->annotation()."\n" );
    return $line;
}

sub comment {
    my ($self, $comment) = @_;
    if ( defined $comment){
        $self->{"COMMENT"} = $comment;
    }
    return $self->{"COMMENT"}; # returns a scalar
}

sub serialize_comment{
    my $self = shift;
    my $line = "COMMENT ".$self->comment()."\n";
    return $line;
}

sub annotation_data {
    my ($self, $annotation_data) = @_;
    if ( defined $annotation_data){
        $self->{"ANNOTATION_DATA"} = $annotation_data;
    }
    return $self->{"ANNOTATION_DATA"}; # returns an array reference to an array of hashes
}

sub serialize_annotation_data{
    my $self = shift;
    my $line = "";
    my $anno_data_nr = 0;
    foreach my $href ( @{ $self->annotation_data() }  ) {
        $anno_data_nr++;
        $line =  "ANNOTATION_DATA:".$anno_data_nr;
        $line .= $self->_href_to_string( $href );
    }
    return $line;
}

sub reactivity {
    my ($self, $reactivity) = @_;
    if ( defined $reactivity){
        $self->{"REACTIVITY"} = $reactivity;
    }
    return $self->{"REACTIVITY"}; # returns an array reference to an array of arrays
}

sub serialize_reactivity{
    my $self = shift;
    my $line = "";
    my $reactivity_nr = 0;
    foreach my $href ( @{ $self->reactivity() }  ) {
        $reactivity_nr++;
        $line =  "REACTIVITY:".$reactivity_nr;
        $line .= $self->_href_to_string( $href );
    }
    return $line;
}

sub reactivity_error {
    my ($self, $reactivity_error) = @_;
    if ( defined $reactivity_error){
        $self->{"REACTIVITY_ERROR"} = $reactivity_error;
    }
    return $self->{"REACTIVITY_ERROR"}; # returns an array reference to an array of arrays
}

sub serialize_reactivity_error{
    my $self = shift;
    my $line = "";
    my $reactivity_error_nr = 0;
    foreach my $href ( @{ $self->reactivity_error() }  ) {
        $reactivity_error_nr++;
        $line =  "REACTIVITY_ERROR:".$reactivity_error_nr;
        $line .= $self->_href_to_string( $href );
    }
    return $line;
}

sub seqpos_reactivity_map {
    my $self = shift;
    return $self->{"SEQPOS_REACTIVITY_MAP"};
}

sub _href_to_string{
    my ($self, $href) = @_;
    my $line = "";
    foreach my $key ( keys( %{ $href } ) ) {
        $line = " ".$key.":".$href->{$key};
    }
    return $line;
}

1;
