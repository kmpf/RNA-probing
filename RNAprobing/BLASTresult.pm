package RNAprobing::BLASTresult;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

sub new {
    my $classname = shift;
    my $filename = shift;
    my $self = {};
    my $logger = get_logger();

    &filename($self);
    &nr_of_entries($self);
    &query_id($self);
    &subject_id($self);
    &percent_identity($self);
    &alignment_length($self);
    &mismatches($self);
    &gap_openings($self);
    &query_start($self);
    &query_end($self);
    &subject_start($self);
    &subject_end($self);
    &e_value($self);
    &bit_score($self);

    bless( $self, $classname );
    $self->read_file($filename);
    return $self;
}

sub read_file {
    my ($self, $filename) = @_;
    my $logger = get_logger();
    my $results_found = 0;
    my @query_ids = ();
    my @subject_ids = ();
    my @percent_identities = ();
    my @alignment_lengths = ();
    my @mismatches = ();
    my @gap_openings = ();
    my @query_starts = ();
    my @query_ends = ();
    my @subject_starts = ();
    my @subject_ends = ();
    my @e_values = ();
    my @bitscores = ();

    open( my $blast_result , "<", $self->filename($filename) ) or
        die "Couldn't open file $self->filename(). Error: $!";
    while (my $line = <$blast_result>) {
        chomp( $line );
        my $string = '^#';
        next if ($line =~ /$string/);
        my @columns = split(/\s+/, $line);
        $logger->debug(scalar(@columns));
        if ( scalar(@columns) == 12 ) {
            $results_found++;
            push( @query_ids, $columns[0] );
            push( @subject_ids, $columns[1] );
            push( @percent_identities, $columns[2] );
            push( @alignment_lengths, $columns[3] );
            push( @mismatches, $columns[4] );
            push( @gap_openings, $columns[5] );
            push( @query_starts, $columns[6] );
            push( @query_ends, $columns[7] );
            push( @subject_starts, $columns[8] );
            push( @subject_ends, $columns[9] );
            push( @e_value, $columns[10] );
            push( @bit_score, $columns[11] );
        } else {
            $logger->error("There are more than 12 columns in:\n$line");
            exit 1;
        }
    }
    close( $blast_result );

    $self->nr_of_entries()
    $self->query_id(\@query_ids);
    $self->subject_id(\@subject_ids);
    $self->percent_identity(\@percent_identities);
    $self->alignment_length(\@alignment_lengths);
    $self->mismatches(\@mismatches);
    $self->gap_openings(\@gap_openings);
    $self->query_start(\@query_starts);
    $self->query_end(\@query_ends);
    $self->subject_start(\@subject_starts);
    $self->subject_end(\@subject_ends);
    $self->e_value(\@e_value);
    $self->bit_score(\@bit_score);
}

###############################################################
##
##  Subroutine section
##
###############################################################


###############################################################
##
##  Getter/Setter subroutines
##
###############################################################

sub filename {
    my ($self, $filename) = @_;
    my $method_key = "FILENAME";
    if ( defined $filename){
        $self->{$method_key} = $filename;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub nr_of_entries {
    my ($self, $nr_of_entries) = @_;
    my $method_key = "NR_OF_ENTRIES";
    if ( defined $nr_of_entries){
        $self->{$method_key} = $nr_of_entries;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub query_id {
    my ($self, $query_id) = @_;
    my $method_key = "QUERY_ID";
    if ( defined $query_id){
        $self->{$method_key} = $query_id;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub subject_id {
    my ($self, $subject_id) = @_;
    my $method_key = "SUBJECT_ID";
    if ( defined $subject_id){
        $self->{$method_key} = $subject_id;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub percent_identity {
    my ($self, $percent_identity) = @_;
    my $method_key = "PERCENT_IDENTITY";
    if ( defined $percent_identity){
        $self->{$method_key} = $percent_identity;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub alignment_length {
    my ($self, $alignment_length) = @_;
    my $method_key = "ALIGNMENT_LENGTH";
    if ( defined $alignment_length){
        $self->{$method_key} = $alignment_length;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub mismatches {
    my ($self, $mismatches) = @_;
    my $method_key = "MISMATCHES";
    if ( defined $mismatches){
        $self->{$method_key} = $mismatches;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub gap_openings {
    my ($self, $gap_openings) = @_;
    my $method_key = "GAP_OPENINGS";
    if ( defined $gap_openings){
        $self->{$method_key} = $gap_openings;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub query_start {
    my ($self, $query_start) = @_;
    my $method_key = "QUERY_START";
    if ( defined $query_start){
        $self->{$method_key} = $query_start;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub query_end {
    my ($self, $query_end) = @_;
    my $method_key = "QUERY_END";
    if ( defined $query_end){
        $self->{$method_key} = $query_end;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub subject_start {
    my ($self, $subject_start) = @_;
    my $method_key = "SUBJECT_START";
    if ( defined $subject_start) {
        $self->{$method_key} = $subject_start;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub subject_end {
    my ($self, $subject_end) = @_;
    my $method_key = "SUBJECT_END";
    if ( defined $subject_end){
        $self->{$method_key} = $subject_end;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub e_value {
    my ($self, $e_value) = @_;
    my $method_key = "E-VALUE";
    if ( defined $e_value){
        $self->{$method_key} = $e_value;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

sub bit_score {
    my ($self, $bit_score) = @_;
    my $method_key = "BIT-SCORE";
    if ( defined $bit_score){
        $self->{$method_key} = $bit_score;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns an array reference
}

1;
