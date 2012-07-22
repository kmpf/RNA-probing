package ProbingRNA::BLASTresult;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

sub new {
    my $self = {};
    my $classname = shift;
    my $logger = get_logger();

    $self->{"FILENAME"} = shift;
    $self->{"NR_OF_ENTRIES"}= 0;
    $self->{"QUERY_ID"} = [];
    $self->{"SUBJECT_ID"} = [];
    $self->{"PERCENT_IDENTITY"} = [];
    $self->{"ALIGNMENT_LENGTH"} = [];
    $self->{"MISMATCHES"} = [];
    $self->{"GAP_OPENINGS"} = [];
    $self->{"QUERY_START"} = [];
    $self->{"QUERY_END"} = [];
    $self->{"SUBJECT_START"} = [];
    $self->{"SUBJECT_END"} = [];
    $self->{"E-VALUE"} = [];
    $self->{"BIT-SCORE"} = [];
    
    my $results_found = 0;
    open ( BLAST_RESULT , "<", $self->{FILENAME}) or die "Couldn't open file $self->{FILENAME}. Error: $!";;
    while (my $line = <BLAST_RESULT>) {
        chomp( $line );
        my $string = '^#';
        next if ($line =~ /$string/);
        my @columns = split(/\s+/, $line);
        $logger->debug(scalar(@columns));
        if ( scalar(@columns) == 12 ) {
            $results_found++;
            $self->{"NR_OF_ENTRIES"}++;
            push( @{ $self->{"QUERY_ID"} }, $columns[0] );
            push( @{ $self->{"SUBJECT_ID"} }, $columns[1] );
            push( @{ $self->{"PERCENT_IDENTITY"} }, $columns[2] );
            push( @{ $self->{"ALIGNMENT_LENGTH"} }, $columns[3] );
            push( @{ $self->{"MISMATCHES"} }, $columns[4] );
            push( @{ $self->{"GAP_OPENINGS"} }, $columns[5] );
            push( @{ $self->{"QUERY_START"} }, $columns[6] );
            push( @{ $self->{"QUERY_END"} }, $columns[7] );
            push( @{ $self->{"SUBJECT_START"} }, $columns[8] );
            push( @{ $self->{"SUBJECT_END"} }, $columns[9] );
            push( @{ $self->{"E-VALUE"} }, $columns[10] );
            push( @{ $self->{"BIT-SCORE"} }, $columns[11] );
        } else {
            $logger->error("There are more than 12 columns in:\n$line");
            exit 1;
        }
        
    }
    close(BLAST_RESULT);

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
##  Getter/Setter subroutines 
##
###############################################################

sub filename {
    my $self = shift;
    return $self->{"FILENAME"};
}

sub nr_of_entries {
    my $self = shift;
    return $self->{"NR_OF_ENTRIES"};
}

sub query_id {
    my $self = shift;
    return $self->{"QUERY_ID"};
}

sub subject_id {
    my $self = shift;
    return $self->{"SUBJECT_ID"};
}

sub percent_identity {
    my $self = shift;
    return $self->{"PERCENT_IDENTITY"}; 
}

sub alignment_length {
    my $self = shift;
    return $self->{"ALIGNMENT_LENGTH"}; 
}

sub mismatches {
    my $self = shift;
    return $self->{"MISMATCHES"}; 
}

sub gap_openings {
    my $self = shift;
    return $self->{"GAP_OPENINGS"}; 
}

sub query_start {
    my $self = shift;
    return $self->{"QUERY_START"} ; 
}

sub query_end {
    my $self = shift;
    return $self->{"QUERY_END"}; 
}

sub subject_start {
    my $self = shift;
    return $self->{"SUBJECT_START"}; 
}

sub subject_end {
    my $self = shift;
    return $self->{"SUBJECT_END"}; 
}

sub e_value {
    my $self = shift;
    return $self->{"E-VALUE"};
}

sub bit_score {
    my $self = shift;
    return $self->{"BIT-SCORE"};
}

1;
