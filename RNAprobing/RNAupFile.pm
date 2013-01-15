package RNAprobing::RNAupFile;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

sub new {
    my $classname = shift;
    my $filename = shift;
    my $self = {};
    my $logger = get_logger();

    &filename();
    $self->{"RNAUP_COMMAND"} = "";
    $self->{"NUCLEOTIDE_POSITIONS"} = [];
    $self->{"RNAUP_VALUES"} = [];

    my %probabilities = ();

    open ( RNAUP_FILE , "<$self->{FILENAME}") or die "Couldn't open file $self->{FILENAME}. Error: $!";;
    my $comment = '^#|^\s*$';
    my $rnaup_cmd = '^# RNAup';
    my @columns = ();
    while (my $line = <RNAUP_FILE>) {
        chomp( $line );
        $line =~ s/^\s+//g;
        if ( $line =~ /$rnaup_cmd/ ) {
            $line =~ s/^#\s*//g;
            $self->{"RNAUP_COMMAND"} = $line;
            next;
        } elsif ( $line =~ /$comment/ ) {
            next;
        }
        @columns = split(/\s+/, $line);
        $probabilities{$columns[0]} = $columns[1];
    }
    
    my @nucleotide_position = sort {$a <=> $b} keys(%probabilities);
    my @rnaup_values = ();
    foreach my $position ( @nucleotide_position ) {
        push(@rnaup_values, $probabilities{$position} );
    }
    

    $self->{"NUCLEOTIDE_POSITIONS"} = \@nucleotide_position;
    $self->{"RNAUP_VALUES"} = \@rnaup_values;

    close(RNAUP_FILE);
    
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

sub nucleotide_positions {
    my $self = shift;
    return $self->{"NUCLEOTIDE_POSITIONS"}; # returns a scalar
}

sub rnaup_values {
    my $self = shift;
    return $self->{"RNAUP_VALUES"}; # returns a scalar
}

1;
