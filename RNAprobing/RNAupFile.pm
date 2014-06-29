package RNAprobing::RNAupFile;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

sub new {
    my ($classname, $filename) = @_;
    my $self = {};
    my $logger = get_logger();

    &filename($self, $filename);
    &nucleotide_positions($self);
    &rnaup_values($self);
    &rnaup_command($self);

    bless( $self, $classname );
    if ( defined $filename ) {
        if ( -e $filename ) {
            $self->read_file($filename);
        } else {
            die "File ".$filename." does not exist.";
        }
    }
    return $self;
}

###############################################################
##
##  Subroutine section
##
###############################################################

################################################################################
##
##  Subroutine that reads in a RDAT file
##
################################################################################

sub read_file {
    my ($self, $filename) = @_;
    print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
    my $logger = get_logger();
    $self->filename( $filename );
 
    open ( my $rnaup_file , "<", $self->filename() ) or die "Couldn't open file $self->filename(). Error: $!";;
    my $comment = '^#|^\s*$';
    my $rnaup_cmd = '^# RNAup';
    my @columns = ();
    my %probabilities = ();
    while (my $line = <$rnaup_file>) {
        chomp( $line );
        $line =~ s/^\s+//g;
        if ( $line =~ /$rnaup_cmd/ ) {
            $line =~ s/^#\s*//g;
            $self->rnaup_command($line);
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
    $self->nucleotide_positions( \@nucleotide_position );
    $self->rnaup_values( \@rnaup_values );

    close($rnaup_file);

    return $self;
}
###############################################################
##
##  Getter/Setter subroutines 
##
###############################################################

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

sub nucleotide_positions {
    my ($self, $nuc_pos) = @_;
    my $method_key = "NUCLEOTIDE_POSITONS";
    if (defined $nuc_pos){
        $self->{$method_key} = $nuc_pos;
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a array ref
}

sub rnaup_values {
    my ($self, $rnaup_val) = @_;
    my $method_key = "RNAUP_VALUES";
    if (defined $rnaup_val){
        $self->{$method_key} = $rnaup_val;
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a array ref
}

sub rnaup_command {
    my ($self, $rnaup_command) = @_;
    my $method_key = "RNAUP_COMMAND";
    if (defined $rnaup_command){
        $self->{$method_key} = $rnaup_command;
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

1;
