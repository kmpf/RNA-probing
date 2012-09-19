package ProbingRNA::OFFFile;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

sub new {
    my $self = {};
    my $classname = shift;
    my $logger = get_logger();

    $self->{"FILENAME"} = shift;
    $self->{"FASTA_HEADER"} = "";
    $self->{"FASTA_ID"} = "";
    $self->{"SEQUENCE"} = "";
    $self->{"STRUCTURE"} = "";
    $self->{"COLUMN_SIZES"} = [];
    $self->{"NR_OF_EDGES"} = 0;
    $self->{"EDGES"} = [];

    $self->read_file() if ( $self->{"FILENAME"} );
    
    bless( $self, $classname );
    return $self;
}

sub read_file {
    my $line_number = 0;
    
    open ( my $off_file , "<", $self->{FILENAME}) or die "Couldn't open file $self->{FILENAME}. Error: $!";
    while (my $line = <OFF_FILE>) {
        next if ($line =~ /^##/);
        $line_number++;
        chomp( $line );
        if ($line_number == 1 && $line =~ /^>/) {
            $self->{"FASTA_HEADER"} = $line;
            $logger->debug("Fasta header: ".$self->{"FASTA_HEADER"});
            $line =~ s/>\s*(\S+).*/$1/g;
            $self->{"FASTA_ID"} = $line;
            $logger->debug("Fasta ID: ".$self->{"FASTA_ID"});
            next;
        } elsif ( $line_number == 1 ) {
            $logger->error("First line isn't a fasta header. It should start with '>'. Something is wrong.");
            exit 1;
        }
        if ($line_number == 2 && $line =~ /^[ACGUacgu]*$/) {
            $self->{"SEQUENCE"} = $line;
            $logger->debug("Sequence: ".$self->{"SEQUENCE"});
            next;
        } elsif ( $line_number == 2 ) {
            $logger->error("Sequence seems not to be an RNA sequence. ONLY those characters are allowed: A,C,G,U,a,c,g,u");
            $logger->error("Sequence: ".$line);
            exit 1;
        }
        if ($line_number == 3 && $line =~ /^[\.\(\)\[\]\{\}<>]*$/) {
            $self->{"STRUCTURE"} = $line;
            $logger->debug("Dot-bracket string: ".$self->{"STRUCTURE"});
            next;
        } elsif ( $line_number == 3 ) {
            $logger->error("Dot-bracket string contains other characters as one of those: .()[]{}<>");
            exit 1;
        }
        if ($line_number == 4 &&  $line =~ /^#\s(\d+;)*\d+/) {
            $line =~ s/^#\s//g;
            push( @{ $self->{"COLUMN_SIZES"} }, split(";", $line) );
            next;
        } elsif ($line_number == 4) {
            $logger->error("Line $line_number should start with '# ' and afterwards semicolon-separated values of the column sizes.");
        }
        if ( $line_number >= 5 && $line =~ /^#\s/ ) {
            $line =~ s/^#\s//g;
            my @split = split(/\s+/, $line) ;
            push(@{ $self->{"EDGES"} }, [@split] );
            $self->{"NR_OF_EDGES"}++;
        } elsif ( $line_number >= 5 && $line =~ /^[^#]/) {
            $logger->error("Line $line_number should start with '# '.");
        }
        
    }
    close($off_file);
}

sub write_file {

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
    return $self->{"FILENAME"}; # returns a scalar
}

sub fasta_header {
    my $self = shift;
    return $self->{"FASTA_HEADER"}; # returns a scalar
}

sub fasta_id {
    my $self = shift;
    return $self->{"FASTA_ID"}; # returns a scalar
}

sub sequence {
    my $self = shift;
    return $self->{"SEQUENCE"}; # returns a scalar
}

sub structure {
    my $self = shift;
    return $self->{"STRUCTURE"}; # returns a scalar (hopefully in dot-bracket notation)
}

sub column_sizes {
    my $self = shift;
    return $self->{"COLUMN_SIZES"};
}

sub nr_of_edges {
    my $self = shift;
    return $self->{"NR_OF_EDGES"};
}

sub edges {
    my $self = shift;
    return $self->{"EDGES"};
}

1;
