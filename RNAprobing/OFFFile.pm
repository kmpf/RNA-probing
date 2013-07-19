package RNAprobing::OFFFile;

use strict;
use warnings;
use Log::Log4perl qw(get_logger :levels);

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

sub new {
    my ($classname, $filename) = @_;
    my $self = {};
    my $logger = get_logger();

    &filename($self, $filename);
    &fasta_header($self);
    &fasta_id($self);
    &sequence($self);
    &structure($self);
    &column_sizes($self);
    &nr_of_edges($self);
    &edges($self);
 
    bless( $self, $classname );
    $self->read_file($filename) if ( defined $filename && -e $filename );
    return $self;
}

sub read_file {
    my ($self, $filename) = @_;
    my $logger = get_logger();
    my $line_number = 0;
    die "No file given to RNAprobing::OFFFile::read_file!" 
        unless (defined $filename );
    open ( my $off_file , "<", $self->filename()) or 
        die "Couldn't open file $self->filename(). Error: $!";
    my @edges = ();
    my $nr_of_edges = "";
    while (my $line = <$off_file>) {
        next if ($line =~ /^##/);
        $line_number++;
        chomp( $line );
        if ($line_number == 1 && $line =~ /^>/) {
            $self->fasta_header( $line );
            $logger->debug("Fasta header: ".$self->fasta_header());
            $line =~ s/>\s*(\S+).*/$1/g;
            $self->fasta_id( $line );
            $logger->debug("Fasta ID: ".$self->fasta_id() );
            next;
        } elsif ( $line_number == 1 ) {
            $logger->error("First line isn't a fasta header."
                ." It should start with '>'. Something is wrong.");
            exit 1;
        }
        if ($line_number == 2 && $line =~ /^[ACGUacgu]*$/) {
            $self->sequence( $line );
            $logger->debug("Sequence: ".$self->sequence() );
            next;
        } elsif ( $line_number == 2 ) {
            $logger->error("Sequence seems not to be an RNA sequence."
                ." ONLY those characters are allowed: A,C,G,U,a,c,g,u\n"
                ."Sequence: ".$line);
            exit 1;
        }
        if ($line_number == 3 && $line =~ /^[\.\(\)\[\]\{\}<>]*$/) {
            $self->structure( $line );
            $logger->debug("Dot-bracket string: ".$self->structure() );
            next;
        } elsif ( $line_number == 3 ) {
            $logger->error("Dot-bracket string contains other characters as one of those:\n"
                .".()[]{}<>");
            exit 1;
        }
        if ($line_number == 4 &&  $line =~ /^#\s(\d+;)*\d+/) {
            $line =~ s/^#\s//g;
            my @col_size = split(";", $line);
            $self->column_sizes(\@col_size);
            next;
        } elsif ($line_number == 4) {
            $logger->error("Line $line_number should start with '# ' followed "
                ."by the semicolon-separated values of the column sizes.");
        }
        if ( $line_number >= 5 && $line =~ /^#\s/ ) {
            $line =~ s/^#\s//g;
            push( @edges, [split(/\s+/, $line)] );
#            push(@{ $self->edges() }, [@split] );
            $nr_of_edges++;
        } elsif ( $line_number >= 5 && $line =~ /^[^#]/) {
            $logger->error("Line $line_number should start with '# '.");
        }
        
    }
    $self->edges(\@edges);
    $self->nr_of_edges($nr_of_edges);
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
    my ($self, $filename) = @_;
    my $method_key = "FILENAME";
    if ( defined $filename ){
        $self->{$method_key} = $filename;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub fasta_header {
    my ($self, $fasta_header) = @_;
    my $method_key = "FASTA_HEADER";
    if ( defined $fasta_header ){
        $self->{$method_key} = $fasta_header;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub fasta_id {
    my ($self, $fasta_id) = @_;
    my $method_key = "FASTA_ID";
    if ( defined $fasta_id ){
        $self->{$method_key} = $fasta_id;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub sequence {
    my ($self, $sequence) = @_;
    my $method_key = "SEQUENCE";
    if ( defined $sequence ){
        $self->{$method_key} = $sequence;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub structure {
    my ($self, $structure) = @_;
    my $method_key = "STRUCTURE";
    if ( defined $structure ){
        $self->{$method_key} = $structure;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar (hopefully in dot-bracket notation)
}

sub column_sizes {
    my ($self, $column_sizes) = @_;
    my $method_key = "COLUMN_SIZES";
    if ( defined $column_sizes ){
        $self->{$method_key} = $column_sizes;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub nr_of_edges {
    my ($self, $nr_of_edges) = @_;
    my $method_key = "NR_OF_EDGES";
    if ( defined $nr_of_edges ){
        $self->{$method_key} = $nr_of_edges;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }

    return $self->{$method_key};
}

sub edges {
    my ($self, $edges) = @_;
    my $method_key = "EDGES";
    if ( defined $edges ){
        $self->{$method_key} = $edges;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub sequence_one_indexed_map {
    my ($self) = @_;
    my $method_key = "SEQUENCE_ONE_INDEXED_MAP";
    return $self->{$method_key} if ( defined $self->{$method_key} );
    if ( defined $self->sequence() ){
	    my @sequence = split(//, $self->sequence());
	    # 1 is used because sequence is 1-indexed
	    my $startpos = 1;
	    my $endpos = scalar(@sequence);
	    $logger->info("Sequence length = ".scalar(@sequence));
	    $logger->info("1-indexed start position = ".$startpos);
	    $logger->info("1-indexed end position = ".$endpos);
	    my $seq_oi_map = {};
	    my $seq_index = 0;
    	for (my $i = $startpos; $i <= $endpos; $i++ ) {
    	    $seq_oi_map->{$i} = $sequence[$seq_index];
    	    $seq_index++;
    	}
        $self->{$method_key} = $seq_oi_map;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = {};
    }
    return $self->{$method_key}; # returns a hash reference
}


1;
