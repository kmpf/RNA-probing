package RNAprobing::RDATFile::Annotation;

use Log::Log4perl qw(get_logger :levels);


sub new {
    my ($classname, $line) = @_;
    my $self = {};

    my @chemicals = (); # array for all chemicals that have been used
    my @processings = (); # array for all processings that have been done
    my $logger = get_logger();

    my @words = ();
    @words = split(/\s+/, $line) if (defined $line);
    for (my $i = 0; $i < $#words; $i++) {
        my $word = $words[$i];
        &experiment_type( $self, $1 )
            if ($word =~ /^experimentType:(\S+)/ && splice(@words, $i, 1));
        &modifier( $self, $1 )
            if ($word =~ /^modifier:(\S+)/ && splice(@words, $i, 1) );
        &mutation( $self, $1 )
            if ($word =~ /^mutation:(\S+)/ && splice(@words, $i, 1) );
        &processings($self, $1) 
            if ($word =~ /^processing:(\S+)/ && splice(@words, $i, 1) );
        &temperature( $self, $1 )
            if ($word =~ /^temperature:(\S+)/ && splice(@words, $i, 1) );
        &chemicals( $self, $1, $2) 
            if ($word =~ /^chemical:(\S+):(\S+)/ && splice(@words, $i, 1) );

    }
    $logger->error("Undigestible rest: ". join(" ",@words) ."\n") if ( scalar(@words) != 0 );
    &chemicals( $self, \@chemicals ) if ( $#chemicals > 0 );
    &processings( $self, \@processings ) if ( $#processings > 0 );

    bless $self, $classname;
}

sub serialize_annotation {
    my ($self) = (@_);
    my $annotation = " ";
    $annotation .= $self->serialize_experiment_type();
    $annotation .= $self->serialize_modifier();
    $annotation .= $self->serialize_mutation();
    $annotation .= $self->serialize_temperature();
    $annotation .= $self->serialize_chemicals();
    $annotation .= $self->serialize_processings();

}

sub experiment_type {
    my ($self, $experiment_type) = @_;
    my $method_key = "EXPERIMENT_TYPE";

    if ( defined $experiment_type){
        push( @{$self->{$method_key}}, $experiment_type );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_experiment_type {
    my $self = shift;
    my $line = "";
    foreach my $experiment_type ( @{$self->experiment_type()} ){
        $line .= "experimentType:".$experimenmt_type." ";
    }
    return $line;
}

sub modifier {
    my ($self, $modifier) = @_;
    my $method_key = "MODIFIER";

    if ( defined $modifier){
        push( @{$self->{$method_key}}, $modifier );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_modifier {
    my $self = shift;
    my $line = "";
    foreach my $modifier ( @{$self->modifier()} ){
        $line .= "modifier:".$modifier." ";
    }
    return $line;
}

sub mutation {
    my ($self, $mutation) = @_;
    my $method_key = "MUTATION";

    if ( defined $mutation){
        push( @{$self->{$method_key}}, $mutation );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_mutation {
    my $self = shift;
    my $line = "";
    foreach my $mutation ( @{$self->mutation()} ){
        $line .= "mutation:".$mutation." ";
    }
    return $line;
}

sub temperature {
    my ($self, $temperature) = @_;
    my $method_key = "TEMPERATURE";

    if ( defined $temperature){
        push( @{$self->{$method_key}}, $temperature );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_temperature {
    my $self = shift;
    my $line = "";
    foreach my $temperature ( @{$self->temperature()} ){
        $line .= "temperature:".$temperature." ";
    }
    return $line;
}

sub chemicals {
    my ($self, $chemical, $conc) = @_;
    my $method_key = "CHEMICALS";

    if ( defined $chemical && defined $conc){
        push( @{$self->{$method_key}}, ( $chemical => $conc) );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a reference to an array of hashes
}

sub serialize_chemicals {
    my $self = shift;
    my $line = "";
#    if ( defined $self->chemicals() ) {
        foreach my $chemical ( @{$self->chemicals()} ) {
            $line .= "chemical:".$_." ";
        }
#    }
    return $line;
}

sub processings {
    my ($self, $processing) = @_;
    my $method_key = "PROCESSINGS";

    if ( defined $processing){
        push( @{$self->{$method_key}}, $processing );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub serialize_processings {
    my $self = shift;
    my $line = "";
    foreach my $processing ( @{$self->processings()} ) {
        $line .= "chemical:".$processing." ";
    }
    return $line;
}


################################################################################
##
##  Subroutine that parses ANNOTATION and ANNOTATION_DATA lines
##
################################################################################

sub _parse_annotation {
    my ($self, $line) = @_;

}


1;
