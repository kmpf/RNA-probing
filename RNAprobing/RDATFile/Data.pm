package RNAprobing::RDATFile::Data;

use Log::Log4perl qw(get_logger :levels);

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

sub new {
    my ($classname, $name, $lines) = @_;
    my $self = {};

    &indices($self);
    &annotation_data($self);
    &reactivity($self);
    &reactivity_error($self);
    &xsel_refine($self);
    &seqpos($self);
    &trace($self);
    &reads($self);

    bless $self, $classname;

    $self->read_data($lines) if ( defined $lines );

    return $self;
}

sub read_data {
    my ($self, $lines) = @_;
    die "No Input for RNAprobing::RDATFile::Name->create_name()" 
        unless ( defined $lines );
    foreach my $line ( split(/\n/, $lines) ) {
        chomp( $line );
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        my @words = split(/\s+/, $line, 2);
        my @split = split( /\s+/, $words[1]);

        # Data section lines
        if ($words[0] =~ /^ANNOTATION_DATA:(\d+)$/) {
            $self->annotation_data( $1, $words[1] );
        }
        if ($words[0] =~ /^REACTIVITY:(\d+)$/) {
            $self->reactivity( $1, [@split] );
        }
        if ($words[0] =~ /^REACTIVITY_ERROR:(\d+)$/) {
            $self->reactivity_error( $1, [@split] );
        }
        if ($words[0] =~ /^XSEL_REFINE:(\d+)$/) {
            $self->xsel_refine( $1, [@split] );
        }
	if ($words[0] =~ /^SEQPOS:(\d+)$/) {
            $self->seqpos( $1, [@split] );
        }
        if ($words[0] =~ /^TRACE:(\d+)$/) {
            $self->trace( $1, [@split] );
        }
        if ($words[0] =~ /^READS:(\d+)$/) {
            $self->reads( $1, [@split] );
        }
    }
}

sub serialize_data {
    my ($self) = @_;
    my $data_content;
    foreach my $index ( $self->indices() ){
	$data_content .= $self->serialize_annotation_data($self).
	    $self->serialize_reactivity($self).
	    $self->serialize_reactivity_error($self).
	    $self->serialize_xsel_refine($self).
	    $self->serialize_seqpos($self).
	    $self->serialize_trace($self).
	    $self->serialize_reads($self);
    }
    return $data_content;
}


# Construct section

sub indices {
    my ($self, $index) = @_;
    my $method_key = "INDICES";
    
    # if the index is defined and not present add it to array and sort it afterwards 
    if ( defined $index ){
	if ( !( $index ~~ @{$self->{$method_key}} ) ) {
	    push ( @{$self->{$method_key}}, $index);
	    $self->{$method_key} = sort { $a <=> $b } @{$self->{$method_key}};
	}
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a array reference
}


sub annotation_data {
    my ($self,  $index, $anno_data) = @_;
    my $method_key = "Annotation_DATA";
    return if ( !(defined $index) );
    if ( defined $anno_data ){
        $self->{$method_key}->{$index} = $anno_data;
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_annotation_data {
    my ($self, $index) = @_;
    return if ( !(defined $index) );
    if ( $self->annotation_data($index) ){
        my $line = "ANNOTATION_DATA:".$index."\t".$self->annotation_data($index)."\n";
        return $line;
    } else {
	return "";
    }
}

sub reactivity {
    my ($self,  $index, $reactivity) = @_;
    my $method_key = "REACTIVITY";
    return if ( !(defined $index) );
    if ( defined $reactivity ){
        $self->{$method_key}->{$index} = $reactivity;
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = [];
    }
    return $self->{$method_key}->{$index}; # returns an array reference 
}

sub serialize_reactivity {
    my ($self, $index) = @_;
    return if ( !(defined $index) );
    if ( @{$self->reactivity($index)} ) {
	my $line = "REACTIVITY:".$index."\t";
	$line .= join( "\t", @{$self->reactivity($index)} )."\n";   
    } else {
        $logger->error("Missing mandatory REACTIVITY entry.");
        $logger->error("RDAT output is erroneous.");
    }
    return $line;
}

sub scaled_reactivity {
    my ($self, $reactivity) = @_;
    my $method_key = "SCALED_REACTIVITY";

    if ( defined $reactivity ) {
        foreach my $aref ( @{$reactivity} ) {
            my $max_reac = max( @{$aref} );
            my $min_reac = min( @{$aref} );
            my $reactivity_span = $max_reac - $min_reac;
            my @scaled_reac = ();
            # scale the entries in @{$aref}
            foreach my $reac_value (@{$aref}) {
                push( @scaled_reac, 
                      ($reac_value - $min_reac) / $reactivity_span );
                
            }
            push( @{$self->{$method_key}}, \@scaled_reac);
        }
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference to an array of
                                 # arrays or just an empty array
}

sub reactivity_error {
    my ($self, $index, $reactivity_error) = @_;
    my $method_key = "REACTIVITY_ERROR";
    return if ( !(defined $index) );
    if ( defined $reactivity_error ){
        $self->{$method_key}->{$index} = $reactivity_error;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = [];
    }
    return $self->{$method_key}->{$index}; # returns an array reference
}

sub serialize_reactivity_error{
    my ($self, $index) = @_;
    my $line = "";
    return if ( !(defined $index) );
    $line =  "REACTIVITY_ERROR:".$index."\t";
    $line .=  join( "\t", @{$self->reactivity_error($index)} )."\n";
    return $line;
}

sub xsel_refine {
    my ($self, $index, $xsel_refine) = @_;
    my $method_key = "XSEL_REFINE";
    return if ( !(defined $index) );
    if ( defined $xsel_refine ){
        $self->{$method_key}->{$index} = $xsel_refine;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = [];
    }
    return $self->{$method_key}->{$index}; # returns an array reference
}

sub serialize_xsel_refine{
    my ($self, $index) = @_;
    my $line = "";
    return if ( !(defined $index) );
    $line =  "XSEL_REFINE:".$index." ";
    $line .=  join( "\t", @{$self->xsel_refine($index)} )."\n";
    return $line;
}

sub trace {
    my ($self, $index, $trace) = @_;
    my $method_key = "TRACE";
    return if ( !(defined $index) );
    if ( defined $trace ){
        $self->{$method_key}->{$index} = $trace;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = [];
    }
    return $self->{$method_key}->{$index}; # returns an array reference
}

sub serialize_trace {
    my ($self, $index) = @_;
    my $line = "";
    return if ( !(defined $index) );
    $line =  "TRACE:".$index."\t";
    $line .=  join( "\t", @{$self->trace($index)} )."\n";
    return $line;
}

sub reads {
    my ($self, $index, $reads) = @_;
    my $method_key = "READS";
    return if ( !(defined $index) );
    if ( defined $reads ){
        $self->{$method_key}->{$index} = $reads;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key}->{$index} = [];
    }
    return $self->{$method_key}->{$index}; # returns an array reference
}

sub serialize_reads{
    my ($self, $index) = @_;
    my $line = "";
    return if ( !(defined $index) );
    $line =  "READS:".$index."\t";
    $line .=  join( "\t", @{$self->reads($index)} )."\n";
    return $line;
}

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
