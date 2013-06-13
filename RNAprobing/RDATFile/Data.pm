package RNAprobing::RDATFile::Data;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

sub new {
    my ($classname, $lines) = @_;
    my $self = {};

    &indices($self);
    &annotation_data($self);
    &reactivity($self);
    &reactivity_error($self);
    &xsel_refine($self);
    &seqpos($self);
    &trace($self);
    &reads($self);

    &scaled_reactivity($self);
    &seqpos_reactivity_map($self);
    &seqpos_scaled_reactivity_map($self);
    bless $self, $classname;

    $self->read_data($lines) if ( defined $lines );

    return $self;
}

sub read_data {
    my ($self, $lines) = @_;
    return if ( !(defined $lines) );
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
	    $self->scaled_reactivity($1, $self->reactivity($1) );
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
    my ($anno_data, $reac, $reac_error, $xsel_refine, $seqpos, $trace, $reads);
    $anno_data = $reac = $reac_error = $xsel_refine =
	$seqpos = $trace = $reads = "";
    foreach my $index ( @{$self->indices()} ){
	$anno_data   .= $self->serialize_annotation_data($index);
	$reac        .= $self->serialize_reactivity($index);
	$reac_error  .= $self->serialize_reactivity_error($index);
	$xsel_refine .= $self->serialize_xsel_refine($index);
	$seqpos      .= $self->serialize_seqpos($index);
	$trace       .= $self->serialize_trace($index);
	$reads       .= $self->serialize_reads($index);
    }
    return $anno_data.$reac.$reac_error.$xsel_refine.$seqpos.$trace.$reads;
}


# Construct section

sub indices {
    my ($self, $index) = @_;
    my $method_key = "INDICES";
    
    # if the index is defined and not present in the index array
    # add it to the array and sort the array afterwards 
    if ( defined $index && !( $index ~~ @{$self->{$method_key}} ) ) {
	push ( @{$self->{$method_key}}, $index);
	@{$self->{$method_key}} = sort { $a <=> $b } @{$self->{$method_key}};
    } elsif ( !(defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a array reference
}


sub annotation_data {
    my ($self,  $index, $anno_data) = @_;
    my $method_key = "ANNOTATION_DATA";
    if ( defined $index ){
	$self->indices( $index );
	if ( defined $anno_data ){
	    $self->{$method_key}->{$index} = $anno_data;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = "";
	}
	return $self->{$method_key}->{$index}; # returns a scalar
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_annotation_data {
    my ($self, $index) = @_;
    if ( defined $index && $self->annotation_data($index) ){
        return "ANNOTATION_DATA:".$index."\t".$self->annotation_data($index)."\n";
    } else {
	return "";
    }
}

sub reactivity {
    my ($self,  $index, $reactivity) = @_;
    my $method_key = "REACTIVITY";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$reactivity} ){
	    $self->{$method_key}->{$index} = $reactivity;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_reactivity {
    my ($self, $index) = @_;
    if ( defined $index && @{$self->reactivity($index)} ) {
	return "REACTIVITY:".$index."\t".
	    join( "\t", @{$self->reactivity($index)} )."\n";
    }
    $logger->error("Missing mandatory REACTIVITY entry.");
    $logger->error("RDAT output is erroneous.");
    return "";
}

sub scaled_reactivity {
    my ($self, $index, $reactivity) = @_;
    my $method_key = "SCALED_REACTIVITY";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$reactivity} ) {
	    my $max_reac = max( @{$reactivity} );
	    my $min_reac = min( @{$reactivity} );
	    my $reactivity_span = $max_reac - $min_reac;
	    my @scaled_reac = ();
	    # scale the entries in @{$aref}
	    foreach my $reac_value ( @{$reactivity} ) {
		push( @scaled_reac, 
		      ($reac_value - $min_reac) / $reactivity_span );
	    }
	    $self->{$method_key}->{$index} = \@scaled_reac;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub reactivity_error {
    my ($self, $index, $reactivity_error) = @_;
    my $method_key = "REACTIVITY_ERROR";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$reactivity_error} ){
	    $self->{$method_key}->{$index} = $reactivity_error;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_reactivity_error{
    my ($self, $index) = @_;
    if ( defined $index && @{$self->reactivity_error($index)} ){
	return  "REACTIVITY_ERROR:".$index."\t".
	    join( "\t", @{$self->reactivity_error($index)} )."\n";
    }
    return "";
}

sub xsel_refine {
    my ($self, $index, $xsel_refine) = @_;
    my $method_key = "XSEL_REFINE";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$xsel_refine} ){
	    $self->{$method_key}->{$index} = $xsel_refine;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_xsel_refine{
    my ($self, $index) = @_;
    if ( defined $index && @{$self->xsel_refine($index)} ) {
	return "XSEL_REFINE:".$index."\t".
	    join( "\t", @{$self->xsel_refine($index)} )."\n";
    }
    return "";
}

sub seqpos {
    my ($self, $index, $seqpos) = @_;
    my $method_key = "SEQPOS";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$seqpos} ){
	    $self->{$method_key}->{$index} = $seqpos;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_seqpos {
    my ($self, $index) = @_;
    my $line = "";
    if ( defined $index && @{$self->seqpos($index)} ) {
	return "SEQPOS:".$index."\t".
	    join( "\t", @{$self->seqpos($index)} )."\n";
    }
    return "";
}

sub trace {
    my ($self, $index, $trace) = @_;
    my $method_key = "TRACE";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$trace} ){
	    $self->{$method_key}->{$index} = $trace;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_trace {
    my ($self, $index) = @_;
    if ( defined $index && @{$self->trace($index)}) {
	return "TRACE:".$index."\t".
	     join( "\t", @{$self->trace($index)} )."\n";
    }
    return "";
}

sub reads {
    my ($self, $index, $reads) = @_;
    my $method_key = "READS";
    if ( defined $index ){
	$self->indices( $index );
	if ( @{$reads} ){
	    $self->{$method_key}->{$index} = $reads;
	}  elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub serialize_reads{
    my ($self, $index) = @_;
    if ( defined $index && @{$self->reads($index)} ) {
	return "READS:".$index."\t".
	    join( "\t", @{$self->reads($index)} )."\n";
    }
    return "";
}

sub seqpos_reactivity_map {
    my ($self, $index, $seqpos) = @_;
    my $method_key = "SEQPOS_REACTIVITY_MAP";
    if ( defined $index ) {
	$self->indices( $index );
	$seqpos = $self->seqpos($index) if ( @{$self->seqpos($index)} );
	if ( scalar(@{$seqpos}) == scalar(@{$self->reactivity($index)}) ) {
	    $self->{$method_key}->{$index} =
		$self->_create_seqpos_reactivity_hash($seqpos, 
						      $self->reactivity($index));
	} elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
}

sub seqpos_scaled_reactivity_map {
    my ($self, $index, $seqpos) = @_;
    my $method_key = "SEQPOS_SCALED_REACTIVITY_MAP";
    if ( defined $index ) {
	$self->indices( $index );
	$seqpos = $self->seqpos($index) if ( @{$self->seqpos($index)} );
	if ( scalar(@{$seqpos}) == scalar(@{$self->scaled_reactivity($index)}) ) {
	    $self->{$method_key}->{$index} =
		$self->_create_seqpos_reactivity_hash($seqpos, 
						      $self->scaled_reactivity($index));
	} elsif ( !(defined $self->{$method_key}->{$index}) ) {
	    $self->{$method_key}->{$index} = [];
	}
	return $self->{$method_key}->{$index}; # returns an array reference
    } else {
	if ( !(defined $self->{$method_key}) ) {
	    $self->{$method_key} = {};
	}
	return $self->{$method_key};
    }
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
            $logger->error(
		"SEQPOS and REACTIVITY have unequal number of entries.");
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
