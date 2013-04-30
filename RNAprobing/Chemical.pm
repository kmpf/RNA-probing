package RNAprobing::Chemical;

use strict;
use warnings;
use Carp;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

sub new {
    my ($classname, $filename) = @_;
    my $self = {};

    &filename($self, $filename);
    &probe_name($self);
    &probe_reac($self);
    &probe_seq($self);
    &probe_str($self);
    &probe_cut($self);

    bless( $self, $classname );

    $self->read_reactivity_file($filename) 
        if (defined $filename && -e $filename );
    return $self;
}

################################################################################
################################################################################
##
##  Subroutine section
##
################################################################################
################################################################################

################################################################################
##
##  Subroutine that reads in a reactivity file
##
################################################################################

sub read_reactivity_file {
    my ($self, $filename) = @_;
    $self->filename($filename);
    open ( my $reac_file , "<", $self->filename() ) or croak "Couldn't open file $self->filename(). Error: $!";
    my $line_nr = 0;
    while (my $line = <$reac_file>) {
        chomp( $line );
        # remove line comments and empty lines
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        $line =~ s/#.*$//g; # remove in-line comments
        $line =~ s/\s+$//;
        $line_nr++;
        $self->probe_name($line) if ($line_nr == 1 && 1 == ($line =~ s/^>//g) );
        $self->probe_reac($line) if ($line_nr == 2 && $line <= 1);
        $self->probe_seq($line) if ($line_nr == 3); # könnte man noch auf zugelassene Nucleotide prüfen
        $self->probe_str($line) if ($line_nr == 4); # könnte man noch auf zugelassene Strukturelemente prüfen
        # get modification position and reset counter for next block
        if ($line_nr == 5) {$self->probe_cut($line); $line_nr = 1;}
        
        
    }
}

################################################################################
##
##  Getter/Setter subroutines 
##
################################################################################

sub filename{
    my ($self, $filename) = @_;
    my $method_key = "FILE_NAME";
    if ( defined $filename){
        $self->{$method_key} = $filename;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar
}

sub probe_name{
    my ($self, $probe_name) = @_;
    my $method_key = "PROBE_NAME";
    if ( defined $probe_name){
        $self->{$method_key} = $probe_name;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "";
    }
    return $self->{$method_key}; # returns a scalar

}

sub probe_reac{
    my ($self, $probe_reac) = @_;
    my $method_key = "PROBE_REACTIVITY";
    if ( defined $probe_reac){
        push @{$self->{$method_key}}, $probe_reac;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference

}

sub probe_seq{
    my ($self, $probe_seq) = @_;
    my $method_key = "PROBE_SEQENCE";

    if ( defined $probe_seq){
        push( @{$self->{$method_key}}, $probe_seq );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub probe_str{
    my ($self, $probe_str) = @_;
    my $method_key = "PROBE_SEC_STRUCTURE";
    if ( defined $probe_str){
        push( $self->{$method_key}, $probe_str );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

sub probe_cut{
    my ($self, $probe_cut) = @_;
    my $method_key = "PROBE_CUT";
    if ( defined $probe_cut){
        push( @{$self->{$method_key}}, $probe_cut );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns an array reference
}

1;
