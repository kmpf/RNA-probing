package RNAprobing::Chemical;

use strict;
use warnings;
use Carp;
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;
use Scalar::Util::Numeric;

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

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

sub read_file {
    my ($self, $filename) = @_;
    $self->filename($filename);
    open ( my $reac_file , "<", $self->filename() ) or 
        die("Couldn't open file ".$self->filename().". Error: $!");
    my $line_nr = 0;
    while (my $line = <$reac_file>) {
        chomp( $line );
        # remove line comments and empty lines
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        $line =~ s/#.*$//g; # remove in-line comments
        $line =~ s/\s+$//; # remove trailing white spaces
        my ($prob_reac, @prob_seq, @prob_str);
        $line_nr++;
        if ($line_nr == 1 && 1 == ($line =~ s/^>//g) ) {
            $self->probe_name($line) ;
        } elsif ($line_nr == 2) {
            $prob_reac = $line;
        } elsif ($line_nr == 3) {
            @prob_seq = split("_", $line);
        } elsif ($line_nr == 4) {
            @prob_str = split("_", $line);
        }
        # get modification position and reset counter for next block
        if ($line_nr == 5) {
            if (Scalar::Util::Numeric::isnum($prob_reac)) {
                $self->probe_reac($prob_reac);
            } else {
                $logger->error("[".$self->filename()."] Given reactivity ".
                               $prob_reac." is not a number.");
                exit 1;
            }
            if (scalar(@prob_seq) == scalar(@prob_str) ) {
                for (my $i = 0; $i < scalar(@prob_seq); $i++) {
                    if ( length($prob_seq[$i]) != length($prob_str[$i]) ) {
                        $logger->error("[".$self->filename()."] ".$prob_str[$i].
                            " and ".$prob_seq[$i]." have unequal length.");
                        exit 1;
                    }
                }
                $self->probe_seq(\@prob_seq);
                $self->probe_str(\@prob_str) ;
            } else {
                $logger->error("[".$self->filename()."] Sequence ".
                               join("_", @prob_seq)." and structure ".join("_", @prob_str).
                               " have unequal number of blocks.");
                exit 1;                               
            }
            $self->probe_cut($line); 
            $line_nr = 1;
        }
                
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
        $probe_name =~ s/\s//g;
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
    return $self->{$method_key}; # returns a refernce to an array of arrays
}

sub probe_str{
    my ($self, $probe_str) = @_;
    my $method_key = "PROBE_SEC_STRUCTURE";
    if ( defined $probe_str){
        push( $self->{$method_key}, $probe_str );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key};  # returns a refernce to an array of arrays
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
