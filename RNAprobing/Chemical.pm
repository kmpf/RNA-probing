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
    &probe_mod($self);
    &probe_mod_pos($self);

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
    my ($prob_reac, @prob_seq, @prob_str, @prob_mod);
    while (my $line = <$reac_file>) {
        chomp( $line );
        # remove line comments and empty lines
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        $line =~ s/#.*$//g; # remove in-line comments
        $line =~ s/\s+$//; # remove trailing white spaces
        $line_nr++;
        if ($line_nr == 1 && 1 == ($line =~ s/^>//g) ) {
            $self->probe_name($line) ;
        } elsif ($line_nr == 2) {
            $prob_reac = $line;
            if (Scalar::Util::Numeric::isnum($prob_reac)) {
                $self->probe_reac($prob_reac);
            } else {
                $logger->error("[".$self->filename()."] Given reactivity ".
                               $prob_reac." is not a number.");
                exit 1;
            }
        } elsif ($line_nr == 3) {
            # get sequence specific info for probe
            if ($line =~ /^[ACGTUMRWSYKVHDBN_]*$/) {
                @prob_seq = split("_", $line);
            } else {
                $logger->error("[".$self->filename()."] ".$line.
                               " contains characters other than: ".
                               "A, C, G, T, U, M, R, W, S, Y, K, V, ".
                               "H, D, B, N, and _");
                exit 1;
            }
        } elsif ($line_nr == 4) {
            # get structure specific info for probe
            if ($line =~ /^[BHIMPU\(\)_]*$/) {
                @prob_str = split("_", $line);
            } else {
                $logger->error("[".$self->filename()."] ".$line.
                               " contains characters other than: ".
                               "B, H, I, M, P, U, (, ), and _");
                exit 1;
            }
        } elsif ($line_nr == 5) {
            # get modification postion info for probe
            if ($line =~ /^[\.|_]*$/) {
                @prob_mod = split("_", $line);
            } else {
                $logger->error("[".$self->filename()."] ".$line.
                               " contains characters other than: ".
                               "., |, and _");
                exit 1;
            }
            # now that we have read all infos check their validity
            # check that sequence, structure, and modification blocks defined in
            # chemical file have same number of entries ... 
            if ( grep { $_ == scalar(@prob_seq) } 
                 (scalar(@prob_str), scalar(@prob_mod)) ){

                my $entry_nr = scalar(@prob_seq);
                # ... and each entry has to be as long as the according entry in 
                # the other arrays ...
                for ( my $i = 0; $i < $entry_nr; $i++) {
                    # next if everything is fine
                    if ( grep { $_ == length($prob_seq[$i])}
                         (length($prob_str[$i]), length($prob_mod[$i])) ) {
                        next;
                    } else {
                        # ... otherwise complain and die ...
                        $logger->error("[".$self->filename()."]\n"
                                       .$prob_str[$i]."\n"
                                       .$prob_seq[$i]."\n"
                                       .$prob_mod[$i]."\n"
                                       ."must have equal length.");
                        exit 1;
                    }
                }
            } else {
                # ... or complain and die
                $logger->error("[".$self->filename()."]\n"
                               .join("_", @prob_seq)."\n"
                               .join("_", @prob_str)."\n"
                               .join("_", @prob_mod)."\n"
                               ."Each line has to contain equal number of '_'.");
                exit 1;                               
            }
            $self->probe_seq(@prob_seq);
            $self->probe_str(@prob_str);
            $self->probe_mod(@prob_mod);
            $self->probe_mod_pos(@prob_mod);
            # Reached end of a probing definition, so prepare for next one
            $line_nr = 1;
            undef @prob_seq;
            undef @prob_str;
            undef @prob_mod;
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
    my ($self, @probe_seq) = @_;
    my $method_key = "PROBE_SEQENCE";

    if ( @probe_seq){
        push( @{$self->{$method_key}}, \@probe_seq );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a reference to an array of arrays
}

sub probe_str{
    my ($self, @probe_str) = @_;
    my $method_key = "PROBE_SEC_STRUCTURE";
    if ( @probe_str){
        push( $self->{$method_key}, \@probe_str );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key};  # returns a reference to an array of arrays
}

sub probe_mod{
    my ($self, @probe_mod) = @_;
    my $method_key = "PROBE_MOD";
    if ( @probe_mod){
        push( @{$self->{$method_key}}, \@probe_mod );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a reference to an array of arrays
}

# Creates an array of arrays containing the positions of the modification point
# set in line #5 per block
sub probe_mod_pos {
    my ($self, @probe_mod) = @_;
    my $method_key = "PROBE_MOD_POS";
    if ( @probe_mod){
        my @matches;
        foreach ( @probe_mod ){
            my @ret;
            while ($_ =~ /(?=\|)/g) {
                push(@ret,  $-[0]);
            }
            push(@matches, \@ret);
        }
        push( @{$self->{$method_key}}, \@matches );
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = [];
    }
    return $self->{$method_key}; # returns a reference to an array of arrays
}

1;
