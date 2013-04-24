package RNAprobing::RDATFile;

use strict;
use warnings;
use Carp;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;
require RNAprobing::RDATFile::Annotation;
require RNAprobing::RDATFile::Name;

print "Logging is not initialized!" unless (Log::Log4perl::initialized() );
my $logger = get_logger();

sub new {
    my ($classname, $filename) = @_;
    my $self = {};

    &filename($self, $filename);
    &rdat_version($self);
    &name($self);

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
##  Subroutine that reads in a RDAT file
##
################################################################################

sub read_file {
    my ($self, $filename) = @_;
    $self->filename( $filename );
    
    open ( my $rdat_file , "<", $self->filename() ) or croak "Couldn't open file $self->filename(). Error: $!";

    my ($name, $lines) = "";
    
    while (my $line = <$rdat_file>) {
        chomp( $line );
        my $comment = '^#|^\s*$';
        next if ($line =~ /$comment/ );
        my @words = split(/\s+/, $line, 2);
        my @split = split( /\s+/, $words[1]);
        $self->rdat_version( $words[1] ) if ($words[0] =~ /^RDAT_VERSION$/);
        if ( $name eq "" && $words[0] =~ /^NAME$/ ){
            $name = $words[1];
            next;
        } elsif ( $name ne "" && $words[0] !~ /^NAME$/ ){
            $lines .= $line."\n";
        } elsif ( $name ne "" && $words[0] =~ /^NAME$/ ){
            
        }        
    }
    close($rdat_file);

    $self->seqpos_reactivity_map( $self->seqpos(), $self->reactivity() );
    $self->scaled_reactivity( $self->reactivity() );

    if ( $self->rdat_version() == 0.24 ) {
        $logger->info("Correct version.");
    } else {
        $logger->error("Incorrect Version. There may occur errors while parsing this file.");
    }

    return $self;
}

################################################################################
##
##  Write subroutine
##
################################################################################

sub write_file {
    my ($self, $filename) = @_;
    $self->filename($filename) if (defined $filename);
    my $file_content = $self->serialize_rdat_version();
    foreach my $key ( keys %{$self->name()} ) {
        $file_content .= $self->name()->{$key}->serialize_construct();
    }
    open( my $rdat_file, ">", $self->filename() ) or croak( "Couldn't open file". $self->filename() );
    print($rdat_file  $file_content);
    close($rdat_file);
    $logger->info($file_content);
}

################################################################################
##
##  Getter/Setter subroutines 
##
################################################################################

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

sub rdat_version {
    my ($self, $rdat_version) = @_;
    my $method_key = "RDAT_VERSION";
    
    if ( defined $rdat_version){
        $self->{$method_key} = $rdat_version;
    } elsif ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = "0.24";
    }
    return $self->{$method_key}; # returns a scalar
}

sub serialize_rdat_version {
    my $self = shift;
    my $line = "RDAT_VERSION ".$self->rdat_version()."\n";
    return $line;
}

sub name {
    my ($self, $name, $lines) = @_;
    my $method_key = "NAME";

    # name and lines means populating the object
    if ( defined $name && defined $lines ) {
        $self->{$method_key}->{$name} = 
            RNAprobing::RDATFile::Name->new($name, $lines);
            return $self->{$method_key}->{$name};
    # just name means either ...
    } elsif ( defined $name ){
        # give me the object with this name ...
        if ( exists $self->{$method_key}->{$name} ) {
            return $self->{$method_key}->{$name};
        # or create a new one with the given name
        } else {
            $self->{$method_key}->{$name} = 
                RNAprobing::RDATFile::Name->new($name);
            return $self->{$method_key}->{$name};
        }
    } 
    
    # this is just for creation of the key when creating a new RDAT object
    if ( !( defined $self->{$method_key}) ) {
        $self->{$method_key} = {};
       return $self->{$method_key}; # returns a hash reference
    }
    # if nothing else is asked for just return the hash
    return $self->{$method_key}; # returns a hash reference
}

sub serialize_name {
    my $self = shift;
    my $line = "NAME ".$self->name()."\n";
    return $line;
}


################################################################################
#
#                       RNAprobing::Annotation
#
################################################################################

#package RNAprobing::Annotation;


1;
