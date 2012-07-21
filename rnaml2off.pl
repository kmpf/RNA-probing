#!/usr/bin/env perl

## Loading modules and initializing variables ##
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Path::Class;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my $dumper = 0;
my $filelist = '';
my $nota = 'LW';
my $help = '';
my $man = 0;
my $verbose = 0;
my @files = ();
GetOptions(
    "dumper" => \$dumper,
    "file|f=s" => \$filelist,
    "notation|n:s" => \$nota,
    "verbose|v+" => \$verbose,
    "help|h" => \$help,
    "man|m" => \$man) or pod2usage(-verbose => 1) && exit;

# Configuration file for the logger
my $log4perl_conf = file(dirname(__FILE__), "rnaml2off.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger = &configureLogger($verbose);

pod2usage(-verbose => 1) && exit if ( $help );
pod2usage(-verbose => 1) && exit if ( !$filelist );
pod2usage(-verbose => 2) && exit if ( $man );



if ($filelist){
    $filelist =~ s/^\s+//;
    $logger->debug($filelist);
    @files = split(/\s+/,$filelist);
    $logger->debug(Dumper(@files));
} else {
    exit;
}





## Check for correct notation selection ##
my $notation = "";
if ($nota eq "LW") {
    $logger->info("Leontis-Westhof notation selected.");
    $notation = "LW";
}
if ($nota eq "BP" ) {
    $logger->info("Base pair notation selected.");
    $notation = "BP";
}
if ($nota eq "LW+" ) {
    $logger->info("Leontis-Westhof notation with nucleotide/amino acid interactions selected.");
    $notation = "LW+";
}
if ($nota eq "LG" ) {
    $logger->info("Lee-Guttel notation selected.");
    $notation = "LG";
}
else {
    $notation = "LW";
}

# Check if files are readable
my @rnamlFiles = ();
foreach (@files){
    $logger->info("Perform file checks for file $_");
    if ( -f $_){
        $logger->info("$_ is a plain file");
        if ( -r $_) {
            $logger->info("$_ can be accessed");
            push(@rnamlFiles, $_);
        } else {
            $logger->error("Can not access $_");
        }

    } else {
        $logger->error("$_ is not a plain file");
        $logger->error("$_ is a directory") if ( -d $_);
    }
}

###############################################################
##
##              Application section
## 
###############################################################

# create XML::Simple object
my $xml = new XML::Simple;

# parse XML files

foreach my $rnaml_file ( @rnamlFiles ){
    ## $data - is the XML::Simple object that holds the representation of the RNAML file
    my $data = $xml->XMLin($rnaml_file, ForceArray => 1);


    foreach my $key ( keys %{$data->{'molecule'}}){
        ## $off_file = file name for the output file (off = our file format)
        $rnaml_file =~ s/-rna1.xml$//g;
        my $ndb_id = $rnaml_file;
        $logger->debug("NDB ID of $rnaml_file is $ndb_id.");
        my $fasta_file = "$ndb_id.$key.fa";
        my $off_file = "$ndb_id.$key.off";

        ## try to open the file $off_file
        open(OURFILE, ">$off_file") or die("Can't open $off_file.");
        $logger->debug("Create output file $off_file.");

        ## try to open the file $fasta_file
        open(FASTAFILE, ">$fasta_file") or die("Can't open $fasta_file.");
        $logger->debug("Create output file $fasta_file.");

        ## $ndb_id - holds the id of the molecule
        print OURFILE "> $ndb_id.$key\n";
        print FASTAFILE "> $ndb_id.$key\n";
        $logger->debug("Header: > $ndb_id");

        ## $sequence - holds the RNA sequence
        my $sequence = join("", @{$data->{'molecule'}->{$key}->{'sequence'}[0]->{'seq-data'}});
        $sequence = &trimString($sequence);
        print OURFILE "$sequence\n";
        print FASTAFILE "$sequence\n";
        $logger->debug("Sequence: $sequence");
        
        ## @dotBracket - holds all Watson-Crick base pairs in Dot-Bracket form
        my @dotBracket;

        # fill @dotBracket with length($sequence) dots
        for (my $i  = 0; $i < length($sequence); $i++){
            $dotBracket[$i] = '.' ;
        }
        $logger->debug(@dotBracket);

        ## $basePairs - is an array reference on an array which contains all information about all base pairs
        my $basePairs = $data->{'molecule'}->{$key}->{'structure'}[0]->{'model'}->{'?'}->{'str-annotation'}->[0]->{'base-pair'};

        if ($dumper) {
            $logger->level($DEBUG);
            $logger->debug(Dumper($basePairs));
        }

        ## %bp - is a hash holding every base-base contact as array under a unique hash-key
        my %bp = ();

        ## every 'base pair' is inserted into %bp
        ## and in case of standard Watson-Crick base pairs the @dotBracket array is extended
        if ($notation eq "LW"){
            my @colSize;
            my %fivePedge_bpLW = ();
            my @bracket_stack = (0);
            foreach (@{$basePairs}){
                # read in all base pairs and store them in %bp
                my @bpLW = &LWnotation($_);
                my $id = &makeId(\@bpLW);
                $bp{$id} = \@bpLW;
            }
            my @bpStart5p = sort keyCmp keys(%bp);
            foreach my $key (@bpStart5p){
                my @bpLW = @{ $bp{$key} };
                ## insert opening and closing bracket for every standard Watson-Crick base pair
                @dotBracket = &insertBrackets($bpLW[0], $bpLW[1], $bpLW[3], $bpLW[4], \@dotBracket, \@bracket_stack);
                @colSize = &colSizeUpdate($bpLW[0], $bpLW[1], $bpLW[2], \@colSize);
            }
            ## $dbString - is the output string of the Dot-Bracket notation
            my $dbLine = join("", @dotBracket);
            print OURFILE "$dbLine\n";
            $logger->debug($dbLine);
            
            ## $colSizeLine - assemble the line that holds the column sizes 
            my $colSizeLine = "# ".length($notation) . ";";
            $colSizeLine .= join(";", @colSize);
            print OURFILE "$colSizeLine\n";
            $logger->debug($colSizeLine);

            # sort the keys with a specific subroutine
            foreach my $key (@bpStart5p){
                my $line = "# ".$notation;
                $line .= &makeColString($bp{$key}, \@colSize);
                print OURFILE "$line\n";
                $logger->debug("$line");
            }
        }
        if ($notation eq "BP"){
            
        }
        if ($notation eq "LW+"){
            
        }
        if ($notation eq "LG"){
            
        }
        # clear %bp or the next
        %bp = ();
        close(OURFILE);
        close(FASTAFILE);
    }
}


###############################################################
##
##              Subroutine section
##
###############################################################

###############################################################
#
# &configureLogger($verbosityLevel)
# - Configures and initialzes the Logger
# - $verbosityLevel = scalar value that sets log level
# -- 0 => $WARN
# -- 1 => $INFO
# -- 2 => $DEBUG
# 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger = get_logger("rnaml2off");
    $logger->info("Verbosity level: $verbose");
    SELECT:{
        if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
        if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
        if ($verbose >= 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
        else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
    }
    return $logger;
}

###############################################################
#
# Invocation - &trimString($string)
# remove all white spaces and newlines from $string
#
###############################################################

sub trimString{
    my $string = $_[0];
    $string =~s/[\n,\s]//g;
    return $string;
};

###############################################################
#
# Invocation - &clearHash()
# guess what? Never used but nice-to-have
#
###############################################################

sub clearHash{
    my %hash = %{$_[0]};
    for (keys %hash){
        delete $hash{$_};
    }
    return %hash;
};

###############################################################
#
# Invocation - &insertBrackets($pos5p, $pos3p, $edge5p, $edge3p, \@dotBracket, \@bracket_stack)
# inserts opening and closing Bracket in Dot-Bracket string
#
###############################################################

sub insertBrackets{
    my $pos5P = $_[0];
    my $pos3P = $_[1];
    my $edge5P = $_[2];
    my $edge3P = $_[3];
    my @dotBracket = @{$_[4]};
    my $bracket_stack = $_[5];
    my $bracket_type = 0;
    my %opening_brackets = ( 0 => '(', 1 => '[', 2 => '{', 3 => '<' );
    my %closing_brackets = ( 0 => ')', 1 => ']', 2 => '}', 3 => '>' );
    $logger->debug("$pos5P $pos3P $edge5P $edge3P");
    $logger->debug("$dotBracket[$pos5P-1] $dotBracket[$pos3P-1]");
    if ($edge5P eq 'W' && $edge3P eq 'W'){
        for (my $i = 0; $i < scalar(@{$bracket_stack}); $i++ ) {
            $logger->debug("Bracket type is $i");
            if ( $pos5P < ${$bracket_stack}[$i] && $pos3P < ${$bracket_stack}[$i] || $pos5P > ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i] ) {
                ${$bracket_stack}[$i] = $pos3P;
                $bracket_type = $i;
                last;
            } elsif ( $pos5P < ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i] && $i + 1 == scalar(@{$bracket_stack}) ) {
                push(@{$bracket_stack}, $pos3P);
                $bracket_type = $i;
                last;
            } elsif ($pos5P < ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i]) {
                $logger->info("Detected pseudoknot crossing ".($i + 1)." edges.");
            } else {
               $logger->error("Encountered a problem with the pseudoknot detection!");
            }
        }
        if ( $bracket_type > 3 ){
            $logger->error("More than four entangled pseudoknots! HELP!");
        } else {
            $dotBracket[$pos5P-1] = $opening_brackets{$bracket_type};
            $dotBracket[$pos3P-1] = $closing_brackets{$bracket_type};
        }
    }
    my $s = join("", @dotBracket);
    $logger->debug($s);
    return @dotBracket;
};

###############################################################
#
# Invocation - internally done via sort()
# sort routine for the hash keys used within %bp
#
###############################################################

sub keyCmp{
    my @a = split(" ", $a);
    my @b = split(" ", $b);
    if ($a[0] < $b[0]){
        return -1;
    }
    if ($a[0] > $b[0]){
        return 1;
    }
    else {
        if ($a[1] < $b[1]){
            return -1;
        }
        if ($a[1] > $b[1]){
            return 1;
        }
        else {
            return $a[3] cmp $b[3];
        }
    }
};

###############################################################
#
# Invocation - &LWnotation($_)
#
#   Returns an array with five elements describing a base pair:
#       0: sequence position of the 5' base
#       1: sequence position of the 3' base
#       2: three letter string describing bond orientation,
#          5' edge and 3' edge
#       3: edge type of the 5' base
#       4: edge type of the 3' base
#
###############################################################

sub LWnotation{
    my $pos3p = $_->{'base-id-3p'}[0]->{'base-id'}[0]->{'position'}[0];
    my $pos5p = $_->{'base-id-5p'}[0]->{'base-id'}[0]->{'position'}[0];
    my $edge3p = &ndb2lw($_->{'edge-3p'}[0]);
    my $edge5p = &ndb2lw($_->{'edge-5p'}[0]);
    my @bpLW = ();
    $bpLW[0] = $pos5p;
    $bpLW[1] = $pos3p;
    $bpLW[2] = $_->{'bond-orientation'}[0].$edge5p.$edge3p;
    $bpLW[3] = $edge5p;
    $bpLW[4] = $edge3p;
    return @bpLW;
}

###############################################################
#
# Invocation - &makeId($_)
# - ID is something like "1 10 cWW"
###############################################################

sub makeId{
    my @bpLW = @{$_[0]};
    my $id;
    foreach my $i (@bpLW){
        $id .= $i." ";
    }
    chop($id);
    $logger->debug("ID: $id");
    return $id;
}

###############################################################
#
# Invocation - &colSizeUpdate(\$bpLW[0], $bpLW[1], $bpLW[2], \@colSize)
#
###############################################################

sub colSizeUpdate{
    my $pos5p = shift;
    my $pos3p = shift;
    my $edge_notation = shift;
    my @bpLW = ($pos5p, $pos3p, $edge_notation);
    my $colSize = shift ;
    for (my $i = 0; $i < scalar(@bpLW); $i++){
        if ( ! defined $colSize->[$i] ){
            $colSize->[$i] = 0;
        }
        if ( length($bpLW[$i]) + 2 > $colSize->[$i] ){
            # the addition of 2 results in a more readable output
            $colSize->[$i] = length($bpLW[$i]) + 2;
        }
    }
    return @{$colSize};
}

###############################################################
#
# Invocation - &makeColString(\@{bp{$key}}, \@colSize, $i)
#
###############################################################

sub makeColString{
    my $bpLW = shift;
    my $colSize = shift;
    my $logger = get_logger();
    my $string;
    for (my $i = 0; $i < scalar(@{$colSize}); $i++) { 
        my $diff = $colSize->[$i] - length($bpLW->[$i]);
        $string = " " x $diff;
        $string .= $bpLW->[$i];
    }
    return $string

}

###############################################################
#
#               Transformation from RNAML to $notation
#               --------------------------------------
#
# Invocation - &ndb2lw($char)
# $char should contain exactly one character
# transform RNAML non-standard into standard Leontis-Westhof Notation
#
###############################################################

sub ndb2lw{
    SELECT:{
        # W = Watson-Crick Edge
        if ($_[0] eq '+' || $_[0] eq '-' ){ return "W"; last SELECT; }
        # N = Non-identified Edge
        if ($_[0] eq '.' || $_[0] eq '?' || $_[0] eq 'X' ){ return "N"; last SELECT; }
        # M = Triplets and higher multiplets
        if ($_[0] eq '+' ){ return "M"; last SELECT; }
        # T = Tertiary interactions
        if ($_[0] eq '!' ){ return "T"; last SELECT; }
        else { return $_[0] }
    }
};


__END__

=head1 SYNOPSIS

 rnaml2off.pl converts a bunch of .rnaml files into .off files.

=head1 DESCRIPTION

 Converts .rnaml files into .off files. The .off files look like normal dot-bracket notation files, but also contain information about Leontis-Westhof style interactions and tertiary interactions between nucleotides.

 Switches that don't define a value can be done in long or short form.
 eg:
   rnaml2off.pl --man
   rnaml2off.pl -m
 This perl script relies on the Log::Log4perl, Path::Class and XML::Simple modules. Please make sure the according packages are around on your system.

=head1 ARGUMENTS





=cut

