#!/usr/bin/env perl
#===============================================================================
#
#         FILE: rnaml2off.pl
#
#        USAGE: ./rnaml2off.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Christoph Kaempf (CK), kaempf@bioinf.uni-leipzig.de
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 05.08.2012 15:28:44
#     REVISION: ---
#===============================================================================

## Loading modules and initializing variables ##
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Find;
use File::Spec qw(rel2abs);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Path::Class;
my $module_dir = dirname(__FILE__);
push(@INC, $module_dir); 

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my $dumper = 0;
my $nota = 'LW';
my $help = '';
my $man = 0;
my $verbose = 0;
my @files = ();
my @directories = ();
GetOptions(
    "dumper" => \$dumper,
    "file|f=s" => \@files,
    "directory|d=s" => \@directories,
    "notation|n:s" => \$nota,
    "verbose|v+" => \$verbose,
    "help|h" => \$help,
    "man|m" => \$man) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if ( $help );
pod2usage(-verbose => 2) && exit if ( $man );
pod2usage(-verbose => 1) && exit if ( scalar(@files) == 0 && scalar(@directories) == 0 );

# allow coma seperated values as input
@files = split(/,/,join(',',@files));
@directories = split(/,/,join(',',@directories));

###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
my $log4perl_conf = file(dirname(__FILE__), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
$logger->info("++++ ".__FILE__." has been started. ++++");

## Lookup @files and @directories for rna1.xml files and insert the found in @rnaml_files
##  - find all rna1.xml files in the @directories given and add them to @files
if ( scalar(@directories ) != 0 ){
    my @checked_directories;
    foreach (@directories) {
        if ( -d $_ ) {
            $logger->info("$_ is a directory");
            push( @checked_directories, $_ );
        } else {
            $logger->info("$_ isn't a directory");
        }
    }
    $logger->error("No valid directory given.") && exit 0 
	if ( scalar(@checked_directories) == 0 && scalar(@files) == 0 );
    $logger->info("Looking for rna1.xml files in directories:\n".
		  join("\n", @checked_directories));
    find(\&wanted, @checked_directories);
}

##  - check found files and add them to @rnaml_files if they passed the checks
my $rnaml_files = &checkFiles(@files);

## Check for correct notation selection ##
my $notation = "";
if ($nota eq "LW") {
    $logger->info("Leontis-Westhof notation selected.");
    $notation = "LW";
} elsif ($nota eq "BP" ) {
    $logger->info("Base pair notation selected.");
    $notation = "BP";
} elsif ($nota eq "LW+" ) {
    $logger->info("Leontis-Westhof notation with nucleotide/amino acid interactions selected.");
    $notation = "LW+";
} elsif ($nota eq "LG" ) {
    $logger->info("Lee-Guttel notation selected.");
    $notation = "LG";
} else {
    $logger->info("Notation set to Leontis-Westhof.");
    $notation = "LW";
}

###############################################################
##
##              Application section
## 
###############################################################

# create XML::Simple object
my $xml = new XML::Simple;

# parse XML files

foreach my $rnaml_file ( @{$rnaml_files} ){
    ## $data - is the XML::Simple object that holds the representation of the RNAML file
    my $data = $xml->XMLin($rnaml_file, ForceArray => 1);


    my($filename, $directories, $suffix) = fileparse($rnaml_file);
    ## $off_files_list name of the file that holds a list of all .off files in a certain directory
    my $off_files_list = $directories ."off_files.list";
    ## $fasta_files_list name of the file that holds a list of all .fasta files in a certain directory
    my $fasta_files_list = $directories ."fasta_files.list";

    foreach my $key ( keys %{$data->{'molecule'}}){
        ## $off_file = file name for the output file (off = our file format)
        $filename =~ s/-rna1\.xml$//g;
        my $ndb_id = $filename;
        $logger->debug("NDB ID of $rnaml_file is $ndb_id.");
        $logger->debug("Key of actual molecule is $key.");

        my $fasta_file = $directories . $ndb_id .".". $key .".fa";
        my $off_file = $directories . $ndb_id .".". $key .".off";

        ## try to open the file $off_file
        open(my $off_fh, ">", $off_file) or die("Can't open $off_file.");
        $logger->debug("Create output file $off_file.");
        
        ## append $off_file to $off_files_list
        open( my $off_list_fh ,">>", $off_files_list) or die("Can't open $off_files_list.");
        print $off_list_fh $off_file ."\n";
        close($off_list_fh);

        ## try to open the file $fasta_file
        open(my $fasta_out_fh, ">", $fasta_file) or die("Can't open $fasta_file.");
        $logger->debug("Create output file $fasta_file.");

        ## append $fasta_file to $fasta_files_list
        open( my $fasta_list_fh ,">>", $fasta_files_list) or die("Can't open $fasta_files_list.");
        print $fasta_list_fh $fasta_file ."\n";
        close($fasta_list_fh);

        ## $ndb_id - holds the id of the molecule
        print $off_fh  "> $ndb_id.$key\n" ;
        print $fasta_out_fh "> $ndb_id.$key\n";
        $logger->debug("Header: > $ndb_id");

        ## $sequence - holds the RNA sequence
        my $sequence = join("", @{$data->{'molecule'}->{$key}->{'sequence'}[0]->{'seq-data'}});
        $sequence = &trimString($sequence);
        print $off_fh "$sequence\n";
        print $fasta_out_fh "$sequence\n";
        $logger->debug("Sequence: $sequence");
        
        ## @dotBracket - holds all Watson-Crick base pairs in Dot-Bracket form
        my @dotBracket;

        # fill @dotBracket with length($sequence) dots
        @dotBracket = split('', '.' x length($sequence) );
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
            my @old_bpLW = undef;
            foreach my $key (@bpStart5p){
                my @bpLW = @{ $bp{$key} };
                @old_bpLW = @bpLW unless ( @old_bpLW );
                ## insert opening and closing bracket for every standard 
                ## Watson-Crick base pair
                @dotBracket = &insertBrackets(\@bpLW, \@old_bpLW, \@dotBracket, \@bracket_stack, $filename, $sequence );
                @colSize = &colSizeUpdate($bpLW[0], $bpLW[1], $bpLW[2], \@colSize);
                @old_bpLW = @bpLW;
                
            }
            ## $dbString - is the output string of the Dot-Bracket notation
            my $dbLine = join("", @dotBracket);
            print $off_fh "$dbLine\n";
            $logger->debug($dbLine);
            
            ## $colSizeLine - assemble the line that holds the column sizes
            if ( scalar(@colSize) > 0 ) {
                my $colSizeLine = "# ". length($notation) .";".join(";", @colSize);
                print $off_fh $colSizeLine."\n" ;
                $logger->debug($colSizeLine);
            }

            # sort the keys with a specific subroutine
            foreach my $key (@bpStart5p){
                my $line = "# ".$notation;
                $line .= &makeColString($bp{$key}, \@colSize);
                print $off_fh "$line\n";
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
        close($off_fh);
        close($fasta_out_fh);
    }
}


###############################################################
##
##              Subroutine section
##
###############################################################

###############################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################

sub checkFiles {
    my @testfiles = @_;
    my @checkedfiles = ();
    my $logger = get_logger("RNAprobing");

    foreach (@testfiles) {
        $logger->info("Perform file checks for file $_");
        if ( -f $_){
            $logger->info("$_ is a plain file");
            if ( -r $_) {
                my $abs_path = File::Spec->rel2abs( $_ );
                $logger->info("$_ can be accessed");
                push(@checkedfiles, $abs_path );
                $logger->debug("Absolute path: $abs_path");
            } else {
                $logger->error("Can not access $_");
            }
        } else {
            $logger->error("$_ is not a plain file");
            $logger->error("$_ is a directory") if ( -d $_);
        }
    }
    return \@checkedfiles;
}

###############################################################################
##              
##  Invocation - \&wanted
##      - Subroutine used by File::Find
##      - Superugly this usage
##
###############################################################################

sub wanted {
    my $logger = get_logger("RNAprobing");
    if ( $_ =~ /rna1\.xml$/ ) {
        $logger->info("$_ is a rna1.xml file");
        # @files is a global variable
        push(@files, $File::Find::name);
    } else {
        $logger->info("$_ is not a rna1.xml file");
    }
}

###############################################################
#
# Invocation - &configureLogger($verbosityLevel)
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
    my $logger = get_logger("RNAprobing");
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
    my ($pos5P, $pos3P, $lw_not, $edge5P, $edge3P) = @{$_[0]};
#    my $pos3P = $_[1];
#    my $edge5P = $_[2];
#    my $edge3P = $_[3];
    my ($old_pos5P, $old_pos3P, $old_lw_not, $old_edge5P, $old_edge3P) = @{$_[1]};
    my @dotBracket = @{$_[2]};
    my $dbn = join("", @dotBracket);
    my $bracket_stack = $_[3];
    my $filename = $_[4];
    my @sequence = split("", $_[5]);

    my $bracket_type = 0;
    my %opening_brackets = ( 0 => '(', 1 => '[', 2 => '{', 3 => '<' );
    my %closing_brackets = ( 0 => ')', 1 => ']', 2 => '}', 3 => '>' );
    
    $logger->debug("$pos5P $pos3P $edge5P $edge3P");
    $logger->debug("$dotBracket[$pos5P-1] $dotBracket[$pos3P-1]");
    if ($edge5P eq 'W' && $edge3P eq 'W' ||
           $edge5P eq 'E' && $edge3P eq 'E' && $sequence[$pos5P-1] eq 'U' && $sequence[$pos3P-1] eq 'G' ||
           $edge5P eq 'E' && $edge3P eq 'E' && $sequence[$pos5P-1] eq 'G' && $sequence[$pos3P-1] eq 'U') {
        for (my $i = 0; $i < scalar(@{$bracket_stack}); $i++ ) {
            $logger->debug("Bracket type is $i");
            $logger->debug("Bracket stack contains ".scalar(@{$bracket_stack})." element(s).");
            if ( $pos5P < ${$bracket_stack}[$i] && $pos3P < ${$bracket_stack}[$i] || $pos5P > ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i] ) {
                ${$bracket_stack}[$i] = $pos3P;
                $bracket_type = $i;
                $logger->debug(${$bracket_stack}[$i]);
                last;
            } elsif ( $pos5P < ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i] && $i + 1 == scalar(@{$bracket_stack}) ) {
                push(@{$bracket_stack}, $pos3P);
                $logger->info("Detected pseudoknot.");
                $bracket_type = $i + 1;
                last;
            } elsif ($pos5P < ${$bracket_stack}[$i] && $pos3P > ${$bracket_stack}[$i]) {
                $logger->info("Detected pseudoknot.");
                $bracket_type = $i + 1;
            } else {
                $logger->error($filename.": Encountered a problem with the ".
                               "pseudoknot detection!");
                $logger->error("$filename : edge $sequence[$pos5P-1]".
                               "$pos5P-$sequence[$pos3P-1]$pos3P probably ".
                               "collides with edge $sequence[$old_pos5P-1]".
                               "$old_pos5P-$sequence[$old_pos3P-1]$old_pos3P");
            }
        }
        if ( $bracket_type > 3 ){
            $logger->error($filename.": More than four entangled pseudoknots! HELP!");
        } else {
            $dotBracket[$pos5P-1] = $opening_brackets{$bracket_type};
            $dotBracket[$pos3P-1] = $closing_brackets{$bracket_type};
        }
    }
    $dbn = join("", @dotBracket);
    $logger->debug($dbn);
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
    my $edge3p = uc(&ndb2lw($_->{'edge-3p'}[0]));
    my $edge5p = uc(&ndb2lw($_->{'edge-5p'}[0]));
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
        if ( $colSize->[$i] < length($bpLW[$i]) + 2 ){
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
    my $logger = get_logger("RNAprobing");
    my $string;
    for (my $i = 0; $i < scalar(@{$colSize}); $i++) { 
        my $diff = $colSize->[$i] - length($bpLW->[$i]);
        $string .= " " x $diff;
        $string .= $bpLW->[$i];
    }
    $logger->debug("This is edge: $string");
    return $string;
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
        # E = Watson-Crick Edge involved in non-standard base pair
        if ($_[0] eq 'W') { return "E"; last SELECT; }
        # N = Non-identified Edge
        if ($_[0] eq '.' || $_[0] eq '?' || $_[0] eq 'X' ){ return "N"; last SELECT; }
        # M = Triplets and higher multiplets
#        if ($_[0] eq '+' ){ return "M"; last SELECT; }
        # T = Tertiary interactions
        if ($_[0] eq '!' ){ return "T"; last SELECT; }
        else { return $_[0] }
    }
};


__END__

=head1 SYNOPSIS

 rnaml2off.pl converts a bunch of rna1.xml files into .off and .fa files.

=head1 DESCRIPTION

 Converts rna1.xml files into .off files. The .off files look like normal dot-bracket notation files, but also contain information about Leontis-Westhof style interactions and tertiary interactions between nucleotides.

 Switches that don't define a value can be done in long or short form.
 eg:
   rnaml2off.pl --man
   rnaml2off.pl -m
 This perl script relies on the Log::Log4perl, Path::Class and XML::Simple modules. Please make sure the according packages are around on your system.

=head1 OPTIONS

=over 8

=item B<-h, --help>

Display help message

=item B< -m, --man>

Display man page

=item B<--dumper>

Switch for debugging purpose. If set all found base pair informations are logged.

=item B<-f, --file>

RNAML file that is transformed into OFF file. This flag can be given multiple times.

=item B<-d, --directory>

Directory recursively searched for RNAML files ending on rna1.xml. 

=item B<-n, --notation>



=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (highest level if used 3 or more times) 

=cut
