#!/usr/bin/env perl
#===============================================================================
#
#         FILE: rdat2fasta.pl
#
#        USAGE: ./rdat2fasta.pl  
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
use utf8;
use Data::Dumper;
use File::Basename;
use File::Find;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;
my $module_dir = dirname(__FILE__);
push(@INC, $module_dir);

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my @directories = ();
my $help = 0;
my $man = 0;
my $verbose = 0;
my $to_dna = 0;
GetOptions(
    "file|f=s" => \@files,
    "directory|d=s" => \@directories,
    "toDNA|t" => \$to_dna,
    "help|h" => \$help,
    "man|m" => \$man,
    "verbose|v+" => \$verbose);

pod2usage(-verbose => 1) && exit if ( $help );
pod2usage(-verbose => 2) && exit if ( $man );
pod2usage({-verbose => 1, -message => "Use this script like this:\n"}) &&
    exit if ( $help || scalar(@files) == 0 && scalar(@directories) == 0 );
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

###############################################################################
#
# Load own perl modules
#
###############################################################################

require RNAprobing::RDATFile;


## Lookup @files and @directories for rna1.xml files and insert the found in @rnamlFiles
##  - find all rna1.xml files in the @directories given and add them to @files
if ( scalar(@directories ) != 0 ){
    my @checked_directories = ();
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
    $logger->info("Looking for rdat files in directories:\n".
		  join("\n", @checked_directories));
   find(\&wanted, @checked_directories);
}

##  - check found files and add them to @rdat_files if they passed the checks
my $rdat_files = &checkFiles(\@files);

my %replace = (U => "T", u => "t");
my $regex = join( "|", keys %replace);
$regex = qr/$regex/;

foreach my $rdat_file ( @{ $rdat_files } ) {
    my($filename, $directories, $suffix) = fileparse($rdat_file);
    my $fasta_file = $directories.$filename;
    $fasta_file =~ s/\.rdat$//g;
    my $fasta_file_short = $fasta_file.".short";
    my $rdat_object = ();
    $rdat_object = RNAprobing::RDATFile->new($rdat_file);
    $rdat_object->write_file($rdat_file.".test");
    # I want to sequences as output
    # 1. $sequence: the complete one for the BLAT mapping
    # 2. $short_seq: the short one only giving the subsequence for
    #    which probing data is available
    my $sequence = $rdat_object->sequence();
    my $short_sequence = "";
    my $seqpos = $rdat_object->seqpos();

    foreach my $pos (@{$seqpos}) {
	if (defined $rdat_object->offset_sequence_map()->{$pos}) {
	    $short_sequence .= $rdat_object->offset_sequence_map()->{$pos};
	}
    }

    if ($to_dna) {
        $sequence =~ s/($regex)/$replace{$1}/g;
        $short_sequence =~ s/($regex)/$replace{$1}/g;
        $fasta_file .= ".dna.fa";
	$fasta_file_short .= ".dna.fa";
    } else {
        $fasta_file .= ".fa";
	$fasta_file_short .= ".fa";
    }

    my $fasta_id = $filename;
    $fasta_id =~ s/\_.*$//g;
    open(my $fasta_fh, ">", $fasta_file) or die("Can't open $fasta_file.");
    print $fasta_fh "> ".$fasta_id."\n";
    print $fasta_fh $sequence."\n";
    close($fasta_fh);

    open(my $fasta_short_fh, ">", $fasta_file_short) or
	die("Can't open $fasta_file_short.");
    print $fasta_short_fh "> ".$fasta_id."\n";
    print $fasta_short_fh $short_sequence."\n";
    close($fasta_short_fh);    
}

my $j = 0;
foreach my $i ( @{$rdat_files}) {
    $j++;
    $logger->info($j.". processed RDAT file: ".$i);
}
###############################################################
##              
##              Subroutine section
##
###############################################################

###############################################################
##              
##  Invocation - \&wanted
##      - Subroutine used by File::Find
##
###############################################################

sub wanted {
    my $logger = get_logger("RNAprobing");
    $logger->info("$_ is a directory") if ( -d $_ );
    $logger->info("$_ isn't a directory") unless ( -d $_ );
    if ( $_ =~ /\.rdat$/ ) {
        $logger->info($File::Find::name." is a .rdat file");
        # @files is a global variable
        push(@files, $File::Find::name);
    } else {
        $logger->info("$_ is not a .rdat file");
    }
}

###############################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $ERROR
## -- 1 => $WARN
## -- 2 => $INFO
## -- >2 => $DEBUG
## 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    $logger->info("Verbosity level: $verbose");
    SELECT:{
	    if ($verbose == 0){$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
	    if ($verbose == 1){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    else { $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
    }
    return $logger;
}

###############################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################

sub checkFiles {
    my $testfiles = shift;
    my @checkedfiles = ();
    my $logger = get_logger("RNAprobing");
    foreach my $testfile (@{ $testfiles}) {
        $logger->info("Perform file checks for file $testfile");
        if ( -f $testfile ){
            $logger->info("$testfile is a plain file");
            if ( -r $testfile) {
                my $abs_path = File::Spec->rel2abs( $testfile );
                $logger->info("$testfile can be accessed");
                push(@checkedfiles, $abs_path );
                $logger->debug("Absolute path: $abs_path");
            } else {
                $logger->error("Can not access $testfile");
            }
        } else {
            $logger->error("$_ is not a plain file");
            $logger->error("$_ is a directory") if ( -d $testfile);
        }
    }
    return \@checkedfiles;
}

__END__


=head1 NAME

rdat2fasta.pl - Creating .fasta files fom .rdat files

=head1 SYNOPSIS

rdat2fasta.pl -f=</path/to/file> -d=</path/to/rdat-directory/> -v -t

=head1 DESCRIPTION

B<rdat2fasta.pl> takes RDAT files as input or searches given directories for them and converts them to FASTA files. The FASTA files are created within the same directory as the original RDAT files.

=head1 OPTIONS

=over 4

=item -f, --file=</path/to/file.rdat>

RDAT file(s) to be converted to FASTA files. Option can be given multiple times.

=item -d, --directory=</path/to/rdat-directory/>

The given directories will be searched for "*.rdat" files which will be converted. Option can be given multiple times.

=item -t, --toDNA

If set RNA sequences are converted into DNA sequences. The resulting file ends on ".dna.fa".

=item -v, --verbose

Increases verbosity level. Option can be given multiple times.

=item -h, --help

Prints the help page.

=item -m, --man

Prints the manual page.

=back

=cut

