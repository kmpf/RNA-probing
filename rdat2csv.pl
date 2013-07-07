#!/usr/bin/env perl
#===============================================================================
#
#         FILE: rdat2csv.pl
#
#        USAGE: ./rdat2csv.pl  
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

foreach my $rdat_file ( @{ $rdat_files } ) {
    # generate output file name
    my($filename, $directories, $suffix) = fileparse($rdat_file);
    my $csv_file = $directories.$filename;
    $csv_file =~ s/\.rdat$/.csv/g;
    # read in RDAT file and process information
    my $rdat_object = RNAprobing::RDATFile->new($rdat_file);
    my @sequence = split(//, $rdat_object->sequence());
    my $startpos = $rdat_object->seq_startpos();
    $logger->debug("Startpos: ".$startpos." nt: ".$rdat_object->offset_sequence_map()->{$startpos});
    my $endpos = $rdat_object->seq_endpos();
    $logger->debug("Endpos: ".$endpos." nt: ".$rdat_object->offset_sequence_map()->{$endpos});
    my $header = "Position\tNucleotide";
    foreach my $index ( @{$rdat_object->data()->indices()} ) {
	$header .= "\texpReac$index\tscaReac$index";
    }
    $header .= "\n";
    my $lines = "";
    for (my $i = $startpos; $i <= $endpos; $i++ ) {

    	$lines .= $i."\t".$rdat_object->offset_sequence_map()->{$i};
	    foreach my $index ( @{$rdat_object->data()->indices()} ) {
	        if (defined $rdat_object->seqpos_reactivity_map($index)->{$i}){
	         	$lines .= "\t".$rdat_object->seqpos_reactivity_map($index)->{$i};
	        } else {
    	    	$lines .= "\t"
	        }
	        if ( defined $rdat_object->seqpos_scaled_reactivity_map($index)->{$i}){
	    	$lines .= "\t"
	    	    .$rdat_object->seqpos_scaled_reactivity_map($index)->{$i};
	        } else {
	    	$lines .= "\t";
	        }
	    }
	    $lines .= "\n";
    }
    $logger->info($lines);
    open(my $csv_fh, ">", $csv_file) or die("Can't open $csv_file.");
    print $csv_fh $header;
    print $csv_fh $lines;
    close($csv_fh);
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

