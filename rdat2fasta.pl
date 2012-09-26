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
use File::Basename;
use File::Find;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;
use Pod::Usage;
my $module_dir = dirname(__FILE__);
$module_dir =~ s/scripts$/RNAprobing/g;
push(@INC, $module_dir);
require RNAprobing::RDATFile;

## Configure Getopt::Long ##
Getopt::Long::Configure ("bundling");
my @files = ();
my @directories = ();
my $help = 0;
my $verbose = 0;
my $to_dna = 0;
GetOptions(
	"file|f=s" => \@files,
    "directory|d=s" => \@directories,
    "toDNA|t" => \$to_dna,
	"verbose|v+" => \$verbose);

if ( $help || scalar(@files) == 0 && scalar(@directories) == 0 ){
    pod2usage( { -verbose => 1,
                 -message => "Use this script like this:\n"});
}
###############################################################################
#                 
# Logger initiation  
#                 
###############################################################################
my $this_file = __FILE__;
$this_file =~ s/scripts/RNAprobing/g;
my $log4perl_conf = file(dirname($this_file), "RNAprobing.log.conf");

# Apply configuration to the logger
Log::Log4perl->init("$log4perl_conf");

# Get the logger
my $logger_name = "RNAprobing";
my $logger = &configureLogger($verbose, $logger_name);
print("++++ ".__FILE__." has been started. ++++\n");
$logger->info("++++ ".__FILE__." has been started. ++++");


## Lookup @files and @directories for rna1.xml files and insert the found in @rnamlFiles
##  - find all rna1.xml files in the @directories given and add them to @files
find( \&wanted, @directories);
$logger->info("\@files after find(): ". join(";", @files));

##  - check found files and add them to @rdat_files if they passed the checks
my $rdat_files = &checkFiles(\@files);

my %replace = (U => "T", u => "t");
my $regex = join( "|", keys %replace);
$regex = qr/$regex/;

foreach my $rdat_file ( @{ $rdat_files } ) {
    my($filename, $directories, $suffix) = fileparse($rdat_file);
    my $fasta_file = $directories.$filename;
    $fasta_file =~ s/\.rdat$//g;
    my $rdat_object = RNAprobing::RDATFile->new($rdat_file);
    my $sequence = $rdat_object->sequence();
    if ($to_dna) {
        $sequence =~ s/($regex)/$replace{$1}/g;
        $fasta_file .= ".dna.fa";
    } else {
        $fasta_file .= ".fa";
    }
    open(my $fasta_fh, ">", $fasta_file) or die("Can't open $fasta_file.");
    print $fasta_fh ">".$rdat_object->name()."\n";
    print $fasta_fh $sequence."\n";
    close($fasta_fh);
}

$logger->info(join(";", @{ $rdat_files}));

###############################################################
##              
##              Subroutine section
##
###############################################################

###############################################################################
##              
##  Invocation - \&wanted
##      - Subroutine used by File::Find
##
###############################################################################

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
## -- 0 => $WARN
## -- 1 => $INFO
## -- 2 => $DEBUG
## 
###############################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    $logger->info("Verbosity level: $verbose");
    SELECT:{
	    if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
	    else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
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

rdat2fasta.pl -f=</path/to/file> -v -v -v -t


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

Prints this help page.

=back


=back
=head1 DESCRIPTION
B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut

