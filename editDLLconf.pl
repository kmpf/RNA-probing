#!/usr/bin/env perl
#===============================================================================
#
#         FILE: editDLconf.pl
#
#        USAGE: ./editDLconf.pl  
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
use File::Basename;
use File::Copy;
use File::Find;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
use Path::Class;

###############################################################################
#
# Options section
#
################################################################################
my @directories_to_search = ();
my $delete = "";
my $insert = "";
my $replace = "";
my $replace_by = "";
my $find_pattern = '\.conf';
my $verbose = 2;
GetOptions(
    "directory|d=s" => \@directories_to_search,
    "insert|i=s" => \$insert,
    "replace|r=s" => \$replace,
    "replaceBy|rb=s" => \$replace_by,
    "delete|d=s" => \$delete,
    "findPattern|f=s" => \$find_pattern,
    "verbose|v+" => \$verbose);

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
$logger->info("++++ ".__FILE__." has been started. ++++");

###############################################################################
#                 
# Application section
#                 
###############################################################################

$logger->error("With what do you want \"$replace\" to be replaced? Use --replaceBy") 
    if ($replace ne "" && $replace_by eq "");
$logger->error("What do you want to be replaced by \"$replace_by\"? Use --replace") 
    if ($replace eq "" && $replace_by ne "");

my @date = localtime;
my $year = $date[5] + 1900;
my $month = sprintf( "%02d", $date[4] + 1);
my $day = sprintf("%02d", $date[3]);

if ( scalar(@directories_to_search) != 0 ) {
    find(\&wanted, @directories_to_search);
} else {
    $logger->error("No directory given to search in. Use --directory");
}


###############################################################################
##              
##              Subroutine section
##
###############################################################################

sub wanted {
    my $logger = get_logger();
    my $file = $_;
    my $backup_file = "$file.$date[5]$date[4]$date[3].old";
    if ( ! -B $file ) {
        if ( $file =~ m/$find_pattern/ ) {
            chomp($_);
            $logger->info("Found pattern \"$find_pattern\" in $file.");
            if (move( $file, $backup_file )) {
                $logger->info("Created a backup file $backup_file of $file");
            } else {
                $logger->info("Couldn't move $file to $backup_file: $!");
                next;
            }
            open(my $conf_file, ">", $file) or die "Couldn't write to $file";
            open(my $backup_conf_file, "<", $backup_file) or die "Couldn't read $backup_file";

            while (<$backup_conf_file>) {
                if ( ($delete ne "") && ($_ =~ m/$delete/) ) {
                    
                } elsif ( ( $_ =~ m/$replace/ ) && ($replace ne "") && ($replace_by ne "") ) {
                    print($conf_file "$replace_by\n");
                } else {
                    print($conf_file "$_\n");
                }
            }
            open(my $insert_file, "<", $insert) or die "Couldn't read from $file";
            my $insert_text = "";
            while (<$insert_file>) {
                $insert_text .= $_;
            }
            print($conf_file $insert_text) unless ( $insert_text eq "" );
            close($conf_file);
            close($backup_conf_file);
        } else {
            $logger->info("$file does not match pattern:\"$find_pattern\"");
        }
    } else {
        $logger->info("File $file is not a text file.");
    }
}

###############################################################################
##
## &configureLogger($verbosityLevel)
## - Configures and initialzes the Logger
## - $verbosityLevel = scalar value that sets log level
## -- 0 => $WARN
## -- 1 => $INFO
## -- 2 => $DEBUG
## 
###############################################################################

sub configureLogger{
    ## Configure the logger ##
    my $verbose = shift;
    my $logger_name = shift;
    my $logger = get_logger($logger_name);
    SELECT:{
	    if ($verbose == 0){ $logger->level($WARN) ; $logger->debug("Log level is WARN") ; last SELECT; }
	    if ($verbose == 1){ $logger->level($INFO) ; $logger->debug("Log level is INFO") ; last SELECT; }
	    if ($verbose == 2){ $logger->level($DEBUG); $logger->debug("Log level is DEBUG") ;  last SELECT; }
	    else {$logger->level($ERROR); $logger->debug("Log level is ERROR") ;  last SELECT; }
    }
    return $logger;
}

