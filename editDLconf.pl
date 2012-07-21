#!/usr/bin/env perl

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
my $find_pattern = "";
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

my $log4perl_conf = file(dirname(__FILE__), "RNAprobing.log.conf");

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
            if (open(CONF_FILE, ">", $file)) {
                $logger->info("Opened $file to write into");
            } else {
                $logger->info("Couldn't write to $file");
                next;
            }
            if (open(BACKUP_FILE, "<", $backup_file)) {
                $logger->info("Opened $backup_file to read from");
            } else {
                $logger->info("Couldn't read from $backup_file");
                next;
            }
            while (<BACKUP_FILE>) {
                if ( ($delete ne "") && ($_ =~ m/$delete/) ) {
                    
                } elsif ( ( $_ =~ m/$replace/ ) && ($replace ne "") && ($replace_by ne "") ) {
                    print(CONF_FILE "$replace_by\n");
                } else {
                    print(CONF_FILE "$_\n");
                }
            }
            print(CONF_FILE "$insert\n") unless ( $insert eq "" );
            close(CONF_FILE);
            close(BACKUP_FILE);
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

