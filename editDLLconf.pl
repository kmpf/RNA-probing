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
use Pod::Usage;

###############################################################################
#
# Options section
#
################################################################################
my @confs_to_combine = ();
my @directories_to_search = ();
my $delete = "";
my $help = 0;
my $insert = "";
my $man = 0;
my $model = "";
my $replace = "";
my $replace_by = "";
my $template = "";
my $find_pattern = '\.conf';
my $verbose = 0;
GetOptions(
    "combine|c=s" => \@confs_to_combine,
    "delete|del=s" => \$delete,
    "directory|dir=s" => \@directories_to_search,
    "findPattern|f=s" => \$find_pattern,
    "help|h" => \$help,
    "insert|i=s" => \$insert,
    "man|m" => \$man,
    "model=s" => \$model,
    "replace|r=s" => \$replace,
    "replaceBy|rb=s" => \$replace_by,
    "template|t=s" => \$template,
    "verbose|v+" => \$verbose);

if ( $help ) {
    pod2usage( -verbose => 1 ) && exit;
} elsif ( $man ) {
    pod2usage( -verbose => 2 ) && exit;
} elsif ( ($replace eq "" || $replace_by eq "") && 
    (scalar(@confs_to_combine) == 0 || $template eq "" || $model eq "") ){
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
$logger->info("++++ ".__FILE__." has been started. ++++");

###############################################################################
#                 
# Application section
#                 
###############################################################################

if ( $replace ne "" && $replace_by ne "") {  
    if ( scalar(@directories_to_search) != 0 ) {
        find(\&wanted, @directories_to_search);
    } 
}


if ( scalar(@confs_to_combine) != 0 && $template ne "" && $model ne "" ) {

# Check files for accessibility
    my @conf_files = ();
    foreach my $file (@confs_to_combine) {
        push(@conf_files, &checkFiles($file) );
    }
    $template = &checkFiles($template);
# If we do only have one file there is nothing to do
    if (scalar(@conf_files) == 0){
        $logger->error("Only one file given to -c , --combine. Nothing to do.");
        exit 0;
    }

# Extract the relevant info from each file and ...
    my @lp_pos = ();
    my @lp_neg = ();
    for ( my $i = 1; $i <= scalar(@conf_files); $i++ ) {
        open(my $conf_file, "<", $conf_files[$i-1]) or 
            die "Couldn't read $conf_files[$i-1].";
        while(<$conf_file>){
            chomp($_);
            if ( $_ =~ m/^lp\.positiveExamples/ ){
                $_ =~ s/^lp\.positiveExamples.*\{//g;
                $_ =~ s/}.*//g;
                my @pos_from_conf = split(',', $_);
                push(@lp_pos, @pos_from_conf);
            }
            if ( $_ =~ m/^lp\.negativeExamples/ ){
                $_ =~ s/^lp\.negativeExamples.*\{//g;
                $_ =~ s/}.*//g;
                my @neg_from_conf = split(',', $_);
                push(@lp_neg, @neg_from_conf);
            }
        }
        close($conf_file);
    }

# ... combine them according to the template
    my $combined = $template;
    $combined =~ s/conf$/combined.conf/;
    open(my $template_file, "<", $template) or die "Couldn't read $template.";
    open(my $combined_file, ">", $combined) or die "Couldn't write $combined.";
    while(<$template_file>){
        my $line = "";
        if ( $_ =~ m/^ks\d*\.fileName/ ){
            $_ =~ m/(.*\")(.*)(\".*)/;
            print($combined_file $1."./".$model.$3."\n");
        } elsif ( $_ =~ m{^// Nr\. of positives: \d*}){
            print($combined_file "// Nr. of positives: ".scalar(@lp_pos)."\n");
        } elsif ($_ =~ m/^lp\.positiveExamples/){
            $line = join(",", @lp_pos);
            print($combined_file "lp.positiveExamples = {".$line."}\n");
        } elsif ( $_ =~ m{^// Nr\. of negatives: \d*}){
            print($combined_file "// Nr. of negatives: ".scalar(@lp_neg)."\n");
        } elsif ($_ =~ m/^lp\.negativeExamples/){
            $line = join(",", @lp_neg);
            print($combined_file "lp.negativeExamples = {". $line ."}\n");
        } else {
            print($combined_file $_);
        }
    }
    close($template_file);
    close($combined_file);
}

###############################################################################
##              
##              Subroutine section
##
###############################################################################

sub wanted {
    my $logger = get_logger("RNAprobing");
    my $file = &checkFiles($_);
    # assemble filename of file to be created
    my @date = localtime;
    my $year = $date[5] + 1900;
    my $month = sprintf( "%02d", $date[4] + 1);
    my $day = sprintf("%02d", $date[3]);
    my $modified_file = "$file.mod.$date[5]$date[4]$date[3].conf";
    if ( ! -B $file ) {
        # see if the file is one we should modify
        if ( $file =~ m/$find_pattern/ ) {
            chomp($_);
            $logger->info("Found pattern \"$find_pattern\" in $file.");
            open(my $conf_file, ">", $modified_file) or die "Couldn't write to $modified_file.";
            open(my $orig_conf_file, "<", $file) or die "Couldn't read $file.";

            # Read in original file and check if we need to delete or replace
            # the current line
            while (<$orig_conf_file>) {
                if ( ($delete ne "") && ($_ =~ m/$delete/) ) {
                    
                } elsif ( ( $_ =~ m/$replace/ ) && ($replace ne "") && ($replace_by ne "") ) {
                    print($conf_file "$replace_by\n");
                } else {
                    print($conf_file "$_\n");
                }
            }
            # Insert specified info at end of file
            open(my $insert_file, "<", $insert) or die "Couldn't read from $insert";
            my $insert_text = "";
            while (<$insert_file>) {
                $insert_text .= $_;
            }
            close($insert_file);
            print($conf_file $insert_text) unless ( $insert_text eq "" );
            close($conf_file);
            close($orig_conf_file);
        } else {
            $logger->info("$file does not match pattern:\"$find_pattern\"");
        }
    } else {
        $logger->info("File $file is not a text file.");
    }
}

###############################################################################
##
## &checkFiles(@filesToBeChecked)
## - Performs file checks and returns an array with all succesfully checked files
## - @filesToBeChecked = array of files to be checked
##
###############################################################################

sub checkFiles {
    my $test_file = shift;
    my $checked_file = "";
    my $logger = get_logger("RNAprobing");
    # Check if files are readable
    if ( -r $test_file){
        $checked_file = $test_file;
        $logger->info("$test_file is readable.");
    } else {
        $logger->error("$test_file is not readable.");
        exit;
    }
    return $checked_file;
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
	    if ($verbose == 0){ 
            $logger->level($WARN) ; 
            $logger->debug("Log level is WARN") ;
            last SELECT;
        }
	    if ($verbose == 1){ 
            $logger->level($INFO) ; 
            $logger->debug("Log level is INFO") ;
            last SELECT;
        }
	    if ($verbose == 2){ 
            $logger->level($DEBUG); 
            $logger->debug("Log level is DEBUG") ;
            last SELECT;
        }
	    else {
            $logger->level($ERROR);
            $logger->debug("Log level is ERROR") ;
            last SELECT;
        }
    }
    return $logger;
}

__END__

=head1 NAME

editDLLconf.pl - 

=head1 SYNOPSIS

editDLLconf.pl [options] 

=head1 DESCRIPTION


=head1 OPTIONS

=over 8

=item B<-h, --help>       

Display help message

=item B<-m, --man>

Display whole man page


    "delete=s" => \$delete,

=item B<-d, --directory>

 \@directories_to_search,

=item B<-f, --findPattern>



=item B<-i,--insert>

Name of file which contains the information to be added to the files found given -d and -f.

=item B<-r, --replace>



=item B<-rb, --replaceBy>



=item B<-v, --verbose>

Increases the verbosity level. Can be used multiple times (verbosest if used 3 or more times) 

=cut
