#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: blat.pl
#
#        USAGE: ./blat.pl
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

## Loading modules and initializing variables ##
use strict;
use warnings;
use utf8;

`blat database.fa query.fa -q=rna -out=blast8 out.psl`;

#system("blat", "-q=rna", "-minIdentity=98", "-out=blast8", "database.fa", "query.fa", "out.blast_out");
