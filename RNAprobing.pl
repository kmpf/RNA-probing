#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: RNAprobing.pl
#
#        USAGE: ./RNAprobing.pl  
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
#      CREATED: 12.09.2012 12:53:10
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Data::Dumper;
use RNA;

my $seq = "AAAATTTTAGCTTAATATAAAGCATCTGATTTGCGTTCAGAAGATGTGAGTGTGTCTTGCAAATTTTA";
# (my $struct,my $mfe) = RNA::fold($seq);  #predict mfe structure of $seq
# RNA::PS_rna_plot($seq, $struct, "rna.ps");  # write PS plot to rna.ps
#RNA::PS_dot_plot($seq, "dot.ps");          # write dot plot to dot.ps


################################################################################
#
#             Stochastic sampling of RNA secondary structures
#
################################################################################

# compute partition function and pair pobabilities
my $gfe = RNA::pf_fold($seq);
print($seq."\n" );

my $structure;
for (my $i = 0; $i < 1; $i++){
    $structure = RNA::pbacktrack($seq);
    print($structure."\n");
}


################################################################################
#
#           Convert Dot-Bracket into structure description 
#
################################################################################

my @db = split('', $structure);

my (@dots, @opening_br, @closing_br) = ();

# at first we expect every nucleotide to be unpaired "U"
my @struc_dec = ("U") x length($structure);
print (join("",@struc_dec)."\n");
print(length($structure)."\n");

# fill arrays with positions
for (my $i = 0; $i < scalar(@db); $i++) {
    push(@dots, $i) if ( $db[$i] eq "." );
    if ( $db[$i] eq "(" ) {
        push(@opening_br, $i);
        $struc_dec[$i] = "P";        
    }

    if ( $db[$i] eq ")" ) {
        $struc_dec[$i] = "P";
    }

    # closing bracket found, so lets classify enclosed unpaired nucleotides
    # if there are any
    if ( $db[$i] eq ")" && $dots[$#dots] > $opening_br[$#opening_br] ) {
        print("Closing bracket found at $i \n");
        # declare enclosing bracket positions
        my $op_br_pos = pop(@opening_br);
        my $cl_br_pos = $i;

        my @enclosed_dots = ();
        # count the number of stems enclosed by the unpaired nucleotides
        my $enclosed_stems = 0;

        # collect all enclosed nucleotides in @enclosed_dots
        while ( $dots[$#dots] > $op_br_pos ) {
            push(@enclosed_dots, pop(@dots));
            # stop if you reached the last enclosed nucleotide
            last if ($dots[$#dots] < $op_br_pos);
            # if found non-continous numbers we found an enclosed stem
            if ( $dots[$#dots]+1 < $enclosed_dots[$#enclosed_dots]) {
                $enclosed_stems++;
            }
        }
       print("Enclosed unpaired nucleotides: ".join(",", @enclosed_dots)."\n");

        # Hairpin or bulge detected
        if ( $enclosed_stems == 0 ) {
            # Hairpin detected if the enclosed dots reach from opening to
            # closing bracket
            if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                $cl_br_pos - 1 == $enclosed_dots[0] ) {            
                foreach my $dot_pos (@enclosed_dots) { $struc_dec[$dot_pos] = "H" }
                print("Hairpin found.\n");
                print(join("",@struc_dec) ."\n");
            }
            # bulge detected if the enclosed dots just touch one bracket
            elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {
                foreach (@enclosed_dots) { $struc_dec[$_] = "B" }
                print("Bulge found.\n");
            }
        }
        # Multi loop with two included stems or an interior loop detected
        elsif ( $enclosed_stems == 1 ) {
            # Interior loop detected if the enclosed dots reach from opening to
            # closing bracket
            if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                $cl_br_pos - 1 == $enclosed_dots[0] ) {
                foreach (@enclosed_dots) { $struc_dec[$_] = "I" }
                print("Interior loop found.\n");
            }
            # Multi loop with two stems detected if the enclosed dots just
            # touch one bracket
            elsif ( $op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] ||
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {
                foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
                print("Multi loop found.\n");
            }
        }
        # Multi loop detected if more than one stem is enclosed
        elsif( $enclosed_stems > 1 ) {
                foreach (@enclosed_dots) { $struc_dec[$_] = "M" }
                print("Multi loop found.\n");
        }
    }
}


print($structure."\n");
print(join("",@struc_dec) ."\n");

################################################################################
#
#           Probe sequence and secondary structure
#
################################################################################





