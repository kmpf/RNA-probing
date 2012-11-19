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

my $sample_size = 3;
my $seq = "ACAATTTTAGCTTAATATAAAGCATCTGATTTGCGTTCAGAAGATGTGAGTGTGTCT";
my @probing_profile = (0) x length($seq);


# definition of probing reagent
my $prob_seq = "NNN"; # specific sequence
my $prob_str = "HHH"; # specific sec. structure
my $prob_cut = "_|_"; # modification point



my @structures = &stochastic_sampling($seq, $sample_size);

my @structure_description = &db2sd(@structures);

@probing_profile = &simulate_probing(\@structure_description, \@probing_profile,
                                     $seq, $prob_seq, $prob_str, $prob_cut);

print("=== Results ===\n");
foreach (@structures) {
    print("$_\n");
}

print(join("", @probing_profile)."\n");

# (my $struct,my $mfe) = RNA::fold($seq);  #predict mfe structure of $seq
# RNA::PS_rna_plot($seq, $struct, "rna.ps");  # write PS plot to rna.ps
#RNA::PS_dot_plot($seq, "dot.ps");          # write dot plot to dot.ps


################################################################################
################################################################################
##
##                              Subroutines
##
################################################################################
################################################################################


################################################################################
#
#             Stochastic sampling of RNA secondary structures
#
################################################################################

sub stochastic_sampling {
    my ($seq, $sample_size) = @_;
    # compute partition function and pair pobabilities
    my $gfe = RNA::pf_fold($seq);
    my @structures;

    for (my $i = 0; $i < $sample_size; $i++){
        push( @structures, RNA::pbacktrack($seq) );
    }

# @structures = ();
# $structures[0] = "...(((..(((.((...))))).(((...)))..((..((...))..))..)))...";

    return @structures;
}

################################################################################
#
#           Convert Dot-Bracket into structure description 
#
################################################################################

sub db2sd {

my (@structures) = @_;
my @structure_description = ();
foreach (@structures) {

    my @db = split('', $_);

    my (@dots, @opening_br, @closing_br) = ();

    # at first we expect every nucleotide to be unpaired "U"
    my @struc_dec = ("U") x length($_);

    # fill arrays with positions
    for (my $i = 0; $i < $#db; $i++) {
        push(@dots, $i) if ( $db[$i] eq "." );
        push(@opening_br, $i) if ( $db[$i] eq "(" );

        # closing bracket found, so lets classify enclosed unpaired nucleotides
        # if there are any
        if ( $db[$i] eq ")" && $dots[$#dots] > $opening_br[$#opening_br] ) {

            # declare enclosing bracket positions
            my $op_br_pos = pop(@opening_br);
            my $cl_br_pos = $i;

            # Declare paired nucleotides
            $struc_dec[$op_br_pos] = "P";
            $struc_dec[$cl_br_pos] = "P";


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
           print("Enclosed unpaired nucleotides: "
                    .join(",", @enclosed_dots)."\n");

            # Hairpin or bulge detected
            if ( $enclosed_stems == 0 ) {
                # Hairpin detected if the enclosed dots reach from opening to
                # closing bracket
                if ($op_br_pos + 1 == $enclosed_dots[$#enclosed_dots] &&
                    $cl_br_pos - 1 == $enclosed_dots[0] ) {            
                    foreach (@enclosed_dots) { $struc_dec[$_] = "H" }
                    print("Hairpin found.\n");
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
                # Interior loop detected if the enclosed dots reach from opening
                # to closing bracket
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
        # No enclosed unpaired nucleotides so lets declare the base pair
        elsif ( $db[$i] eq ")" ) {
            $struc_dec[$i] = "P";
            my $rel_open_br = pop(@opening_br);
            $struc_dec[$rel_open_br] = "P";
        }
    }
    my $str_desc = join("",@struc_dec);
    print("$_\n$str_desc\n");
    push(@structure_description, join("",@struc_dec));
}
return @structure_description;

}

################################################################################
#
#           Probe sequence and secondary structure
#
################################################################################

sub simulate_probing {
    my ($str, $probing_profile, $seq, $prob_seq, $prob_str, $prob_cut) = @_;

    foreach my $str_desc (@{$str}) {
        # Find the modification/cut point
        $prob_cut =~ /\|/;
        my $cut_pos = $-[0];

        # create regex from sequence
        $prob_seq =~ s/N/[ACGTU]/g;

        my @seq_matches = &match_all_positions($prob_seq, $seq, $cut_pos);
        my @str_matches = &match_all_positions($prob_str, $str_desc, $cut_pos);


        if ($#seq_matches > $#str_matches) {
            my $more_matches = join(",",@seq_matches);
            foreach (@str_matches) {
                ${probing_profile}[$_]++ if ($more_matches =~ /,$_,/);
            }
        }
        else {
             my $more_matches = join(",",@str_matches);
            foreach (@seq_matches) {
                ${probing_profile}[$_]++  if ($more_matches =~ /,$_,/);
            }
        }

#        print(join("",@{$probing_profile})."\n");
#        print("$seq\n");
    }
    return @{$probing_profile};
}

sub match_all_positions {
    my ($regex, $string, $cut_pos) = @_;
    my @ret;
    while ($string =~ /(?=$regex)/g) {
        push(@ret,  $-[0] + $cut_pos);
    }
    return @ret
}


