#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: tagcleaner
#   Date: 2011-10-18
#   Version: 0.12
#
#   Usage:
#      tagcleaner [options]
#
#      Try 'tagcleaner -h' for more information.
#
#    Purpose: TagCleaner can be used to automatically detect and efficiently
#              remove the tag sequences (e.g. WTA primer sequences) from your
#              sequences.
#
#===============================================================================

use strict;
use warnings;

use integer; #force integer arithmetic for bit-parallel sequence matching

#use Data::Dumper; ###
use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile); #for output file(s)
use File::Path;
use Fcntl qw(:flock SEEK_END :DEFAULT); #for log file

$| = 1; # Do not buffer output

my $LINE_WIDTH  = 60;
my $TAG5        = '';
my $TAG5rev     = '';
my $TAG3        = '';
my $TAG3rev     = '';
my $MM5         = 0;
my $MM3         = 0;
my $MMFTF35     = 0;
my $MMFTF5R     = 0;
my $MMFTF3R     = 0;
my $ALPHABET    = 'ACGTN';
my $ALPHABET_IUPAC = $ALPHABET.'RYMKWSBDHV';
my $MATRIX      = ();
my $MATRIX_DEFAULT = {A => {A => 1},
                     C => {        C => 1},
                     G => {                G => 1},
                     T => {                        T => 1},
                     R => {A => 1,         G => 1,         R => 1},
                     Y => {        C => 1,         T => 1,         Y => 1},
                     M => {A => 1, C => 1,                                 M => 1},
                     K => {                G => 1, T => 1,                         K => 1},
                     W => {A => 1,                 T => 1,                                 W => 1},
                     S => {        C => 1, G => 1,                                                 S => 1},
                     B => {        C => 1, G => 1, T => 1,         Y => 1,         K => 1,         S => 1, B => 1},
                     D => {A => 1,         G => 1, T => 1, R => 1,                 K => 1, W => 1,                 D => 1},
                     H => {A => 1, C => 1,         T => 1,         Y => 1, M => 1,         W => 1,                         H => 1},
                     V => {A => 1, C => 1, G => 1,         R => 1,         M => 1,                 S => 1,                         V => 1},
                     N => {A => 1, C => 1, G => 1, T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1}};
my $MATRIX_SUBSET = {A => {A => 1,                         R => 1,         M => 1,         W => 1,                 D => 1, H => 1, V => 1, N => 1},
                     C => {        C => 1,                         Y => 1, M => 1,                 S => 1, B => 1,         H => 1, V => 1, N => 1},
                     G => {                G => 1,         R => 1,                 K => 1,         S => 1, B => 1, D => 1,         V => 1, N => 1},
                     T => {                        T => 1,         Y => 1,         K => 1, W => 1,         B => 1, D => 1, H => 1,         N => 1},
                     R => {A => 1,         G => 1,         R => 1,                                                 D => 1,         V => 1, N => 1},
                     Y => {        C => 1,         T => 1,         Y => 1,                                 B => 1,         H => 1,         N => 1},
                     M => {A => 1, C => 1,                                 M => 1,                                         H => 1, V => 1, N => 1},
                     K => {                G => 1, T => 1,                         K => 1,                 B => 1, D => 1,                 N => 1},
                     W => {A => 1,                 T => 1,                                 W => 1,                 D => 1, H => 1,         N => 1},
                     S => {        C => 1, G => 1,                                                 S => 1, B => 1,                 V => 1, N => 1},
                     B => {        C => 1, G => 1, T => 1,         Y => 1,         K => 1,         S => 1, B => 1,                         N => 1},
                     D => {A => 1,         G => 1, T => 1, R => 1,                 K => 1, W => 1,                 D => 1,                 N => 1},
                     H => {A => 1, C => 1,         T => 1,         Y => 1, M => 1,         W => 1,                         H => 1,         N => 1},
                     V => {A => 1, C => 1, G => 1,         R => 1,         M => 1,                 S => 1,                         V => 1, N => 1},
                     N => {A => 1, C => 1, G => 1, T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1}};
my $MATRIX_UNION  = {A => {A => 1,                         R => 1,         M => 1,         W => 1,                 D => 1, H => 1, V => 1, N => 1},
                     C => {        C => 1,                         Y => 1, M => 1,                 S => 1, B => 1,         H => 1, V => 1, N => 1},
                     G => {                G => 1,         R => 1,                 K => 1,         S => 1, B => 1, D => 1,         V => 1, N => 1},
                     T => {                        T => 1,         Y => 1,         K => 1, W => 1,         B => 1, D => 1, H => 1,         N => 1},
                     R => {A => 1,         G => 1,         R => 1,         M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     Y => {        C => 1,         T => 1,         Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     M => {A => 1, C => 1,                 R => 1, Y => 1, M => 1,         W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     K => {                G => 1, T => 1, R => 1, Y => 1,         K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     W => {A => 1,                 T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1,         B => 1, D => 1, H => 1, V => 1, N => 1},
                     S => {        C => 1, G => 1,         R => 1, Y => 1, M => 1, K => 1,         S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     B => {        C => 1, G => 1, T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     D => {A => 1,         G => 1, T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     H => {A => 1, C => 1,         T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     V => {A => 1, C => 1, G => 1,         R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1},
                     N => {A => 1, C => 1, G => 1, T => 1, R => 1, Y => 1, M => 1, K => 1, W => 1, S => 1, B => 1, D => 1, H => 1, V => 1, N => 1}};
my $MATRIX_EXACT  = {A => {A => 1},
                     C => {        C => 1},
                     G => {                G => 1},
                     T => {                        T => 1},
                     R => {                                R => 1},
                     Y => {                                        Y => 1},
                     M => {                                                M => 1},
                     K => {                                                        K => 1},
                     W => {                                                                W => 1},
                     S => {                                                                        S => 1},
                     B => {                                                                                B => 1},
                     D => {                                                                                        D => 1},
                     H => {                                                                                                H => 1},
                     V => {                                                                                                        V => 1},
                     N => {                                                                                                                N => 1}};
my $INIT_ADD    = 10;
my $CONTINUE    = 0;
my $STATS       = 0;
my $SPLIT5R     = 0;
my $SPLIT3R     = 0;
my $SPLIT35     = 0;
my $PREDICT     = 0;
my $KMER_LEN    = 5;
my $TRIM_WITHIN = 0;
my $RANGE_THRESHOLD  = 40;
my $MEDIAN_THRESHOLD = 10;
my $FILTERED    = 0;
my $INFO        = 0;
my $NOMATCH     = 0;
my $WEB         = 0;
my $WEBID       = 0;
my $MINLEN      = 0;
my $MAX_TAG_LEN = 32; #depending on the system, use up to 64 for 64 bit systems
my $MAX_WEB_LEN = 50; #length of sequence used for frequency plots
my $WHAT        = 'standalone';
my $VERSION     = '0.12';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print "TagCleaner-$WHAT $VERSION\n"; exit; },
            'fastq=s',
            'fasta=s',
            'qual=s',
            'tag5=s',
            'tag3=s',
            'trim_within=i',
            '64',
            'matrix=s',
            'cont',
            'split:i',
            'split5r:i',
            'split3r:i',
            'splitall:i',
            'out_format=i',
            'out=s',
            'line_width',
            'mm5=i',
            'mm3=i',
            'stats',
            'predict',
            'log:s',
            'filtered',
            'nomatch=i',
            'info',
            'minlen=i',
            'web=s' #for web version use
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

TagCleaner - For cleaner sequences.

=head1 VERSION

TagCleaner standalone 0.12

=head1 SYNOPSIS

perl tagcleaner.pl [-h] [-help] [-version] [-man] [-verbose] [-64] [-fastq input_fastq_file] [-fasta input_fasta_file] [-qual input_quality_file] [-out_format int_value] [-out filename_prefix] [-trim_within bases_from_ends] [-cont] [-split number_of_mismatches] [-split5r number_of_mismatches] [-split3r number_of_mismatches] [-splitall number_of_mismatches] [-tag5 5_prime_tag_sequence] [-tag3 3_prime_tag_sequence] [-mm5 number_of_mismatches] [-mm3 number_of_mismatches] [-stats] [-predict] [-log optional_file_name] [-filtered] [-nomatch int_value] [-info] [-minlen minimum_read_length] [-matrix matrix_type]

=head1 DESCRIPTION

The TagCleaner tool can be used to automatically detect and efficiently remove tag sequences (e.g. WTA tags) from genomic and metagenomic datasets. Go to http://tagcleaner.sourceforge.net/ for access to a user-friendly web-interface, the documentation and the latest updates.

Tag sequence trimming should be performed before quality trimming and sequence dereplication. The trimming of low-quality bases at the ends might truncate the tag sequence and reduce the ability to recognize the remainder of the tag sequence. In those cases, large parts of the tag sequences might still remain for further analysis and data processing steps. Dereplication before trimming may miss duplicated sequences due to variations in the tag sequences that will be trimmed off later and would therefore require an additional dereplication step after the trimming.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<-64>

This option should only be used with caution. The default value for the maximum tag length is set to 32 assuming a 32 bit system. If you have a 64 bit system and need to work with tags longer than 32 bp, you can use this option to tell the program to assume a 64 bit system. If you set 64 for a 32 bit system, it will produce incorrect results.

=item B<***** INPUT OPTIONS *****>

=item B<-fastq> <file>

Input file in FASTQ format that contains the sequence and quality data.

=item B<-fasta> <file>

Input file in FASTA format that contains the sequence data.

=item B<-qual> <file>

Input file in QUAL format that contains the quality data.

=item B<-matrix> <matrix_type>

Use union, subset or exact matrix instead of default matrix when matching tag characters to the sequence characters. The default matrix requires that the ambiguous characters of the sequence represent bases that are a subset of the ambigous characters of the tag sequence (e.g. R = AG in the sequence is a subset of D = AGT and is marked with an X in the matrix). The subset matrix requires that either the tag character or the sequence character is a subset of the other. The union matrix requires that at least one base is common in tag and sequence (e.g. R = AG and M = AC have base A in common and are marked with an X in the matrix). The exact matrix allows to trim ambiguous bases without considering the ambiguity code (e.g. N in tag only trims N, not ACGT). Optional options are either "subset", "union" or "exact".

Default matrix:

      A C G T R Y M K W S B D H V N (tag)
    A X       X   X   X     X X X X
    C   X       X X     X X   X X X
    G     X   X     X   X X X   X X
    T       X   X   X X   X X X   X
    R         X             X   X X
    Y           X         X   X   X
    M             X           X X X
    K               X     X X     X
    W                 X     X X   X
    S                   X X     X X
    B                     X       X
    D                       X     X
    H                         X   X
    V                           X X
    N                             X
    (sequence)

Subset matrix:

      A C G T R Y M K W S B D H V N (tag)
    A X       X   X   X     X X X X
    C   X       X X     X X   X X X
    G     X   X     X   X X X   X X
    T       X   X   X X   X X X   X
    R X   X   X             X   X X
    Y   X   X   X         X   X   X
    M X X         X           X X X
    K     X X       X     X X     X
    W X     X         X     X X   X
    S   X X             X X     X X
    B   X X X   X   X   X X       X
    D X   X X X     X X     X     X
    H X X   X   X X   X       X   X
    V X X X   X   X     X       X X
    N X X X X X X X X X X X X X X X
    (sequence)

Union matrix:

      A C G T R Y M K W S B D H V N (tag)
    A X       X   X   X     X X X X
    C   X       X X     X X   X X X
    G     X   X     X   X X X   X X
    T       X   X   X X   X X X   X
    R X   X   X   X X X X X X X X X
    Y   X   X   X X X X X X X X X X
    M X X     X X X   X X X X X X X
    K     X X X X   X X X X X X X X
    W X     X X X X X X   X X X X X
    S   X X   X X X X   X X X X X X
    B   X X X X X X X X X X X X X X
    D X   X X X X X X X X X X X X X
    H X X   X X X X X X X X X X X X
    V X X X   X X X X X X X X X X X
    N X X X X X X X X X X X X X X X
    (sequence)

Exact matrix:

      A C G T R Y M K W S B D H V N (tag)
    A X
    C   X
    G     X
    T       X
    R         X
    Y           X
    M             X
    K               X
    W                 X
    S                   X
    B                     X
    D                       X
    H                         X
    V                           X
    N                             X
    (sequence)

=item B<***** OUTPUT OPTIONS *****>

=item B<-out> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with an additional "_tagcleaner_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option.

Example: use "file_passed" to generate the output file file_passed.fasta (fasta output) in the current directory

=item B<-out_format> <integer>

To change the output format, use one of the following options. If not defined, the output format will be the same as the input format.

1 (FASTA only), 2 (FASTA and QUAL) or 3 (FASTQ)

=item B<-line_width> <integer>

Sequence characters per line. Use 0 if you want each sequence in a single line. Use 80 for line breaks every 80 characters. Note that this option only applies to FASTA output files, since FASTQ files store sequences without additional line breaks. [default: 60]

=item B<-stats>

Prints the number of tag sequences matching for different numbers of mismatches. In combination with -split, the number of sequences with fragment-to-fragment concatenations is printed as well. Cannot be used in combination with -predict and will not perform any trimming. The output values are separated by tabs with the header line: "#Param Mismatches_or_Splits Number_of_Sequences Percentage Percentage_Sum". Cannot be used in combination with -predict and require -tag5 or -tag3.

=item B<-predict>

Use this option to have TagCleaner predict the tag sequences. It will attempt to predict the tag at either or both sites, if possible. The algorithm implemented for the tag prediction assumes the randomness of a typical metagenome. Datasets that do not contain random sequences from organisms in an environment, but rather contain, for example, 16S data may cause incorrect detection of the tag sequences. However, the tag sequences will most likely be over-predicted and can be redefined by the user prior to data processing. The tag sequence prediction uses filtered base frequencies instead of raw base frequencies. This allows a more accurate prediction as it accounts for incomplete and shifted tag sequences. The output values are separated by tabs with the header line: "#Param Tag_Sequence Tag_Length Percent_Explained". If no tags are reported, then no tags could be identified in the data set. Cannot be used in combination with -tag3 or -tag5 or -stats. When using this option, no trimming will be performed.

=item B<-log> <file>

Log file to keep track of parameters, errors, etc. The log file name is optional. If no file name is given, the log file name will be "inputname.log". If the log file already exists, new content will be added to the file.

=item B<-info>

This option will provide the trimming and splitting information in the header line after the sequence identifier. The following information is given and separated by a single space: initial length, length after trimming, 5'-end trimming position, 3'-end trimming position, number of mismatches at 5'-end, number of mismatches at 3'-end and number of sequences after splitting. In case of a splitting event, the number of mismatches at the 5'- and 3'-end will be the number of mismatches of the concatenated tags.

=item B<-nomatch> <integer>

This option allows to filter sequences that do not match the tag sequence at the ends or do not contain inner tags within the maximum number of allowed mismatches. The following options allow to filter reads B<not matching> the tag at:

(1) 5'-end, (2) 3'-end, (3) either end - requires both ends to match, (4) both ends - requires either end to match, or (5) reads B<splitting> into two or more reads.

=item B<-minlen> <integer>

Minimum read length of trimming and/or splitting.

=item B<-filtered>

Output the sequences that would be filtered out instead of the sequences passing the filters. This includes sequences that e.g. are tag sequence repeats or do not fulfill the -nomatch option.

=item B<***** TRIM / SPLIT OPTIONS *****>

=item B<-tag5> <string>

Tag sequence at 5'-end. Use option -predict if unknown.

=item B<-tag3> <string>

Tag sequence at 3'-end. Use option -predict if unknown.

=item B<-mm5> <integer>

Maximum number of allowed mismatches at the 5'-end. [default: 0]

=item B<-mm3> <integer>

Maximum number of allowed mismatches at the 3'-end. The independent definition for the 5'- and 3'-end of the reads accounts for the differences in tag sequences due to the limitations of the sequencing method used to generate the datasets. The 3'-end will in most cases show a lower number of matching tag sequences with low number of mismatches due to incomplete or missing tags at the ends of incompletely sequenced fragments. [default: 0]

=item B<-trim_within> <integer>

The sequence of the tag could occur not only at the sequence end, but also at any other position of the sequence. To assure that only tags are trimmed, the tag sequences can be defined to occur only at the ends allowing a certain number of variable bases. The default value for -trim_within for a tag sequence of at least 10 bp is 1.5 times the tag sequence length. [default: 1.5x tag length]

Example: Use -tag3 NNNNNNNNNCCAAACACACCCAACACACCAC -trim_within 60 to trim B<ATCCATTTCCCAAACACACCCAACACACCAC>AAAAAAAAAAAAAAAACAAACAACACC

=item B<-cont>

Trim tag sequences continuously. This is helpful if you have sequence with tag sequence repeats or sequence that are concatenated tag sequences. Note that a high number of allowed mismatches and continuous trimming can cause over-trimming. Use more than 20% mismatches with continuous trimming only with caution.

=item B<-split> <integer>

This features removes tag contaminations inside the sequences and splits fragment-to-fragment concatenations into separate sequences. The optional integer value specifies the maximum number of allowed mismatches for the internal (concatenated) tag sequence(s). This feature should be used with caution for inputs with only a 5' or 3' tag sequence (likely splits too many false positive that naturally occur for single tags compared to much longer concatenated 5' and 3' tags). The number of mismatches is used as maximum value for -stats. This option will cause a decrease in speed. [default: 0]

sequence1-tag3-tag5-sequence2

=item B<-split5r> <integer>

This feature is similar to -split, but instead of search for 3'-5'-tag repeats, it will search for 5'-tag repeats (2 or more). This option only applies if both -tag3 and -tag5 are specified. To split at a single 5'-tag, run the program without -tag3 and with -split. [default: 0]

sequence1-tag5-tag5-sequence2

=item B<-split3r> <integer>

This feature is similar to -split, but instead of search for 3'-5'-tag repeats, it will search for 3'-3'-tag repeats (2 or more). This option only applies if both -tag3 and -tag5 are specified. To split at a single 3'-tag, run the program without -tag5 and with -split. [default: 0]

sequence1-tag3-tag3-sequence2

=item B<-splitall> <integer>

This feature is for convenience only and applies the integer value to all split options (-split, -split3r and -split5r). This option only applies if both -tag3 and -tag5 are specified and will overwrite all other split options with the given integer value. [default: 0]

Combining all split options can split things like sequence1-tag3-tag3-tag5-tag5-tag5-sequence2 into sequence1 and sequence2.


=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >> so that I can make TagCleaner better.

=head1 COPYRIGHT

Copyright (C) 2010-2011  Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

my ($file1,$command);

#check for the system option if 32 or 64 bit should be assumed
if(exists $params{64}) {
    $command .= ' -64';
    $MAX_TAG_LEN = 64;
}

#check for iupac option
if(exists $params{matrix}) {
    $command .= ' -mattrix '.$params{matrix}||'';
    if($params{matrix} eq 'union') {
        $MATRIX = $MATRIX_UNION;
    } elsif($params{matrix} eq 'subset') {
        $MATRIX = $MATRIX_SUBSET;
    } elsif($params{matrix} eq 'exact') {
        $MATRIX = $MATRIX_EXACT;
    } else {
        &printError("unknow matrix type \"".$params{matrix}."\". Please use either \"union\", \"subset\" or \"exact\" as matrix type");
    }
} else {
    $MATRIX = $MATRIX_DEFAULT;
}

#Check if output dir is defined and if not, if it can be created
if(exists $params{out}) {
    $command .= ' -out '.$params{out};
    ($params{out_dir},$params{out_file}) = &getPathAndFile($params{out});
    delete($params{out_dir}) unless($params{out_dir});
}
if(exists $params{out_dir}) {
    mkpath($params{out_dir}, { verbose => 0, mode => 0711, error => \my $err });
    if(@$err) {
        &printError("problems while trying to create directory \"".$params{'out_dir'}."\"");
    }
} else {
    $params{out_dir} = cwd();
}
$params{out_dir} = cwd().'/'.$params{out_dir} unless($params{out_dir} =~ /^\//);
$params{out_dir} .= '/' unless($params{out_dir} =~ /\/$/);

#Check if input file exists and check if file format is correct
if(exists $params{fasta} && exists $params{fastq}) {
    &printError('fasta and fastq cannot be used together');
} elsif(exists $params{fasta}) {
    $command .= ' -fasta '.$params{fasta};
    if(-e $params{fasta}) {
        $file1 = $params{fasta};
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fasta') {
            &printError('input file for -fasta is in '.uc($format).' format not in FASTA format');
        }
    } else {
        &printError("could not find input file \"".$params{fasta}."\"");
    }
} elsif(exists $params{fastq}) {
    $command .= ' -fastq '.$params{fastq};
    if(-e $params{fastq}) {
        $file1 = $params{fastq};
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fastq') {
            &printError('input file for -fastq is in '.uc($format).' format not in FASTQ format');
        }
    } else {
        &printError("could not find input file \"".$params{fastq}."\"");
    }
} else {
    &printError("you did not specify an input file containing the query sequences");
}
if(exists $params{fastq} && exists $params{qual}) {
    &printError('fastq and qual cannot be used together');
} elsif(exists $params{qual}) {
    $command .= ' -qual '.$params{qual};
    if(-e $params{qual}) {
        #check for file format
        my $format = &checkFileFormat($params{qual});
        unless($format eq 'qual') {
            &printError('input file for -qual is in '.uc($format).' format not in QUAL format');
        }
    } else {
        &printError("could not find input file \"".$params{qual}."\"");
    }
}

#check if anything todo
unless( exists $params{tag5} ||
        exists $params{tag3} ||
        exists $params{stats} ||
        exists $params{predict}
        ) {
    &printError('nothing to do with input data');
}

#check if output format is possible
if(exists $params{out_format}) {
    $command .= ' -out_format '.$params{out_format};
    if($params{out_format} =~ /\D/) {
        &printError('output format option has to be an integer value');
    } elsif($params{out_format} == 2 || $params{out_format} == 3) {
        unless(exists $params{fastq} || exists $params{qual}) {
            &printError('cannot use this output format option without providing quality data as input');
        }
    } elsif($params{out_format} != 1) {
        &printError('output format option not available');
    }
} else {
    if(exists $params{fastq}) {
        $params{out_format} = 3;
    } elsif(exists $params{fasta} && exists $params{qual}) {
        $params{out_format} = 2;
    } else {
        $params{out_format} = 1;
    }
}

#check if stats and predict
if(exists $params{stats} && exists $params{predict}) {
    &printError('-stats and -predict cannot be used together');
}

#check stats and prevent out of files for stats
if(exists $params{stats}) {
    $command .= ' -stats';
    $params{out} = 'null';
    $STATS = 1;
} else {
    $STATS = 0;
}

#check if tag sequence prediction and prevent out of files for stats
if(exists $params{predict}) {
    $command .= ' -predict';
    $params{out} = 'null';
    $PREDICT = 1;
} else {
    $PREDICT = 0;
}

#check if output needs to be generated
if(exists $params{out} && $params{out} eq 'null' && !$STATS && !$PREDICT) {
    &printError('no output selected (set to null)');
}

#check tag sequences
if(exists $params{tag5}) {
    $command .= ' -tag5 '.$params{tag5};
    if(&checkTagSequence($params{tag5})) {
        $TAG5 = $params{tag5};
        $TAG5rev = &rev($TAG5);
    } else {
        &printError("the tag sequence \"".$params{tag5}."\" is not in valid format: ");
    }
    #check if tag length is within possible length restricted by the architecture (32 or 64 for 32bit or 64bit systems)
    if(length($TAG5) > $MAX_TAG_LEN) {
        &printError("The tag sequence is too long. The maximum allowed length currently set for your system is $MAX_TAG_LEN. If you are using a 64 bit system and you see 32 as maximum length, you can change this value to 64 using the -system option");
    }
}
if(exists $params{tag3}) {
    $command .= ' -tag3 '.$params{tag3};
    if(&checkTagSequence($params{tag3})) {
        $TAG3 = $params{tag3};
        $TAG3rev = &rev($TAG3);
    } else {
        &printError("the tag sequence \"".$params{tag3}."\" is not in valid format: ");
    }
    #check if tag length is within possible length restricted by the architecture (32 or 64 for 32bit or 64bit systems)
    if(length($TAG3rev) > $MAX_TAG_LEN) {
        &printError("The tag sequence is too long. The maximum allowed length currently set for your system is $MAX_TAG_LEN. If you are using a 64 bit system and you see 32 as maximum length, you can change this value to 64 using the -system option");
    }
}
unless($TAG3 || $TAG5 || $PREDICT) {
    &printError("no tag sequence specified");
}

#checking number of mismatches
if(exists $params{mm5}) {
    $command .= ' -mm5 '.$params{mm5};
    if($params{mm5} =~ m/^\d+$/) {
        if($params{mm5} < 0) {
            &printError("The number of allowed mismatches has to be a positive number");
        } else {
            $MM5 = $params{mm5};
        }
    } else {
        &printError("The number of allowed mismatches has to be an integer number");
    }
}
if(exists $params{mm3}) {
    $command .= ' -mm3 '.$params{mm3};
    if($params{mm3} =~ m/^\d+$/) {
        if($params{mm3} < 0) {
            &printError("The number of allowed mismatches has to be a positive number");
        } else {
            $MM3 = $params{mm3};
        }
    } else {
        &printError("The number of allowed mismatches has to be an integer number");
    }
}

#check for read length value
if(exists $params{minlen}) {
    $command .= ' -minlen '.$params{minlen};
    if($params{minlen} =~ m/^\d+$/) {
        if($params{minlen} < 1) {
            &printError("The minimum read length after trimming and/or splitting has to be a positive number greater than 0");
        } else {
            $MINLEN = $params{minlen};
        }
    } else {
        &printError("The minimum read length after trimming and/or splitting has to be an integer number");
    }
}

#checking if trim_within is valid
if(exists $params{trim_within}) {
    $command .= ' -trim_within '.$params{trim_within};
    if($params{trim_within} =~ m/^\d+$/) {
        if($params{trim_within} < length($TAG5) || $params{trim_within} < length($TAG3)) {
             &printError("The -trim_within value has to be at least the number of bases in the tag sequence");
        } else {
            $TRIM_WITHIN = $params{trim_within};
        }
    } else {
        &printError("The -trim_within value has to be a positive integer number");
    }
}
#check for continuously trimming
if(exists $params{cont}) {
    $command .= ' -cont';
    $CONTINUE = 1;
} else {
    $CONTINUE = 0;
}

#check for filtered
if(exists $params{filtered}) {
    $command .= ' -filtered';
    $FILTERED = 1;
} else {
    $FILTERED = 0;
}

#check for info
if(exists $params{info}) {
    $command .= ' -info';
    $INFO = 1;
} else {
    $INFO = 0;
}

#check for web predict
if(exists $params{web}) {
    $command .= ' -web '.$params{web};
    $WEB = 1;
    $WEBID = $params{web};
} else {
    $WEB = 0;
}

#check if splitall
if(exists $params{splitall}) {
    if($params{splitall} >= 0) {
        $params{split} = $params{splitall};
        $params{split5r} = $params{splitall};
        $params{split3r} = $params{splitall};
    } else {
        &printError("the number of mismatches for tag sequences inside the sequence cannot be negative");
    }
}

#check if split
if(exists $params{split}) {
    $command .= ' -split '.$params{split};
    $SPLIT35 = 1;
    if($params{split} >= 0) {
        $MMFTF35 = $params{split};
    } else {
        &printError("the number of mismatches for tag sequences inside the sequence cannot be negative");
    }
} else {
    $SPLIT35 = 0;
}
if(exists $params{split5r}) {
    $command .= ' -split5r '.$params{split5r};
    $SPLIT5R = 1;
    if($params{split5r} >= 0) {
        $MMFTF5R = $params{split5r};
    } else {
        &printError("the number of mismatches for tag sequences inside the sequence cannot be negative");
    }
} else {
    $SPLIT5R = 0;
}
if(exists $params{split3r}) {
    $command .= ' -split3r '.$params{split3r};
    $SPLIT3R = 1;
    if($params{split3r} >= 0) {
        $MMFTF3R = $params{split3r};
    } else {
        &printError("the number of mismatches for tag sequences inside the sequence cannot be negative");
    }
} else {
    $SPLIT3R = 0;
}

#check for no match
if(exists $params{nomatch}) {
    $command .= ' -nomatch '.$params{nomatch};
    if(($params{nomatch} == 1 && $TAG5) || ($params{nomatch} == 2 && $TAG3) || ($params{nomatch} == 3 && $TAG5 && $TAG3) || ($params{nomatch} == 4 && $TAG5 && $TAG3) || ($params{nomatch} == 5 && $TAG5 && $TAG3)) {
        $NOMATCH = $params{nomatch};
    } else {
        &printError("the value specified for -nomatch is invalid. The valid values are either of 1, 2, 3, 4 or 5 and require the respective tag(s) to be defined");
    }
} else {
    $NOMATCH = 0;
}

#set remaining parameters
if($params{out_format} == 3) {
    $LINE_WIDTH = 0;
} elsif(exists $params{line_width}) {
    $command .= ' -line_width '.$params{line_width};
    $LINE_WIDTH = $params{line_width};
}

if(exists $params{verbose}) {
    $command .= ' -verbose';
}

#Check for log file
if(exists $params{log}) {
    $command .= ' -log'.(length($params{log}) ? ' '.$params{log} : '');
    unless($params{log}) {
        $params{log} = join("__",$file1||'nonamegiven').'.log';
    }
    $params{log} = cwd().'/'.$params{log} unless($params{log} =~ /^\//);
    &printLog("Executing TagCleaner with command: \"perl tagcleaner.pl".$command."\"");
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

my $filename = $file1;
while($filename =~ /[\w\d]+\.[\w\d]+/) {
    $filename =~ s/\.[\w\d]+$//;
}

#create filehandles for the output data
my ($fhgood,$fhgood2);
my ($filenamegood,$filenamegood2);
my ($nogood);
$nogood = 0;
if(exists $params{out}) {
    if($params{out} eq 'null') {
        $nogood = 1;
    } else {
        open($fhgood,">".$params{out}.'.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
        $filenamegood = $params{out}.'.fast'.($params{out_format} == 3 ? 'q' : 'a');
    }
} else {
    $fhgood = File::Temp->new( TEMPLATE => $filename.'_tagcleaner_XXXX',
                               SUFFIX => '.fast'.($params{out_format} == 3 ? 'q' : 'a'),
                               UNLINK => 0);
    $filenamegood = $fhgood->filename;
}
if($filenamegood) {
    &printLog("Write results to file: \"".$filenamegood."\"");
    &addToIdFile($WEBID,'generate_file',$filenamegood) if($WEB);
}
if($params{out_format} == 2) {
    if(exists $params{out}) {
        unless($nogood) {
            open($fhgood2,">".$params{out}.'.qual') or &printError('cannot open output file');
            $filenamegood2 = $params{out}.'.qual';
        }
    } else {
        $fhgood2 = File::Temp->new( TEMPLATE => $filename.'_tagcleaner_XXXX',
                                    SUFFIX => '.qual',
                                    UNLINK => 0);
        $filenamegood2 = $fhgood2->filename;
    }
}
if($filenamegood2) {
    &printLog("Write quaity results to file: \"".$filenamegood2."\"");
}

my $numlines = 0;
my ($progress,$counter,$part);
$progress = $counter = $part = 1;
if(exists $params{verbose} || $WEB) {
    print STDERR "Estimate size of input data for status report (this might take a while for large files)\n" unless($WEB);
    $numlines = &getLineNumber($file1);
    print STDERR "\tdone\n" unless($WEB);

    #for progress bar
    $progress = 0;
    $counter = 1;
    $part = int($numlines/100);
}

#pre-process tag sequences into Peq hash
my $Peq5 = &getPeq($TAG5,$ALPHABET_IUPAC);
my $Peq3 = &getPeq($TAG3rev,$ALPHABET_IUPAC);
my $Peq5rev = &getPeq($TAG5rev,$ALPHABET_IUPAC);
my $Peq3rev = &getPeq($TAG3,$ALPHABET_IUPAC);

#number of Ns at the ends
my @ns;
(my $ntmp = $TAG5) =~ s/^(N*)[^N]*(N*)$/$1$2/;
push(@ns,length($ntmp));
($ntmp = $TAG3) =~ s/^(N*)[^N]*(N*)$/$1$2/;
push(@ns,length($ntmp));

#for approx matching algo
my $m5 = length($TAG5);
my $m3 = length($TAG3rev);

#length of substr where to first search for matching pattern
my $initlen5 = ($TRIM_WITHIN ? $TRIM_WITHIN : $m5+&max($INIT_ADD,($m5/2)));
my $initlen3 = ($TRIM_WITHIN ? $TRIM_WITHIN : $m3+&max($INIT_ADD,($m3/2)));

#minimum number of bases to move if too many mismatches
my $minmove5 = &max($INIT_ADD,($m5/2));
my $minmove3 = &max($INIT_ADD,($m3/2));

#define Pv0 = Pv used for each new sequence
my $Pv05 = 2**$m5-1;
my $Pv03 = 2**$m3-1;

#get sp values for loop in patsearch
my @sp;
push(@sp,0) unless(!$TAG5);
push(@sp,1) unless(!$TAG3rev);
my $sps = (scalar(@sp) > 1 ? 1 : 0);


#parse input data
print STDERR "Parse and process input data\n" if(exists $params{verbose});
print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
&printLog("Parse and process input data: \"$file1\"".($params{qual} ? " and \"".$params{qual}."\"" : ''));
my ($tsvfile,$seqid,$header,$seq,$qual,$count,$numbases,$numseqs,$numbasesout,$numseqsout,$length);
$numseqs = $numseqsout = $numbases = $numbasesout = 0;
$count = 0;
$qual = '';
$seq = '';
#stats data
my %stats;
#prediction data
my (%kmers,%basefreqs);
#parse data
my $seqcount = 0;

open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file1 |") or &printError("Could not open file $file1: $!");
open(FILE2,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < ".$params{qual}." |") or &printError("Could not open file ".$params{qual}.": $!") if(exists $params{qual});
if(exists $params{fastq}) {
    while(<FILE>) {
        chomp();
        if($count == 0 && /^\@(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                $numseqs++;
                $numbases += $length;
                if($PREDICT) {
                    &getKmersAndBaseFreqs(\%kmers,\%basefreqs,$seq,$length,$KMER_LEN,&max($MAX_WEB_LEN,$MAX_TAG_LEN));
                } else {
                    #identify tag
                    my ($trims,$minscores,$splitstatus) = &searchForTag([$m5,$m3],[$Peq5,$Peq3,$Peq5rev,$Peq3rev],[$Pv05,$Pv03],[$initlen5,$initlen3],[$MM5,$MM3],[$minmove5,$minmove3],$seq,$length,\@sp,\@ns);
                    #process sequence
                    if($STATS) {
                        #calc summary stats
                        foreach my $k (keys %$splitstatus) {
                            $stats{$k}->{$splitstatus->{$k}}++;
                        }
                        $stats{tag5}->{$minscores->[0]}++ if($TAG5);
                        $stats{tag3}->{$minscores->[($TAG5 ? 1 : 0)]}++ if($TAG3rev);
                    } else {
                        #process data
                        &processData($seqid,$seq,$qual,$header,$trims,$length);
                    }
                }
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
            $qual = '';
        } elsif($count == 1) {
            $seq = $_;
        } elsif($count == 3) {
            $qual = $_;
            $count = -1;
        }
        $count++;
        #progress bar stuff
        $counter++;
        if($counter > $part) {
            $counter = 1;
            $progress++;
            $progress = 99 if($progress > 99);
            print STDERR "\r\tstatus: ".$progress." \%" if(exists $params{verbose});
            &addToIdFile($WEBID,'status_tc',$progress) if($WEB);
        }
    }
    #add last one
    $length = length($seq);
    if($seq && $length) {
        #remove space and dash from sequences
        $seq =~ s/[\s\-]//g;
        $numseqs++;
        $numbases += $length;
        if($PREDICT) {
            &getKmersAndBaseFreqs(\%kmers,\%basefreqs,$seq,$length,$KMER_LEN,&max($MAX_WEB_LEN,$MAX_TAG_LEN));
        } else {
            #identify tag
            my ($trims,$minscores,$splitstatus) = &searchForTag([$m5,$m3],[$Peq5,$Peq3,$Peq5rev,$Peq3rev],[$Pv05,$Pv03],[$initlen5,$initlen3],[$MM5,$MM3],[$minmove5,$minmove3],$seq,$length,\@sp,\@ns);
            #process sequence
            if($STATS) {
                #calc summary stats
                foreach my $k (keys %$splitstatus) {
                    $stats{$k}->{$splitstatus->{$k}}++;
                }
                $stats{tag5}->{$minscores->[0]}++ if($TAG5);
                $stats{tag3}->{$minscores->[($TAG5 ? 1 : 0)]}++ if($TAG3rev);
            } else {
                #process data
                &processData($seqid,$seq,$qual,$header,$trims,$length);
            }
        }
    }
} elsif(exists $params{fasta}) {
    while(<FILE>) {
        chomp();
        if(/^>(\S+)\s*(.*)$/) {
            $length = length($seq);
            if($seq && $length) {
                #get qual data if provided
                if(exists $params{qual} && !$PREDICT) {
                    while(<FILE2>) {
                        chomp();
                        last if(/^>/ && $qual);
                        next if(/^>/);
                        $qual .= $_.' ';
                    }
                }
                #remove space and dash from sequences
                $seq =~ s/[\s\-]//g;
                #print to tsv file
                $numseqs++;
                $numbases += $length;
                if($PREDICT) {
                    &getKmersAndBaseFreqs(\%kmers,\%basefreqs,$seq,$length,$KMER_LEN,&max($MAX_WEB_LEN,$MAX_TAG_LEN));
                } else {
                    #identify tag
                    my ($trims,$minscores,$splitstatus) = &searchForTag([$m5,$m3],[$Peq5,$Peq3,$Peq5rev,$Peq3rev],[$Pv05,$Pv03],[$initlen5,$initlen3],[$MM5,$MM3],[$minmove5,$minmove3],$seq,$length,\@sp,\@ns);
                    #process sequence
                    if(exists $params{stats}) {
                        #calc summary stats
                        foreach my $k (keys %$splitstatus) {
                            $stats{$k}->{$splitstatus->{$k}}++;
                        }
                        $stats{tag5}->{$minscores->[0]}++ if($TAG5);
                        $stats{tag3}->{$minscores->[($TAG5 ? 1 : 0)]}++ if($TAG3rev);
                    } else {
                        #process data
                        &processData($seqid,$seq,&convertQualNumsToAsciiString($qual),$header,$trims,$length);
                    }
                }
                $qual = '';
            }
            $seqid = $1;
            $header = $2 || '';
            $seq = '';
        } else {
            $seq .= $_;
        }
        #progress bar stuff
        $counter++;
        if($counter > $part) {
            $counter = 1;
            $progress++;
            $progress = 99 if($progress > 99);
            print STDERR "\r\tstatus: ".$progress." \%" if(exists $params{verbose});
            &addToIdFile($WEBID,'status_tc',$progress) if($WEB);
        }
    }
    #add last one
    $length = length($seq);
    if($seq && $length) {
        #get qual data if provided
        if(exists $params{qual} && !$PREDICT) {
            while(<FILE2>) {
                chomp();
                last if(/^>/ && $qual);
                next if(/^>/);
                $qual .= $_.' ';
            }
        }
        #remove space and dash from sequences
        $seq =~ s/[\s\-]//g;
        #print to tsv file
        $numseqs++;
        $numbases += $length;
        if($PREDICT) {
            &getKmersAndBaseFreqs(\%kmers,\%basefreqs,$seq,$length,$KMER_LEN,&max($MAX_WEB_LEN,$MAX_TAG_LEN));
        } else {
            #identify tag
            my ($trims,$minscores,$splitstatus) = &searchForTag([$m5,$m3],[$Peq5,$Peq3,$Peq5rev,$Peq3rev],[$Pv05,$Pv03],[$initlen5,$initlen3],[$MM5,$MM3],[$minmove5,$minmove3],$seq,$length,\@sp,\@ns);
            #process sequence
            if(exists $params{stats}) {
                #calc summary stats
                foreach my $k (keys %$splitstatus) {
                    $stats{$k}->{$splitstatus->{$k}}++;
                }
                $stats{tag5}->{$minscores->[0]}++ if($TAG5);
                $stats{tag3}->{$minscores->[($TAG5 ? 1 : 0)]}++ if($TAG3rev);
            } else {
                #process data
                &processData($seqid,$seq,&convertQualNumsToAsciiString($qual),$header,$trims,$length);
            }
        }
        $qual = '';
    }
}

print STDERR "\r\tdone          \n" if(exists $params{verbose});

&printLog("Finished processing input data: \"$file1\"".($params{qual} ? " and \"".$params{qual}."\"" : ''));
if($WEB) {
    &addToIdFile($WEBID,'status_tc',100);
    &addToIdFile($WEBID,'generate_num',$numseqsout);
}

close(FILE);
close(FILE2) if(exists $params{qual});

if(!$STATS && !$PREDICT) {
    print STDERR "Clean up empty files\n" if(exists $params{verbose} && $nogood);
    #close filehandles
    close($fhgood) unless($nogood);
    if($params{out_format} == 2) {
        close($fhgood2) unless($nogood);
    }
}

if(exists $params{verbose} && !$STATS && !$PREDICT) {
    { no integer;
        print STDERR "Input/Output stats:\n";
        print STDERR "\tInput sequences: ".&addCommas($numseqs)."\n";
        print STDERR "\tInput bases: ".&addCommas($numbases)."\n";
        print STDERR "\tInput mean length: ".sprintf("%.2f",$numbases/$numseqs)."\n";
        print STDERR "\tOutput sequences: ".&addCommas($numseqsout)."\n";
        print STDERR "\tOutput bases: ".&addCommas($numbasesout)."\n";
        print STDERR "\tOutput mean length: ".sprintf("%.2f",$numbasesout/$numseqsout)."\n";
        &printLog("Input sequences: ".&addCommas($numseqs));
        &printLog("Input bases: ".&addCommas($numbases));
        &printLog("Input mean length: ".sprintf("%.2f",$numbases/$numseqs));
        &printLog("Output sequences: ".&addCommas($numseqsout));
        &printLog("Output bases: ".&addCommas($numbasesout));
        &printLog("Output mean length: ".sprintf("%.2f",$numbasesout/$numseqsout));
    }
}

if($STATS) {
    if($WEB) {
        &addToIdFile($WEBID,'status_tc',100);
        &addToIdFile($WEBID,'numseqs',$numseqs);
        &addToIdFile($WEBID,'numbases',$numbases);
    } else {
        print STDOUT '#'.join("\t", qw(Param Number_of_Mismatches_or_Splits Number_of_Sequences Percentage Percentage_Sum))."\n";
    }
    { no integer;
        foreach my $k (sort keys %stats) {
            my $sum = 0;
            foreach my $n (sort {$a <=> $b} keys %{$stats{$k}}) {
                $sum += $stats{$k}->{$n};
                if($WEB) {
                    &addToIdFile($WEBID,'tag'.($k eq 'tag5' ? 0 : 1).'mm'.$n,join(",",$stats{$k}->{$n},sprintf("%.2f",100/$numseqs*$stats{$k}->{$n})));
                } else {
                    print STDOUT join("\t",$k,$n,$stats{$k}->{$n},sprintf("%.2f",100/$numseqs*$stats{$k}->{$n}),sprintf("%.2f",100/$numseqs*$sum))."\n";
                }
            }
        }
    } #no integer
}

if($PREDICT) {
    #filter and shift kmers
    my ($kmersum,$kmershift) = &filterAndShiftKmers(\%kmers);

    #filter and shift base frequencies
    my ($basefreqsshifted,$basefreqsraw) = &filterAndShiftBaseFreqs($kmershift,$kmersum,\%basefreqs,$numseqs);

    { no integer;

    #predict tag sequence for ends with at least 10% in aligned kmers
    my %tags;
    foreach my $sp (keys %$basefreqsshifted) {
        if((100/$numseqs*$kmersum->{$sp}) > 10) {
            $tags{$sp} = &predictTag($basefreqsshifted,$sp);
        }
    }

    #print results
    if($WEB) {
        &addToIdFile($WEBID,'status_tc',100);
        my (%vars,%varsf,@tmpvars,@tmpvarsf);
        foreach my $sp (keys %$basefreqsraw) {
            %vars = %varsf = ();
            foreach my $i (0..($MAX_WEB_LEN-1)) {
                @tmpvars = @tmpvarsf = ();
                foreach my $base ('A','C','G','T','N') {
                    #raw
                    if(exists $basefreqsraw->{$sp}->{$i}->{$base} && $basefreqsraw->{$sp}->{$i}->{$base} > 0) {
                        push(@{$vars{L}},$basefreqsraw->{$sp}->{$i}->{$base});
                        push(@tmpvars,$basefreqsraw->{$sp}->{$i}->{$base}) unless($base eq 'N');
                    } else {
                        push(@{$vars{L}},0);
                        push(@tmpvars,0) unless($base eq 'N');
                    }
                    #filtered
                    if(exists $basefreqsshifted->{$sp}->{$i}->{$base} && $basefreqsshifted->{$sp}->{$i}->{$base} > 0) {
                        push(@{$varsf{L}},$basefreqsshifted->{$sp}->{$i}->{$base});
                        push(@tmpvarsf,$basefreqsshifted->{$sp}->{$i}->{$base}) unless($base eq 'N');
                    } else {
                        push(@{$varsf{L}},0);
                        push(@tmpvarsf,0) unless($base eq 'N');
                    }
                }
                #raw
                push(@{$vars{R}},&range(\@tmpvars));
                push(@{$vars{M}},&median(\@tmpvars));
                #filtered
                push(@{$varsf{R}},&range(\@tmpvarsf));
                push(@{$varsf{M}},&median(\@tmpvarsf));
            }
            #raw
            foreach(keys %vars) {
#                print join("\t",'p'.$_.$sp,join(",",@{$vars{$_}}))."\n";
                &addToIdFile($WEBID,'p'.$_.$sp,join(",",@{$vars{$_}}));
            }
            #filtered
            foreach(keys %varsf) {
#                print join("\t",'pf'.$_.$sp,join(",",@{$varsf{$_}}))."\n";
                &addToIdFile($WEBID,'pf'.$_.$sp,join(",",@{$varsf{$_}}));
            }
#            print join("\t",'tag'.($sp ? 3 : 5),$tags{$sp})."\n";
            &addToIdFile($WEBID,'tag'.$sp,$tags{$sp});
            &addToIdFile($WEBID,'explained'.$sp,sprintf("%.2f",(100/$numseqs*$kmersum->{$sp})));
        }
    } else {
        print '#'.join("\t",qw(Param Tag_Sequence Tag_Length Percentage_Explained))."\n";
        foreach my $sp (keys %tags) {
            $tags{$sp} =~ s/^E+//;
            $tags{$sp} =~ s/E+$//;
            print join("\t",'tag'.($sp ? 3 : 5),$tags{$sp},length($tags{$sp}),sprintf("%.2f",(100/$numseqs*$kmersum->{$sp})))."\n" if($tags{$sp});
        }
    }

    } #no integer
}

#
###########################################################
## APPROX MATCHING STUFF
###########################################################
#
#see the project website or publication for more information on the bit-parallel algorithm
sub getPeq {
    my ($tag,$alphabet) = @_;

    #check if tag sequence is defined
    return 0 unless($tag || $MATRIX);

    #generate Peq for probe string according to Meyers
    #sequence left-most character will be represented in
    #least significant bit (LSB) and right-most character in
    #most significant bit (MSB)
    #reverse sequence to allow left-shift the bit-vector,
    #then OR in value for current character
    #
    #this could overflow if the pattern is longer than 32/64 letters

    my @rev_seq = split(//,scalar reverse $tag);

    #initialize Peq from alphabet
    my @tmp = split(//,$alphabet);
    my %Peq;
    $Peq{$_} = 0 foreach(@tmp);
    #add zero value for non-iupac chars
    $Peq{0} = 0;

    foreach my $char (keys %Peq) {
	foreach my $base (@rev_seq) {
            $Peq{$char} <<= 1;
            if(exists $MATRIX->{$base}->{$char}) { #include wild cards (iupac chars)
                $Peq{$char} |= 1;
            } # else OR with 0 => do nothing
	}
    }

    return \%Peq;
}

sub searchForTag {
    my ($ms,$Peqs,$Pv0s,$initlens,$mms,$minmoves,$seq,$length,$sp,$ns) = @_;

    my ($minpos,$minscore,$subseq,@trimpos,@trimscore,$minposI,$minscoreI,@minscores,@starts,@innertags,$tmplength,@trims,@seqs,%splitstatus);
    $seq     = uc($seq);
    @seqs = ($seq,&rev($seq));
    @trimpos = (0)x2;
    @starts  = (0)x2;
    @trimscore = (0)x2;

    #take substring of initial length from both ends and search for tag sequence
    foreach my $i (@$sp) {
	if($length < $initlens->[$i]) { #too short --> mismatches in any case
	    ($minpos,$minscore) = &matchImperfect($ms->[$i],$Peqs->[$i],$Pv0s->[$i],$seqs[$i]);
            if($STATS) {
                push(@minscores,$minscore);
            } elsif($minscore <= $mms->[$i]) {
                $trimpos[$i] = $minpos;
                $trimscore[$i] = $minscore;
            }
	} else { #long enough to get substring
	    my $again = 1; #check if match directly after match
	    my $start = 0;
	    while($again == 1) {
		$tmplength = $length-$start;
		$tmplength = ($initlens->[$i] < $tmplength ? $initlens->[$i] : $tmplength);
		$subseq = substr($seqs[$i],($start > $ns->[$i] ? $start-$ns->[$i] : $start),$tmplength+($start > $ns->[$i] ? $ns->[$i] : 0));
		#check for match at ends
		($minpos,$minscore) = &matchImperfect($ms->[$i],$Peqs->[$i],$Pv0s->[$i],$subseq);
		#save start for matching in rest
                if($minscore > $mms->[$i]) {
                    $again = 0;
                    $starts[$i] += ($minmoves->[$i] + $minscore);
                    if($STATS) {
                        push(@minscores,$minscore);
                    }
                } else {
                    $start += $minpos-($start > $ns->[$i] ? $ns->[$i] : 0);
		    $starts[$i] = $start;
                    if($STATS) {
                        $again = 0;
                        push(@minscores,$minscore);
                    } else {
                        $trimpos[$i] = $start;
                        $trimscore[$i] = $minscore;
                        if($CONTINUE && $start < $length) {
                            $again = 1;
                        } else {
                            $again = 0;
                        }
                    }
                }
	    }
	}
    }

    #search for inner tag sequences
    my $num5 = 0;
    my $num3 = 0;
    my $num5rev = 0;
    my $num3rev = 0;
    if(scalar(@$sp) == 2 && $length > $initlens->[0] && $starts[0] < $length && $starts[1] < $length) { # 5' AND 3' tag
        my ($minposs5,$minscores5,$minposs3,$minscores3,$minposs5rev,$minscores5rev,$minposs3rev,$minscores3rev);
        if($SPLIT35 || $SPLIT5R) {
            #search for 5' tag matches to get their 3' ends
            $subseq = substr($seqs[0],$starts[0],($length-$starts[0]));
            ($minposs5,$minscores5) = &matchImperfect2($ms->[0],$Peqs->[0],$Pv0s->[0],$subseq,&max($MMFTF35,$MMFTF5R),$ns->[0]);
            $num5 = scalar(@$minposs5);
        }
        if(($SPLIT35 && $num5) || $SPLIT3R) {
            #search for 3' tag matches to get their 5' ends
            $subseq = substr($seqs[1],$starts[1],($length-$starts[1]));
            ($minposs3,$minscores3) = &matchImperfect2($ms->[1],$Peqs->[1],$Pv0s->[1],$subseq,&max($MMFTF35,$MMFTF3R),$ns->[1]);
            $num3 = scalar(@$minposs3);
        }
        if($SPLIT5R && $num5) {
            #search for 5' tag matches to get their 5' ends
            $subseq = substr($seqs[1],0,($length-$starts[0]));
            ($minposs5rev,$minscores5rev) = &matchImperfect2($ms->[0],$Peqs->[2],$Pv0s->[0],$subseq,$MMFTF5R,$ns->[0]);
            @$minposs5rev = reverse @$minposs5rev;
            @$minscores5rev = reverse @$minscores5rev;
            $num5rev = scalar(@$minposs5rev);
        }

        if($SPLIT3R && $num3) {
            #search for 3' tag matches to get their 3' ends
            $subseq = substr($seqs[0],0,($length-$starts[1]));
            ($minposs3rev,$minscores3rev) = &matchImperfect2($ms->[1],$Peqs->[3],$Pv0s->[1],$subseq,$MMFTF3R,$ns->[1]);
            @$minposs3rev = reverse @$minposs3rev;
            @$minscores3rev = reverse @$minscores3rev;
            $num3rev = scalar(@$minposs3rev);
        }

        #combine results
        #3-3-3-3-3-5-5-5-5-5
        my %used;
        if($SPLIT3R && $num3 && $num3rev) {
            foreach my $i (0..$num3rev-1) {
#                next if(exists $used{3}->{$i});
                foreach my $j (0..$num3-1) {
                    next if($i == $j || exists $used{3}->{$j});
                    # (1) 3'-tag after 3'-tag
                    # (2+3) 3'-end to 3'-end of tags is within tag lengths +- the combined number of allowed mismatches
                    # (4) the combined number of mismatches is less or equal to the allowed maximum
                    if($minposs3rev->[$i] > $length-$minposs3->[$j]-$starts[1] &&
                       ($minposs3rev->[$i] - ($length-$minposs3->[$j]-$starts[1])) <= ($ms->[1]*2+$MMFTF3R) &&
                       ($minposs3rev->[$i] - ($length-$minposs3->[$j]-$starts[1])) >= ($ms->[1]*2-$MMFTF3R) &&
                       $minscores3->[$j]+$minscores3rev->[$i] <= $MMFTF3R) {
                        push(@innertags,[33,($length-$minposs3->[$j]-$starts[1]),$minposs3rev->[$i],($minscores3->[$j]+$minscores3rev->[$i])]);
#                        $used{3}->{$i} = 1;
                        $used{3}->{$j} = 1;
                    }
                }
            }
        }
        if($SPLIT35 && $num5 && $num3) {
            foreach my $i (0..$num5-1) {
                foreach my $j (0..$num3-1) {
                    next if(exists $used{'3'}->{$j});
                    # (1) 5'-tag after 3'-tag
                    # (2+3) 3'-end to 5'-end of tags is within tag lengths +- the combined number of allowed mismatches
                    # (4) the combined number of mismatches is less or equal to the allowed maximum
                    if($starts[0]+$minposs5->[$i] > $length-$minposs3->[$j]-$starts[1] &&
                       ($starts[0]+$minposs5->[$i] - ($length-$minposs3->[$j]-$starts[1])) <= ($ms->[0]+$ms->[1]+$MMFTF35) &&
                       ($starts[0]+$minposs5->[$i] - ($length-$minposs3->[$j]-$starts[1])) >= ($ms->[0]+$ms->[1]-$MMFTF35) &&
                       $minscores5->[$i]+$minscores3->[$j] <= $MMFTF35) {
                        push(@innertags,[35,($length-$minposs3->[$j]-$starts[1]),($starts[0]+$minposs5->[$i]),($minscores5->[$i]+$minscores3->[$j])]);
                        $used{5}->{$i} = 1;
#                        $used{3}->{$j} = 1;
                    }
                }
            }
        }
        if($SPLIT5R && $num5 && $num5rev) {
            foreach my $i (0..$num5-1) {
                next if(exists $used{5}->{$i});
                foreach my $j (0..$num5rev-1) {
                    next if($i == $j || exists $used{5}->{$j});
                    # (1) 5'-tag after 5'-tag
                    # (2+3) 5'-end to 5'-end of tags is within tag lengths +- the combined number of allowed mismatches
                    # (4) the combined number of mismatches is less or equal to the allowed maximum
                    if($starts[0]+$minposs5->[$i] > $length-$minposs5rev->[$j] &&
                       ($starts[0]+$minposs5->[$i] - ($length-$minposs5rev->[$j])) <= ($ms->[0]*2+$MMFTF5R) &&
                       ($starts[0]+$minposs5->[$i] - ($length-$minposs5rev->[$j])) >= ($ms->[0]*2-$MMFTF5R) &&
                       $minscores5->[$i]+$minscores5rev->[$j] <= $MMFTF5R) {
                        push(@innertags,[55,($length-$minposs5rev->[$j]),($starts[0]+$minposs5->[$i]),($minscores5->[$i]+$minscores5rev->[$j])]);
#                        $used{5}->{$i} = 1;
                        $used{5}->{$j} = 1;
                    }
                }
            }
        }
    } elsif($SPLIT35 && scalar(@$sp) == 1 && $length > $initlens->[$sp->[0]] && $starts[$sp->[0]] < $length) { # 5' OR 3' tag
        #search for tag matches
        $subseq = substr($seqs[$sp->[0]],$starts[$sp->[0]]);
        my ($minposs,$minscores) = &matchImperfect2($ms->[$sp->[0]],$Peqs->[$sp->[0]],$Pv0s->[$sp->[0]],$subseq,$MMFTF35,$ns->[$sp->[0]]);
        my $tmp = scalar(@$minposs);
        if($tmp) {
            #check for tag matches
            foreach my $i (0..$tmp-1) {
                push(@innertags,[($minposs->[$i]+$starts[$sp->[0]]),$minscores->[$i]]);
            }
        }
    }

    #get [start,length] touple from end and internal tags
    if($TAG5 && $TAG3) { #tag 5 and tag 3
        if(scalar(@innertags)) {
            @innertags = sort {$a->[1] <=> $b->[1]} @innertags;
            my $tmpstart = $trimpos[0];
            my $tmpscore = $trimscore[0];
            while(@innertags) {
                my $tmp = shift(@innertags);
                if($tmp->[1] > $tmpstart) {
                    push(@trims,[$tmpstart,($tmp->[1]-$tmpstart),$tmpscore,$tmp->[3]]);
                    $splitstatus{'split'.($tmp->[0] == 35 ? '' : $tmp->[0])}++;
                }
                $tmpstart = $tmp->[2];
                $tmpscore = $tmp->[3];
            }
            #last one
            push(@trims,[$tmpstart,($length-$tmpstart-$trimpos[1]),$tmpscore,$trimscore[1]]);
        } else {
            push(@trims,[$trimpos[0],($length-$trimpos[0]-$trimpos[1]),$trimscore[0],$trimscore[1]]);
        }
    } elsif($TAG5) { #tag 5 only
        if(scalar(@innertags)) {
            my $tmpstart = $trimpos[0];
            my $tmpscore = $trimscore[0];
            while(@innertags) {
                my $tmp = shift(@innertags);
                push(@trims,[$tmpstart,($tmp->[0]-$tmpstart-$ms->[0]),$tmpscore,0]);
                $tmpstart = $tmp->[0];
                $tmpscore = $tmp->[1];
                $splitstatus{'split'}++;
            }
            #last one
            push(@trims,[$tmpstart,($length-$tmpstart),$tmpscore,$trimscore[1]]);
        } else {
            push(@trims,[$trimpos[0],($length-$trimpos[0]),$trimscore[0],$trimscore[1]]);
        }
    } elsif($TAG3) { #tag 3 only
        if(scalar(@innertags)) {
            my $tmpstart = 0;
            my $tmpscore = 0;
            while(@innertags) {
                my $tmp = pop(@innertags);
                push(@trims,[$tmpstart,($length-$tmp->[0])]);
                $tmpstart = $length-$tmp->[0]+$ms->[1];
                $tmpscore = $tmp->[1];
                $splitstatus{'split'}++;
            }
            #last one
            push(@trims,[$tmpstart,($length-$tmpstart-$trimpos[1]),$tmpscore,$trimscore[1]]);
        } else {
            push(@trims,[0,($length-$trimpos[1]),$trimscore[0],$trimscore[1]]);
        }
    }

    #remove 0-length touples
    my $index = 0;
    while($index < scalar(@trims)) {
        if($trims[$index]->[1] <= 0) {
            splice(@trims,$index,1);
        } else {
            $index++;
        }
    }

    return (\@trims,\@minscores,\%splitstatus);
}

sub matchImperfect {
    my ($m,$Peq,$Pv0,$seq) = @_;
    my $score = $m;
    my $minscore = $score;
    my $minpos = 0;
    my $count = 0;
    my $Pv = $Pv0;
    my $Mv = 0;
    my ($Eq, $Xv, $Xh, $Ph, $Mh);
    foreach(split(//,$seq)) {
	$Eq = $Peq->{$_}||$Peq->{0};
	$Xv = $Eq | $Mv;
	$Xh = ((($Eq & $Pv) + $Pv) ^ $Pv) | $Eq;

	$Ph = $Mv | ~($Xh | $Pv);
	$Mh = $Pv & $Xh;

	if($Ph & 1 << ($m-1)) {
	    $score++;
	} elsif($Mh & 1 << ($m-1)) {
	    $score--;
	}

	$count++;
	if($score < $minscore) {
	    $minscore = $score;
	    $minpos = $count; # location counts begin at zero
	    last if($minscore == 0);
	}

	$Ph <<= 1;
	$Pv = ($Mh << 1) | ~($Xv | $Ph);
	$Mv = $Ph & $Xv;
    }

    return ($minpos,$minscore);
}

sub matchImperfect2 {
    my ($m,$Peq,$Pv0,$seq,$max,$n) = @_;
    my $score = $m;
    my (@minscores,@minposs,$count);
    $count = 0;
    my $Pv = $Pv0;
    my $Mv = 0;
    my $counting = 0;
    my ($Eq, $Xv, $Xh, $Ph, $Mh);
    foreach(split(//,$seq)) {
	$Eq = $Peq->{$_}||$Peq->{0};
	$Xv = $Eq | $Mv;
	$Xh = ((($Eq & $Pv) + $Pv) ^ $Pv) | $Eq;

	$Ph = $Mv | ~($Xh | $Pv);
	$Mh = $Pv & $Xh;

	if($Ph & 1 << ($m-1)) {
	    $score++;
	} elsif($Mh & 1 << ($m-1)) {
	    $score--;
	}

	$counting++;
	unless($score > $max) {
	    if(!$count || ($minposs[$count-1]+$m-2-$n) <= $counting) {
		push(@minscores,$score);
		push(@minposs,$counting); # location counts begin at zero
		$count++;
	    } elsif($count && $minscores[$count-1] > $score) {
		    $minscores[$count-1] = $score;
		    $minposs[$count-1] = $counting;
	    }
	}

	$Ph <<= 1;
	$Pv = ($Mh << 1) | ~($Xv | $Ph);
	$Mv = $Ph & $Xv;
    }

    return (\@minposs,\@minscores);
}

#
###########################################################
## PREDICT TAG
###########################################################
#

#add the kmers and base frequencies to the respective hash
sub getKmersAndBaseFreqs {
    my ($kmers,$basefreqs,$seq,$length,$kmerlen,$taglen) = @_;
    my $i;
    unless($length < $kmerlen) {
        my $kmer5 = substr($seq,0,$kmerlen);
        my $kmer3 = substr($seq,($length-$kmerlen));
        $kmers{0}->{$kmer5}++;
        $kmers{1}->{$kmer3}++;

        if($length < $taglen) {
	    $i = 0;
            my @bases = split(//,$seq);
	    foreach(@bases) {
		$basefreqs->{0}->{$kmer5}->{$i}->{$_}++;
		$i++;
	    }
	    $i = 0;
	    foreach(@bases) {
		$basefreqs->{1}->{$kmer3}->{$i}->{$_}++;
		$i++;
	    }
	} else {
	    $i = 0;
	    foreach(split(//,substr($seq,0,$taglen))) {
		$basefreqs->{0}->{$kmer5}->{$i}->{$_}++;
		$i++;
	    }
	    $i = 0;
	    foreach(split(//,substr($seq,($length-$taglen)))) {
		$basefreqs->{1}->{$kmer3}->{$i}->{$_}++;
		$i++;
	    }
	}
    }
}

#get the frequency of possible tags by shifting kmers by max 2 positions when aligned
sub filterAndShiftKmers {
    my ($kmers) = @_;

    #find most abundant kmer counts
    my $percentone = $numseqs/100;
    my $percentten = $numseqs/10;
    my %most;
    foreach my $sp (keys %$kmers) {
	$most{$sp}->{max} = 0;
	foreach(keys %{$kmers->{$sp}}) {
	    next if($_ eq 'A'x5 || $_ eq 'T'x5 || $_ eq 'C'x5 || $_ eq 'G'x5 || $_ eq 'N'x5);
	    if($kmers->{$sp}->{$_} >= $percentten) {
		$most{$sp}->{ten}++;
	    } elsif($kmers->{$sp}->{$_} >= $percentone) {
		$most{$sp}->{one}++;
	    }
	    #get max count
	    $most{$sp}->{max} = $kmers->{$sp}->{$_} if($most{$sp}->{max} < $kmers->{$sp}->{$_});
	}
    }

    #filter kmers by frequency - threshold of >10% occurrence -> max of 9 different kmers or more if there is non with >10% occurrence
    my $numseqssub = $numseqs/10;
    my $onecount = 2;
    foreach my $sp (keys %$kmers) {
	foreach(keys %{$kmers->{$sp}}) {
	    if(exists $most{$sp}->{ten} && $most{$sp}->{ten} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} <= $percentten);
	    } elsif(exists $most{$sp}->{one} && $most{$sp}->{one} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} <= $percentone);
	    } else {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} != $most{$sp}->{max});
	    }
	}
    }

    my (%kmersum,%kmershift);
    foreach my $sp (sort {$b <=> $a} keys %$kmers) { #5' before 3'
	#if more than one kmer in array, test if shifted by max 2 positions
	my $numkmer = scalar(keys %{$kmers->{$sp}});

	if($numkmer > 1) {
	    my @matrix;
	    my @kmersort = sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}};
	    foreach my $i (0..($numkmer-2)) {
		foreach my $j (($i+1)..($numkmer-1)) {
		    $matrix[$i]->[$j-($i+1)] = &align2seqs($kmersort[$j],$kmersort[$i]);
		}
	    }
	    my $countgood = 0;
	    foreach my $i (0..($numkmer-2)) {
		if(!defined @{$matrix[0]->[$i]}) { #not matching
		    my $count = 0;
		    foreach my $j (1..($numkmer-2)) {
			$count++;
			last if(defined $matrix[$j]->[$i-$j]); #found shift using other kmers
		    }
		    if($count < ($numkmer-1) && $i > 0) {
			my $sum = 0;
			my $sign;
			foreach my $j (0..$count) {
			    next unless(defined $matrix[$j] && defined $matrix[$j]->[$i-1]); #fix: 08/2010
			    if(defined $sign) {
				if(defined $matrix[$j]->[$i-1]->[0] && (($sign < 0 && $matrix[$j]->[$i-1]->[0] < 0) || ($sign > 0 && $matrix[$j]->[$i-1]->[0] > 0))) { #fix: 02/2011
				    $sum += $matrix[$j]->[$i-1]->[0];
				} elsif(($sign < 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] < 0)) || $sign > 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] > 0)) {
				    $sum += $matrix[$j]->[$i-1]->[1];
				}
			    } else {
				$sum += $matrix[$j]->[$i-1]->[0];
			    }
			    $sign = ($matrix[$j]->[$i-1]->[0] < 0 ? -1 : 1) if(defined $matrix[$j]->[$i-1]->[0]); #fix: 02/2011
			}
			$matrix[0]->[$i] = [$sum] if(defined $sign); #fix: 08/2010
		    }
		}
		if(!defined @{$matrix[0]->[$i]}) {
		    last;
		} else {
		    $countgood++;
		}
	    }
	    if($countgood) {
		my $min;
		if($sp == 3) { #3' prime end, 5 for 5' end
		    #find maximum shift to right (pos value)
		    $min = -100;
		    foreach my $i (0..($countgood-1)) {
			$min = ($min > $matrix[0]->[$i]->[0] ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min > 0) {
			$min = -$min;
		    } else {
			$min = 0;
		    }
		} else {
		    #find maximum shift to left (neg value)
		    $min = 100;
		    foreach my $i (0..($countgood-1)) {
			$min = ($min < $matrix[0]->[$i]->[0] ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min < 0) {
			$min = abs($min);
		    } else {
			$min = 0;
		    }
		}
		$kmershift{$sp}->{$kmersort[0]} = $min;
		$kmersum{$sp} += $kmers->{$sp}->{$kmersort[0]};
		foreach my $i (0..($countgood-1)) {
		    $kmershift{$sp}->{$kmersort[$i+1]} = $matrix[0]->[$i]->[0]+$min;
		    $kmersum{$sp} += $kmers->{$sp}->{$kmersort[$i+1]};
		}
	    } else {
		my $tmp = (sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}})[0];
		$kmershift{$sp}->{$tmp} = 0;
		$kmersum{$sp} += $kmers->{$sp}->{$tmp};
	    }
	} elsif($numkmer == 1) {
	    my $tmp = (keys %{$kmers->{$sp}})[0];
	    $kmershift{$sp}->{$tmp} = 0;
	    $kmersum{$sp} += $kmers->{$sp}->{$tmp};
	}
    }

   return (\%kmersum,\%kmershift);
}

sub filterAndShiftBaseFreqs {
    my ($kmershift,$kmersum,$basefreqs,$numseqs) = @_;
    my (%freqs,%freqsraw);
    foreach my $sp (keys %$basefreqs) {
	foreach my $kmer (keys %{$basefreqs->{$sp}}) {
            foreach my $i (keys %{$basefreqs->{$sp}->{$kmer}}) {
                foreach my $base (keys %{$basefreqs->{$sp}->{$kmer}->{$i}}) {
                    $freqsraw{$sp}->{$i}->{$base} += $basefreqs->{$sp}->{$kmer}->{$i}->{$base};
                    if(exists $kmershift->{$sp}->{$kmer}) {
                        $freqs{$sp}->{$i+$kmershift->{$sp}->{$kmer}}->{$base} += $basefreqs->{$sp}->{$kmer}->{$i}->{$base};
                    }
                }
            }
	}
    }
    %$basefreqs = () unless($WEB);

    { no integer;
        #get percentage values
        foreach my $sp (keys %freqs) {
            foreach my $i (keys %{$freqs{$sp}}) {
                foreach my $base (keys %{$freqs{$sp}->{$i}}) {
                    $freqs{$sp}->{$i}->{$base} = $freqs{$sp}->{$i}->{$base} * 100 / $kmersum->{$sp};
                }
            }
        }
        #get percentage values raw
        foreach my $sp (keys %freqsraw) {
            foreach my $i (keys %{$freqsraw{$sp}}) {
                foreach my $base (keys %{$freqsraw{$sp}->{$i}}) {
                    $freqsraw{$sp}->{$i}->{$base} = $freqsraw{$sp}->{$i}->{$base} * 100 / $numseqs;
                }
            }
        }
    }

    return (\%freqs,\%freqsraw);
}

sub predictTag {
    my ($freqs,$sp) = @_;

    #get base frequency values
    my (@vals,@bases);
    my @index = 0..(&max($MAX_TAG_LEN,$MAX_WEB_LEN)-1);
    @index = reverse @index if($sp);
    foreach my $i (@index) {
	my $max = 0; #only keep base with maximum frequency
	my $maxF = 0;
	foreach my $base ('A','C','G','T') {
	    if(exists $freqs->{$sp}->{$i}->{$base} && $freqs->{$sp}->{$i}->{$base} > 0) {
		push(@{$vals[$i]},$freqs->{$sp}->{$i}->{$base});
		if($freqs->{$sp}->{$i}->{$base} > $max) {
		    $bases[$i] = $base;
		    $max = $freqs->{$sp}->{$i}->{$base};
		}
	    } else {
		push(@{$vals[$i]},0);
	    }
	}
    }
    @bases = reverse @bases if($sp);

    my ($tag,$len,$tmp);

    #get range values
    my @ranges;
    foreach(@vals) {
	push(@ranges,&range($_));
    }

    #get median values
    my @medians;
    foreach(@vals) {
	push(@medians,&median($_));
    }

    #check for differnce between range and median
    my $diff = 0;
    foreach my $i (@index[0..9]) {
	$diff += ($ranges[$i+($sp ? -1: 1)]-$medians[$i]);
    }
    $diff /= 10;

    #filtered frequencies
    #find two jumps in range values for WTA tags
    $len = 0;
    foreach my $i (reverse(@index)) {
	if($ranges[$i] <= $medians[$i]+5) {
	    $len++;
	} else {
	    last;
	}
    }
    $len = &max($MAX_TAG_LEN,$MAX_WEB_LEN)-$len-1 unless($sp);

    #get tag sequence
    $tmp = $tag = '';
    $count = 0;
    foreach my $i (@index) {
	if(($i == 0 || $ranges[$i] > $medians[$i]*3 || ($medians[$i] == 0 && $ranges[$i] > 0)) && $tmp ne 'E' && $tmp ne 'N' && scalar(@bases)) { #fix: 02/2011
	    $tag .= shift(@bases);
	} elsif((($sp && $i >= $len) || (!$sp && $i <= $len)) && $tmp ne 'E' && $count < 15) {
	    $tag .= 'N';
	    $tmp = 'N';
	    $count++;
	} else {
	    $tag .= 'E';
	    $tmp = 'E';
	}
    }
    $tag = substr($tag,0,$MAX_TAG_LEN) if(length($tag) > $MAX_TAG_LEN);
    $tag = scalar reverse $tag if($sp);

    #raw frequencies
#    #find two jumps in range values for WTA tags
#    my @changes;
#    foreach my $i (0..&max($MAX_TAG_LEN,$MAX_WEB_LEN)-2) {
#	$changes[$i] = [abs($ranges[$i+1]-$ranges[$i]),$i];
#    }
#    my @sort = sort {$b->[0] <=> $a->[0]} @changes;

#    #get tag sequence
#    $len = $sort[1]->[1];
#    my $tmp = '';
#    my $ncount = 0;
#    foreach my $i (@index) {
#	if((($diff >= 25 && $ranges[$i] >= $RANGE_THRESHOLD) || ($diff < 25 && $ranges[$i] >= $medians[$i])) && $tmp ne 'E' && $tmp ne 'N') {
#	    $tag .= shift(@bases);
#	} elsif((($diff < 25 && $ranges[$i]*2 >= $medians[$i]) || ($diff >= 25 && (($sp && $i > $len) || (!$sp && $i <= $len)))) && $tmp ne 'E' && $ncount < 15) {
#	    $tag .= 'N';
#	    $tmp = 'N';
#	    $ncount++;
#	} else {
#	    $tag .= 'E';
#	    $tmp = 'E';
#	}
#    }
#    $tag = scalar reverse $tag if($sp);

    return $tag;
}


#
################################################################################
## MISC FUNCTIONS
################################################################################
#

sub printError {
    my $msg = shift;
    print STDERR "ERROR: ".$msg.".\n\nTry \'perl tagcleaner.pl -h\' for more information.\nExit program.\n";
    &printLog("ERROR: ".$msg.". Exit program.\n");
    exit(0);
}

sub printWarning {
    my $msg = shift;
    print STDERR "WARNING: ".$msg.".\n";
    &printLog("WARNING: ".$msg.".\n");
}

sub printLog {
    my $msg = shift;
    if(exists $params{log}) {
        my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
        open(FH, ">>", $params{log}) or die "ERROR: Can't open file ".$params{log}.": $! \n";
        flock(FH, LOCK_EX) or die "ERROR: Cannot lock file ".$params{log}.": $! \n";
        print FH "[tagcleaner-".$WHAT."-$VERSION] [$time] $msg\n";
        flock(FH, LOCK_UN) or die "ERROR: cannot unlock ".$params{log}.": $! \n";
        close(FH);
    }
}

sub getPathAndFile {
    my $str = shift;
    my ($path,$file);
    if($str =~ m/^(.+\/)([^\/]+)$/) {
        $path = $1;
        $file = $2;
    } else {
        $file = $str;
    }
    return ($path,$file);
}

sub getLineNumber {
    my $file = shift;
    my $lines = 0;
    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    $lines += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);
    close(FILE);
    return $lines;
}

sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    open(FILE,"perl -p -e 's/\r/\n/g;s/\n\n/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && /^[ACGTNacgtn]+/) {
            $fasta = 2;
        } elsif($qual == 1 && /^\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && /^[ACGTNacgtn]+/) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/);
        }
    }

    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

sub checkTagSequence {
	my $tag = shift;
	$tag = uc($tag);
	if($tag =~ m/^[ACGTNMRWSYKVHDB]+$/) {
		return 1;
	} else {
		return 0;
	}
}

sub rev {
    my $seq = shift;
    $seq = scalar reverse $seq;
    return $seq;
}

sub max {
    my ($a,$b) = @_;
    return ($a < $b ? $b : $a);
}

sub processData {
    my ($seqid,$seq,$qual,$header,$trims,$length) = @_;
    my ($seqn,$qualn,$seqidn,$filter);
    my $idcount = 1;
    my $num = scalar(@$trims);

    #check if read is to be filtered
    $filter = 0;
    if(($NOMATCH == 5 && $num > 1) ||
       ($NOMATCH == 1 && $trims->[0]->[0] == 0) ||
       ($NOMATCH == 2 && ($trims->[$num-1]->[0]+$trims->[$num-1]->[1] == $length)) ||
       ($NOMATCH == 3 && ($trims->[0]->[0] == 0 || ($trims->[$num-1]->[0]+$trims->[$num-1]->[1] == $length))) ||
       ($NOMATCH == 4 && $trims->[0]->[0] == 0 && ($trims->[$num-1]->[0]+$trims->[$num-1]->[1] == $length))) {
        $filter = 1;
    }

    if($MINLEN) {
        my $test = 0;
        foreach my $t (@$trims) {
            if($t->[1] >= $MINLEN) {
                $test = 1;
                last;
            }
        }
        $filter = 1 unless($test);
    }

    if($filter && $FILTERED) {
        $trims = [[0,$length,0,0]];
    } elsif(($filter && !$FILTERED) || (!$filter && $FILTERED)) {
        return 0;
    }

    foreach my $t (@$trims) {
        $seqidn = $seqid;
        if(scalar(@$trims) > 1) {
            #add count to original id
            $seqidn .= '_'.$idcount++;
        }
        #trim if necessary
        $seqn = substr($seq,$t->[0],$t->[1]);
        $qualn = substr($qual,$t->[0],$t->[1]) if(defined $qual && length($qual));

        #change header if necessary
        if($INFO && !$FILTERED) {
            $header = join(" ",$length,length($seqn),$t->[0],($length-$t->[0]-$t->[1]),$t->[2],$t->[3],$num);
        }

        $numseqsout++;
        $numbasesout += length($seqn);

        #set line length
        if($LINE_WIDTH && $params{out_format} != 3) {
            $seqn =~ s/(.{$LINE_WIDTH})/$1\n/g;
            $seqn =~ s/\n$//;
        }
        #write data
        if($params{out_format} == 1) { #FASTA
            print $fhgood '>'.$seqidn.($header ? ' '.$header : '')."\n";
            print $fhgood $seqn."\n";
        } elsif($params{out_format} == 3) { # FASTQ
            print $fhgood '@'.$seqidn.($header ? ' '.$header : '')."\n";
            print $fhgood $seqn."\n";
            print $fhgood '+'.$seqidn.($header ? ' '.$header : '')."\n";
            print $fhgood $qualn."\n";
        } elsif($params{out_format} == 2) { #FASTA+QUAL
            print $fhgood '>'.$seqidn.($header ? ' '.$header : '')."\n";
            print $fhgood $seqn."\n";
            print $fhgood2 '>'.$seqidn.($header ? ' '.$header : '')."\n";
            print $fhgood2 &convertQualArrayToString(&convertQualAsciiToNums($qualn),$LINE_WIDTH)."\n";
        }
    }
}

sub convertQualNumsToAscii {
    my $qual = shift;
    my @ascii;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	push(@ascii,&convertQualNumToAscii($_));
    }
    return \@ascii;
}

sub convertQualNumsToAsciiString {
    my $qual = shift;
    my $ascii;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	$ascii .= &convertQualNumToAscii($_);
    }
    return $ascii;
}

sub convertQualAsciiToNums {
    my $qual = shift;
    my @nums;
    my @ascii = split('',$qual);
    foreach(@ascii) {
	push(@nums,&convertQualAsciiToNum($_));
    }
    return \@nums;
}

sub convertQualNumToAscii {
    my $qual = shift;
    return chr(($_<=93? $_ : 93) + 33);
}

sub convertQualAsciiToNum {
    my $qual = shift;
    return (ord($qual) - 33);
}

sub convertQualArrayToString {
    my ($nums,$linelen) = @_;
    $linelen = 80 unless($linelen);
    my $str;
    my $count = 0;
    foreach my $n (@$nums) {
        $str .= ($n < 10 ? ' '.$n : $n).' ';
        if(++$count > $linelen) {
            $count = 0;
            $str =~ s/\s$//;
            $str .= "\n";
        }
    }
    $str =~ s/[\s\n]$//;
    return $str;
}

sub addCommas {
    my $num = shift;
    return unless(defined $num);
    return $num if($num < 1000);
    $num = scalar reverse $num;
    $num =~ s/(\d{3})/$1\,/g;
    $num =~ s/\,$//;
    $num = scalar reverse $num;
    return $num;
}

sub align2seqs {
    my ($seq1,$seq2) = @_;
    my @shift;

    #get number of shifted positions
    if(substr($seq1,0,4) eq substr($seq2,1,4)) { #shift right by 1
	push(@shift,1);
    } elsif(substr($seq1,0,3) eq substr($seq2,2,3)) { #shift right by 2
	push(@shift,2);
    }
    if(substr($seq1,1,4) eq substr($seq2,0,4)) { #shift left by 1
	push(@shift,-1);
    } elsif(substr($seq1,2,3) eq substr($seq2,0,3)) { #shift left by 2
	push(@shift,-2);
    }

    return \@shift;
}

sub range {
    my $vals = shift;
    my @sort = sort {$a <=> $b} @$vals;
    my $num = scalar(@sort);
    my $range = abs($sort[$num-1]-$sort[0]);
    $range = 95 if($range > 95);
    return $range;
}

sub median {
    my $vals = shift;
    my @sort = sort {$a <=> $b} @$vals;
#    pop(@sort);
    my $num = scalar(@sort);
    if($num%2) {
	return $sort[$num/2];
    } else {
	return (($sort[$num/2]+$sort[$num/2-1])/2);
    }
}

sub addToIdFile {
    my ($idfile,$type,$value) = @_;
    my (@out,@lineargs,$check);
    $check = 0;
    sysopen(FH, $idfile, O_RDWR|O_CREAT) or &printError("Can't open file $idfile: $!");
    flock(FH, LOCK_EX) or &printError("Can't write-lock file $idfile: $!");
    while(<FH>) {
	chomp();
	@lineargs = ();
	@lineargs = split(/\t/);
	next unless(scalar(@lineargs));
	#
	if($lineargs[0] eq $type) {
	    $lineargs[1] = $value;
	    $check = 1;
	}
	push(@out,join("\t",@lineargs)."\n");
    }
    unless($check) {
	push(@out,join("\t",$type,$value||'')."\n");
    }
    seek(FH, 0, 0) or &printError("Can't rewind file $idfile: $!");
    truncate(FH, 0) or &printError("Can't truncate file $idfile: $!");
    foreach(@out) {
	print FH $_;
    }
    close(FH);
    return 1;
}






#