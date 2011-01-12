#!/usr/bin/perl -w
use strict;
use Carp;
=pod
Thomas Keane, Fri May 02 12:12:20 BST 2008 @508 /Internet Time/
Pathogen Genomics, Sanger Institute, UK

This is a script that takes SSAHA2 cigar format output of reads vs. a reference genome
and creates Artemis output coverage plots from ONLY the uniquely aligned reads

Code Licence:
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
=cut

#check arguments
if( @ARGV != 2 )
{
	print "Usage: cigar2Coverage.pl ssaha_output_file reference_genome_fasta\n";
	exit;
}

my $cigar_file = $ARGV[ 0 ];
my $fasta_file = $ARGV[ 1 ];

if( ! -f $cigar_file || -z $cigar_file ){print "Cannot find ssaha output file or file is empty - please check\n";exit;}
if( ! -f $fasta_file || -z $fasta_file ){print "Cannot find reference fasta file or file is empty - please check\n";exit;}

print "Reading reference fasta file....\n";
my $ref = createHash( $fasta_file );
my %referenceFasta = %{ $ref };

my %coverageArrays_forward;
my %coverageArrays_reverse;
foreach( keys( %referenceFasta ) )
{
	my $key = $_;
	#strip any spaces from the reference headers (ssaha will not include them in the cigar output)
	if( $_ =~ /\s+/ )
	{
		my @s = split( /\s+/, $_ );
		$referenceFasta{ $s[ 0 ] } = $referenceFasta{ $_ };
		delete( $referenceFasta{ $_ } );
		$key = $s[ 0 ];
	}
	
	#create the empty coverage arrays
	my @temp_fwd;
	my @temp_rev;
	my $contigLength = length( ( (split( /\n/, $referenceFasta{ $key } ) )[ 1 ] ) );
	for( my $i = 0; $i < $contigLength; $i ++ )
	{
		push( @temp_fwd, 0 );
		push( @temp_rev, 0 );
	}
	$coverageArrays_forward{ $key } = \@temp_fwd;
	$coverageArrays_reverse{ $key } = \@temp_rev;
}

print "Hashing ssaha cigar output...\n";

if( $cigar_file =~ /\.gz$/ )
{
	open( SSAHA_OUTPUT, "gunzip -c $cigar_file |" ) or die "Cant open ssaha2 output file: $!\n";
}
else
{
	open( SSAHA_OUTPUT, $cigar_file ) or die "Cant open ssaha2 output file: $!\n";
}

my $processed = 0;
my %ssaha_read_matches;
my $ssaha_current_read = "";
my $ssaha_current_read_ref = [];
while( <SSAHA_OUTPUT> )
{
	chomp;
	
	if( $_ =~ /^cigar:/ )
	{
		my $line = $_;
		my @ssaha_splits = split( /\s+/, $line );
		shift( @ssaha_splits );
		$line = join(' ',@ssaha_splits);
		
		if( $ssaha_current_read eq "" )
		{
			$ssaha_current_read = $ssaha_splits[ 0 ];
			$ssaha_current_read_ref = [];
			push( @{ $ssaha_current_read_ref }, $line );
		}
		elsif( $ssaha_current_read ne $ssaha_splits[ 0 ] )
		{
			#at a changover point between two different read matches
			if( @{ $ssaha_current_read_ref } == 1 )
			{
				#read has a single hit so include it
				my @hits = @{ $ssaha_current_read_ref };
				my @splits = split( /\s+/, $hits[ 0 ] );
				
				my $start;
				my $finish;
				if( $splits[ 5 ] > $splits[ 6 ] )
				{
					$start = $splits[ 6 ] - 1;
					$finish = $splits[ 5 ] - 1;
				}
				else
				{
					$start = $splits[ 5 ] - 1;
					$finish = $splits[ 6 ] - 1;
				}
				
				if( $splits[ 3 ] eq '+' )
				{
					for( my $i = $start; $i <= $finish ; $i ++ )
					{
						$coverageArrays_forward{ $splits[ 4 ] }[ $i ] ++;
					}
				}
				else
				{
					for( my $i = $start; $i <= $finish ; $i ++ )
					{
						$coverageArrays_reverse{ $splits[ 4 ] }[ $i ] ++;
					}
				}
			}
			
			$processed ++;
			print "Processed $processed reads...\n" unless $processed % 10000 != 0;
			
			$ssaha_current_read_ref = undef;
			$ssaha_current_read_ref = [];
			
			$ssaha_current_read = $ssaha_splits[ 0 ];
			push( @{ $ssaha_current_read_ref } , $line );
		}
		else
		{
			push( @{ $ssaha_current_read_ref } , $line );
		}
	}
}

#put the entries for the last read in the hash table
if( @{ $ssaha_current_read_ref } == 1 )
{
	#read has a single hit so include it
	my @hits = @{ $ssaha_current_read_ref };
	my @splits = split( /\s+/, $hits[ 0 ] );
	
	my $start;
	my $finish;
	if( $splits[ 5 ] > $splits[ 6 ] )
	{
		$start = $splits[ 6 ] - 1;
		$finish = $splits[ 5 ] - 1;
	}
	else
	{
		$start = $splits[ 5 ] - 1;
		$finish = $splits[ 6 ] - 1;
	}
	
	if( $splits[ 3 ] eq '+' )
	{
		for( my $i = $start; $i <= $finish ; $i ++ )
		{
			$coverageArrays_forward{ $splits[ 4 ] }[ $i ] ++;
		}
	}
	else
	{
		for( my $i = $start; $i <= $finish ; $i ++ )
		{
			$coverageArrays_reverse{ $splits[ 4 ] }[ $i ] ++;
		}
	}
}
close( SSAHA_OUTPUT );

print "Creating plot files....\n";

my $count = 0;
foreach( keys( %coverageArrays_forward ) )
{
	open( OUT, ">".$_.".plot" ) or die "Cannot create output file\n";
	
	my @t = @{ $coverageArrays_forward{ $_ } };
	my @t1 = @{ $coverageArrays_reverse{ $_ } };
	
	for( my $i = 0; $i < @t; $i ++ )
	{
		print OUT $t[ $i ].' '.$t1[ $i ]."\n";
	}
	close( OUT );
	
	$count ++;
}

#a subroutine to take ssaha cigar output and group the records according to reads
sub hash_ssaha_cigar_output
{
	if(  @_ != 1 )
	{
		print "Usage: hash_ssaha_cigar_output ssaha_output_file";
		exit;
	}
	my $ssaha_output = shift;
	
	if( $ssaha_output =~ /\.gz$/ )
	{
		open( SSAHA_OUTPUT, "gunzip -c $ssaha_output |" ) or die "Cant open ssaha2 output file: $!\n";
	}
	else
	{
		open( SSAHA_OUTPUT, $ssaha_output ) or die "Cant open ssaha2 output file: $!\n";
	}
	
	my %ssaha_read_matches;
	my $ssaha_current_read = "";
	my $ssaha_current_read_ref = [];
	while( <SSAHA_OUTPUT> )
	{
		chomp;
		
		if( $_ =~ /^cigar:/ )
		{
			my $line = substr( $_, 7 );
		my @ssaha_splits = split( /\s+/, $line );
		if( $ssaha_current_read eq "" )
		{
			$ssaha_current_read = $ssaha_splits[ 0 ];
			$ssaha_current_read_ref = [];
			push( @{ $ssaha_current_read_ref }, $line );
			#push( @ssaha_current_read_matches, $_ );
		}
		elsif( $ssaha_current_read ne $ssaha_splits[ 0 ] )
		{
			#at a changover point between two different read matches
			my $tmp_ref = \@{ $ssaha_current_read_ref };
			$ssaha_read_matches{ $ssaha_current_read } = $tmp_ref;
			
			$ssaha_current_read_ref = undef;
			$ssaha_current_read_ref = [];
			
			$ssaha_current_read = $ssaha_splits[ 0 ];
			push( @{ $ssaha_current_read_ref } , $line );
		}
		else
		{
			push( @{ $ssaha_current_read_ref } , $line );
		}
		}
	}
	#put the entries for the last read in the hash table
	$ssaha_read_matches{ $ssaha_current_read } = $ssaha_current_read_ref;
	close( SSAHA_OUTPUT );
	
	return \%ssaha_read_matches;
}

#a subroutine to create a hash table from a fastq OR fasta file
#keys are the read/contig names
#entries are the sequences and quality values
sub createHash
{
	croak "Usage: createHash file_name" unless @_ == 1; 
	my $file = shift;
	
	if( ! (-f $file) ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	my $read_name = <READS>; #first readname in first line of file
	chomp( $read_name );
	my $read = "$read_name\n";
	$read_name = substr( $read_name, 1, length( $read_name ) - 1 ); #remove the @/> sign
	my $line;
	while( <READS> ) #read down file until hit start of next read
	{
		chomp;
		$line = $_;
		if( $line =~ /^[\@|>].*/ ) #if hit start of next read
		{
			if( defined $reads_file{ $read_name } )
			{
				#check the sequence entry is same as one that is already defined
				if( $reads_file{ $read_name } eq $read )
				{
					print "WARNING: Read entry with identical sequence already exists in hash: $read_name\n";
				}
				else
				{
					print "WARNING: Read entry with different sequence already exists in hash: $read_name\n";
				}
			}
			else
			{
				#print "$read\n" unless $read_name ne "big961b12.p1k";
				$reads_file{ $read_name } = $read; #enter into the hash
			}
			
			$read_name = substr( $line, 1, length( $line ) - 1 ); #remove the @/> sign - next read name is in the line variable
			$read = "$line\n";
		}
		else
		{
			if( $line =~ /^\+/ )
			{
				$read = $read."\n".$line; #add to info for the current read
			}
			else
			{
				$read = $read.$line; #add to info for the current read
			}
		}
	}
	close( READS );
	
	#enter final value into the hash
	$reads_file{ $read_name } = $read; #enter into the hash
	
	#my @keys = keys %reads_file;
	#my $size = @keys;
	#print "Fastq file has $size reads\n";
	
	return \%reads_file;
}
