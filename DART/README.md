```
#!/usr/bin/perl -w

$file_seqmap_a = "/home/rchavez/Langebio/Ruri/DART/RSRR1_RSSR2/seqmap/SNPs_a_AGPv4_3mm.tsv";
$file_seqmap_b = "/home/rchavez/Langebio/Ruri/DART/RSRR1_RSSR2/seqmap/SNPs_b_AGPv4_5mm.tsv";
$file_dart = "/home/rchavez/Langebio/Ruri/DART/RSRR1_RSSR2/Report-SNP-DMz18-185_RSRR1-RSRR2.csv";
$file_out = "/home/rchavez/Langebio/Ruri/DART/RSRR1_RSSR2/calls_03.csv";
$delim = ",";
$sep = "___";

# get column ids as plate_row_col_num, which are spread across the first six rows
$count = 0;
%col_to_pos = ();
$r = undef;
open ($r, "<", $file_dart) or die $!;
while ($line = <$r>)
{
	chomp $line;
	next if ($line eq "");
	#$line_1 = $line if ($count == 0); # DMz
	#$line_2 = $line if ($count == 1); # DART plate
	$line_3 = $line if ($count == 2); # client plate
	$line_4 = $line if ($count == 3); # row
	$line_5 = $line if ($count == 4); # col
	$line_6 = $line if ($count == 5); # num_sample
	$count++;
	if ($count == 6)
	{
		@a_3 = split($delim, $line_3); # client plate
		@a_4 = split($delim, $line_4); # row
		@a_5 = split($delim, $line_5); # col
		@a_6 = split($delim, $line_6); # num_sample
		for ($i = 16; $i <= $#a_3; $i++)
		{
			undef $id;
			$id = join($sep, ($a_3[$i], $a_4[$i], $a_5[$i], $a_6[$i]));
			$col_to_pos{$i} = $id;
		}
	}
}
close $r;

# Ruri's plate_row_col_num to genotype table, with B73s and PTs modified to B73_1, _2, etc. because we'll store this in a hash later
$file_gen = "/home/rchavez/Langebio/Ruri/DART/RSRR1_RSSR2/RSRR1_RSRR2_Sample_Mapping.csv";
%pos_to_gen = ();
$r = undef;
open ($r, "<", $file_gen) or die $!;
<$r>; # header
while ($line = <$r>)
{
	chomp $line;
	next if ($line eq "");
	$line =~ s/\"//g;
	@a = split(",", $line);
	undef $id;
	$id = join($sep, ($a[1], $a[2], $a[3], $a[4]));
	$pos_to_gen{$id} = $a[5];
}
close $r;

# get alignment coordinates of a-sequences
%hash_a = ();
$r = undef;
open ($r, "<", $file_seqmap_a) or die $!;
<$r>; # header
while ($line = <$r>)
{
	chomp $line;
	next if ($line eq "");
	@a = split("\t", $line); # seqmap output is 0:trans_id 1:trans_coord 2:target_seq 3:probe_id 4:probe_seq 5:num_mismatch 6:strand
	#next if (length($a[3]) > 21);
	undef $chr_raw;
	$chr_raw = $a[0];
	next if ($chr_raw !~ /^\d/);
	$chr_raw =~ /^(.+?)\sdna:/;
	undef $chr;
	$chr = $1;
	undef $coord;
	$coord = $chr.$sep.$a[1]; # use $sep string to concatenate chromosome and start position
	$hash_a{$a[3]}{$coord}{$a[5]} = 1;
}
close $r;

# get alignment coordinates of b-sequences
%hash_b = ();
$r = undef;
open ($r, "<", $file_seqmap_b) or die $!;
<$r>; # header
while ($line = <$r>)
{
	chomp $line;
	next if ($line eq "");
	@a = split("\t", $line); # seqmap output is 0:trans_id 1:trans_coord 2:target_seq 3:probe_id 4:probe_seq 5:num_mismatch 6:strand
	#next if (length($a[3]) > 21);
	undef $chr_raw;
	$chr_raw = $a[0];
	next if ($chr_raw !~ /^\d/);
	$chr_raw =~ /^(.+?)\sdna:/;
	undef $chr;
	$chr = $1;
	undef $coord;
	$coord = $chr.$sep.$a[1]; # use $sep string to concatenate chromosome and start position
	$hash_b{$a[3]}{$coord}{$a[5]} = 1;
}
close $r;

# let's parse the DART table
%calls = %snps = ();
$count = 0;
$r = undef;
open ($r, "<", $file_dart) or die $!;
while ($line = <$r>)
{
	chomp $line;
	next if ($line eq "");
	next if ($line !~ /^\d/); # skip all headers, they start with *
	#print "count=".$count."\n".$line."\n";
	#
	$line_a = $line if ($count == 0);
	$line_b = $line if ($count == 1);
	$count++;
	#
	if ($count == 2) # we have a sequence pair, let's parse it
	{
		$count = 0;
		@a = split($delim, $line_a);
		@b = split($delim, $line_b);
		#print "lasta:".$a[$#a]."\tlastb:".$b[$#b]."\n";
		undef $alleleid_a;
		$alleleid_a = $a[0];
		undef $id_a;
		$alleleid_a =~ /^(\d+?)\|/; # get the numeric id so we're sure we have a pair
		$id_a = $1;
		undef $alleleid_b; # id_b_snp is the full id
		$alleleid_b = $b[0];
		undef $id_b; # id_b is the numeric id
		$alleleid_b =~ /^(\d+?)\|/; # get the numeric id so we're sure we have a pair
		$id_b = $1;
		die "ids not eq!\n".$line_a."\n".$line_b."\n" if ($id_a ne $id_b); # numeric ids of the pair should be the same, otherwise something strange is going on
		#
		# get array of alignment coordinates for each id_a
		@c_a = ();
		@c_a = keys(%{$hash_a{$id_a}});
		if (scalar(@c_a) == 0) # if no coordinates then id_a did not align
		{
			#print "id_a:".$id_a." did not align!\n";
			next;
		}
		elsif (scalar(@c_a) > 1) # if aligned more than once we don't want it
		{
			#print "id_a:".$id_a." aligned more than once!\n";
			next;
		}
		# get the alignment coordinate
		undef $coord;
		$coord = $c_a[0]; # after our ifs @c_a contains only one value
		#
		# get the mismatch value for id_a
		@mms = ();
		@mms = keys(%{$hash_a{$id_a}{$coord}});
		die "mms!\n" if (scalar(@mms) == 0); # there should be at least one mm for each coordinate, if not there's a problem
		die "mms!\n" if (scalar(@mms) > 1); # since we only keep id_a that aligns only once there should only be one mm value
		if ($mms[0] != 0) # aligned with mismatches; we don't want it
		{
			#print "id_a:".$id_a." aligned with more than 0 mismatches!\n";
			next;
		}
		#
		# now id_b;  by definition if we skip non-aligned with same id_a coordinate and aligned more than once then all is left is aligned once at id_a coordinates
		if (!exists($hash_b{$alleleid_b}{$coord})) # id_b did not align at id_a coordinate; probably has > 5 mismatches and failed to align; consider as not in genome
		{
			#print "id_a:".$alleleid_a."\tid_b:".$alleleid_b."\t".$coord." !exists in b!\n";
			next;
		}
		elsif (scalar(keys(%{$hash_b{$alleleid_b}})) > 1) # aligned more than once, we don't want it
		{
			#print "alleleid_b:".$alleleid_b." aligned more than once!\n";
			next;
		}
		#
		# get chr and start of alignment; add snp position
		# snp pos is col 3
		# make id as snp_chr_startsnp
		@chr_start = split($sep, $coord); # coordinates were concatenated with $sep
		undef $chr;
		undef $start;
		$chr = $chr_start[0];
		$start = $chr_start[1];
		undef $start_snp;
		$start_snp = $start + $a[3] + 1; # +1 because snp positions in table are 0-indexed, but our alignments are 1-indexed
		undef $snp;
		$snp = $b[2];
		$snp =~ /\d{1,2}\:(\w)\>(\w)/;
		undef $base_a;
		undef $base_b;
		$base_a = $1;
		$base_b = $2;
		undef $snp_name;
		$snp_name = $alleleid_a.$sep.$chr.$sep.$start_snp;
		#print $snp_name."\n";
		if (exists($snps{$snp_name})) # snps should be unique, otherwise something strange is going on
		{
			print "already parsed ".$snps{$snp_name}." !\n".$line."\n";
			die;
		}
		$snps{$snp_name} = $line;
		#
		# get calls for each individual; samples start at col 16, end at col 145, that is the last @a (or @b) index
		for ($i = 16; $i <= $#a; $i++)
		{
			undef $pos;
			$pos = $col_to_pos{$i}; # get plate_row_col_num according to the column headers
			die "pos error at ".$i."\n" if (!defined($pos));
			undef $genotype;
			$genotype = $pos_to_gen{$pos}; # get the genotype according to Ruri's table
			die "genotype error at ".$pos."\n" if (!defined($genotype));
			die "genotype-snp pair already parsed!\n" if (exists($calls{$snp_name}{$genotype}));
			if ($a[$i] eq "-" || $b[$i] eq "-") # undefined call
			{
				$calls{$snp_name}{$genotype} = "--";
			}
			elsif ($a[$i] == 1 && $b[$i] == 0) # this would be a B73, sensu AGPv4.39
			{
				$calls{$snp_name}{$genotype} = $base_a.$base_a;
			}
			elsif ($a[$i] == 0 && $b[$i] == 1) # SNP
			{
				$calls{$snp_name}{$genotype} = $base_b.$base_b;
			}
			elsif ($a[$i] == 1 && $b[$i] == 1) # B73 and SNP
			{
				$calls{$snp_name}{$genotype} = $base_a.$base_b;
			}
			else # anything else means something strange is going on
			{
				die "wierd calls!\n".$i."\n".$line_a."\n".$line_b."\n"; # fails here with a .tsv input table; don't know why parsing at line ends fails... why??? I suspect hyphen characters
			}
		}
		#die;
	}
}
close $r;
#print join($delim, sort(keys(%{$calls{$snp_name}})))."\n";
#print scalar(keys(%calls))."\n";
#die;

# write calls
$w = undef;
open ($w, ">", $file_out) or die $!;
print $w "alleleid,chr,start".$delim.join($delim, sort(keys(%{$calls{$snp_name}})))."\n"; # the array keys(%{$calls{$snp_name}}) has the individuals, which are our column names
foreach $snp_name (sort(keys(%calls))) # sort is not needed, but... I mean... why not?
{
	@a = split($sep, $snp_name);
	print $w join($delim, @a);
	foreach $genotype (sort(keys(%{$calls{$snp_name}})))
	{
		print $w $delim.$calls{$snp_name}{$genotype};
	}
	print $w "\n";
}
close $w;
```
