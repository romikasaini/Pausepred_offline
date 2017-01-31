##TODO: Pausepred is designed to predict ribosomal pauses using sorted BAM alignment files.
##Author: Romika Kumari
use strict;
use Bio::DB::Fasta;
if(@ARGV<8){print "Input is <BAM_file><window_size><foldchange for pause><reference fasta file><read_length_min><read_length_max><coverage><upstream sequence><downstream sequence>\n";}
my $bam_file=$ARGV[0];
my $window = $ARGV[1];
my $foldchange=$ARGV[2]; 
my $fasta_file=$ARGV[3];
my $readlength_min=$ARGV[4];
my $readlength_max=$ARGV[5];
my $cov=$ARGV[6];
my $US_seq= $ARGV[7];
my $DS_seq= $ARGV[8];

my $db = Bio::DB::Fasta->new($fasta_file);

open (FH, ">$bam_file-pausepred_output.csv") or die "$!";


open F1 ,"samtools view $bam_file |";

print FH "gene_name\tcoordinate_position\tnumber_of_reads_mapped\tPause_score\tcoverage(%)\t50_upstream_seq\t50_downstream_seq(including pause position)\n";
my ($win_start,$win_end,@values,$seq_id, %id_sort_chk, %type,@overlap_values);		## Window start and window end; array to store values of current window; sequence ID (CHR/GEne id)
my @out_file;
my %freq_track;

while(<F1>)
{
	chomp;
	next if(/^(\@)/);
	my @array=  split(/\s+/);
##Conditions to create new windows
	if($array[2]=~/[a-zA-Z]/ && length($array[9]) >=$ARGV[4] && length($array[9]) <= $ARGV[5])
	{
		
		if($seq_id ne $array[2])
			{			## if New gene or chromosome is found
				if($seq_id)
					{
						
						process(\@values,\%type,$win_start,$win_end,$seq_id,\@out_file,\%freq_track) if (scalar @values) >1;
					}
		
				$seq_id = $array[2];($win_start,$win_end)=(1,$window); 
				if(exists $id_sort_chk{$array[2]}){ warn"WARNING!!: Looks like unsorted SAM input. (first seq_id, then coordinates)\n";}else{$id_sort_chk{$array[2]}=1;}
			}
		elsif( $array[3]>$win_end )			## if the current cordinate falls out of the current wondow
			{
				process(\@values,\%type,$win_start,$win_end,$seq_id,\@out_file,\%freq_track)if (scalar @values) >1;
				x:	
				$win_start=($win_end+1)-($window*.5);
				$win_end=($win_start+$window)-1;
				undef @values;		# values vector with only 1 value need to be reset
				push(@values,@overlap_values);
				undef @overlap_values;
				goto x if( $array[3]>$win_start and $array[3]>$win_end ) ## if the current cordinate(45) falls doesnt belong to the new window(28-34); simply skip the window
				
			}
if($array[3]>=$win_start+750 && $array[3]<=$win_end)
	{	
	push @overlap_values, $array[3];	
	}
	push @values, $array[3];			## after defining windows keep tracj of coordinates// ## gene name and coordinates will be kept contast by condition above
	
	}
}

process(\@values,\%type,$win_start,$win_end,$seq_id,\@out_file,\%freq_track) if (scalar @values) >1;

sub process
{
my ($values,$type,$win_start,$win_end,$seq_id,$output_arr,$freq_track)=@_;
$type->{$seq_id}->{$_}++ foreach (@$values);
my @occu = keys %{$type->{$seq_id}};	##this array contais the value that how many read positions are present in particular gene
my @sites = values %{$type->{$seq_id}}; ## this array contains the value how many times one read position is mapped.
my $sum = eval join '+', @sites;	##$sum will tell the total number of reads mapped				
my $average = $sum/$window;
my $coverage=(scalar(@occu)/$window)*100;	
	foreach (sort keys %{$type->{$seq_id}})
	{	
		my $seq_down = $db->seq($seq_id, $_ => $_+$DS_seq);
		my $seq_up = $db->seq($seq_id, $_-$US_seq => $_ -1);
		my $pause_score=$type->{$seq_id}->{$_}/$average;
		if ($type->{$seq_id}->{$_} >= $average*$foldchange && $coverage>=$cov)
		{		
			$output_arr->[scalar @{$output_arr}]= [$seq_id, $_, $type->{$seq_id}->{$_},$pause_score,$coverage,$seq_up,$seq_down];
			$freq_track->{$seq_id.$_}++;
			
		}
     
    		
   
		
	}
undef (@$values);		## clean the array for next window
undef (%$type);			## clean the hash for next window;
}
my(@duplicates, @nondup);

foreach(@out_file)
{

  push @duplicates, $_ if $freq_track{$_->[0].$_->[1]}>1;
  push @nondup, $_ if $freq_track{$_->[0].$_->[1]}==1;
}

my (@uniq,@final_uniq);
my @test = sort {my $avg_pause=(($a->[3]+$b->[3])/2);my $avg_cov=(($a->[4]+$b->[4])/2); push @uniq,[$a->[0],$a->[1],$a->[2],$avg_pause,$avg_cov,$a->[5],$a->[6]] if $b->[1] eq $a->[1]; $b->[1] cmp $a->[1]} @duplicates;

push @final_uniq,@uniq,@nondup;

################################
################################ Zscore calculation
my (@lines,@zscore_values);
my $total_sum;
my %zscore_hash;
my $window;
my @sorted = sort { $b->[4] <=> $a->[4] } @final_uniq;
if(scalar @sorted >=300)
{
$window=300;
}
else
{
$window=scalar(@sorted)
}
my $no_windows=((scalar(@sorted))/$window)*2;
my$count = 1;
my $i = 1;
my $total=0;
my $beginat =1;
my $endat =$window;
while ($count <= $no_windows)
{
my @values;
my $total=0;
for(my $j=(($beginat)-1); $j <= ($endat-1); $j++)
{
push @values,[$sorted[$j]];
$total+=$sorted[$j]->[3];
}
$count = $count+1;
my $mean=($total/$window);
my $sqtotal = 0;
foreach my $k(@values){ foreach my $m (@$k){$sqtotal += ($mean-$m->[3]) ** 2;}}
my $variance=($sqtotal / ($window));

my $std = ($sqtotal / ($window)) ** 0.5;

foreach my $k(@values){foreach my $m (@$k){my $zscore=(($m->[3]-$mean)/$std); push @zscore_values,[$m->[0],$m->[1],$m->[2],$m->[3],$m->[4],$m->[5],$m->[6],$zscore];}}
$beginat = ($endat)-($window*.5);
$endat = ($beginat+$window);
}

my (%h,%uniq_zscore_values);
foreach (@zscore_values)
{
if(! $h{$_->[0].$_->[1]}++) {$uniq_zscore_values{$_->[0].$_->[1]}=$_; }
else{$uniq_zscore_values{$_->[0].$_->[1]}->[7]= ($uniq_zscore_values{$_->[0].$_->[1]}->[7]+$_->[7])/$h{$_->[0].$_->[1]};}
}
undef %h;


foreach (sort keys %uniq_zscore_values)
{
if($_ ne '')
{
print FH join(',',@{$uniq_zscore_values{$_}}),"\n";
}
}
print "Output has been written to file $bam_file-pausepred_output.csv\n";