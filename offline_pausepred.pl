##TODO: Pausepred is designed to predict ribosomal pauses using sorted BAM alignment files.
##Author: Romika Kumari
use strict;
use Bio::DB::Fasta;
if(@ARGV<9){die "\nInputs are <BAM_file><window_size><foldchange for pause><reference fasta file><comma separated read lengths><window coverage><upstream sequence><downstream sequence><comma separated offset_value>\n\n";}
my $bam_file=$ARGV[0];
my $window = $ARGV[1];
my $foldchange=$ARGV[2]; 
my $fasta_file=$ARGV[3];
my $readlength=$ARGV[4];
my $cov=$ARGV[5];
my $US_seq= $ARGV[6];
my $DS_seq= $ARGV[7];
my $offset= $ARGV[8];
print "Please enter output file name\n";
chomp (my $outfile=<STDIN>);
my $db = Bio::DB::Fasta->new($fasta_file);

open (FH, ">$outfile") or die "$!";


open F1 ,"samtools view $bam_file |";
print FH "Command used: perl offline_pausepred.pl $bam_file $window $foldchange $fasta_file $readlength $cov $US_seq $DS_seq $offset \n";
print FH "gene_name,coordinate_position,number_of_reads_mapped,Pause_score,coverage(%),50_upstream_seq,50_downstream_seq(including_pause_position),Z-score\n";
my ($win_start,$win_end,@values,$seq_id, %id_sort_chk, %type,@overlap_values);		## Window start and window end; array to store values of current window; sequence ID (CHR/GEne id)
my @out_file;
my %freq_track;


while(<F1>)
{
	chomp;
	next if(/^(\@)/);
	my @array=  split(/\s+/);
	my @read_len=split/,/,$readlength;
my @off_set=split/,/,$offset;

	

my $pos;
##Conditions to create new windows
	if( $array[2]=~/[a-zA-Z]/ && length($array[9]) >=@read_len[0] && length($array[9]) <= @read_len[-1])
	{
		##to add offsets to the positions
		for(my $i=0;$i<=scalar @read_len;$i++){
   if ($off_set[$i] >=0 && $array[1] eq 0 && $read_len[$i]==length($array[9]))     {$pos=($array[3]+$off_set[$i]);}  ##5' offset for forward strand
                if ($off_set[$i] >=0 && $array[1] eq 16 && $read_len[$i]==length($array[9])){$pos=(($array[3]+length($array[9])-1)-$off_set[$i]);} ##5'offset for reverse strand
                if ($off_set[$i]<0 && $array[1] eq 0 && $read_len[$i]==length($array[9])) {$pos=(($array[3]+length($array[9])-1)+$off_set[$i]);} ##3'offset for forward strand
                if ($off_set[$i]<0 && $array[1] eq 16 && $read_len[$i]==length($array[9])){$pos=($array[3]-$off_set[$i]);} ##3' offset for reverse strand
                                }	
		if($seq_id ne $array[2])
			{			## if New gene or chromosome is found
				if($seq_id)
					{
						
						process(\@values,\%type,$win_start,$win_end,$seq_id,\@out_file,\%freq_track) if (scalar @values) >1;
						undef @overlap_values;
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
if($array[3]>=($win_start+$win_end)/2 && $array[3]<=$win_end)
	{	
	push @overlap_values,$pos;
	}
	push @values,$pos;
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
			my $pause_score2=sprintf("%.2f",$pause_score);
			$output_arr->[scalar @{$output_arr}]= [$seq_id, $_, $type->{$seq_id}->{$_},$pause_score2,$coverage,$seq_up,$seq_down];
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
if(scalar @sorted ==0)
{
print "No pauses predicted for these parameters. Please change parameters and try again. Thanks.\n";
exit;
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
print "Output has been written to file $outfile\n";
