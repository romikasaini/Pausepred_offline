#!/usr/bin/perl
#Rfeet is designed to create ribosome profiles
##Command: perl offline_rfeet.pl perl plot_inputfile_optional.pl example1_sorted.bam example_sequence.fa chr:3347-4347 example2_sorted.bam
use strict;
use List::Util qw( reduce );
use List::Util qw( min max );
use Bio::DB::Fasta;
use Statistics::R;

if(@ARGV<5){print "Inputs are <BAM alignment file><reference fasta file><gene list><comma separated list of offsets><comma separated list of read lengths>\n"; print "\n"; die "If you have second input file, then above mentioned inputs will be followed by <BAM alignment file2><comma separated offsets for second file><read lengths for second file>\n";}
my $bam_file=$ARGV[0];
my $fasta_file=$ARGV[1];
my $gene_list=$ARGV[2];
my $offset1=$ARGV[3];
my $readlen1=$ARGV[4];
my $bam_file_rnaseq=$ARGV[5];
my $offset2=$ARGV[6];
my $readlen2=$ARGV[7];
print "Please enter plot type(i.e absolute or normalized)\n";
chomp (my $plottype=<STDIN>);
print "Please enter legend to be shown for first file\n";
chomp (my $legend1 =<STDIN>);
print "Please enter legend to be shown for second file. Press enter if only one file is used\n";
chomp (my $legend2 =<STDIN>);
my $cov_plot;
if($plottype eq 'absolute' && defined $bam_file_rnaseq)
{
print "Please write 'coverage' if you want second file to be a coverage plot. Press enter if only one file is used or if you dont want a coverage plot\n";
chomp ($cov_plot=<STDIN>);
}
#print "Please enter offset value.\n";
#chomp (my $offset=<STDIN>);

my ($plot_strt,$plot_end,$frame_loop);
my $gene_name;
my @gene_list_array=split/\,/,$gene_list;

print "Command used :perl offline_rfeet $bam_file\t$fasta_file\t$gene_list\t$bam_file_rnaseq\t$plottype\t$legend1\t$legend2\t$cov_plot\t$offset1\t$offset2\n";
print "Output files will be created shortly...\n";
foreach my $gene_name(@gene_list_array)
{

if($gene_name=~/chr/)
{
my @array=  split/[:-]/,$gene_name;
$gene_name=$array[0];
$plot_strt=$array[1];
$plot_end=$array[2];
$frame_loop=5;
}

my %type;
my %new_filled_hash;
my %typerna;
my %new_filled_hashrna;
my (@names,@namesrna,@scores,@data,@names_rnaseq,@scores_rnaseq,@data_rnaseq,@scores_fwd,@scores_rev,@scores_fwdrna,@scores_revrna);


my $db = Bio::DB::Fasta->new("$fasta_file");
my $seq = $db->seq($gene_name, $plot_strt => $plot_end);
my $rcseq= & reverse($seq);

my $length = length $seq;

if($gene_name!~/chr/)
{
$plot_strt=1;
$plot_end=$length;
$frame_loop=2;
}


system ("samtools index $bam_file");
open PLOTCOD ,"samtools view $bam_file '$gene_name:$plot_strt-$plot_end' | awk '{if(\$4>=$plot_strt && \$4<=$plot_end) print \$0}' |";
if(defined $bam_file_rnaseq)
{
system ("samtools index $bam_file_rnaseq");
open RNASEQPLOT, "samtools view  -b $bam_file_rnaseq '$gene_name:$plot_strt-$plot_end' | genomeCoverageBed -d -split -ibam stdin | awk '{if(\$2>=$plot_strt && \$2<=$plot_end) print \$0}' |" or die "Cant write rnaseq file";

open ABSRNA , "samtools view $bam_file_rnaseq '$gene_name:$plot_strt-$plot_end' | awk '{if(\$4>=$plot_strt && \$4<=$plot_end) print \$0}' |";
}
my ($profile, $frame_plot);
my %unique_orientation ;
my %unique_orientationrna;
my @read_len1=split/,/,$readlen1;
my @off_set1=split/,/,$offset1;


while(<PLOTCOD>)
{
	print PLOTCOD;
	chomp;
	next if(/^(\@)/);
	my @array=  split(/\s+/);
      my $pos;
        for(my $i=0;$i<=scalar @read_len1;$i++) ##loop to add different offsets to different read lengths        
                {

                if ($off_set1[$i] >=0 && $array[1] eq 0 && $read_len1[$i]==length($array[9])){ $pos=$array[3]+$off_set1[$i];}
                if ($off_set1[$i] >=0 && $array[1] eq 16 && $read_len1[$i]==length($array[9])){ $pos=(($array[3]+length($array[9])-1)-$off_set1[$i]);}
                if ($off_set1[$i]<0 && $array[1] eq 0 && $read_len1[$i]==length($array[9])){($pos=($array[3]+length($array[9])-1)+$off_set1[$i]);}
                if ($off_set1[$i]<0 && $array[1] eq 16 && $read_len1[$i]==length($array[9])){$pos=($array[3]-$off_set1[$i]);}
                }


	$unique_orientation{$array[3]}=$array[1];
	      if (defined $pos)
        {
        $type{$pos}++;
        }

}

my $new_filled_hash = fill_hash(\%type,5);
my @read_len2=split/,/,$readlen2;
my @off_set2=split/,/,$offset2;	


while(<ABSRNA>)
{
        print ABSRNA;
        chomp;
        next if(/^(\@)/);
        my @array=  split(/\s+/);
             my $pos2;                                                                                                                                                                           
        for(my $i=0;$i<=scalar @read_len2;$i++) ##loop to add different offsets to different read lengths                                                                                         
                {                                                                                                                                                                                
                                                                                                                                                                                                 
                if ($off_set2[$i] >=0 && $array[1] eq 0 && $read_len2[$i]==length($array[9])){ $pos2=$array[3]+$off_set2[$i];}                                                                   
		if ($off_set2[$i] >=0 && $array[1] eq 16 && $read_len2[$i]==length($array[9])){ $pos2=(($array[3]+length($array[9])-1)-$off_set2[$i]);}
                if ($off_set2[$i]<0 && $array[1] eq 0 && $read_len2[$i]==length($array[9])){($pos2=($array[3]+length($array[9])-1)+$off_set2[$i]);}                                       
                if ($off_set2[$i]<0 && $array[1] eq 16 && $read_len2[$i]==length($array[9])){$pos2=($array[3]-$off_set2[$i]);}                                                                    
                }   
        $unique_orientationrna{$array[3]}=$array[1];
         if (defined $pos2)                                                                                                                                                                      
        {                                                                                                                                                                                     
        $typerna{$pos2}++;                                                                                                                                                                      
        } 

}

my $new_filled_hashrna = fill_hash(\%typerna,5);
       

foreach my$key (sort {$a<=>$b} keys %$new_filled_hash)
        {
 next unless(exists $new_filled_hashrna->{$key} and exists $new_filled_hash->{$key});

                if($unique_orientation{$key} eq 0) {push @scores_fwd, $new_filled_hash->{$key}; push @scores_rev,0; push @names,$key; }
                elsif($unique_orientation{$key} eq 16) {push @scores_rev, $new_filled_hash->{$key};push @scores_fwd,0; push @names,$key; }
                elsif(! exists ($unique_orientation{$key})){push @scores_fwd, $new_filled_hash->{$key};push @scores_rev,0; push @names,$key;}
                if ($unique_orientationrna{$key} eq 0) {push @scores_fwdrna, $new_filled_hashrna->{$key}; push @scores_revrna,0; push @namesrna,$key; }
                elsif($unique_orientationrna{$key} eq 16) {push @scores_revrna, $new_filled_hashrna->{$key};push @scores_fwdrna,0; push @namesrna,$key;}
                elsif(! exists ($unique_orientationrna{$key})){push @scores_fwdrna, $new_filled_hashrna->{$key};push @scores_revrna,0; push @namesrna,$key;}

}


while(<RNASEQPLOT>)
{
	print RNASEQPLOT;
	chomp;
        next if(/^(\@)/);
        my @array=  split(/\s+/);
	push @names_rnaseq, $array[1];
	push @scores_rnaseq,$array[2];
}


my $max_ribo_rev = max @scores_rev;
my $max_ribo_fwd = max @scores_fwd;
my $max_ribo=[ $max_ribo_rev => $max_ribo_fwd ] -> [ $max_ribo_rev <= $max_ribo_fwd ];
my $max_rnaseq = max @scores_rnaseq;
my $highestvalue = [ $max_ribo => $max_rnaseq ] -> [ $max_ribo <= $max_rnaseq ]; 
						my $R = Statistics::R->new();
		$R->set( 'pos_riboseq', \@names);
		$R->set( 'pos_rnaseq', \@names_rnaseq);
		$R->set('scores_rnaseq',\@scores_rnaseq);
		$R->set('scores_fwd',\@scores_fwd);
		$R->set('scores_rev',\@scores_rev);
		$R->set('legend1',$legend1);
                $R->set('legend2',$legend2);
                $R->set('scores_fwdrna',\@scores_fwdrna);
                $R->set('scores_revrna',\@scores_revrna);

sub fill_hash
{

	my $type=shift;
	my $interval = shift || 5;
	my @l=sort {$b<=>$a} keys %$type;
	my ($min, $max)=($l[-1], $l[0]);
	
	my %hash2;		
	for (my $i=$plot_strt; $i<=$plot_end; $i++)
	{		
		$hash2{$i}=0;
		
	}	
	$hash2{$_}=$type->{$_} foreach @l;
	
return \%hash2;
}

close PLOTCOD;
close RNASEQPLOT; 



my @graphs;
my $mygraph_frame;


foreach my $frame_c(0..$frame_loop)##three frames as three loops 
{
	my $plus_minus='+';
	if($frame_c>2){$seq=$rcseq; $plus_minus='-';}
	my(@data_frame,@scores_frame,@names_frame);		##single hash %frame contains data from the current $frame_c
	my %frame;
	my $j=0;

		for(my $i=0; $i < $length; $i++)
			{
					$frame{$i}=0 if !(exists($frame{$i}));	
					next if ($i<$j or $i< $frame_c);
					$j=$i;
					my $frame=substr ($seq, $j, 3);
					if($frame=~m/ATG/){$frame{$j}=1;}
					elsif($frame=~m/TAG|TAA|TGA/){$frame{$j}=2;}
					else{$frame{$j}=0;}
					$j=$i+3;
			}

	my (@k, @c1, @c2);
if($frame_c<3)
	{
        foreach my$key (sort {$a<=>$b} keys %frame)
                                {
               
                                my$value =  $frame{$key};
                                        push( @k, $key);

                                        if($value == 1)
                                        {
                                                push( @c1, 1);
                                                push( @c2, 0);
                                        }
                                        elsif($value == 2)
                                        {
                                                push( @c1, 0);
                                                push( @c2, 2);
                                        }
                                        elsif($value == 0)
                                        {
                                                push( @c1, 0);
                                                push( @c2, 0);

                                        }
                                }
	}
else
	{
        foreach my$key (sort {$b<=>$a} keys %frame)
                                {
                           
                                my$value =  $frame{$key};
                                        push( @k, $key);

                                        if($value == 1)
                                        {
                                                push( @c1, 1);
                                                push( @c2, 0);
                                        }
                                        elsif($value == 2)
                                        {
                                                push( @c1, 0);
                                                push( @c2, 2);
                                        }
                                        elsif($value == 0)
                                        {
                                                push( @c1, 0);
                                                push( @c2, 0);

                                        }
                                }
	}

if($frame_c==0){$R->set( 'f0_pos', \@k); $R->set('f0_start',\@c1); $R->set('f0_end',\@c2); }
if($frame_c==1){$R->set( 'f1_pos', \@k); $R->set('f1_start',\@c1); $R->set('f1_end',\@c2); }
if($frame_c==2){$R->set( 'f2_pos', \@k); $R->set('f2_start',\@c1); $R->set('f2_end',\@c2); }							
if($frame_c==3){$R->set( 'f3_pos', \@k); $R->set('f3_start',\@c1); $R->set('f3_end',\@c2); }
if($frame_c==4){$R->set( 'f4_pos', \@k); $R->set('f4_start',\@c1); $R->set('f4_end',\@c2); }
if($frame_c==5){$R->set( 'f5_pos', \@k); $R->set('f5_start',\@c1); $R->set('f5_end',\@c2); }
							
								

}
############
############ Absolute plot
if($plottype eq 'absolute')
{
if($frame_loop==2)
{
$R->run( q`head(f0_pos)` );
$R->run( qq`png("$gene_name-$plot_strt-$plot_end.png",res = 500, pointsize =4, width = 2000, height = 2000)` );
		$R->set('xmin',$plot_strt);
		$R->set('xmax',$plot_end);
		$R->set('name',$gene_name);
		$R->run(q`ymax1<-max(scores_rnaseq,scores_fwd)`);
                $R->run(q`ymax2<-max(scores_fwd,scores_fwdrna)`);
		$R->run(q`scores_fwdrna[which(scores_fwdrna==0)]=NA`);
		$R->run(q`scores_fwd[which(scores_fwd==0)]=NA`);
		$R->run( q`layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6), nrow = 6, ncol = 2, byrow = TRUE),heights = c(14,1,1,1,1))`);
		$R->run( q`par(mar=c(8,8,8,8),mgp=c(5,1,0))`);
		if($cov_plot eq 'coverage')
		{
		$R->run( q`plot(pos_rnaseq,scores_rnaseq,xlim=c(xmin,xmax), ylim=c(0,ymax1),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='gray',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
			if($bam_file_rnaseq=~/bam/i)
                                {
				$R->run( q`legend("topright", c(legend1,legend2), cex=2, col=c("red","gray"), lwd=c(2,2,2),bty="n")`);
				}
			else
				{
				$R->run( q`legend("topright", c(legend1), cex=2, col=c("red"), lwd=c(2,2,2),bty="n")`);
				}
		$R->run( q`par(new=TRUE)` );
                $R->run( q`plot(pos_riboseq,scores_fwd,xlim=c(xmin,xmax),ylim=c(0,ymax1),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=1, green=0, blue=0, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );

		}
		else
		{
		$R->run( q`plot(pos_riboseq,scores_fwdrna,xlim=c(xmin,xmax),ylim=c(0,ymax2),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='black',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(new=TRUE)` );
			if($bam_file_rnaseq=~/bam/i)
				{
				$R->run( q`legend("topright", c(legend1,legend2), cex=2, col=c("red","black"), lwd=c(2,2,2),bty="n")`);
				}
			else
				{
				$R->run( q`legend("topright", c(legend1), cex=2, col=c("red"), lwd=c(2,2,2),bty="n")`);
				}
		$R->run( q`par(new=TRUE)` );
                $R->run( q`plot(pos_riboseq,scores_fwd,xlim=c(xmin,xmax),ylim=c(0,ymax2),type='h',xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=1, green=0, blue=0, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		}		
		$R->run( q`par(mar=c(0,8,1,4))`);
		$R->run(q`f0_start[which(f0_start==0)]=NA`);
                $R->run(q`f1_start[which(f1_start==0)]=NA`);
		$R->run(q`f2_start[which(f2_start==0)]=NA`);
		$R->run(q`f0_end[which(f0_end==0)]=NA`);
                $R->run(q`f1_end[which(f1_end==0)]=NA`);
                $R->run(q`f2_end[which(f2_end==0)]=NA`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab =NA, ylab =NA,col='white',cex.lab=2,axes=FALSE)` );
                $R->run( q`legend("topleft", c("Start","Stop"), cex=1.5, col=c("green","red"), lwd=c(2,2),horiz=TRUE,bty="n")`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab ="+1",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f0_pos,f0_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = "+1",col="red",cex.lab=2)` );
		$R->run( q`plot(f1_pos,f1_start,ylim=c(0,2),type='h',xaxt='n',yaxt='n',xlab = NA, ylab = "+2",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f1_pos,f1_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab ="+2",col="red",cex.lab=2)` );
		$R->run( q`plot(f2_pos,f2_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab ="Open Reading Frames", ylab ="+3",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f2_pos,f2_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab ="Open Reading Frames", ylab ="+3",col="red",cex.lab=2)` );
		$R->run( q`mtext(side=1, text="Open Reading Frames", line=3,cex=1.5)`);
		$R->run( q`dev.off()` );
		$R->run(q`scores_fwdrna[is.na(scores_fwdrna)] <- 0`);
                $R->run(q`scores_fwd[is.na(scores_fwd)] <- 0`);
		 if(defined $bam_file_rnaseq)
                {
		$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd,reads_mapped_file2=scores_fwdrna), file=paste(name,"-",xmin,"-",xmax,"-",".csv",sep=""))`);  
		}
		else{$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd), file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);}
}

if($frame_loop==5)
{
$R->run( q`head(f0_pos)` );
$R->run( qq`png("$gene_name-$plot_strt-$plot_end.png",res = 500, pointsize =4, width = 2000, height = 2000)` );
		$R->set('xmin',$plot_strt);
		$R->set('xmax',$plot_end);
		 $R->set('name',$gene_name);
		$R->run(q`ymax<-max(scores_rnaseq,scores_fwdrna,scores_revrna,scores_fwd,scores_rev)`);
		$R->run(q`scores_fwdrna[which(scores_fwdrna==0)]=NA`);
		$R->run(q`scores_fwd[which(scores_fwd==0)]=NA`);
		$R->run(q`scores_revrna[which(scores_revrna==0)]=NA`);
		$R->run(q`scores_rev[which(scores_rev==0)]=NA`);
		$R->run( q`layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9), nrow = 9, ncol = 2, byrow = TRUE),heights = c(14,1,1,1,1,1,1,1))`);
		$R->run( q`par(mar=c(8,8,4,4),mgp=c(5,1,0))`);
		if($cov_plot eq 'coverage')
		{
		$R->run( q`plot(pos_rnaseq,scores_rnaseq,xlim=c(xmin,xmax), ylim=c(0,ymax),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='gray',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
			if($bam_file_rnaseq=~/bam/i)
				{	
				$R->run( q`legend("topright", c(paste(legend1, " (", "foward",")", sep=""),paste(legend1, " (", "reverse",")", sep=""),legend2), cex=2, col=c("red","blue","gray"), lwd=c(2,2),bty="n")`);
				}
			else
				{
				$R->run( q`legend("topright", c(paste(legend1, " (", "foward",")", sep=""),paste(legend1, " (", "reverse",")", sep="")), cex=2, col=c("red","blue"), lwd=c(2,2),bty="n")`);
				}
		}
		else
		{
		$R->run( q`plot(pos_riboseq,scores_fwdrna,xlim=c(xmin,xmax), ylim=c(0,ymax),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='black',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_revrna,xlim=c(xmin,xmax), ylim=c(0,ymax),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='chartreuse4',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
				if($bam_file_rnaseq=~/bam/i)
                                {
				$R->run( q`legend("topright", c(paste(legend1, " (", "foward",")", sep=""),paste(legend1, " (", "reverse",")", sep=""),paste(legend2, " (", "foward",")", sep=""),paste(legend2, " (", "reverse",")", sep="")), cex=2, col=c("red","blue","black","chartreuse4"), lwd=c(2,2),bty="n")`);
				}
				else
				{
				$R->run( q`legend("topright", c(paste(legend1, " (", "foward",")", sep=""),paste(legend1, " (", "reverse",")",sep="")), cex=2, col=c("red","blue"), lwd=c(2,2),bty="n")`);
				}
		}
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_fwd,xlim=c(xmin,xmax),ylim=c(0,ymax),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=1, green=0, blue=0, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_rev,xlim=c(xmin,xmax),ylim=c(0,ymax),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=0, green=0, blue=1, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(mar=c(0,8,1,4))`);
		$R->run(q`f0_start[which(f0_start==0)]=NA`);
                $R->run(q`f1_start[which(f1_start==0)]=NA`);
                $R->run(q`f2_start[which(f2_start==0)]=NA`);
                $R->run(q`f0_end[which(f0_end==0)]=NA`);
                $R->run(q`f1_end[which(f1_end==0)]=NA`);
                $R->run(q`f2_end[which(f2_end==0)]=NA`);
                $R->run(q`f3_start[which(f3_start==0)]=NA`);
                $R->run(q`f4_start[which(f4_start==0)]=NA`);
                $R->run(q`f5_start[which(f5_start==0)]=NA`);
                $R->run(q`f3_end[which(f3_end==0)]=NA`);
                $R->run(q`f4_end[which(f4_end==0)]=NA`);
                $R->run(q`f5_end[which(f5_end==0)]=NA`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab =NA, ylab =NA,col='white',cex.lab=2,axes=FALSE)` );
		$R->run( q`legend("topleft", c("Start","Stop"), cex=1.5, col=c("green","red"), lwd=c(2,2),horiz=TRUE,bty="n")`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),xlim=(range(f0_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="+1",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f0_pos,f0_end,ylim=c(0,2),xlim=(range(f0_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="+1",col="red",cex.lab=2)` );
		$R->run( q`plot(f1_pos,f1_start,ylim=c(0,2),xlim=(range(f1_pos)),type='h',xaxt='n',yaxt='n',xlab =NA, ylab ="+2",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f1_pos,f1_end,ylim=c(0,2),xlim=(range(f1_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="+2",col="red",cex.lab=2)` );
		$R->run( q`plot(f2_pos,f2_start,ylim=c(0,2),xlim=(range(f2_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="+3",col='green',cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f2_pos,f2_end,ylim=c(0,2),xlim=(range(f2_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="+3",col="red",cex.lab=2)` );
		$R->run( q`plot(f3_pos,f3_start,ylim=c(0,2),xlim=rev(range(f3_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="-1",col="green",cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f3_pos,f3_end,ylim=c(0,2),xlim=rev(range(f3_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="-1",col="red",cex.lab=2)` );
		$R->run( q`plot(f4_pos,f4_start,ylim=c(0,2),xlim=rev(range(f4_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="-2",col="green",cex.lab=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f4_pos,f4_end,ylim=c(0,2),xlim=rev(range(f4_pos)),type='h', xaxt='n',yaxt='n',xlab =NA, ylab ="-2",col="red",cex.lab=2)` );
		$R->run( q`plot(f5_pos,f5_start,ylim=c(0,2),xlim=rev(range(f5_pos)),type='h', xaxt='n',yaxt='n',xlab ="Open Reading Frames", ylab ="-3",col="green",cex.lab=2)` );
		$R->run( q`par(new=TRUE)`);
		$R->run( q`plot(f5_pos,f5_end,ylim=c(0,2),xlim=rev(range(f0_pos)),type='h', xaxt='n',yaxt='n',xlab ="Open Reading Frames", ylab ="-3",col="red",cex.lab=2)` );
		$R->run( q`mtext(side=1, text="Open Reading Frames", line=3,cex=1.5)`);
		$R->run( q`dev.off()` ); 
		$R->run(q`scores_fwdrna[is.na(scores_fwdrna)] <- 0`);
                $R->run(q`scores_fwd[is.na(scores_fwd)] <- 0`);
                $R->run(q`scores_revrna[is.na(scores_revrna)] <- 0`);
                $R->run(q`scores_rev[is.na(scores_rev)] <- 0`);
 		if(defined $bam_file_rnaseq)
                {
		$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd+scores_rev,reads_mapped_file2=scores_fwdrna+scores_revrna),file=paste(name,"-",xmin,"-",xmax,"-",".csv",sep=""))`);
		}
		else{$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd+scores_rev),file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);}

}
}
######################################
###################################### Normalized plot
if($plottype eq 'normalized')
{
if($frame_loop==2)
{
$R->run( q`head(f0_pos)` );
$R->run( qq`png("$gene_name-$plot_strt-$plot_end.png",res = 500, pointsize =4, width = 2000, height = 2000)` );
		$R->set('xmin',$plot_strt);
		$R->set('xmax',$plot_end);
		$R->set('ymax',$highestvalue);
		$R->set('name',$gene_name);
		$R->run( q`layout(matrix(c(1,1,2,2,3,3,4,4,5,5), nrow = 5, ncol = 2, byrow = TRUE),heights = c(14,1,1,1,1))`);
		$R->run(q`mini<-min(sum(scores_fwd),sum(scores_fwdrna))`);
		$R->run(q`file1_total<-sum(scores_fwd)`);
                $R->run(q`file2_total<-sum(scores_fwdrna)`);
                if(defined $bam_file_rnaseq)
                {
                $R->run(q`ymax_local<-max((scores_fwdrna/file2_total),(scores_fwd/file1_total))`);
                }
                else{$R->run(q`ymax_local<-max((scores_fwd/file1_total))`);}
		$R->run(q`scores_fwdrna[which(scores_fwdrna==0)]=NA`);
		$R->run(q`scores_fwd[which(scores_fwd==0)]=NA`);
		$R->run( q`par(mar=c(8,8,4,4),mgp=c(5,1,0))`);
		$R->run( q`plot(pos_riboseq,(scores_fwdrna/file2_total),xlim=c(xmin,xmax), ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='black',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`legend("topright", c("riboseq","rnaseq"), cex=2, col=c("red","black"), lwd=c(2,2,2),bty="n")`);
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,(scores_fwd/file1_total),xlim=c(xmin,xmax),ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=1, green=0, blue=0, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(mar=c(0,8,1,4))`);
		$R->run(q`f0_start[which(f0_start==0)]=NA`);
                $R->run(q`f1_start[which(f1_start==0)]=NA`);
                $R->run(q`f2_start[which(f2_start==0)]=NA`);
                $R->run(q`f0_end[which(f0_end==0)]=NA`);
                $R->run(q`f1_end[which(f1_end==0)]=NA`);
                $R->run(q`f2_end[which(f2_end==0)]=NA`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f0_pos,f0_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f1_pos,f1_start,ylim=c(0,2),type='h',xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f1_pos,f1_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f2_pos,f2_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f2_pos,f2_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`dev.off()` );
		         $R->run(q`scores_fwdrna[is.na(scores_fwdrna)] <- 0`);
                $R->run(q`scores_fwd[is.na(scores_fwd)] <- 0`);  
		if(defined $bam_file_rnaseq)
                {
		$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd,reads_mapped_file2=scores_fwdrna,normalized_file1=(scores_fwd/file1_total),normalized_file2=(scores_fwdrna/file2_total)),file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);
		}
		else{$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd,normalized_file1=(scores_fwd/file1_total)),file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);}
}

if($frame_loop==5)
{
$R->run( q`head(f0_pos)` );
$R->run( qq`png("$gene_name-$plot_strt-$plot_end.png",res = 500, pointsize =4, width = 2000, height = 2000)` );
		$R->set('xmin',$plot_strt);
		$R->set('xmax',$plot_end);
		$R->set('ymax',$highestvalue);
		$R->set('name',$gene_name);
		$R->run(q`file1_total<-sum(scores_fwd+scores_rev)`);
                $R->run(q`file2_total<-sum(scores_fwdrna+scores_revrna)`);
                if(defined $bam_file_rnaseq)
                {
                $R->run(q`ymax_local<-max((scores_fwdrna/file2_total),(scores_fwd/file1_total),(scores_revrna/file2_total),(scores_rev/file1_total))`);
                }
                else{$R->run(q`ymax_local<-max((scores_fwd/file1_total),(scores_rev/file1_total))`);}
		$R->run(q`scores_fwdrna[which(scores_fwdrna==0)]=NA`);
		$R->run(q`scores_fwd[(scores_fwd==0)]=NA`);
		$R->run(q`scores_revrna[which(scores_revrna==0)]=NA`);
		$R->run(q`scores_rev[which(scores_rev==0)]=NA`);
		$R->run( q`layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8), nrow = 8, ncol = 2, byrow = TRUE),heights = c(14,1,1,1,1,1,1,1))`);
		$R->run( q`par(mar=c(8,8,4,4),mgp=c(5,1,0))`);
		$R->run( q`plot(pos_riboseq,scores_fwdrna/file2_total,xlim=c(xmin,xmax), ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='black',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_revrna/file2_total,xlim=c(xmin,xmax), ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col='chartreuse4',bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`legend("topright", c(paste("riboseq", " (", "foward",")", sep=""),paste("riboseq", " (", "reverse",")", sep=""),paste("rnaseq", " (", "foward",")", sep=""),paste("rnaseq", " (", "reverse",")", sep="")), cex=2, col=c("red","blue","black","chartreuse4"), lwd=c(2,2),bty="n")`);
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_fwd/file1_total,xlim=c(xmin,xmax),ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=1, green=0, blue=0, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(pos_riboseq,scores_rev/file1_total,xlim=c(xmin,xmax),ylim=c(0,ymax_local),type='h', xlab="coordinate positions",ylab="no.of reads mapped",col=rgb(red=0, green=0, blue=1, alpha=0.5),bty='l',lwd=1,cex.lab=2,cex.axis=2,cex.main=2)` );
		$R->run( q`par(mar=c(0,8,1,4))`);
		$R->run(q`f0_start[which(f0_start==0)]=NA`);
                $R->run(q`f1_start[which(f1_start==0)]=NA`);
                $R->run(q`f2_start[which(f2_start==0)]=NA`);
                $R->run(q`f0_end[which(f0_end==0)]=NA`);
                $R->run(q`f1_end[which(f1_end==0)]=NA`);
                $R->run(q`f2_end[which(f2_end==0)]=NA`);
		$R->run(q`f3_start[which(f3_start==0)]=NA`);
                $R->run(q`f4_start[which(f4_start==0)]=NA`);
                $R->run(q`f5_start[which(f5_start==0)]=NA`);
                $R->run(q`f3_end[which(f3_end==0)]=NA`);
                $R->run(q`f4_end[which(f4_end==0)]=NA`);
                $R->run(q`f5_end[which(f5_end==0)]=NA`);
		$R->run( q`plot(f0_pos,f0_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f0_pos,f0_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f1_pos,f1_start,ylim=c(0,2),type='h',xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f1_pos,f1_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f2_pos,f2_start,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col='green')` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f2_pos,f2_end,ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f3_pos,f3_start,xlim=rev(range(f3_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="green")` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f3_pos,f3_end,xlim=rev(range(f3_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f4_pos,f4_start,xlim=rev(range(f4_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="green")` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f4_pos,f4_end,xlim=rev(range(f4_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`plot(f5_pos,f5_start,xlim=rev(range(f5_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="green")` );
		$R->run( q`par(new=TRUE)` );
		$R->run( q`plot(f5_pos,f5_end,xlim=rev(range(f5_pos)),ylim=c(0,2),type='h', xaxt='n',yaxt='n',xlab = NA, ylab = NA,col="red")` );
		$R->run( q`dev.off()` ); 
 		$R->run(q`scores_fwdrna[is.na(scores_fwdrna)] <- 0`);
                $R->run(q`scores_fwd[is.na(scores_fwd)] <- 0`);
                $R->run(q`scores_revrna[is.na(scores_revrna)] <- 0`);
                $R->run(q`scores_rev[is.na(scores_rev)] <- 0`);
		if(defined $bam_file_rnaseq)
                {
		$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd+scores_rev,reads_mapped_file2=scores_fwdrna+scores_revrna,normalized_file1=((scores_fwd+scores_rev)/file1_total),normalized_file2=((scores_fwdrna+scores_revrna)/file2_total)),file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);
		}
		else{$R->run(q`write.csv(data.frame(coordinate_position=pos_riboseq, reads_mapped_file1=scores_fwd+scores_rev,normalized_file1=((scores_fwd+scores_rev)/file1_total)),file=paste(name,"-",xmin,"-",xmax,".csv",sep=""))`);}
		
}
}

sub reverse {
        my $seq = shift;
        my $revcomp=$seq;
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return reverse($revcomp);
}
print "Output plot saved to:$gene_name-$plot_strt-$plot_end.png\n";
print "Output csv saved to:$gene_name-$plot_strt-$plot_end.csv\n";
 }
