## 1. Pausepred
PausePred predicts ribosomal pause sites using mapped Ribo-seq data (in sorted BAM format). It can be used as a web application (http://pausepred.ucc.ie) or as a stand-alone tool (github: https://github.com/romikasaini/Pausepred_offline). The PausePred tool is written using Perl scripting language (Perl 5, version 18, subversion 2 (v5.18.2)).

## Inputs required to run the stand-alone version:

BAM_file, window_size, foldchange for pause, reference fasta file, comma separated read lengths, window coverage%, upstream_seq and downstream_seq, comma separated offset values

Inputs should be mentioned in the order given above, delimited by space.

Note: The 5' offset should be specified in positive numbers whereas 3' offset should be specified in negative numbers.

## Example Command:
##Genome alignment example:
>perl offline_pausepred.pl example1-genome.bam 1000 20 example-genome.fa 28,29,30 10 50 50 0,0,0

##Transcriptome alignment example:
>perl offline_pausepred.pl example-transcriptome.bam 1000 20 example-transcriptome.fa 28,29,30 10 50 50 0,0,0

or if there is an annotation file available for transcriptome alignments:

>perl offline_pausepred.pl example-transcriptome.bam 1000 20 example-transcriptome.fa 28,29,30 10 50 50 0,0,0 example_annotation.txt

Note: Please enter equal number of comma separated read lengths and offset values

## Modules and Packages required for implementing the stand-alone version:
1. SAMtools(Follow this link to install http://www.htslib.org/download/)


2. Bio::DB::Fasta

Bio::DB::Fasta is a perl module and can be installed by installing BioPerl-1.6.1. Follow this link to install bioperl http://bioperl.org/INSTALL.html
                                                  
or 

Please follow the instructions given at the bottom of this document ("Perl module installation").


## 2. Rfeet
Rfeet generates a graphical view of the ribosome footprint density across a gene/chromosome length. It can be used as a web application (http://pausepred.ucc.ie/rfeet2.html) or as a standalone tool (github: https://github.com/romikasaini/Pausepred_offline).
Rfeet is written using Perl scripting language(Perl 5, version 18, subversion 2 (v5.18.2 )) and R programming language.

## Inputs required to run the stand-alone version:

##### For single input file 

first_bam_file, Fasta sequence, file, gene/transscript/chr:strt-end, comma_separated_offsets, comma_separated_read_lengths

##### Note: second file is optional
first_bam_file, Fasta sequence, file, gene/transscript/chr:strt-end, comma_separated_offsets, comma_separated_read_lengths, second_bam_file,comma_separated_offsets_2, comma_separated_read_lengths_2

Inputs should be provided in the order mentioned above, delimited by space. Additional inputs such as legend names, plot type etc. will be requested through the console once the script is initiated.

## Example command:
##Genome alignment example
>perl offline_rfeet.pl example1-genome.bam example-genome.fa chr:3347-4347 0,0,0 28,29,30

##Transcriptome alignment example
>perl offline_rfeet.pl example-transcriptome.bam example-transcriptome.fa yaaA 0,0,0 28,29,30

note: Second file is optional and multiple genes/chromosome locations can ploted at a time using a comma separated input as shown in the example given below.
>perl offline_rfeet.pl example1-genome.bam example-genome.fa chr:3347-4347,chr:2347-3346,chr:4348-5346 0,0,0 28,29,30 example2-genome.bam 0,0,0 28,29,30

## Modules/Packages required to run standalone version:
1. SAMtools(Follow this link to install http://www.htslib.org/download/)

Following are perl modules and they can be installed by following the "Perl module installation" steps given below.

2. List::Util
3. Bio::DB::Fasta(which can be installed by installing BioPerl-1.6.1)
4. Statistics::R
5. Bedtools (follow link for installation: http://bedtools.readthedocs.io/en/latest/content/installation.html)

To install Statistics::R package you need to install R in your OS first, since Statistics::R need to find R path to work. 

## Perl module installation:

You can use the CPAN shell to install these modules. Open CPAN shell by running following command.

>perl -MCPAN -e shell

and then modules can be installed by running following commands.

>install Bio::DB::Fasta

>install List::Util

>install Statistics::R

Note: If there is a problem with module installation using above commands, please use force installation (for e.g. "force install Bio::DB::Fasta").

## Citation
If you use this tool, please cite: Kumari et al (2018) RNA (http://rnajournal.cshlp.org/content/early/2018/07/26/rna.065235.117.long)
