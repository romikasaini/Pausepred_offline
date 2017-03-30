## 1. Pausepred
PausePred predicts ribosomal pause sites using mapped Ribo-seq data (in sorted BAM format). It can be used as a web application (http://pausepred.ucc.ie) or as a stand-alone tool (github: https://github.com/romikasaini/Pausepred_offline). The PausePred tool is written using Perl scripting language (Perl 5, version 18, subversion 2 (v5.18.2)).

## Inputs required to run the stand-alone version:

BAM_file, window_size, foldchange for pause, reference fasta file, read_length_min, read_length_max, coverage%, upstream_seq and downstream_seq, offset value

Inputs should be mentioned in the order given above, delimited by space.

## Example Command:
>perl offline_pausepred.pl example1_sorted.bam 1000 10 example_sequence.fa 20 35 10 50 50 12


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
first_bam_file, Fasta sequence, file, gene/transscript/chr:strt-end, second_bam_file

##### Note: second file is optional

Inputs should be provided in the order mentioned above, delimited by space. Additional inputs such as legend names, plot type etc. will be requested through the console once the script is initiated.

## Example command:
>perl offline_rfeet.pl example1_sorted.bam example_sequence.fa chr:3347-4347 example2_sorted.bam

note: Second file is optional and multiple genes/chromosome locations can ploted at a time using a comma separated input as shown in the example given below.
>perl offline_rfeet.pl example1_sorted.bam example_sequence.fa chr:3347-4347,chr:2347-3346,chr:4348-5346 example2_sorted.bam

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
