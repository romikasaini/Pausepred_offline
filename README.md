#1. Pausepred
Pausepred is designed to predict ribosomal pause sites using sorted BAM alignment file as input. It can be used as a web application (http://pausepred.ucc.ie) or as a stand-alone tool (github: https://github.com/romikasaini/Pausepred_offline).
This tool is written using Perl scripting language(Perl 5, version 18, subversion 2 (v5.18.2))

## Inputs required to run the standalone version:

BAM_file, window_size, foldchange for pause, reference fasta file, read_length_min, read_length_max, coverage%, upstream_seq and downstream_seq

Inputs should be mentioned in the order given above, delimited by space.

## Example Command:
>perl offline_pausepred.pl example1_sorted.bam 1000 10 example_sequence.fa 20 35 10 50 50


## Modules and Packages required for implementing the stand-alone version:
1. SAMtools(Follow this link to install http://www.htslib.org/download/)

Following is a perl module.
2. Bio::DB::Fasta

Bio::DB::Fasta be installed by installing BioPerl-1.6.1 follow this link to install bioperl http://bioperl.org/INSTALL.html
                                                  
or 
Please follow the instructions given at the bottom of this document(Perl module installation).


#2. Rfeet
Rfeet is designed to generate ribosome profiles, to get a graphical view of the footprint density across gene/chromosome length. It can be used as a web application (http://pausepred.ucc.ie/rfeet2.html) or as a standalone tool (github: https://github.com/romikasaini/Pausepred_offline).
This tool is written using Perl scripting language(Perl 5, version 18, subversion 2 (v5.18.2 )) and R programming language.

## Inputs required to run the standalone version:
first_bam_file, Fasta sequence, file, gene/transscript/chr:strt-end, second_bam_file

note: second file is optional

Inputs should be provided in the order mentioned above, delimited by space. Additional inputs such as legend names, plot type etc. will be requested through the console once the script is initiated.

##Example command
>perl offline_rfeet.pl perl plot_inputfile_optional.pl example1_sorted.bam example_sequence.fa chr:3347-4347 example2_sorted.bam

note: Second file is optional and multiple genes/chromosome locations can ploted at a time using a comma separated input as shown in the example given below.
perl offline_rfeet.pl perl plot_inputfile_optional.pl example1_sorted.bam example_sequence.fa chr:3347-4347,chr:2347-3346,chr:4348-5346 example2_sorted.bam

## Modules/Packages required to run standalone version
1. SAMtools(Follow this link to install http://www.htslib.org/download/)

Following are perl modules and they can be installed by following the perl module installation given below.
2. List::Util
3. Bio::DB::Fasta(which can be installed by installing BioPerl-1.6.1)
4. Statistics::R
To install this package you need to install R in your OS first, since Statistics::R need to find R path to work. 

##Perl module installation

You can use the CPAN shell to install these modules. Open CPAN shell by running following command.

>perl -MCPAN -e shell

and then Bio::DB::Fasta can be installed by running following command.

>install Bio::DB::Fasta

Please follow same steps to install List::Util and Statistics::R
