#1. Pausepred
Pausepred is designed to predict ribosomal pauses using sorted BAM alignment files. It can be run as a webtool(http://pausepred.ucc.ie) or as a standalone version(github: https://github.com/romikasaini/Pausepred_offline).
It is written in perl language(perl 5, version 18, subversion 2 (v5.18.2))

## Inputs required to run the standalone version

###BAM_file ###window_size ###foldchange for pause ###reference fasta file ###read_length_min ###read_length_max ###coverage% ###upstream_seq ###downstream_seq

## Example Command
perl offline_pausepred.pl example1_sorted.bam 1000 10 example_sequence.fa 20 35 10 50 50


## Modules/Packages required to run standalone version
1. SAMtools(Follow this link to install http://www.htslib.org/download/)
1. Bio::DB::Fasta

Can be installed by installing BioPerl-1.6.1 follow this link to install bioperl http://bioperl.org/INSTALL.html
or 

You can use the CPAN shell to install these modules. For example

perl -MCPAN -e shell

install Bio::DB::Fasta


#2. Rfeet
Rfeet is designed to create ribosome profiles, to get a graphical veiw of the ribosomal density across gene/chr length. It can be run as a webtool(http://pausepred.ucc.ie/rfeet2.html) or as a standalone version(github: https://github.com/romikasaini/Pausepred_offline).
It is written in perl language(perl 5, version 18, subversion 2 (v5.18.2))

## Inputs required to run the standalone version
###first_bam_file ###Fasta sequence ###file ###gene/transscript/chr:strt-end ###second_bam_file

note: second file is optional

##Example command
perl offline_rfeet.pl perl plot_inputfile_optional.pl example1_sorted.bam example_sequence.fa chr:3347-4347 example2_sorted.bam

note: Second file is optional and multiple genes/chr locations can ploted at one time using a comma separated list. Example is given below.
perl offline_rfeet.pl perl plot_inputfile_optional.pl example1_sorted.bam example_sequence.fa chr:3347-4347,chr:2347-3346,chr:4348-5346 example2_sorted.bam

## Modules/Packages required to run standalone version
1. SAMtools(Follow this link to install http://www.htslib.org/download/)
2. List::Util
3. Bio::DB::Fasta(which can be installed by installing BioPerl-1.6.1)
4. Statistics::R
To install this package you need to install R in your OS first, since Statistics::R need to find R path to work fine. 

You can use the CPAN shell to install these modules. For example

perl -MCPAN -e shell

install Statistics::R
