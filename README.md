kraken
======

Getting started - 
Download the package, unzip and run make from the main directory.  Use gcc version 4.2 or higher for this. You can start by running the code with the data provided in the samples directory. This is explained below.


Example - 
The following example performs a mapping of items in the GTF input file "dmel.gtf" from specie "dmel" to "dyak". The data for this example can be found in the sample folder. The output is written to mapped.gtf as default. You can run the command:
../runKraken  -c dere_dyak_dmel.config -s dmel.gtf -S dmel -T dyak

For more information, see the software guide (doc/softwareGuide.pdf)
