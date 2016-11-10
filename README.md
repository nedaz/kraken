kraken
======

- Getting started 
  - Create a git clone by running git clone https://github.com/nedaz/kraken.git
  N.B. Make sure you use git to clone the repository and not use other methods such as svn checkout or download az zip 

- Installation
  - Use gcc version 5 or higher and CMake higher than 3.5
  - Run ./configure in parent directory
  - Run make -C build -j 10
  - Binaries will be located in /bin directory

- Get started by example 
The following example performs a mapping of items in the GTF input file "dmel.gtf" from specie "dmel" to "dyak". The data for this example can be found in the sample folder. The output is written to mapped.gtf as default. You can run the command:
../runKraken  -c dere_dyak_dmel.config -s dmel.gtf -S dmel -T dyak

For more information, see the software guide (doc/softwareGuide.pdf)
