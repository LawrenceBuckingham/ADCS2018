# ADCS2018

Source code and Supplementary Data to accompany paper presented at ADCS2018 in Dunedin, NZ; 11-12 December 2018.

Contents:
*	Source files used to generate data shown in paper. 
*	SP100000 test dataset
*	SP100000 homologs file (defines relevant sequences for each sequence)
*	SP100000 domains (defines extent of instances of each protein family in containing sequences)
*	Sample scripts (cross-validation for supervised method; all-vs-all for unsupervised method)
*	Charts showing precision-recall curves. *** TODO ***

Notes:
*	Make files are provided for g++ 7.3.0 under Cygwin and g++ 7.4.0 on Linux.
*	Linux executables are statically linked, because I had issues getting a suitably modern compiler installed on the workstations where the experiments were executed. Remove the various -static flags from the Linux makefile to create dynamically linked versions.
*	Scripts will have to be updated to suit your configuration. In particular, you will have to alter the number of threads and the directory structure to match how you place the data files and run the experiments. A few variables near the end of the scripts covers this.
*	Occasionally GIT does not cooperate with respect to line-endings in the scripts. I have done what I can to ensure that these are strictly UNIX.

L.B. 2018-12-19
