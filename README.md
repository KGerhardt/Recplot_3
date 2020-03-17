# Recplot_4
Recruitment plot scripts designed to allow for the translation of multiple formats of reads into recruitment plots, to mass-generate recplots early in a project, and to generate interactive recruitment plots that allow for more thorough exploration of your data.

Currently, the recplot scripts are not packages that can be downloaded and loaded into R with library(...). That comes later.

There are 3 scripts that are now involved here:

*recplot_matrix.py*

Usage: Recplot_Prep.py -c [contigs FASTA file] -r [blast, magicblast, or SAM format reads] -f [format (sam, blast)], -o [output file name prefix] --interactive -m [MAGs file]

A python script which generates two files (prefix.lim and prefix.rec) used by the recplot scripts from a set of contigs and reads aligned to those contigs. 

Please note that you should filter these beforehand to have only contigs that you are interested in inside the contigs file, and only reads aligning to those contigs in the reads file.

*Recplot4_Static.R*

Requires these packages:

data.table

ggplot2

enveomics.R

cowplot

(for parallel processing) doParallel

Usage: Rscript/R CMD BATCH Recplot4_Static.R -dir [directory with lim and rec files]* -lib [R library directory]* -t [num threads] -id [pct. ID bin size] -w [genome bin width in bp] -sep [T or F to combine plots into 1 PDF or produce 1 PDF per plot, resp.] -lin_hist [T or F; base pair by pct ID histogram linear or log scale, resp.] -in_grp [Lower pct. ID boundary for a read to be classified as an in-group read]

* denotes a required option

This is a script designed for mass-producing recruitment plots from a set of many lim/rec files in a non-interactive session. Place all of the lim/rec files into a directory and point this script at them using -dir. All of the files in the directory will be transformed into recruitment plots, each in their own PDF if -sep T is used, or in a single PDF if -sep F is used.

If multiple threads are called, the script expects a unix environment. Parallelization occurs over lim/rec files, so each thread will handle one set of lim/rec files at a time. All other options are still available.

*Recplot4_Interactive.R*

Requires these packages:

data.table

ggplot2

enveomics.R

plotly

shiny

This script contains functions for producing interactive recruitment plots. Its contents are to be run from an interactive session of R, and after running the contents of the script, it expects manual usage of the recplot_suite function, where a user gives the prefix for a lim/rec pair in their current R session directory.

Use should be intuitive from there.
