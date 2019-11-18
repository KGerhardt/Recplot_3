# Recplot_4
Further development on enveomics.R's recruitment plots

Currently, this is not a package that can be downloaded and loaded into R with library(...). That comes later.

You will have to ensure that you have the following packages installed before anything else:

data.table
ggplot2
shiny
plotly
enveomics.R

As a result, download the R script Recplot4.R and run the full thing in an interactive session of R (meaning R studio, or R on the command line). Provided that you have the needed libraries, this will create the set of functions needed for the interactive recruitment plots to run, and will not execute any other code.

After running the script, begin an interactive plot from the same session by using the following command:

recplot_suite(prefix = "prefix_of_lim_or_rec_file")

This is the exact same file that one of the original recruitment plots would have used, and is produced from a set of blast-aligned reads in the same way as before using Miguel's script available here: http://enve-omics.ce.gatech.edu/enveomics/docs?t=BlastTab.catsbj.pl

From there, it should be intuitive.

I am developing this on R version 3.6.0. I do not know if older versions of the above packages will support this tool.
