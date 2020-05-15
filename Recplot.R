library(reticulate)
library(ggplot2)
library(shiny)
library(data.table)
library(plotly)
library(cowplot)
library(enveomics.R)
library(shinyBS)

#use_condaenv()
use_python("/usr/local/bin/python")

###########################Development variables

#setwd("/mnt/c/Users/Kenji/Desktop/Recplot4")

#setwd("C:/Users/Kenji/Desktop/Recplot4/recplot_final_build/")

#This certainly makes the script accessible...
source_python("https://raw.githubusercontent.com/KGerhardt/Recplot_4/master/recplot_database_carlos_genes.py")

##############################################
#The python import is a space-efficient, but strcutrually awkward data object
#List of 2 items: list of lists of contig starts, stops, assoc. counts per %ID bin, and %ID bins.
#This function converts the structure into a recplot-ready data.table with appropriate labelling and returns some other key values for building the plots.
pydat_to_recplot_dat <- function(extracted_MAG, contig_names){
  
  id_breaks <- unlist(extracted_MAG[[2]])
  #id_breaks <- paste0(id_breaks - id_width, "-", id_breaks)
  
  extracted_MAG <- extracted_MAG[[1]]
  names(extracted_MAG) = contig_names
  
  ends <- lapply(extracted_MAG, function(x){
    
    return(x[[2]])
    
  })
  
  maximum_table <- data.table(contig = names(ends), length = lapply(ends, max))
  maximum_table[, relative_end := cumsum(length) + 1:nrow(maximum_table) - 1]
  
  pos.max <- maximum_table$relative_end[nrow(maximum_table)]
  
  bp_unit <- c("(bp)", "(Kbp)", "(Mbp)", "(Gbp)")[findInterval(log10(pos.max), c(0,3,6,9,12,Inf))]
  bp_div <- c(1, 1e3, 1e6, 1e9)[findInterval(log10(pos.max), c(0,3,6,9,12,Inf))]
  
  ends <- unlist(ends)
  
  bins <- rbindlist(lapply(extracted_MAG, function(x){
    
    return(data.table(cbind(do.call(rbind, x[[3]]))))
    
  }))
  
  colnames(bins) = as.character(id_breaks)
  
  starts <- lapply(extracted_MAG, function(x){
    
    return(x[[1]])
    
  })
  #needs start, end, contig name for annotation, rel.pos absolute for plotting
  
  bins[, Start := unlist(starts)]
  
  bins[, End := ends]
  
  bins[, seq_pos := seq(1/bp_div, pos.max/bp_div , length.out = nrow(bins))]
  
  bins[, contig := rep(names(starts), times = lengths(starts))]
  
  bins <- melt.data.table(bins, id.vars = c("contig", "Start", "End", "seq_pos"))
  
  colnames(bins)[5:6] = c("Pct_ID_bin", "bp_count")
  
  bins[, Pct_ID_bin := as.numeric(levels(bins$Pct_ID_bin))[bins$Pct_ID_bin]]
  
  #return(bins)
  return(list(bins, bp_unit,bp_div, pos.max))
  
}

#Transform the base data into the static recplot - use it for a print me! function.
#Add in the title as sample : MAG
create_static_plot <- function(base, bp_unit, bp_div, pos_max, in_grp_min, id_break, width, linear, showpeaks, ...){
  
  ends <- base[, max(End), by = contig]
  ends[, V1 := cumsum(V1) - 1 + 1:nrow(ends)]
  
  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
  
  #Lower left panel
  
  p <- ggplot(base, aes(x = seq_pos, y = Pct_ID_bin, fill=log10(bp_count)))+ 
    scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
    ylab("Percent Identity") +
    xlab(paste("Position in Genome", bp_unit)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "none", 
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    geom_raster()
  
  p <- p + annotate("rect", xmin = 0, xmax = pos_max/bp_div, 
                    ymin = in_grp_min, 
                    ymax = 100, fill = "darkblue", alpha = .15)
  
  read_rec_plot <- p + geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40")
  
  base[, group_label := ifelse(base$Pct_ID_bin-id_break >= in_grp_min, "depth.in", "depth.out")]
  setkeyv(base, c("group_label", "seq_pos"))
  
  #upper left panel
  
  depth_data <- base[, sum(bp_count/(End-Start+1)), by = key(base)]
  colnames(depth_data)[3] = "count"
  
  ddSave <- base[, sum(bp_count), by = key(base)]
  
  nil_depth_data <- depth_data[count == 0]
  nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
  
  seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
  
  depth_data$count[depth_data$count == 0] <- NA
  
  seq_depth_chart <- ggplot(depth_data, aes(x = seq_pos, y = count, colour=group_label, group = group_label))+
    geom_step(alpha = 0.75) +
    scale_y_continuous(trans = "log10", labels = scales::scientific) +
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = "none", 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          panel.background = element_blank(), 
          axis.title.x = element_blank(), 
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)) +
    scale_color_manual(values = group.colors) +
    ylab("Depth")
  
  if(nrow(nil_depth_data) > 0){
    seq_depth_chart <- seq_depth_chart + geom_segment(data = nil_depth_data, aes(x = seq_pos, xend = seq_pos, y = 0, yend = seg_upper_bound, color = group_label, group = group_label))
    
  }
  
  #Seq. depth histograms (top right panel)
  {
    depth_data <- depth_data[count > 0]
    seqdepth.lim <- range(c(depth_data$count[depth_data[,group_label == "depth.in"]], depth_data$count[depth_data[,group_label == "depth.out"]])) * c(1/2, 2)
    hist_binwidth <- (log10(seqdepth.lim[2]/2) - log10(seqdepth.lim[1] * 2))/199
    
    depth_data[,count := log10(count)]
    depth_data$group_label <-factor(depth_data$group_label, levels = c("depth.out", "depth.in"))
    depth_data <- depth_data[order(group_label),]
    
    p4 <- ggplot(depth_data, aes(x = count, fill = group_label)) +
      geom_histogram(binwidth = hist_binwidth) +
      scale_fill_manual(values = group.colors) +
      scale_y_continuous(expand=c(0,0)) +
      theme(legend.position = "none", 
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"), 
            panel.background = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.text.x = element_text(colour = "white"),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 14)) +
      xlab(label = element_blank()) +
      ylab(label = element_blank())
    
    if(showpeaks){
    
    enve.recplot2.__peakHist <- function
    ### Internal ancilliary function (see `enve.RecPlot2.Peak`).
    (x, mids, counts=TRUE){
      d.o <- x$param.hat
      if(length(x$log)==0) x$log <- FALSE
      if(x$log){
        d.o$x <- log(mids)
      }else{
        d.o$x <- mids
      }
      prob  <- do.call(paste('d', x$dist, sep=''), d.o)
      if(!counts) return(prob)
      if(length(x$values)>0) return(prob*length(x$values)/sum(prob))
      return(prob*x$n.hat/sum(prob))
    }
    
    h.breaks <- seq(log10(seqdepth.lim[1] * 2), log10(seqdepth.lim[2]/2), 
                    length.out = 200)
    h.mids <- (10^h.breaks[-1] + 10^h.breaks[-length(h.breaks)])/2
    
    min_info <- list()
    #min_info$pos.breaks = c(ends$V1, pos_max/bp_div)
    
    ends[, adjust := c(-1, V1[1:nrow(ends)-1])+1]
    
    base[, contiguous_end := End + ends$adjust[match(contig, ends$contig)]]
    
    min_info$pos.breaks = c(0, base$contiguous_end[match(ddSave$seq_pos[ddSave$group_label == "depth.in"], base$seq_pos)])
    
    min_info$pos.counts.in = ddSave$V1[ddSave$group_label == "depth.in"]
    
    #This is finicky
    
    try({
      peaks <- enve.recplot2.findPeaks(min_info)
      
      dpt <- signif(as.numeric(lapply(peaks, function(x) x$seq.depth)), 2)
      frx <- signif(100 * as.numeric(lapply(peaks,function(x) ifelse(length(x$values) == 0, x$n.hat, length(x$values))/x$n.total)), 2)
      
      
      
      if (peaks[[1]]$err.res < 0) {
        err <- paste(", LL:", signif(peaks[[1]]$err.res,3))
      }  else {
        err <- paste(", err:", signif(as.numeric(lapply(peaks, function(x) x$err.res)), 2))
      }
      labels <- paste(letters[1:length(peaks)], ". ", dpt, "X (", frx, "%", err, ")", sep = "")
      
      peak_counts <- lapply(peaks, enve.recplot2.__peakHist, h.mids)
      
      plot_breaks = h.breaks[-length(h.breaks)]
      
      gg_peak_info <- data.table(plot_breaks = rep(plot_breaks, length(peak_counts)), count = unlist(peak_counts), grp = rep(labels, each = length(plot_breaks)))
      
      p4 <- p4 + geom_line(data = gg_peak_info, aes(x = plot_breaks, y = count, color = grp, group = grp), inherit.aes = F, color = "red", lwd = 1.13)
      
    })
    
    }
    
    p4 <- p4 + coord_flip()
    
    if(showpeaks){
    
    try({
      o_max <- max(table(findInterval(depth_data$count[depth_data$group_label == "depth.in"], h.breaks)))*.75
      x_start = 0
      
      if(length(labels) > 0){
        for(i in labels){
          p4 <- p4 + annotate("text", label = i, y = o_max, x = x_start, size = 3)
          x_start = x_start - .185
        }
      }
    })
      
    }
    
    seq_depth_hist <- p4
    
    rm(p4)
    
  }
  
  #bp counts histogram (bottom right panel)
  {
    
    bp_data <- base[,sum(bp_count), by = Pct_ID_bin]
    
    if(linear == 1){
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0)) + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"), 
              panel.background = element_rect(fill = "#EEF7FA"), 
              axis.text.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 14))+
        ylab("Base Pair Count by ANI")
    } else {
      
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0), trans = "log10") + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"), 
              panel.background = element_rect(fill = "#EEF7FA"), 
              axis.text.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 14))+
        ylab("Base Pair Count by ANI")
    }
    
    bp_count_hist <- p4 + annotate("rect", xmin = in_grp_min, xmax = 100, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15) + coord_flip()
    
    rm(p4)
    
  }
  
  overall_plot <- plot_grid(seq_depth_chart, seq_depth_hist, read_rec_plot, bp_count_hist, align = "hv", ncol = 2, rel_widths = c(2.7, 1), rel_heights = c(1, 2.3))
  
  return(overall_plot)
  
}

gene_pydat_to_recplot_dat_prodigal <- function(prodigal_gene_mess){
  
  contigs <- names(prodigal_gene_mess)
  
  lengths <- unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(length(x[[1]]))
    
  })))
  
  pretty_data <- data.table(contig = rep(contigs, times = lengths))
  
  
  
  
  pretty_data[, gene_name := unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(x[[1]])
    
  }))) ]
  pretty_data[, start := unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(x[[2]])
    
  }))) ]
  pretty_data[, end  := unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(x[[3]])
    
  }))) ]
  pretty_data[, localizer :=  (start+end)/2]
  pretty_data[, strand  := unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(x[[4]])
    
  }))) ]
  pretty_data[, annotation  := unname(unlist(lapply(prodigal_gene_mess, function(x){
    
    return(x[[5]])
    
  }))) ]
  
  
  return(pretty_data)
  
  
}

#This is the GUI function
recplot_landing_page <- function(){
  
  if (interactive()) {
    
    system <- get_sys()
    
    if(system == "Windows"){
      format_choices <- c("Tabular BLAST" = "blast", "SAM" = "sam")
    }else{
      format_choices <- c("Tabular BLAST" = "blast", "SAM" = "sam", "BAM" = "bam")
    }
    
    gene_choices <- c("Prodigal GFF" = "prodigal")
    
    ui <- fluidPage(
      
      tabsetPanel(id = "tabs",
                  tabPanel("Database Management",
                           fluidRow(
                             column(4, 
                                    h2("Create a new Database"),
                                    actionButton('dir', '(1) Choose Directory', icon = icon("folder-open")),
                                    textInput("cur_dir",label = NULL, value = paste("Working in:", getwd())),
                                    br(),
                                    actionButton('contigs', '(2) Select Contigs', icon = icon("file-upload")),
                                    textInput("contig_file",label = NULL, value = "No contigs selected."),
                                    br(),
                                    actionButton('mags', '(3) Select MAGs', icon = icon("file-upload")),
                                    textInput("mag_file",label = NULL, value = "No MAGs file selected."),
                                    br(),
                                    actionButton('reads', '(4) Select Reads', icon = icon("file-upload")),
                                    textInput("read_file",label = NULL, value = "No mapped read file selected."),
                                    
                                    
                                    selectInput('fmt', 'Mapped Read Format', selected = "Tabular BLAST", choices = format_choices),
                                    br(),
                                    textInput("dbname",label = "(5) Name the database", value = "Enter name here."),
                                    br(),
                                    actionButton('db' , "(6) Create database", icon("coins")),
                                    
                                    
                                    #Tooltips
                                    
                                    bsTooltip("dir", "(Optional) Select a working directory. The database and any saved plots will be placed here.", placement = "right"),
                                    bsTooltip("contigs", "Select a FASTA format file containing assembled contig DNA sequences.", placement = "right"),
                                    bsTooltip("mags", "Select a MAGs association file. This should be a 2-column, tab-separated file with the name of each contig in the contigs file in the first column, and the name of the MAG to which each contig belongs in the second column.", placement = "right"),
                                    bsTooltip("reads", "Select a mapped read file. These reads should be mapped to the contigs in the contig file selected above.", placement = "right"),
                                    bsTooltip("fmt", "Select the format of the mapped reads to be added to the new database.", placement = "right"),
                                    bsTooltip("dbname", "Name your database. A .db extension will be added to the end of the name you give it.", placement = "right"),
                                    bsTooltip("db", "Click me after selecting all input files to create your database.", placement = "right")
                             ),
                             
                             column(4, 
                                    h2("Work with an existing database"),
                                    br(),
                                    actionButton('exist_db', 'Select an existing DB', icon = icon("coins")),
                                    textInput("exist_dbname",label = NULL, value = "No DB currently selected"),
                                    br(),
                                    h4("Add More Reads"),
                                    actionButton('add_sample', 'Select another sample to add?', icon = icon("file-upload")),
                                    textInput("add_samp",label = NULL, value = "No new sample to add."),
                                    selectInput('fmt_add', 'Mapped Read Format', selected = "Tabular BLAST", choices = format_choices),
                                    actionButton('new_samp_commit', "Add this sample to the DB", icon = icon("coins")),
                                    br(),
                                    br(),
                                    h4("Add Genes"),
                                    actionButton('add_genes', 'Add genes to database?', icon = icon("file-upload")),
                                    textInput("add_gen", label = NULL, value = "No genes to add."),
                                    selectInput('fmt_gen', 'Gene format', selected = "Prodigal GFF", choices = gene_choices),
                                    actionButton('genes_commit', "Add these genes to the DB", icon = icon("coins")),
                                    
                                    selectInput('task', 'Plot contigs or plot genes?', selected = "Contigs", choices = c("Contigs" = "contigs", "Genes" = "genes")),
                                    
                                    bsTooltip("exist_db", "Select a database previously created with Recruitment Plot.", placement = "right"),
                                    bsTooltip("add_sample", "(Optional) Select another set of reads mapped to the same contigs to be added to the database.", placement = "right"),
                                    bsTooltip("fmt_add", "Select the format of the mapped reads to be added to the existing database.", placement = "right"),
                                    bsTooltip("add_genes", "(Optional) Add genes for existing contigs.", placement = "right"),
                                    bsTooltip("fmt_gen", "Select the format of the genes to be added to the existing database. Currently only Prodigal GFF format is supported.", placement = "right"),
                                    bsTooltip("new_samp_commit", "Once you have selected another set of mapped reads to add and chosen the format, click this to add the sample. The sample will not be added until you do.", placement = "right"),
                                    bsTooltip("genes_commit", "Once you have selected a set of genes to add and chosen the format (currently only Prodigal GFF), click this to add the genes. The genes will not be added until you do.", placement = "right")
                                    
                                    
                             ),
                             column(4,
                                    h2("Build Report"),
                                    textOutput("message")
                             )        
                           )
                  ),
                  
                  tabPanel("Recruitment Plot", 
                           sidebarPanel(
                             width = 3,
                             
                             h4("Choose a MAG"),
                             
                             selectInput("samples", "(1) Select a sample in the database", selected = NULL, choices = NULL),
                             
                             selectInput("mags_in_db", "(2) Select a MAG in the sample", selected = NULL, choices = NULL),
                             
                             h4("Select Bin Resolution"),
                             
                             numericInput("width", "(3) Genome Resolution", min = 75, max = 5000, value = 1000),
                             
                             numericInput("height", "(4) Pct. ID Resolution", min = 0.05, max = 3, value = 0.5),
                             
                             numericInput("low_bound", "(5) Minimum Pct. ID", min = 50, max = 95, value = 70),
                             
                             h4("Load Selected MAG"),
                             
                             actionButton('get_a_mag', '(6) View Selected MAG', icon = icon("jedi-order")),
                             
                             h4("Fine Tuning (Interactive)"),
                             
                             numericInput("in_group_min_stat", "(7) In-Group Pct. ID", min = 50, max = 99.5, value = 90),
                             selectInput("linear_stat", "(8) BP Histogram Scale", choices = c("Linear" = 1, "Logarithmic" = 2), selected = 1),
                             
                             checkboxInput("show_peaks", "(9) Display Depth Peaks?"),
                             
                             textInput("pdf_name", "Name and Save?"),
                             actionButton("print_stat", "Save to PDF", icon = icon("save")),
                             
                             bsTooltip("samples", "This menu contains a list of samples within the database. Select one, and the MAGs field will be populated with the MAGs found in that sample.", placement = "right"),
                             bsTooltip("mags_in_db", "This menu contains the set of MAGs in currently selected sample. Select one, then select resolution parameters.", placement = "right"),
                             bsTooltip("width", "Approximate number of base pairs in each genome window. The Recruitment Plot attempts to normalize bin width for each contig to this size. Lower values = higher resolution, but is slower. Higher values = lower resolution, but is faster.", placement = "right"),
                             bsTooltip("height", "Controls the resolution of percent identity to the reference. Lower values here will result in finer resolution, but will be slower. Hint: The default 0.5% window means a resolution of 1 base pair mismatch per 200 bases; finer resolution is probably uneccessary.", placement = "right"),
                             bsTooltip("low_bound", "Reads mapping below this percent identity will not be included in the current recruitment plot.", placement = "right"),
                             bsTooltip("get_a_mag", "Click this to load the current MAG into the viewer and plot it. Please wait for the plot to appear after clicking this. This loads the data for all tabs.", placement = "right"),
                             bsTooltip("in_group_min_stat", "Controls the lower edge of the shaded region in the recruitment plot's lower panels. Reads mapping at or above this percent identity are regarded as the \"in-group\" for the Recruitment Plot, and are represented by the dark blue lines in the upper panels.", placement = "right"),
                             bsTooltip("linear_stat", "Causes the lower right panel to display base pair counts per percent identity bin in linear scale or log scale.", placement = "right"),
                             bsTooltip("print_stat", "After loading a plot (meaning you should be able to see it), add a name in the associated text box and then click this to print a PDF of the current", placement = "right"),
                             bsTooltip("show_peaks", "Calculate and overlay peaks for the depth of coverage histogram (top right panel)", placement = "right"),
                             
                             bsTooltip("recplot_main", "Bottom left panel: a 2-D histogram of the counts of base pairs falling into a bin defined by position in the genome (x-axis) and percent identity (y-axis) Bins are as wide as the genome resolution parameter                             if viewing contigs, and cover genes & intergenic regions in contiguous chunks if viewing genes. The shaded section is the current in-group, which is the dark blue line on the top two plots Top left panel: Average sequencing depth for each x-axis bin on in the bottom left panel. Dark blue corresponds to depth of coverage for bins in the in-group, and light blue to the out-group. Segments at the bottom of the plot have zero coverage. Top right panel: a histogram of the depths of coverage observed in the corresponding in/out group in the sequencing depths chart (top left). If peaks are selected, they correspond to the estimates of the genome's average sequencing depth. Bottom right panel: a histogram of the number of bases falling into each percent identity bin across the entire genome, displayed in linear or log scale depending on your selection.", trigger = "click", placement = "left")
                             ),
                           mainPanel(id = "recplot_main",
                            #Spacing
                            fluidRow(
                              column(12,
                                     div(style = "height:60px;background-color: white;", "")
                            )
                            ),
                            fluidRow(
                              
                              plotOutput("read_recruitment_plot", height = "850px")
                              
                            )
                            
                            
                           )
                  ),
                  tabPanel("Hover Plot", 
                           sidebarPanel(
                             width = 3,
                             
                             h3("Choose a MAG"),
                             
                             selectInput("samples_interact", "(1) Select a sample in the database", selected = NULL, choices = NULL),
                             
                             selectInput("mags_in_db_interact", "(2) Select a MAG in the sample", selected = NULL, choices = NULL),
                             
                             h3("Select Bin Resolution"),
                             
                             numericInput("width_interact", "(3) Genome Resolution", min = 75, max = 5000, value = 1000),
                             
                             numericInput("height_interact", "(4) Pct. ID Resolution", min = 0.05, max = 3, value = 0.5),
                             
                             numericInput("low_bound_interact", "(5) Minimum Pct. ID", min = 50, max = 95, value = 70),
                             
                             h3("Load Selected MAG"),
                             
                             actionButton('get_a_mag_interact', '(6) View selected MAG', icon = icon("jedi-order")),
                             
                             h3("Fine Tuning (Interactive)"),
                             
                             numericInput("in_group_min_interact", "(7) In-Group Pct. ID", min = 50, max = 99.5, value = 95),
                             
                             bsTooltip("samples_interact", "This menu contains a list of samples within the database. Select one, and the MAGs field will be populated with the MAGs found in that sample.", placement = "right"),
                             bsTooltip("mags_in_db_interact", "This menu contains the set of MAGs in currently selected sample. Select one, then select resolution parameters.", placement = "right"),
                             bsTooltip("width_interact", "Approximate number of base pairs in each genome window. The Recruitment Plot attempts to normalize bin width for each contig to this size. Lower values = higher resolution, but is slower. Higher values = lower resolution, but is faster.", placement = "right"),
                             bsTooltip("height_interact", "Controls the resolution of percent identity to the reference. Lower values here will result in finer resolution, but will be slower. Hint: The default 0.5% window means a resolution of 1 base pair mismatch per 200 bases; finer resolution is probably uneccessary.", placement = "right"),
                             bsTooltip("low_bound_interact", "Reads mapping below this percent identity will not be included in the current recruitment plot.", placement = "right"),
                             bsTooltip("get_a_mag_interact", "Click this to load the current MAG into the viewer and plot it. Please wait for the plot to appear after clicking this.", placement = "right"),
                             bsTooltip("in_group_min_interact", "Controls the lower edge of the shaded region in the recruitment plot's lower panels. Reads mapping at or above this percent identity are regarded as the \"in-group\" for the Recruitment Plot, and are represented by the dark blue lines in the upper panel.", placement = "right")
                           ),
                           mainPanel(
                             column(12, 
                                    plotlyOutput("Plotly_seq_depth", height = "286px"),
                                    plotlyOutput("Plotly_read_rec", height = "574px")
                             )
                           )
                  )
      )
      #End fluid page
    )
    
    
    
    server <- function(input, output, session) {
      
      output$message <- renderText("Welcome to Recruitment Plot! Please select your inputs.")
      
      directory <- getwd()
      reads <- "No file selected. Try again?"
      contigs <- "No file selected. Try again?"
      mags <- "No file selected. Try again?"
      new_samp <- "No new sample. Try again?"
      db <- "No existing database selected. Try again?"
      new_genes <- "No genes selected. Try again?"
      
      plotting_materials <- NA
      
      gene_data <- NA
      
      exist_db <- "No existing database selected. Try again?"
      
      samples_in_db <- "No database selected or built yet."
      
      #Database building
      
      observeEvent(input$dir, {
        
        tryCatch({
          directory <- choose.dir()
        },
        error = function(cond){
          directory <- getwd()
          return(directory)
        })
        
        setwd(directory)
        
        updateTextInput(session, "cur_dir", value = paste("Working in:", directory))
        
      })
      
      observeEvent(input$reads, {
        tryCatch({
          reads <- file.choose()
        },
        error = function(cond){
          reads <- "No file selected. Try again?"
          return(reads)
        })
        updateTextInput(session, "read_file", value = reads)
        
      })
      
      observeEvent(input$contigs, {
        tryCatch({
          contigs <- file.choose()
        },
        error = function(cond){
          contigs <- "No file selected. Try again?"
          return(contigs)
        })
        
        updateTextInput(session, "contig_file", value = contigs)
      })
      
      observeEvent(input$mags, {
        tryCatch({
          mags <- file.choose()
        },
        error = function(cond){
          mags <- "No file selected. Try again?"
          return(mags)
        })
        updateTextInput(session, "mag_file", value = mags)
      })
      
      observeEvent(input$db, {
        
        ready_to_make = TRUE
        
        #check to make sure that inputs are set and exist:
        if(input$read_file == "No mapped read file selected."){
          ready_to_make <- F
        }
        if(input$contig_file == "No contigs selected."){
          ready_to_make <- F
        }
        if(input$mag_file == "No MAGs file selected."){
          ready_to_make <- F
        }
        
        if(!ready_to_make){
          output$message <- renderText("Not all files have been selected. Please select a set of reads, contigs, and a MAGs file.")
          return(NA)
        }
        
        if(!file.exists(input$read_file)){
          ready_to_make <- F
        }
        if(!file.exists(input$contig_file)){
          ready_to_make <- F
        }
        if(!file.exists(input$mag_file)){
          ready_to_make <- F
        }
        
        if(input$dbname == "Enter name here.."){
          output$message <- renderText(paste("The database needs a name. Please give it one!"))
          ready_to_make <- F
        }
        
        if(!ready_to_make){
          output$message <- renderText("Not all files selected seem to be found, or the database hasn't been given a name. Please select files and Enter name here..")
        }else{
          output$message <- renderText(paste("Database in creation. Please wait..."))
          
          sqldb_creation(contigs = input$contig_file, mags = input$mag_file, sample_reads = list(input$read_file), map_format = input$fmt, database = paste0(input$dbname, ".db"))
          
          cat("done!\n")
          
          output$message <- renderText(paste("Database built!"))
          
          updateTextInput(session, "exist_dbname", value = paste0(input$dbname, ".db"))
          
        }
        
      })
      
      observeEvent(input$exist_db,{
        
        tryCatch({
          db <- file.choose()
        },
        error = function(cond){
          db <- "No existing database selected. Try again?"
          return(db)
        })
        
        updateTextInput(session, "exist_dbname", value = db)
        
      })
      
      observeEvent(input$add_sample,{
        
        #Add a don't-do-this if there's not a selected DB
        if(input$exist_dbname == "" | input$exist_dbname == "No existing database selected. Try again?"){
          output$message <- renderText("You have to select an existing database first.")
          return(NA)
        }
        
        tryCatch({
          new_samp <- file.choose()
        },
        error = function(cond){
          new_samp <- "No new sample. Try again?"
          return(new_samp)
        })
        
        updateTextInput(session, "add_samp", value = new_samp)
        
      })
      
      observeEvent(input$new_samp_commit, {
        
        if(input$exist_dbname == "No DB currently selected" | input$exist_dbname == "No existing database selected. Try again?"){
          output$message <- renderText("You cannot add a new sample without making a new database or choosing an existing one first.")
          return(NA)
        }
        
        if(input$add_samp == "No new sample to add." | input$add_samp == "No new sample. Try again?"){
          
          output$message <- renderText("You must choose a sample before committing it to the database.")
          
          return(NA)
          
        }
        
        
        add_sample(input$exist_dbname, list(input$add_samp), input$fmt_add)
        
        cat("done!\n")
        
        output$message <- renderText("Sample added!")
        
        samples_in_db = assess_samples(input$exist_dbname)
        
        labels <- unlist(samples_in_db)
        
        samples_in_db <- unlist(samples_in_db)
        names(samples_in_db) = labels
        
        updateSelectInput(session, "samples", choices = samples_in_db)
        updateSelectInput(session, "samples_interact", choices = samples_in_db)
        
        
      })
      
      observeEvent(input$add_genes,{
        
        #Add a don't-do-this if there's not a selected DB
        if(input$exist_dbname == "" | input$exist_dbname == "No existing database selected. Try again?"){
          output$message <- renderText("You have to select an existing database first.")
          return(NA)
        }
        
        tryCatch({
          new_genes <- file.choose()
        },
        error = function(cond){
          new_genes <- "No genes selected. Try again?"
          return(new_genes)
        })
        
        updateTextInput(session, "add_gen", value = new_genes)
        
      })
      
      observeEvent(input$genes_commit, {
        
        if(input$exist_dbname == "No DB currently selected" | input$exist_dbname == "No existing database selected. Try again?"){
          output$message <- renderText("You cannot add a new sample without making a new database or choosing an existing one first.")
          return(NA)
        }
        
        if(input$add_gen == "No genes to add." | input$add_gen == "No genes selected. Try again?"){
          
          output$message <- renderText("You must choose genes before committing them to the database.")
          
          return(NA)
          
        }
        
        
        add_genes_to_db(input$exist_dbname, input$add_gen, input$fmt_gen)
        
        cat("done!\n")
        
        output$message <- renderText("Genes added!")
        
        
      })
      
      observeEvent(input$exist_dbname, {
        
        if(input$exist_dbname == "No DB currently selected" | input$exist_dbname == "No existing database selected. Try again?"){
          samples_in_db <- c("nothing selected"="Nothing Selected")
        }else{
          
          tryCatch({
            samples_in_db = assess_samples(input$exist_dbname)

          },
          error = function(cond){
            output$message <- renderText(paste("No existing database selected. Try again?"))
            return(NA)
          })
          

          labels <- unlist(samples_in_db)
          
          samples_in_db <- unlist(samples_in_db)
          names(samples_in_db) = labels
          
        }
        
        updateSelectInput(session, "samples", choices = samples_in_db)
        updateSelectInput(session, "samples_interact", choices = samples_in_db)
        
        
      })
      
      #Plot pages
      
      observeEvent(input$samples, {
        
        
        if(input$exist_dbname == "No DB currently selected"){
          
          mags_in_samp <- c("nothing selected"="Nothing Selected")
          
        }else{
          
          if(input$samples == "" | input$samples == "Nothing Selected"){
            
            mags_in_samp <- c("nothing selected"="Nothing Selected")
            
          }else{
            mags_in_samp = unlist(assess_MAGs(input$exist_dbname, input$samples))
            
            labels <- mags_in_samp
            
            mags_in_samp <- unlist(mags_in_samp)
            
            names(mags_in_samp) = labels
            
            #print(mags_in_samp)
          }
        }
        
        updateSelectInput(session, "mags_in_db", choices = mags_in_samp)
        updateSelectInput(session, "mags_in_db_interact", choices = mags_in_samp)
        
      })
      
      observeEvent(input$get_a_mag, {
        if(input$exist_dbname == "No DB currently selected" | input$samples == "Select a sample in the database" | input$mags_in_db == "Select a MAG in the sample"){
          recplot_data <- data.frame(placeholder = 1)
        }else{
          
          if(input$task == "contigs"){
            recplot_data <- extract_MAG_for_R(input$exist_dbname, input$samples, input$mags_in_db, input$width, input$height, input$low_bound)
          }else{
            recplot_data <- extract_genes_MAG_for_R(input$exist_dbname, input$samples, input$mags_in_db, input$height, input$low_bound)
            gene_data <<- gene_pydat_to_recplot_dat_prodigal(recplot_data[[3]])
            
          }
          
          contig_names <- unlist(get_contig_names(input$exist_dbname, input$mags_in_db))
          
          plotting_materials <<- pydat_to_recplot_dat(recplot_data, contig_names)
          
          cat("done!\n")
          
          old_value <- input$in_group_min_stat
          
          updateNumericInput(session, "in_group_min_stat", value = old_value-1)
          updateNumericInput(session, "in_group_min_interact", value = old_value-1)
          
          updateNumericInput(session, "in_group_min_stat", value = old_value)
          updateNumericInput(session, "in_group_min_interact", value = old_value)
          
        }
        
      })
      
      observeEvent(input$get_a_mag_interact, {
        if(input$exist_dbname == "No DB currently selected" | input$samples == "Select a sample in the database" | input$mags_in_db == "Select a MAG in the sample"){
          recplot_data <- data.frame(placeholder = 1)
        }else{
          
          if(input$task == "contigs"){
            recplot_data <- extract_MAG_for_R(input$exist_dbname, input$samples, input$mags_in_db, input$width, input$height, input$low_bound)
          }else{
            recplot_data <- extract_genes_MAG_for_R(input$exist_dbname, input$samples, input$mags_in_db, input$height, input$low_bound)
            gene_data <<- gene_pydat_to_recplot_dat_prodigal(recplot_data[[3]])
          }
          
          contig_names <- unlist(get_contig_names(input$exist_dbname, input$mags_in_db))
          
          plotting_materials <<- pydat_to_recplot_dat(recplot_data, contig_names)
          
          old_value <- input$in_group_min_stat
          
          updateNumericInput(session, "in_group_min_stat", value = old_value-1)
          updateNumericInput(session, "in_group_min_interact", value = old_value-1)
          
          updateNumericInput(session, "in_group_min_stat", value = old_value)
          updateNumericInput(session, "in_group_min_interact", value = old_value)
          
          cat("done!\n")
          
        }
        
      })
      
      observeEvent(input$print_stat, {
        
        if(is.na(plotting_materials)[1]){
          print("Cannot plot an unloaded mag!")
          return(NA)
        }else{
          
          if(input$pdf_name == ""){
            print("I need a name first!")
            return(NA)
          }else{
            
            base <- one_mag()
            
            req(!is.na(base))
            
            bp_unit <- base[[2]]
            bp_div <- base[[3]]
            pos_max <- base[[4]]
            base <- base[[1]]
            
            static_plot <- create_static_plot(base = base,
                                              bp_unit = bp_unit,
                                              bp_div = bp_div,
                                              pos_max = pos_max,
                                              in_grp_min = input$in_group_min_stat,
                                              id_break = input$height,
                                              width = input$width,
                                              linear = input$linear_stat
            )
            
            title <- ggdraw() + 
              draw_label(
                paste("Read Recruitment Plot for sample:", input$samples, "MAG:", input$mags_in_db),
                fontface = 'bold',
                x = 0,
                hjust = 0
              ) +
              theme(
                # add margin on the left of the drawing canvas,
                # so title is aligned with left edge of first plot
                plot.margin = margin(0, 0, 0, 7)
              )
            static_plot <- plot_grid(
              title, static_plot,
              ncol = 1,
              # rel_heights values control vertical title margins
              rel_heights = c(0.1, 1)
            )
            
            pdf(paste0(input$pdf_name, ".pdf"), height = 11, width = 17)
            
            print(static_plot)
            
            dev.off()
            
            
          }
          
        }
        
        
        
      })
      
      one_mag <- reactive({
        
        input$in_group_min_stat
        input$linear_stat
        
        plotting_materials
        
      })
      
      #Static plots
      
      output$read_recruitment_plot <- renderPlot({
        base <- one_mag()
        
        req(!is.na(base))
        
        bp_unit <- base[[2]]
        bp_div <- base[[3]]
        pos_max <- base[[4]]
        base <- base[[1]]
        
        static_plot <- create_static_plot(base = base,
                                          bp_unit = bp_unit,
                                          bp_div = bp_div,
                                          pos_max = pos_max,
                                          in_grp_min = input$in_group_min_stat,
                                          id_break = input$height,
                                          width = input$width,
                                          linear = input$linear_stat,
                                          showpeaks = input$show_peaks
        )
        
        return(static_plot)
      })
      
      #Hover plots
      
      output$Plotly_read_rec<-renderPlotly({
        
        base <- one_mag()
        
        req(!is.na(base))
        
        bp_unit <- base[[2]]
        bp_div <- base[[3]]
        pos_max <- base[[4]]
        base <- base[[1]]
        
        ends <- base[, max(End), by = contig]
        ends[, V1 := cumsum(V1) - 1 + 1:nrow(ends)]
        
        widths <- base$End - base$Start +1
        
        base$bp_count <- base$bp_count*(input$width/widths)
        
        if(input$task == "contigs"){
          p <- ggplot(base, aes(x = seq_pos, y = Pct_ID_bin, fill=log10(bp_count), text = paste0("Contig: ", contig,
                                                                                                 "\nPos. in Contig: ", Start, "-", End,
                                                                                                 "\nNorm. Bin Count: ", round(bp_count))))+ 
            scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
            ylab("Percent Identity") +
            xlab(paste("Position in Genome", bp_unit)) +
            scale_y_continuous(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0, 0)) +
            theme(legend.position = "none", 
                  axis.line = element_line(colour = "black"),
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 14)) +
            annotate("rect", xmin = 0, xmax = pos_max/bp_div, 
                     ymin = input$in_group_min_interact, 
                     ymax = 100, fill = "darkblue", alpha = .15)+
            geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40") +
            geom_raster()
        }else{
          
          base[, gene_id := "Non-Genic"]
          base[, gene_start := "N/A"]
          base[, gene_end := ""]
          base[, strand := "N/A"]
          base[, annotation := "N/A"]
          
          matcher <- base[, list(unique(Start), unique(End)), by = contig]
          
          
          p <- ggplot(base, aes(x = seq_pos, y = Pct_ID_bin, fill=log10(bp_count), text = paste0("Contig: ", contig,
                                                                                                 "\nPos. in Contig: ", Start, "-", End,
                                                                                                 "\nNorm. Bin Count: ", round(bp_count))))+ 
            scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
            ylab("Percent Identity") +
            xlab(paste("Position in Genome", bp_unit)) +
            scale_y_continuous(expand = c(0, 0)) +
            scale_x_continuous(expand = c(0, 0)) +
            theme(legend.position = "none", 
                  axis.line = element_line(colour = "black"),
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 14)) +
            annotate("rect", xmin = 0, xmax = pos_max/bp_div, 
                     ymin = input$in_group_min_interact, 
                     ymax = 100, fill = "darkblue", alpha = .15)+
            geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40") +
            geom_raster()
        }
        
        read_rec_plot <- ggplotly(p, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "grey90") %>% 
          style(hoverinfo = "none", traces = c(1, 2)) 
        
        return(read_rec_plot)
        
      })
      
      output$Plotly_seq_depth <- renderPlotly({
        
        base <- one_mag()
        
        req(!is.na(base))
        
        bp_unit <- base[[2]]
        bp_div <- base[[3]]
        pos_max <- base[[4]]
        base <- base[[1]]
        
        base$group_label <- ifelse(base$Pct_ID_bin-input$height >= input$in_group_min_interact, "depth.in", "depth.out")
        
        setkeyv(base, c("group_label", "seq_pos"))
        
        depth_data <- base[, list(sum(bp_count/(End-Start+1)), unique(Start), unique(End), unique(contig)), by = key(base)]
        colnames(depth_data)[3:6] = c("count", "Start", "End", "contig")
        
        group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
        
        ddSave <- base[, sum(bp_count), by = key(base)]
        
        nil_depth_data <- depth_data[count == 0]
        nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
        
        seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
        
        depth_data$count[depth_data$count == 0] <- NA
        
        if(input$task == "contigs"){
        
        seq_depth_chart <- ggplot(depth_data, aes(x = seq_pos, y = count, colour=group_label, group = group_label, text = paste0("Contig: ", contig,
                                                                                                                                 "\nPos. in Contig: ", Start, "-", End,
                                                                                                                      "\nSeq. Depth: ", round(count))))+
          geom_step(alpha = 0.75) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(expand=c(0,0))+
          theme(legend.position = "none", 
                panel.border = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"), 
                panel.background = element_blank(), 
                axis.title.x = element_blank(), 
                #axis.text.y = element_text(angle = 90, hjust = 0.5),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 14)) +
          scale_color_manual(values = group.colors) +
          ylab("Log 10 Depth")
        }else{
          
          
          depth_data[, gene_id := "Non-Genic"]
          depth_data[, gene_start := "N/A"]
          depth_data[, gene_end := ""]
          depth_data[, strand := "N/A"]
          depth_data[, annotation := ""]
          
          seq_depth_chart <- ggplot(depth_data, aes(x = seq_pos, y = count, colour=group_label, group = group_label, text = paste0("Contig: ", contig,
                                                                                                                                   "\nPos. in Contig: ", Start, "-", End,
                                                                                                                                   "\nSeq. Depth: ", round(count))))+
            geom_step(alpha = 0.75) +
            scale_y_continuous(trans = "log10") +
            scale_x_continuous(expand=c(0,0))+
            theme(legend.position = "none", 
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"), 
                  panel.background = element_blank(), 
                  axis.title.x = element_blank(), 
                  #axis.text.y = element_text(angle = 90, hjust = 0.5),
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 14)) +
            scale_color_manual(values = group.colors) +
            ylab("Log 10 Depth")
          
        }
        
        a <- list(
          range = c(10^min(depth_data$count, na.rm = T), 10^max(depth_data$count, na.rm = T)),
          showticklabels = TRUE,
          exponentformat = "e"
        )
        
        seq_depth_chart <- ggplotly(seq_depth_chart, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "grey90", yaxis = a)
        
        
        return(seq_depth_chart)
        
      })
      
      #Consistency between plot panels
      
      observeEvent(input$tabs, {
        
        if(input$tabs == "Recruitment Plot"){
          if(!is.null(input$samples_interact)){
            updateSelectInput(session, "samples", selected = input$samples_interact)
          }
          if(!is.null(input$mags_in_db_interact)){
            updateSelectInput(session, "mags_in_db", selected = input$mags_in_db_interact)
          }
          updateNumericInput(session, "width", value = input$width_interact)
          updateNumericInput(session, "height", value = input$height_interact)
          updateNumericInput(session, "low_bound", value = input$low_bound_interact)
          updateNumericInput(session, "in_group_min_stat", value = input$in_group_min_interact)
        }
        if(input$tabs == "Hover Plot"){
          if(!is.null(input$samples)){
            updateSelectInput(session, "samples_interact", selected = input$samples)
          }
          if(!is.null(input$mags_in_db)){
            updateSelectInput(session, "mags_in_db_interact", selected = input$mags_in_db)
          }
          updateNumericInput(session, "width_interact", value = input$width)
          updateNumericInput(session, "height_interact", value = input$height)
          updateNumericInput(session, "low_bound_interact", value = input$low_bound)
          updateNumericInput(session, "in_group_min_interact", value = input$in_group_min_stat)
        }
        
      })
      
      
    }
    
    runApp(list(ui = ui, server = server), launch.browser = T)
  }
  
  
}

recplot_landing_page()

#testing peaks
#setwd("C:/Users/Kenji/Desktop/Recplot4/recplot_final_build/")
#source_python("recplot_database_carlos_genes.py")

#db <- "Nouveau.db"

#genes <- "genes.fa"

#samps <- assess_samples(db)[[1]][[1]]

#mags <- assess_MAGs(db, samps)[[1]][[1]]

#contigs <- unlist(get_contig_names(db, mags))

#recplot_raw <- extract_genes_MAG_for_R(db, samps, mags, 0.5, 70)

#recplot_list <- pydat_to_recplot_dat(recplot_raw, contigs)

#base <- recplot_list[[1]]
#bp_unit <- recplot_list[[2]]
#bp_div <- recplot_list[[3]]
#pos_max <- recplot_list[[4]]
#in_grp_min <- 95
#id_break <- 0.5
#width <- 1000
#linear <- 1

#create_static_plot(base, bp_unit, bp_div, pos_max, 95, 0.5, 1000, 1, showpeaks = F)


#parse_prodigal_genes(genes)

#add_genes_to_db(db, genes, "prodigal")

#add_sample(db, list("test_sam.sam"), "sam")


#source_python("recplot_database_carlos_genes.py")
#hide_the_print <- prepare_matrices_genes(db, "GS3.D25", 1000, 0.5, 70)

#a_gene_mag <- extract_genes_MAG_for_R(db, "C:\\Users\\Kenji\\Desktop\\Recplot4\\recplot_final_build\\test_sam _2.sam", "GS3.D25", 0.5, 70)
  
#fmt <- pydat_to_recplot_dat(a_gene_mag)

