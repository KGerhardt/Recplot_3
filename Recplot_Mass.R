#Package check and install.
{
  check <- suppressWarnings(suppressMessages(require(reticulate)))
  
  if(!check){
    install.packages("reticulate")
    library(reticulate)
  }
  
  check <- suppressWarnings(suppressMessages(require(ggplot2)))
  
  if(!check){
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  check <- suppressWarnings(suppressMessages(require(shiny)))
  
  if(!check){
    install.packages("shiny")
    library(shiny)
  }
  
  check <- suppressWarnings(suppressMessages(require(data.table)))
  
  if(!check){
    install.packages("data.table")
    library(data.table)
  }
  
  check <- suppressWarnings(suppressMessages(require(plotly)))
  
  if(!check){
    install.packages("plotly")
    library(plotly)
  }
  
  check <- suppressWarnings(suppressMessages(require(cowplot)))
  
  if(!check){
    install.packages("cowplot")
    library(cowplot)
  }
  
  check <- suppressWarnings(suppressMessages(require(enveomics.R)))
  
  if(!check){
    install.packages("enveomics.R")
    library(enveomics.R)
  }
  
  check <- suppressWarnings(suppressMessages(require(shinyBS)))
  
  if(!check){
    install.packages("shinyBS")
    library(shinyBS)
  }
  
  check <- suppressWarnings(suppressMessages(require(hms)))
  
  if(!check){
    install.packages("hms")
    library(hms)
  }
  
  check <- suppressWarnings(suppressMessages(require(easycsv)))
  
  if(!check){
    install.packages("easycsv")
    library(easycsv)
  }
  
  check <- suppressWarnings(suppressMessages(require(shinyalert)))
  
  if(!check){
    install.packages("shinyalert")
    library(shinyalert)
  }
  
  check <- suppressWarnings(suppressMessages(require(htmlwidgets)))
  
  if(!check){
    install.packages("htmlwidgets")
    library(htmlwidgets)
  }
  
  
}

#Helper functions
{
  #This will download whatever the current python script is. You have to run it before landing page.
  get_python <- function(){
    
    source_python("https://raw.githubusercontent.com/KGerhardt/Recplot_4/master/Recplot.py", envir = globalenv())
    
  }
  #Checks for necessary setup steps and makes the requisite installs as needed.
  prepare_environment <- function(){
    
    cat("Checking for Miniconda and installing if necessary...\n")
    try({
      install_miniconda()
    })
    
    #Checking for first-time use of recplots
    if(!"recruitment_plots" %in% conda_list()$name){
      cat("Creating Miniconda environment: 'recruitment_plots'\n")
      conda_create(envname = "recruitment_plots")
    }
    
    use_miniconda(condaenv = "recruitment_plots", required = T)
    
    if(get_sys() != "Windows"){
      if(py_module_available("pysam")){
        cat("Attempting to install pysam to recruitment_plots... ")
        try({
          py_install(packages = "pysam", envname = "recruitment_plots", pip = T)
          cat("Done!\n")
        }) 
      }else{
        cat("Pysam already installed. You probably shouldn't be seeing this warning. Did you call prepare_environment() twice?\n")
      }
      
      get_python()
      
    }
    
  }
  
  #Prepares the background miniconda env if necessary; otherwise, sets the environment and loads the python script functions
  initiate <- function(){
    cat("Initiating recruitment plot environment. Please wait a moment. ")
    
    tryCatch({
      
      use_miniconda(condaenv = "recruitment_plots", required = T )
      get_python()
      
    }, error = function(cond){
      
      cat("\nPerforming first-time setup. Wait a moment, please.\n")
      
      prepare_environment()
      
      return("Environment prepared.")
      
    }  )
  }
  
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
  
  create_static_plot <- function(base, bp_unit, bp_div, pos_max, in_grp_min, id_break, width, linear, showpeaks, ends, trunc_behavior = "ends", trunc_degree = as.integer(75), ...){
    
    group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
    
    #Sets any starts < trunc degree to trunc_degree
    base <- base[Start < trunc_degree, Start := trunc_degree]
    #Selects the bins at the end of each contig and subtracts trunc degree from it
    base[base[, End > (max(End)-trunc_degree), by = contig]$V1, End := (End - trunc_degree),]
    #If the final bin was too small, removes it.
    base <- base[Start <= End,]
    
    #Allows for count normalization by bin width across all bins
    norm_factor <- min(base$End-base$Start) + 1
    
    #Lower left panel
    {    
      p <- ggplot(base, aes(x = seq_pos, y = Pct_ID_bin, fill=log10((bp_count * (norm_factor/(End-Start+1))))))+ 
        scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
        ylab("Percent Identity") +
        xlab(paste("Position in Genome", bp_unit)) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, pos_max/bp_div), breaks = scales::pretty_breaks(n = 10)) +
        theme(legend.position = "none", 
              axis.line = element_line(colour = "black"),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 14)) +
        geom_raster()
      
      
      p <- p + annotate("rect", xmin = 0, xmax = pos_max/bp_div, 
                        ymin = in_grp_min, 
                        ymax = 100, fill = "darkblue", alpha = .15)
      
      read_rec_plot <- p + geom_vline(xintercept = ends$V1/bp_div[-nrow(ends)], col = "#AAAAAA40")
    }
    
    base[, group_label := ifelse(base$Pct_ID_bin-id_break >= in_grp_min, "depth.in", "depth.out")]
    setkeyv(base, c("group_label", "seq_pos"))
    
    #upper left panel
    {
      depth_data <- base[, sum(bp_count/(End-Start+1), na.rm = T), by = key(base)]
      colnames(depth_data)[3] = "count"
      
      ddSave <- base[, sum(bp_count, na.rm = T), by = key(base)]
      
      nil_depth_data <- depth_data[count == 0]
      nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
      
      seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
      
      depth_data$count[depth_data$count == 0] <- NA
      
      seq_depth_chart <- ggplot(depth_data, aes(x = seq_pos, y = count, colour=group_label, group = group_label))+
        geom_step(alpha = 0.75) +
        scale_y_continuous(trans = "log10", labels = scales::scientific) +
        scale_x_continuous(expand=c(0,0), limits = c(0, pos_max/bp_div), breaks = scales::pretty_breaks(n = 10))+
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
        
        if(nrow(ends) > 1){
          ends[, adjust := c(-1, V1[1:nrow(ends)-1])+1]
        }else{
          ends[, adjust := 0]
        }
        
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
      
      bp_data <- base[,sum(bp_count, na.rm = T), by = Pct_ID_bin]
      
      if(linear == 1){
        
        p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin-(id_break))) +
          geom_step() +
          scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n = 3)) + 
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
          ylab("Base Pair Count by % ID")
        
      } else {
        
        p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin-(id_break))) +
          geom_step() +
          scale_y_continuous(expand = c(0,0), trans = "log10", breaks = scales::log_breaks(n = 4)) + 
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
          ylab("Base Pair Count by % ID")
        
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
    
    pretty_data[, gene_start := unname(unlist(lapply(prodigal_gene_mess, function(x){
      
      return(x[[2]])
      
    }))) ]
    
    pretty_data[, gene_end := unname(unlist(lapply(prodigal_gene_mess, function(x){
      
      return(x[[3]])
      
    }))) ]
    
    pretty_data[, strand  := unname(unlist(lapply(prodigal_gene_mess, function(x){
      
      return(x[[4]])
      
    }))) ]
    pretty_data[, annotation  := unname(unlist(lapply(prodigal_gene_mess, function(x){
      
      return(x[[5]])
      
    }))) ]
    
    
    return(pretty_data)
    
    
  }
  
}

#dev
#options <- c("-database", "C:/Users/Kenji/Desktop/Recplot4/recplot_final_build/mass_rec/with_mags.db", "-res", "low")
#setwd("C:/Users/Kenji/Desktop/Recplot4/recplot_final_build/mass_rec")

#Commandline Args
{
  options <- commandArgs()
  
  print("Options selected")
  
  print(options)
  
  default_check <- integer(0)
  
  database <- which(grepl("-database", options))
  
  #Allow user to supply an excel file of mags/genomes to plot from the DB
  specified_mags <- which(grepl("-plot_these",options))
  #If they don't, assume the user wants to do it all
  do_all <- F
  
  #Default resolution is the medium 1000 BP window, or if genes, then genes and long IGR
  #Low res: 5000 bp, genes only
  #high res: 333 bp, all regions
  resolution <- which(grepl("-res", options))
  
  in_group <- which(grepl("-species_boundary",options))
  
  if(identical(default_check, database)){
    
    print("Cannot proceed without a database. Exiting.")
    quit(save = "no")
        
  }else{
    
    database <- options[database + 1]
  
  }
 
  if(identical(default_check, specified_mags)){
    
    do_all <- T
    
  }else{
    
    default_check <- read_excel(options[specified_mags + 1])
    
  }
  
  if(identical(default_check, resolution)){
    
    resolution <- "med"
    
  }else{
    
    resolution <- options[resolution + 1]
    
  }
  
  if(identical(default_check, in_group)){
    
    in_group <- 95
  
  }else{
    
    in_group <- as.numeric(options[in_group + 1])
    
  }
   
}

multi_rec <- function(plot_all, database, res = "med", selected, group = 95, linear = 2, threads = 1){

if(res == "low"){
  wid <- 3000
  step <- 1
}
if(res == "med"){
  wid <- 1000
  step <- 0.5
}
if(res == "high"){
  wid <- 500
  step <- 0.25
}
if(res == "unreasonable"){
    wid <- 250
    step <- 0.125
} 
   
initiate()

samps <- unlist(assess_samples(database))

mags <- lapply(samps, function(x){
  
  unlist(assess_MAGs(database, x))
  
})

names(mags) = samps

pdf(paste0(database, ".pdf"), height = 11, width = 17)

for(sample in samps){
  
  for(genome in mags[[sample]]){
    
    contigs <- get_contig_names(database, genome)
    base <- pydat_to_recplot_dat(extract_MAG_for_R(database, sample, genome, wid, step, 70), contigs)
    
    bp_unit <- base[[2]]
    bp_div <- base[[3]]
    pos_max <- base[[4]]
    base <- base[[1]]
    
    ending <- base[, max(End), by = contig]
    ending[, V1 := cumsum(V1) - 1 + 1:nrow(ending)]
    
    static_plot <- create_static_plot(base = base,
                                      bp_unit = bp_unit,
                                      bp_div = bp_div,
                                      pos_max = pos_max,
                                      in_grp_min = group,
                                      id_break = step,
                                      width = wid,
                                      linear = linear,
                                      showpeaks= T,
                                      ends = ending
    )
    
    title <- ggdraw() + 
      draw_label(
        paste("Read Recruitment Plot for sample:", sample, "MAG:", genome),
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
    
    print(static_plot)
    

  }
  
}

dev.off()

return(NA)

}

multi_rec(database = database, res = resolution)
