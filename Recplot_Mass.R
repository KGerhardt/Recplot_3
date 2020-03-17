options=commandArgs()

check_for_default <- integer(0)

super_rec <- which(grepl("-input", options))
threads <- which(grepl("-threads", options))
separate <- which(grepl("-combine", options))
lin_hist <- which(grepl("-lin_hist", options))
in_group <- which(grepl("-in_grp", options))

help <- which(grepl("-h", options))
help2 <- which(grepl("-help", options))
help3 <- which(grepl("--help", options))

if(any(!identical(help, check_for_default) | !identical(help2, check_for_default) | !identical(help3, check_for_default))){
  print("Usage Options:")
  
  print("-input [recruitment matrices file to plot]")
  print("-threads [INT: Number of threads to use for parallel. Only works with unix OS. Def. 1]")
  print("-combine; Produce 1 PDF with all recplots if present. Otherwise 1 PDF per MAG")
  print("-lin_hist [T or F. Linear scale on BP histogram (lower right panel) if T, log10 scale if F]")
  print("-in_grp [Decimal: Lower pct. ID boundary for in-group (Dark blue) on recplot. Def. 95]")
  
  quit(save = "no") 
}

if(identical(super_rec, check_for_default)){
  print("You must specify a recruitment matrices file using the -input option.")
  print("This should be the [prefix]_recruitment_matrcies.tsv produced by recplot_matrix.py")
  print("All MAGs in this file will be plotted.")
  print("Exiting program without plotting.")
  quit(save = "no")
}

super_rec <- options[super_rec+ 1]

if(!file.exists(super_rec)){
  print(paste("Recruitment matrices file:", super_rec, "not found.")) 
  print("Are you sure the name is correct, and you're in the correct directory?")
  print("Exiting program without plotting.")
}

print(paste("Input specified as:", super_rec))

library_location <- which(grepl("-lib", options))

if(identical(check_for_default, library_location)){
  print("Using the default location for R libraries during this session. This will likely fail in supercomputer environments.")
  flush.console()
  library_location <- .libPaths()
}else{
  library_location <- options[library_location + 1]
}

if(identical(check_for_default, in_group)){
  print("Defaulting in-group pct ID min to 95.")
  flush.console()
  in_group <- 95
}else{
  in_group <- as.numeric(options[in_group + 1])
}

if(identical(check_for_default, lin_hist)){
  print("Base pair count histogram (lower right panel) will be plotted with linear scale. If you wanted log scale, set this to false (-lin_hist F).")
  flush.console()
  lin_hist <- T
}else{
  lin_hist <- (options[lin_hist + 1])
  lin_hist <- eval(parse(text="lin_hist==T"))
}

if(identical(check_for_default, separate)){
  print("Printing a separate PDF for each MAG. If you want to group MAGs into a single PDF, add -combine flag.")
  flush.console()
  separate <- T
}else{
  print("Printing a single PDF for all MAGs")
  flush.console()
  separate <- F
}

if(identical(check_for_default, threads)){
  print("Defaulting to 1 thread. Multiple threads only supported on unix environments.")
  flush.console()
  threads <- 1
}else{
  threads <- as.numeric(options[threads + 1])
}


suppressMessages(library(enveomics.R, lib.loc = library_location))
suppressMessages(library(data.table, lib.loc = library_location))
suppressMessages(library(ggplot2, lib.loc = library_location))
suppressMessages(library(cowplot, lib.loc = library_location))

if(threads > 1){
  suppressMessages(library(doParallel, lib.loc = library_location))  
}

read_super_rec <- function(file = super_rec){
  
  #Read data
  matrices <- fread(file, sep = "\t")
  
  #Sort by MAG
  setkey(matrices, "MAG_name")
  matrices_names <- unique(matrices$MAG_name)
  
  
  #Essentially gets the data... needs reshaping
  matrices <- matrices[, list(list(.SD)), by = key(matrices)]$V1
  names(matrices) = matrices_names
  
  return(matrices)
  
}

recplot_4_static <- function(matrices_object, in_grp_min = 95, linear = T){
  
  ends <- matrices_object[,max(End), by = Contig_name]
  
  bp_unit <- c("(bp)", "(Kbp)", "(Mbp)", "(Gbp)")[findInterval(log10(max(cumsum(ends$V1))), c(0,3,6,9,12,Inf))]
  bp_div <- c(1, 1e3, 1e6, 1e9)[findInterval(log10(max(cumsum(ends$V1))), c(0,3,6,9,12,Inf))]
  
  ends[, V1 := V1/bp_div]
  
  pos_max <- max(cumsum(ends$V1))
  
  ID_bins <- colnames(matrices_object)[4:ncol(matrices_object)]
  
  id_break <- diff(as.numeric(unlist(strsplit(ID_bins[1], "-"))))
  
  colnames(matrices_object)[4:ncol(matrices_object)] = unlist(strsplit(ID_bins, "-"))[c(F,T)]
  
  matrices_object$seq_pos = seq(1/bp_div, pos_max, length.out = nrow(matrices_object))

  base <- melt.data.table(matrices_object, id.vars = c("Contig_name", "Start", "End", "seq_pos"))
  colnames(base)[5:6] = c("Pct_ID_bin", "bp_count")
  
  base$Pct_ID_bin <- as.numeric(levels(base$Pct_ID_bin))[base$Pct_ID_bin]
  
  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
  
  #Read Rec. Plot (bottom left panel)
  {  
  p <- ggplot(base, aes(x = seq_pos, y = Pct_ID_bin, fill=log10(bp_count)))+ 
    scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
    ylab("Percent Identity") +
    xlab(paste("Position in Genome", bp_unit)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "none", axis.line = element_line(colour = "black")) +
    geom_raster()
  
  p <- p + annotate("rect", xmin = 0, xmax = pos_max, 
                    ymin = in_grp_min, 
                    ymax = 100, fill = "darkblue", alpha = .15)
  
  read_rec_plot <- p + geom_vline(xintercept = cumsum(ends$V1)[-nrow(ends)], col = "#AAAAAA40")
  }
  
  #Seq. Depth chart (top left panel)
  {  
    base$group_label <- ifelse(base$Pct_ID_bin-id_break >= in_grp_min, "depth.in", "depth.out")
    setkeyv(base, c("group_label", "seq_pos"))
    
    
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
      theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
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
      theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), 
            axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(colour = "white"), axis.ticks.x = element_blank()) +
      xlab(label = element_blank()) +
      ylab(label = element_blank())
    

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
    min_info$pos.breaks = c(cumsum(ends$V1), pos_max*bp_div)
    min_info$pos.counts.in = ddSave$count[ddSave$group_label == "depth.in"]
    
    #This is finicky and just doesn't work like half the time.
    
    try({peaks <- enve.recplot2.findPeaks(min_info)
    
    dpt <- signif(as.numeric(lapply(peaks, function(x) x$seq.depth)), 2)
    frx <- signif(100 * as.numeric(lapply(peaks,function(x) ifelse(length(x$values) == 0, x$n.hat, length(x$values))/x$n.total)), 2)
    
    #No bad peaks, please
    peaks <- peaks[frx >= 5]
    dpt <- dpt[frx >= 5]
    frx <- frx[frx >= 5]
    
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
    
    p4 <- p4 + coord_flip()
    
    
    try({
      o_max <- max(table(findInterval(depth_data$count[depth_data$group_label == "depth.in"], h.breaks)))*.6
      x_start = 0
      
      if(length(labels) > 0){
      for(i in labels){
        p4 <- p4 + annotate("text", label = i, y = o_max, x = x_start, size = 3)
        x_start = x_start - .185
      }
      }
    })
    
    seq_depth_hist <- p4
    
    rm(p4)
    
  }
  
  #bp counts histogram (bottom right panel)
  {
    
    bp_data <- base[,sum(bp_count), by = Pct_ID_bin]
    
    if(linear){
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0)) + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "#EEF7FA"), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+
        ylab("Base Pair Count by ANI")
    } else {
      
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = Pct_ID_bin)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0), trans = "log10") + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "#EEF7FA"), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+
        ylab("Base Pair Count by ANI")
    }
    
    bp_count_hist <- p4 + annotate("rect", xmin = in_grp_min, xmax = 100, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15) + coord_flip()
    
    rm(p4)
    
  }
  
  overall_plot <- plot_grid(seq_depth_chart, seq_depth_hist, read_rec_plot, bp_count_hist, align = "hv", ncol = 2, rel_widths = c(3, 1), rel_heights = c(1, 2.5))
  
  return(overall_plot)
  
}

generate_static_plots <- function(super_rec_matrix, minimum_id, linear, thd, separ){
  
  if(separ){
    
    if(thd == 1){
      
      for(i in names(super_rec_matrix)){

        recplot <- recplot_4_static(super_rec_matrix[[i]], in_grp_min = minimum_id, linear = linear)
        
        title <- ggdraw() + 
          draw_label(
            i,
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )
        
        recplot <- plot_grid(
          title, recplot,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1)
        )
        
        pdf(paste0(i, "_recplot.pdf"), width = 16, height = 10)
        print(recplot)
        dev.off()
        
      }
      
    }else{

      cl <- makeCluster(min(thd, detectCores(), length(super_rec_matrix)))
      clusterExport(cl, varlist = c("recplot_4_static", "super_rec_matrix"), envir = environment())
      
      registerDoParallel(cl)
      
      silent <- foreach(i = names(super_rec_matrix), .packages= c("data.table", "ggplot2", "cowplot", "enveomics.R")) %dopar% {
        
        recplot <- recplot_4_static(super_rec_matrix[[i]], in_grp_min = minimum_id, linear = linear)
        
        title <- ggdraw() + 
          draw_label(
            i,
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )
        recplot <- plot_grid(
          title, recplot,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1)
        )
        
        pdf(paste0(i, "_recplot.pdf"), width = 16, height = 10)
        print(recplot)
        dev.off()
        
        return(NA)
        
      }
      
      stopCluster(cl)
      
    }
    
  }else{
    
    if(thd == 1){
      
      pref <- substr(super_rec, 1, nchar(super_rec)-25)
      
      pdf(paste0(pref, "_static_recruitment_plots.pdf"), height = 10, width = 16)
      
      for(i in names(super_rec_matrix)){
        
        recplot <- recplot_4_static(super_rec_matrix[[i]], in_grp_min = minimum_id, linear = linear)
        
        title <- ggdraw() + 
          draw_label(
            i,
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )
        recplot <- plot_grid(
          title, recplot,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1)
        )
        
        print(recplot)
        
      }
      
      dev.off()
      
    }else{
      
      grp_count <- min(thd, detectCores(), length(super_rec_matrix))
      
      cl <- makeCluster(grp_count)
      clusterExport(cl, varlist = c("recplot_4_static", "super_rec_matrix"), envir = environment())
      
      registerDoParallel(cl)
      
      pref <- substr(super_rec, 1, nchar(super_rec)-25)
      
      mag_names <- names(super_rec_matrix)
      
      groups <- (1:length(mag_names)) %/% grp_count
      
      pdf(paste0(pref, "_static_recruitment_plots.pdf"), height = 10, width = 16)
      
      for(j in unique(groups)){
      
      current_names <- mag_names[groups == j]
      
      a_few_plots <- foreach(i = current_names, .packages= c("data.table", "ggplot2", "cowplot", "enveomics.R")) %dopar% {
        
        recplot <- recplot_4_static(super_rec_matrix[[i]], in_grp_min = minimum_id, linear = linear)
        
        title <- ggdraw() + 
          draw_label(
            i,
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )
        recplot <- plot_grid(
          title, recplot,
          ncol = 1,
          # rel_heights values control vertical title margins
          rel_heights = c(0.1, 1)
        )
        
        return(recplot)
        
      }
      
      for(k in a_few_plots){
        print(k)
      }
      
      
      }
      
      stopCluster(cl)
     
      dev.off()
       
    }
    
  }
  
}


#Remove me later
mat <- read_super_rec(file = super_rec)

generate_static_plots(super_rec_matrix = mat, minimum_id = in_group, linear = lin_hist, thd = threads, separ = separate)



