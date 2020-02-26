options=commandArgs()

check_for_default <- integer(0)

directory <- which(grepl("-dir", options))
id_lower <- which(grepl("-lb", options))
id_upper <- which(grepl("-hb", options))
threads <- which(grepl("-t", options))
step <- which(grepl("-id", options))
bin_width <- which(grepl("-w", options))
#I don't yet have the code in my lim/rec script for this to be anywhere but 3rd.
#pct_id_col <- which(grepl("-col", options))
separate <- which(grepl("-sep", options))
lin_hist <- which(grepl("-lin_hist", options))
in_group <- which(grepl("-in_grp", options))

help <- which(grepl("-h", options))
help2 <- which(grepl("-help", options))
help3 <- which(grepl("--help", options))

if(any(!identical(help, check_for_default) | !identical(help2, check_for_default) | !identical(help3, check_for_default))){
  print("Usage Options:")
 
  print("-dir [directory containing lim/rec files to make recplots for]")
  print("-lb [INT: Lower pct. ID boundary for recplot plotting area Def. 70]")
  print("-hb [INT: Uper pct. ID boundary for recplot plotting area. You shouldn't use this. Def. 100]")
  print("-threads [INT: Number of threads to use for parallel. Only works with unix OS. Def. 1]")
  print("-id [Decimal: Width of pct ID bins. Def. 0.5]")
  print("-w [INT: Approx. num. bp. for genome bins. Def. 1000]")
  print("-sep [T or F. Produce one PDF per lim/rec file if T; Produce 1 PDF with all recplots if F. Def. T.]")
  print("-lin_hist [T or F. Linear scale on BP histogram (lower right panel) if T, log10 scale if F]")
  print("-in_grp [Decimal: Lower pct. ID boundary for in-group (Dark blue) on recplot. Def. 95]")
  
  quit(save = "no") 
}

if(identical(directory, check_for_default)){
  print("You must specify a directory using the -dir option.")
  print("This should be the directory which contains all of the .lim and .rec files you wish to generate recplots from.")
  print("Every .lim and .rec file in the directory will have a recplot generated.")
  quit(save = "no")
}

directory <- options[directory + 1]
print(paste("Directory specified as:", directory))

setwd(directory)

library_location <- which(grepl("-lib", options))

if(identical(check_for_default, library_location)){
  print("Using the default location for R libraries during this session. This will likely fail in supercomputer environments.")
  flush.console()
  library_location <- .libPaths()
}else{
  library_location <- options[library_location + 1]
}

if(identical(check_for_default, id_lower)){
  print("Lower pct. ID defaulting to 70")
  flush.console()
  id_lower <- 70
}else{
  print("You really shouldn't be changing this one. Reads with pct. ID values outside the range of [lower bound,upper bound] may break the code.")
  flush.console()
  id_lower <- as.numeric(options[id_lower + 1])
}

if(identical(check_for_default, id_upper)){
  print("Upper pct. ID defaulting to 100.")
  flush.console()
  id_upper <- 100
}else{
  print("You really shouldn't be changing this one. Reads with pct. ID values outside the range of [lower bound,upper bound] may break the code.")
  flush.console()
  id_upper <- as.numeric(options[id_upper + 1])
}

if(identical(check_for_default, step)){
  print("Defaulting pct. ID bin height to 0.5")
  flush.console()
  step <- 0.5
}else{
  step <- as.numeric(options[step + 1])
}

if(identical(check_for_default, bin_width)){
  print("Defaulting approx. genome bin width to 1000 bp.")
  flush.console()
  bin_width <- 1e3
}else{
  bin_width <- as.numeric(options[bin_width + 1])
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
  print("Printing a separate PDF for each lim/rec prefix. If you want to group a directory into a single PDF, set this to false (-sep F).")
  flush.console()
  separate <- T
}else{
  separate <- (options[separate + 1])
  separate <- eval(parse(text="separate==T"))
}

if(identical(check_for_default, threads)){
  print("Defaulting to 1 thread. Multiple threads only supported on unix environments.")
  flush.console()
  threads <- 1
}else{
  threads <- as.numeric(options[threads + 1])
}

recplot_options <- data.frame(option = c("Pct ID Lower",
                                         "Pct ID Upper",
                                         "Pct ID Step",
                                         "Pct ID In-Group Min",
                                         "Approx. Bin Width",
                                         "BP Hist. Style"),
                              choice = as.character(c(id_lower,
                                         id_upper,
                                         step,
                                         in_group,
                                         bin_width,
                                         lin_hist)))

suppressMessages(library(enveomics.R, lib.loc = library_location))
suppressMessages(library(data.table, lib.loc = library_location))
suppressMessages(library(ggplot2, lib.loc = library_location))
suppressMessages(library(cowplot, lib.loc = library_location))

if(threads > 1){
  suppressMessages(library(doParallel, lib.loc = library_location))  
}



prefixes_lim <- list.files(pattern = ".lim")
prefixes_rec <- list.files(pattern = ".rec")

prefixes_lim <- substr(prefixes_lim, 1, nchar(prefixes_lim)-4)
prefixes_rec <- substr(prefixes_rec, 1, nchar(prefixes_rec)-4)

if(!(all(prefixes_lim %in% prefixes_rec) & all(prefixes_rec %in% prefixes_lim))){
  print("Some lim/rec files appear to lack a mate. Check to make sure these files are paired. Only lims/recs with mates will be plotted.")
  prefixes_lim <- prefixes_lim[prefixes_lim %in% prefixes_rec]
  prefixes_rec <- prefixes_rec[prefixes_rec %in% prefixes_lim]
}

cat("\n\n\n")

print("Plotting the following prefixes:")
for(i in prefixes_lim){
  print(i)
  flush.console()
}

print("Using the following options:")
print(recplot_options)
flush.console()

get_counts <- function(start_pos, end_pos, pos.breaks){  
  start_x_bin <- findInterval(start_pos, pos.breaks, left.open = T)
  len_set <- end_pos+1-start_pos
  
  output <- rep(0, length(pos.breaks)-1)
  while_counter <- output
  
  for(i in 1:length(start_pos)){
    start_bin <- start_x_bin[i]
    bin_end <- pos.breaks[start_bin+1]
    ends <- end_pos[i]
    len <- len_set[i]
    
    output[start_bin] <- output[start_bin] + len

    while(ends > bin_end){
      while_counter[start_bin] <- while_counter[start_bin] + 1
      
      len <- ends-floor(bin_end)
      
      output[start_bin] <- output[start_bin] - len + 1
      start_bin <- start_bin + 1
      bin_end <- pos.breaks[start_bin+1]
      
      output[start_bin] <- output[start_bin] + len
      
    }
    
  }
  output <- output-while_counter
  
  return(list(output))
}

get_matrix <- function(rec, pos.breaks, id.breaks, rec.idcol){
  rec$y_bin <- findInterval(unlist(rec[, ..rec.idcol]), id.breaks, left.open = T)
  rec_counts <- rec[, list(get_counts(V1, V2, pos.breaks)), by = list(y_bin)]
  counts <- lapply(1:(length(id.breaks)-1), function(x){return(rep(0, length(pos.breaks)-1))})
  for(i in 1:nrow(rec_counts)){
    counts[[rec_counts$y_bin[i]]] <- rec_counts$V1[[i]]
  }
  counts <- t(do.call(rbind, counts))
}

get_lim_rec <- function(prefix){
  rec <- fread(cmd = paste("grep -v '^#' ", prefix, ".rec", sep = ""), 
               sep = "\t", quote = "")
  lim <- fread(cmd = paste("grep -v '^#' ", prefix, ".lim", sep = ""), 
               sep = "\t", quote = "")
  
  if(any(rec$V2 >= max(lim$V3))){
    print("Oddity in rec file: read aligned at greater than end of last contig.")
    rec <- rec[V2<max(lim$V3),]
  }
  
  return(list(lim, rec))
}

recplot_4_make_frame <- function(lim,
                                 rec,
                                 approx_break_size = 1e3, 
                                 ani_break_size = 0.5, 
                                 ani_bounds = c(70, 100),
                                 ani_col = 3)
{
  
  # By default, we will split the genome into chunks of roughly equal size - here 1000 BP.
  # "Roughly equal" being that boundaries will start and stop at contig boundaries, and will divide a contig into (length(contig) int_div approx_break_size) + 1
  pos.breaks <- c(unlist(lapply(1:nrow(lim), function(i){
    
    #The first value is always lim$V2[i], and the last is lim$V3[i]. This results in the next index always being exactly 1 greater than the max of the previous
    # By removing this last one each time, we retain the break behavior, and close gaps.
    # Last, we append on the end of the sequence, since it gets removed through the above process.
    brk <- (round(seq(lim$V2[i], lim$V3[i], length.out = ((lim$V3[i]-lim$V2[i]+1)%/%approx_break_size)+1)))
    brk <- brk[1:(length(brk)-1)]
    
  })), max(lim$V3))
  
  
  id.breaks <- rev(unique(c(seq(ani_bounds[2], ani_bounds[1], by = -ani_break_size), 70)))
  
  read_rec_mat <- get_matrix(rec=rec, pos.breaks = pos.breaks, id.breaks = id.breaks, rec.idcol = ani_col)
  
  rownames(read_rec_mat) <- pos.breaks[1:(length(pos.breaks)-1)]
  colnames(read_rec_mat) <- id.breaks[2:(length(id.breaks))]
  
  
  result <- data.table(melt(read_rec_mat, keep.rownames = T))
  colnames(result) = c("Pos_In_Genome", "ANI_Range", "Bin_Count")
  
  return(result)
  
}


#A function which takes a prefix and settings as input and retruns a full recruitment plot.
recplot4_static_plot <- function(prefix,
                                 break_size = bin_width, 
                                 id_break = step, 
                                 bounds=c(id_lower, id_upper),
                                 id_col=3,
                                 in_grp_min = in_group,
                                 linear = lin_hist)
{
  
  LR <- get_lim_rec(prefix)
  base <- recplot_4_make_frame(LR[[1]], 
                                    LR[[2]],
                                    approx_break_size = break_size,
                                    ani_break_size = id_break,
                                    ani_bounds = bounds,
                                    ani_col = 3)
  
  bp_unit <- c("(bp)", "(Kbp)", "(Mbp)", "(Gbp)")[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
  bp_div <- c(1, 1e3, 1e6, 1e9)[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
  
  pos_max <- max(LR[[1]]$V3)/bp_div
  
  seq_breaks <- LR[[1]]$V3/bp_div
  
  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
  
  contigs <- data.table(pos =  LR[[1]]$V2, contig = LR[[1]]$V1[findInterval(seq_breaks*bp_div, LR[[1]]$V2)])
  
  
  #Recruitment plot main (bottom left) panel
  {
  widths <- c(base$Pos_In_Genome[base$ANI_Range == base$ANI_Range[1]], pos_max*bp_div)
  widths <- widths[-1]-widths[-length(widths)]
  widths <- rep(widths, nrow(base)/length(widths))
  base$Bin_Count <- base$Bin_Count*(break_size/widths)
  
  p <- ggplot(base,aes(x=Pos_In_Genome/bp_div, y=ANI_Range, fill=log10(Bin_Count)))+ 
    scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
    ylab("Percent Identity") +
    xlab(paste("Position in Genome", bp_unit)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    geom_raster()
  p <- p + annotate("rect", xmin = 0, xmax = pos_max, 
                    ymin = in_grp_min, 
                    ymax = 100, fill = "darkblue", alpha = .15)
  read_rec_plot <- p + geom_vline(xintercept = seq_breaks, col = "#AAAAAA40") 
  }
  
  #Seq. Depth chart (top left panel)
  {  
  base$group_label <- ifelse(base$ANI_Range-id_break >= in_grp_min, "depth.in", "depth.out")
  setkeyv(base, c("group_label", "Pos_In_Genome"))
  
  depth_data <- base[, sum(Bin_Count), by = key(base)]
  colnames(depth_data)[3] = "count"
  
  bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
  bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
  
  ddSave <- depth_data
  
  depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
  
  depth_data[,Pos_In_Genome := Pos_In_Genome/bp_div]
  
  nil_depth_data <- depth_data[count == 0]
  nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
  
  seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
  
  depth_data$count[depth_data$count == 0] <- NA
  
  seq_depth_chart <- ggplot(depth_data, aes(x = Pos_In_Genome, y = count, colour=group_label, group = group_label))+
    geom_step(alpha = 0.75) +
    scale_y_continuous(trans = "log10", labels = scales::scientific) +
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    scale_color_manual(values = group.colors) +
    ylab("Depth")
  
  if(nrow(nil_depth_data) > 0){
    seq_depth_chart <- seq_depth_chart + geom_segment(data = nil_depth_data, aes(x = Pos_In_Genome, xend = Pos_In_Genome, y = 0, yend = seg_upper_bound, color = group_label, group = group_label))
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
    min_info$pos.breaks = c(bin_ends$start, pos_max*bp_div)
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
    for(i in labels){
      p4 <- p4 + annotate("text", label = i, y = o_max, x = x_start, size = 3)
      x_start = x_start - .185
    }
    })
    
    seq_depth_hist <- p4
    
    rm(p4)
    
  }
  
  #bp counts histogram (bottom right panel)
  {
    
    bp_data <- base[,sum(Bin_Count), by = ANI_Range]
    
    if(linear){
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = ANI_Range)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0)) + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "#EEF7FA"), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+
        ylab("Base Pair Count by ANI")
    } else {
      
      
      p4 <- ggplot(data = bp_data, aes(y = V1, x = ANI_Range)) +
        geom_step() +
        scale_y_continuous(expand = c(0,0), trans = "log10") + 
        scale_x_continuous(expand = c(0,0)) +
        theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "#EEF7FA"), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())+
        ylab("Base Pair Count by ANI")
    }
    
    bp_count_hist <- p4 + annotate("rect", xmin = in_grp_min, xmax = 100, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15) + 
      coord_flip()
    
    rm(p4)
    
  }
  
  overall_plot <- plot_grid(seq_depth_chart, seq_depth_hist, read_rec_plot, bp_count_hist, align = "hv", ncol = 2, rel_widths = c(3, 1), rel_heights = c(1, 2.5))
  
  return(overall_plot)
  
}

#The function which controls the flow of the static plot creation.
recplot4_generate_static_plots <- function(prefixes, 
                                           break_size = 1e3, 
                                           id_break = 0.5, 
                                           bounds=c(70, 100),
                                           id_col=3,
                                           threads = 1,
                                           separate_outputs=T,
                                           in_grp_min = 95,
                                           linear = T)
{
  
  threads <- min(threads, length(prefixes))
  
  if(threads == 1){
    
    if(separate_outputs){
      
      for(i in prefixes_lim){
        
        recplot = ""
        
        try({
          
          recplot <- recplot4_static_plot(prefix = i,
                                   break_size = break_size,
                                   id_break = id_break,
                                   bounds = bounds,
                                   id_col=3,
                                   in_grp_min = in_grp_min,
                                   linear = linear)
          
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
        
        })

        if(length(recplot)==9){
          pdf(paste(i, ".pdf"), width = 17, height = 11)
          print(recplot)
          dev.off()
        }else{
          print(paste("There was an error plotting sample", i, "and there will be no recplot for this sample"))
        }

      }
      
    }else{
      
      pdf(paste("recruitment_plots.pdf"), width = 17, height = 11)
      
      for(i in prefixes_lim){
        
        recplot = ""
        
        try({
          
          recplot <- recplot4_static_plot(prefix = i,
                                          break_size = break_size,
                                          id_break = id_break,
                                          bounds = bounds,
                                          id_col=3,
                                          in_grp_min = in_grp_min,
                                          linear = linear)
          
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
          
        })
        
        if(length(recplot)==9){
          print(recplot)
        }else{
          print(paste("There was an error plotting sample", i, "and there will be no recplot for this sample"))
        }
     
        
    }
    
      dev.off() 
      
    }
    
    }else{
    
    print(paste("Operating on:", min(threads, detectCores()), "cores"))
    
    cl <- makeCluster(min(threads, detectCores(), length(prefixes_lim)))
    
    clusterExport(cl=cl, varlist = c("library_location", 
                                     "prefixes_lim",
                                     "recplot4_static_plot",
                                     "get_counts",
                                     "get_lim_rec",
                                     "get_matrix",
                                     "recplot_4_make_frame"))
                                    
    clusterEvalQ(cl=cl, expr = library(enveomics.R, lib.loc = library_location))
    clusterEvalQ(cl=cl, expr = library(data.table, lib.loc = library_location))
    clusterEvalQ(cl=cl, expr = library(ggplot2, lib.loc = library_location))
    clusterEvalQ(cl=cl, expr = library(cowplot, lib.loc = library_location))
    
    clusterExport(cl=cl, varlist = c("break_size",
                                     "id_break",
                                     "id_col",
                                     "bounds",
                                     "in_grp_min",
                                     "linear"), 
                  envir=environment())
    
    if(separate_outputs){
      
      registerDoParallel(cl)
      
      plots <- foreach(i = prefixes_lim) %dopar% {
        
        recplot = ""
        
        recplot <- try({
          
          recplot <- recplot4_static_plot(prefix = i,
                                          break_size = break_size,
                                          id_break = id_break,
                                          bounds = bounds,
                                          id_col=id_col,
                                          in_grp_min = in_grp_min,
                                          linear = linear)
          
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
          
        })
        
        return(recplot)
        
      }
      
      stopCluster(cl)
      
      for(i in 1:length(plots)){
        if(length(plots[[i]])==9){
          pdf(paste0(prefixes_lim[i],".pdf"), width = 17, height = 11)
          print(plots[[i]])
          dev.off()
        }
      }
      
      
    }else{

      
      
      registerDoParallel(cl)
      
      plots <- foreach(i = prefixes_lim) %dopar% {
        
        recplot = ""
        
        recplot <- try({
          
          recplot <- recplot4_static_plot(prefix = i,
                                          break_size = break_size,
                                          id_break = id_break,
                                          bounds = bounds,
                                          id_col=3,
                                          in_grp_min = in_grp_min,
                                          linear = linear)
          
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
          
        })
        
        return(recplot)
        
      }
      
      stopCluster(cl)
      
      pdf(paste("recruitment_plots.pdf"), width = 17, height = 11)
      
      for(i in 1:length(plots)){
        if(length(plots[[i]])==9){
          print(plots[[i]])
        }
      }
      
      dev.off()
      
    }
    
    
    
  }
  
  
  
}



recplot4_generate_static_plots(prefixes_lim,
                               break_size = bin_width,
                               id_break = step,
                               bounds = c(id_lower, id_upper),
                               id_col = 3,
                               threads = threads,
                               separate_outputs = separate,
                               in_grp_min=in_group,
                               linear = lin_hist)




