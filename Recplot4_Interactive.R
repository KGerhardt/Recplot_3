#PACKAGES

library(enveomics.R)
library(data.table)
library(ggplot2)
library(shiny)
library(plotly)

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

handle_gff3 <- function(gff_file, adjust=LR[[1]]){
  
  gff_dat <- fread(cmd=paste("grep -v '#'", gff_file), sep = "\t")
  colnames(gff_dat) = c("contig","prog","feature","start","end","score","strand","phase","annotation")
  
  if(!all(gff_dat$contig %in% adjust$V1)){
    print("Sequence names do not match between lim/rec file and predicted genes file. Check to ensure the names match.")
  }
  
  gff_dat$pad = adjust$V2[match(gff_dat$contig, adjust$V1)]
  
  gff_dat[,start := start+pad-1]
  gff_dat[,end := end+pad-1]
  
  gff_dat$pad = NULL
  
  gff_dat$endpoint = ifelse(gff_dat$strand=="+", 1, -1)
  
  return(gff_dat)
  
}

#Primary interactive function 
recplot_suite <- function(prefix){
  
  LR <- get_lim_rec(prefix)
  
  bp_unit <- c("(bp)", "(Kbp)", "(Mbp)", "(Gbp)")[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
  bp_div <- c(1, 1e3, 1e6, 1e9)[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
  
  pos_max <- max(LR[[1]]$V3)/bp_div
  
  seq_breaks <- LR[[1]]$V3/bp_div

  group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
  
  contigs <- data.table(pos =  LR[[1]]$V2, contig = LR[[1]]$V1[findInterval(seq_breaks*bp_div, LR[[1]]$V2)])
  
    ui <- fluidPage(
      #Controls the inputs
      sidebarPanel(
                   submitButton("Update View", icon("refresh")),
                   numericInput("binSz", "Pos. Bin Width:", 1000, min = 1, max = max(LR[[1]]$V3, na.rm = T), step = 100),
                   numericInput("ani_step", "Pct. Identity Bin Width:", 0.5, min = 0.05, max = 5, step = 0.1),
                   numericInput("in_group_min", "In-Group Pct. Identity:", 95, min = 70, max = 100, step = 1),
                   numericInput("genome_region_min", "Start Pos in Genome.", 0, min = 0, max = pos_max-1),
                   numericInput("genome_region_max", "End Pos in Genome", pos_max, min = 1, max = pos_max),
                   selectInput("BP_hist_log", "Log Scale BP histogram?", choices = list("Linear" = 1, "Log" = 2), selected = 1),
                   selectInput("show_peaks", "Calculate Peaks?", choices = list("No" = 1, "Yes" = 2), selected = 1),
                   width = 2),
      
      mainPanel(
                tabsetPanel(type = "tabs",
                  tabPanel("Overview", column(9, 
                plotOutput("seq_depth_chart", height = "286px"),
                plotOutput("read_rec_plot", height = "574px"
                           )),
                column(3,
                      plotOutput("seq_depth_hist", height = "286px"),
                      plotOutput("bp_hist", height = "574px"))),
                tabPanel("Read Rec Hover", 
                         plotlyOutput("Plotly_read_rec", height = "860px")),
                tabPanel("Seq. Dep. Hover", 
                         plotlyOutput("Plotly_SD", height = "860px"))
                )
                )
                
    )
  
  

  
  server <- function(input,output){
    
    ranges <- reactiveValues(x = NULL, y = NULL)
    
    dat <- reactive({
          recplot_4_make_frame(LR[[1]],
                                 LR[[2]], 
                                 approx_break_size = input$binSz,
                                 ani_break_size = input$ani_step)
    })
    

      
      output$read_rec_plot<-renderPlot({
         base <- dat()
         
         base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
         
         widths <- c(base$Pos_In_Genome[base$ANI_Range == base$ANI_Range[1]], pos_max*bp_div)
         widths <- widths[-1]-widths[-length(widths)]
         widths <- rep(widths, nrow(base)/length(widths))
         base$Bin_Count <- base$Bin_Count*(input$binSz/widths)
         
         p <- ggplot(base,aes(x=Pos_In_Genome/bp_div, y=ANI_Range, fill=log10(Bin_Count)))+ 
              scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
              ylab("Percent Identity") +
              xlab(paste("Position in Genome", bp_unit)) +
              scale_y_continuous(expand = c(0, 0)) +
              scale_x_continuous(expand = c(0, 0)) +
              theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
              geom_raster()
         p <- p + annotate("rect", xmin = input$genome_region_min, xmax = input$genome_region_max, 
                           ymin = input$in_group_min, 
                           ymax = 100, fill = "darkblue", alpha = .15)
         p <- p + geom_vline(xintercept = seq_breaks[seq_breaks >= input$genome_region_min & seq_breaks <=  input$genome_region_max], col = "#AAAAAA40") 
         print(p)
      })
      
      output$seq_depth_chart <- renderPlot({
        base <- dat()
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
                
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        
        depth_data[,Pos_In_Genome := Pos_In_Genome/bp_div]
        
        nil_depth_data <- depth_data[count == 0]
        nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
        
        seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
        
        depth_data$count[depth_data$count == 0] <- NA
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min & Pos_In_Genome <= input$genome_region_max]

        p2 <- ggplot(depth_data, aes(x = Pos_In_Genome, y = count, colour=group_label, group = group_label))+
              geom_step(alpha = 0.75) +
              scale_y_continuous(trans = "log10", labels = scales::scientific) +
              scale_x_continuous(expand=c(0,0))+
              theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
              scale_color_manual(values = group.colors) +
              ylab("Depth")
        
        if(nrow(nil_depth_data) > 0){
          p2 <- p2 + geom_segment(data = nil_depth_data, aes(x = Pos_In_Genome, xend = Pos_In_Genome, y = 0, yend = seg_upper_bound, color = group_label, group = group_label))
        }
        

        print(p2)
        
      })
      
      output$seq_depth_hist <- renderPlot({
        base <- dat()
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
        
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        ddSave <- depth_data
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        
        
        depth_data <- depth_data[count > 0]
        
        seqdepth.lim <- range(c(depth_data$count[depth_data[,group_label == "depth.in"]], depth_data$count[depth_data[,group_label == "depth.out"]])) * c(1/2, 2)
        hist_binwidth <- (log10(seqdepth.lim[2]/2) - log10(seqdepth.lim[1] * 2))/199
        
        
        depth_data[,count := log10(count)]
        depth_data$group_label <-factor(depth_data$group_label, levels = c("depth.out", "depth.in"))
        depth_data <- depth_data[order(group_label),]
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]

        
        p4 <- ggplot(depth_data, aes(x = count, fill = group_label)) +
          geom_histogram(binwidth = hist_binwidth) +
          scale_fill_manual(values = group.colors) +
          scale_y_continuous(expand=c(0,0)) +
          theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), 
                axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(colour = "white"), axis.ticks.x = element_blank()) +
          xlab(label = element_blank()) +
          ylab(label = element_blank())
        
        #Peaks
        if(input$show_peaks == 2){
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
          peaks <- enve.recplot2.findPeaks(min_info)
          
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
          
          
        }
        
        p4 <- p4 + coord_flip()
        
        if(input$show_peaks == 2){
          o_max <- max(table(findInterval(depth_data$count[depth_data$group_label == "depth.in"], h.breaks)))*.6
          x_start = 0
          for(i in labels){
            p4 <- p4 + annotate("text", label = i, y = o_max, x = x_start, size = 5)
            x_start = x_start - .185
          }
          
        }
        
        print(p4)
        
      })
      
      output$bp_hist <- renderPlot({
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        
        
        bp_data <- base[,sum(Bin_Count), by = ANI_Range]
        
        if(input$BP_hist_log == 1){

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
        
        p4 <- p4 + annotate("rect", xmin = input$in_group_min, xmax = 100, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15) + 
              coord_flip()
      
        
        print(p4)
        
      })
      
      #TAB 2
      
      output$Plotly_read_rec <- renderPlotly({
        
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        #Widths cycle on contig change
        widths <- c(base$Pos_In_Genome[base$ANI_Range == base$ANI_Range[1]], pos_max*bp_div)
        widths <- widths[-1]-widths[-length(widths)]
        widths <- rep(widths, nrow(base)/length(widths))
        base$Bin_Count <- round(base$Bin_Count*(input$binSz/widths), 0)
        
        p <- ggplot(base, aes(x=Pos_In_Genome/bp_div, y=ANI_Range, fill=log10(Bin_Count), text = paste0("Contig: ", contigs$contig[(findInterval(base$Pos_In_Genome, contigs$pos))], 
                                                                                                        "\nPos. in Contig: ", (Pos_In_Genome + 1 - contigs$pos[(findInterval(base$Pos_In_Genome, contigs$pos))]), "-", (Pos_In_Genome - contigs$pos[(findInterval(base$Pos_In_Genome, contigs$pos))])+widths, 
                                                                                                        "\nNorm. Bin Count: ", base$Bin_Count, "\nPct. ID Range: ", base$ANI_Range-input$ani_step, "-", base$ANI_Range)))+ 
          scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
          ylab("Percent Identity") +
          xlab(paste("Position in Genome", bp_unit)) +
          scale_y_continuous(expand = c(0, 0)) +
          scale_x_continuous(expand = c(0, 0)) +
          theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
          annotate("rect", xmin = input$genome_region_min, xmax = input$genome_region_max, 
                     ymin = input$in_group_min, 
                     ymax = 100, fill = "darkblue", alpha = .15) + 
          geom_vline(xintercept = seq_breaks[seq_breaks >= input$genome_region_min & seq_breaks <=  input$genome_region_max], col = "#AAAAAA40") +
          geom_raster()

        p <- ggplotly(p, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "#EEF7FA") %>% 
          style(hoverinfo = "none", traces = c(1, 2)) 
        
        print(p)
        
      })
      
      #TAB 3
      
      output$Plotly_SD <- renderPlotly({
        
        base <- dat()
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
        
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        depth_data$width <- bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        depth_data[,Pos_In_Genome := Pos_In_Genome/bp_div]
        depth_data[,width := width/bp_div]
        
        depth_data <- depth_data[depth_data$count > 0,]
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min & Pos_In_Genome <= input$genome_region_max]
        
      
        p2 <- ggplot(depth_data, aes(x = Pos_In_Genome, y = count, colour=group_label, group = group_label,
                                     text = paste0("Contig: ", contigs$contig[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))], 
                                                  "\nPos. in Contig: ", (depth_data$Pos_In_Genome*bp_div + 1 - contigs$pos[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))]), "-",(depth_data$Pos_In_Genome*bp_div - contigs$pos[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))])+(depth_data$width*bp_div), 
                                                  "\nAvg. depth on interval: ", round(depth_data$count, 2))))+
          geom_step(alpha = 0.75) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(expand=c(0,0))+
          theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_text("Pos. In Genome")) +
          scale_color_manual(values = group.colors) +
          ylab("Depth")
          
        
        p2 <- ggplotly(p2, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "white")
        
        print(p2)
        
      })
      

      
  
  }
  
  shinyApp(ui, server)
  
}

recplot_suite_genes_summary <- function(prefix, genes_file){
    
    LR <- get_lim_rec(prefix)
    
    bp_unit <- c("(bp)", "(Kbp)", "(Mbp)", "(Gbp)")[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
    bp_div <- c(1, 1e3, 1e6, 1e9)[findInterval(log10(max(LR[[1]]$V3, na.rm = T)), c(0,3,6,9,12,Inf))]
    
    pos_max <- max(LR[[1]]$V3)/bp_div
    
    seq_breaks <- LR[[1]]$V3/bp_div
    
    group.colors <- c(depth.in = "darkblue", depth.out = "lightblue", depth.in.nil = "darkblue", depth.out.nil = "lightblue")
    
    contigs <- data.table(pos =  LR[[1]]$V2, contig = LR[[1]]$V1[findInterval(seq_breaks*bp_div, LR[[1]]$V2)])
    
    gene_data <- handle_gff3(genes_file, LR[[1]])
    gene_data[,start := start/bp_div]
    gene_data[,end := end/bp_div]
    
    ui <- fluidPage(
      #Controls the inputs
      sidebarPanel(
        submitButton("Update View", icon("refresh")),
        numericInput("binSz", "Pos. Bin Width:", 1000, min = 1, max = max(LR[[1]]$V3, na.rm = T), step = 100),
        numericInput("ani_step", "Pct. Identity Bin Width:", 0.5, min = 0.05, max = 5, step = 0.1),
        numericInput("in_group_min", "In-Group Pct. Identity:", 95, min = 70, max = 100, step = 1),
        numericInput("genome_region_min", "Start Pos in Genome.", 0, min = 0, max = pos_max-1),
        numericInput("genome_region_max", "End Pos in Genome", pos_max, min = 1, max = pos_max),
        selectInput("BP_hist_log", "Log Scale BP histogram?", choices = list("Linear" = 1, "Log" = 2), selected = 1),
        selectInput("show_peaks", "Calculate Peaks?", choices = list("No" = 1, "Yes" = 2), selected = 1),
        width = 2),
      
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Summary Report", tableOutput("sum_rep")), 
                    tabPanel("Overview", column(9, 
                                                plotOutput("seq_depth_chart", height = "286px"),
                                                plotOutput("read_rec_plot", height = "574px"
                                                )),
                             column(3,
                                    plotOutput("seq_depth_hist", height = "286px"),
                                    plotOutput("bp_hist", height = "574px"))),
                    tabPanel("Read Rec Hover", 
                             plotlyOutput("Plotly_read_rec", height = "860px")),
                    tabPanel("Seq. Dep. Hover", 
                             plotlyOutput("Plotly_SD", height = "860px"))
        )
      )
      
    )
    
    
    
    
    server <- function(input,output){
      
      ranges <- reactiveValues(x = NULL, y = NULL)
      
      dat <- reactive({
        recplot_4_make_frame(LR[[1]],
                             LR[[2]], 
                             approx_break_size = input$binSz,
                             ani_break_size = input$ani_step)
      })
      
      
      output$sum_rep <- renderTable({
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        tot_reads <- nrow(LR[[2]])
        in_reads <- sum(LR[[2]]$V3 > input$in_group_min)
        in_grp_bp <- sum((LR[[2]]$V2[LR[[2]]$V3 > input$in_group_min] - LR[[2]]$V1[LR[[2]]$V3 > input$in_group_min]))
        percent_in <- round(in_grp_bp/sum((LR[[2]]$V2-LR[[2]]$V1)), 4)*100
        ani <- round(sum((LR[[2]]$V2[LR[[2]]$V3 > input$in_group_min] - LR[[2]]$V1[LR[[2]]$V3 > input$in_group_min])*LR[[2]]$V3[LR[[2]]$V3 > input$in_group_min])/in_grp_bp, 2)
        
        
        stats = data.table(Summary = c("Total Reads",
                                       "Reads in In-Group",
                                       "Total In-Group BP",
                                       "In-Group BP Pct.",
                                       "ANIr"),
                           Result = c(tot_reads,
                                      in_reads,
                                      in_grp_bp,
                                      percent_in,
                                      ani))
        
        stats
        
        
      })
      
      
      output$read_rec_plot<-renderPlot({
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        widths <- c(base$Pos_In_Genome[base$ANI_Range == base$ANI_Range[1]], pos_max*bp_div)
        widths <- widths[-1]-widths[-length(widths)]
        widths <- rep(widths, nrow(base)/length(widths))
        base$Bin_Count <- base$Bin_Count*(input$binSz/widths)
        
        p <- ggplot(base,aes(x=Pos_In_Genome/bp_div, y=ANI_Range, fill=log10(Bin_Count)))+ 
          scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
          ylab("Percent Identity") +
          xlab(paste("Position in Genome", bp_unit)) +
          scale_y_continuous(expand = c(0, 0)) +
          scale_x_continuous(expand = c(0, 0)) +
          theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
          geom_raster()
        p <- p + annotate("rect", xmin = input$genome_region_min, xmax = input$genome_region_max, 
                          ymin = input$in_group_min, 
                          ymax = 100, fill = "darkblue", alpha = .15)
        p <- p + geom_vline(xintercept = seq_breaks[seq_breaks >= input$genome_region_min & seq_breaks <=  input$genome_region_max], col = "#AAAAAA40") 
        print(p)
      })
      
      output$seq_depth_chart <- renderPlot({
        base <- dat()
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
        
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        
        depth_data[,Pos_In_Genome := Pos_In_Genome/bp_div]
        
        nil_depth_data <- depth_data[count == 0]
        nil_depth_data$group_label <- ifelse(nil_depth_data$group_label == "depth.in", "depth.in.nil", "depth.out.nil")
        
        seg_upper_bound <- min(depth_data$count[depth_data$count > 0])
        
        depth_data$count[depth_data$count == 0] <- NA
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min & Pos_In_Genome <= input$genome_region_max]
        
        p2 <- ggplot(depth_data, aes(x = Pos_In_Genome, y = count, colour=group_label, group = group_label))+
          geom_step(alpha = 0.75) +
          scale_y_continuous(trans = "log10", labels = scales::scientific) +
          scale_x_continuous(expand=c(0,0))+
          theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
          scale_color_manual(values = group.colors) +
          ylab("Depth")
        
        if(nrow(nil_depth_data) > 0){
          p2 <- p2 + geom_segment(data = nil_depth_data, aes(x = Pos_In_Genome, xend = Pos_In_Genome, y = 0, yend = seg_upper_bound, color = group_label, group = group_label))
        }
        
        
        print(p2)
        
      })
      
      output$seq_depth_hist <- renderPlot({
        base <- dat()
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
        
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        ddSave <- depth_data
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        
        
        depth_data <- depth_data[count > 0]
        
        seqdepth.lim <- range(c(depth_data$count[depth_data[,group_label == "depth.in"]], depth_data$count[depth_data[,group_label == "depth.out"]])) * c(1/2, 2)
        hist_binwidth <- (log10(seqdepth.lim[2]/2) - log10(seqdepth.lim[1] * 2))/199
        
        
        depth_data[,count := log10(count)]
        depth_data$group_label <-factor(depth_data$group_label, levels = c("depth.out", "depth.in"))
        depth_data <- depth_data[order(group_label),]
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        
        p4 <- ggplot(depth_data, aes(x = count, fill = group_label)) +
          geom_histogram(binwidth = hist_binwidth) +
          scale_fill_manual(values = group.colors) +
          scale_y_continuous(expand=c(0,0)) +
          theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), 
                axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(colour = "white"), axis.ticks.x = element_blank()) +
          xlab(label = element_blank()) +
          ylab(label = element_blank())
        
        #Peaks
        if(input$show_peaks == 2){
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
          peaks <- enve.recplot2.findPeaks(min_info)
          
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
          
          
        }
        
        p4 <- p4 + coord_flip()
        
        if(input$show_peaks == 2){
          o_max <- max(table(findInterval(depth_data$count[depth_data$group_label == "depth.in"], h.breaks)))*.6
          x_start = 0
          for(i in labels){
            p4 <- p4 + annotate("text", label = i, y = o_max, x = x_start, size = 5)
            x_start = x_start - .185
          }
          
        }
        
        print(p4)
        
      })
      
      output$bp_hist <- renderPlot({
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        
        
        bp_data <- base[,sum(Bin_Count), by = ANI_Range]
        
        if(input$BP_hist_log == 1){
          
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
        
        p4 <- p4 + annotate("rect", xmin = input$in_group_min, xmax = 100, ymin = 0, ymax = Inf, fill = "darkblue", alpha = .15) + 
          coord_flip()
        
        
        print(p4)
        
      })
      
      #TAB 2
      
      output$Plotly_read_rec <- renderPlotly({
        
        base <- dat()
        
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        #Widths cycle on contig change
        widths <- c(base$Pos_In_Genome[base$ANI_Range == base$ANI_Range[1]], pos_max*bp_div)
        widths <- widths[-1]-widths[-length(widths)]
        widths <- rep(widths, nrow(base)/length(widths))
        base$Bin_Count <- round(base$Bin_Count*(input$binSz/widths), 0)
        
        p <- ggplot(base, aes(x=Pos_In_Genome/bp_div, y=ANI_Range, fill=log10(Bin_Count), text = paste0("Contig: ", contigs$contig[(findInterval(base$Pos_In_Genome, contigs$pos))], 
                                                                                                        "\nPos. in Contig: ", (Pos_In_Genome + 1 - contigs$pos[(findInterval(base$Pos_In_Genome, contigs$pos))]), "-", (Pos_In_Genome - contigs$pos[(findInterval(base$Pos_In_Genome, contigs$pos))])+widths, 
                                                                                                        "\nNorm. Bin Count: ", base$Bin_Count, "\nPct. ID Range: ", base$ANI_Range-input$ani_step, "-", base$ANI_Range)))+ 
          scale_fill_gradient(low = "white", high = "black",  na.value = "#EEF7FA")+
          ylab("Percent Identity") +
          xlab(paste("Position in Genome", bp_unit)) +
          scale_y_continuous(expand = c(0, 0)) +
          scale_x_continuous(expand = c(0, 0)) +
          theme(legend.position = "none", axis.line = element_line(colour = "black"), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
          annotate("rect", xmin = input$genome_region_min, xmax = input$genome_region_max, 
                   ymin = input$in_group_min, 
                   ymax = 100, fill = "darkblue", alpha = .15) + 
          geom_vline(xintercept = seq_breaks[seq_breaks >= input$genome_region_min & seq_breaks <=  input$genome_region_max], col = "#AAAAAA40") +
          geom_raster()
        
        p <- ggplotly(p, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "#EEF7FA") %>% 
          style(hoverinfo = "none", traces = c(1, 2))
        
        p2 <- ggplot(gene_data, aes(x = (start+end)/2, y = endpoint, text = annotation)) +
          geom_col(width=((gene_data$end-gene_data$start)/2)-1/bp_div, color = "grey75")+
          theme(rect = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.5)) +
          geom_hline(yintercept = 0, size = 1)
        
        p2 <- ggplotly(p2, dynamicTicks = T, tooltip = c("text")) %>%
          style(hoverinfo = "none", traces = 2)
        
        fullPlot <- subplot(p, p2, nrows=2, shareX = T, heights = c(.8, .2))
        
        print(fullPlot)
        
      })
      
      #TAB 3
      
      output$Plotly_SD <- renderPlotly({
        
        base <- dat()
        base <- base[Pos_In_Genome >= input$genome_region_min*bp_div & Pos_In_Genome <= input$genome_region_max*bp_div]
        
        base$group_label <- ifelse(base$ANI_Range-input$ani_step >= input$in_group_min, "depth.in", "depth.out")
        setkeyv(base, c("group_label", "Pos_In_Genome"))
        
        depth_data <- base[, sum(Bin_Count), by = key(base)]
        colnames(depth_data)[3] = "count"
        
        bin_ends <- data.table(start = unique(depth_data$Pos_In_Genome))
        bin_ends$width <- c(bin_ends$start[2:nrow(bin_ends)], pos_max * bp_div)+1-bin_ends$start
        
        depth_data$count <- depth_data$count/bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        depth_data$width <- bin_ends$width[match(depth_data$Pos_In_Genome, bin_ends$start)]
        
        depth_data[,Pos_In_Genome := Pos_In_Genome/bp_div]
        depth_data[,width := width/bp_div]
        
        depth_data <- depth_data[depth_data$count > 0,]
        
        depth_data <- depth_data[Pos_In_Genome >= input$genome_region_min & Pos_In_Genome <= input$genome_region_max]
        
        
        p2 <- ggplot(depth_data, aes(x = Pos_In_Genome, y = count, colour=group_label, group = group_label,
                                     text = paste0("Contig: ", contigs$contig[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))], 
                                                   "\nPos. in Contig: ", (depth_data$Pos_In_Genome*bp_div + 1 - contigs$pos[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))]), "-",(depth_data$Pos_In_Genome*bp_div - contigs$pos[(findInterval(depth_data$Pos_In_Genome*bp_div, contigs$pos))])+(depth_data$width*bp_div), 
                                                   "\nAvg. depth on interval: ", round(depth_data$count, 2))))+
          geom_step(alpha = 0.75) +
          scale_y_continuous(trans = "log10") +
          scale_x_continuous(expand=c(0,0))+
          theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), axis.title.x = element_text("Pos. In Genome")) +
          scale_color_manual(values = group.colors) +
          ylab("Depth")
        
        
        p2 <- ggplotly(p2, dynamicTicks = T, tooltip = c("text")) %>% 
          layout(plot_bgcolor = "white")
        
        p3 <- ggplot(gene_data, aes(x = (start+end)/2, y = endpoint, text = annotation)) +
          geom_col(width=((gene_data$end-gene_data$start)/2)-1/bp_div, color = "grey75")+
          theme(rect = element_blank(), axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0), limits = c(-1.5, 1.5)) +
          geom_hline(yintercept = 0, size = 1)
        
        p3 <- ggplotly(p3, dynamicTicks = T, tooltip = c("text")) %>%
          style(hoverinfo = "none", traces = 2)
        
        fullPlot <- subplot(p2, p3, nrows=2, shareX = T, heights = c(.8, .2))
        
        print(fullPlot)
        
      })
      
      
      
      
    }
    
    shinyApp(ui, server)
    
  }
  
