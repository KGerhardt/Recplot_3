#Options
options=commandArgs()

check_for_default <- integer(0)

directory <- which(grepl("-dir", options))
threads <- which(grepl("-t", options))
library_location <- which(grepl("-lib", options))

help <- which(grepl("-h", options))
help2 <- which(grepl("-help", options))
help3 <- which(grepl("--help", options))


if(any(!identical(help, check_for_default) | !identical(help2, check_for_default) | !identical(help3, check_for_default))){
  print("Usage Options:")
  
  print("-dir [directory containing lim/rec files to make recplot stats for]")
  print("-threads [INT: Number of threads to use for parallel. Only works with unix OS. Def. 1]")
  print("-lib [R library location] Only data.table is required; data.table and doParallel if threads are used.")
  
  quit(save = "no") 
}

directory <- options[directory + 1]

if(identical(directory, character(0))){
  directory <- getwd()
}

print(paste("Directory specified as:", directory))

setwd(directory)

if(identical(check_for_default, library_location)){
  print("Using the default location for R libraries during this session. This will likely fail in supercomputer environments.")
  flush.console()
  library_location <- .libPaths()
}else{
  library_location <- options[library_location + 1]
}

if(identical(check_for_default, threads)){
  print("Defaulting to 1 thread. Multiple threads only supported on unix environments.")
  flush.console()
  threads <- 1
}else{
  threads <- as.numeric(options[threads + 1])
}


#PACKAGES

library(data.table, lib.loc = library_location)

get_lim_rec <- function(prefix){
  rec <- fread(cmd = paste("grep -v '^#' ", prefix, ".rec", sep = ""), 
               sep = "\t", quote = "")
  lim <- fread(cmd = paste("grep -v '^#' ", prefix, ".lim", sep = ""), 
               sep = "\t", quote = "")
  return(list(lim, rec))
}

recplot_4_statistics <- function(prefixes, threads = 1){
  
  if(threads == 1){
    
    silent <- lapply(prefixes, function(i){
      
      LR <- get_lim_rec(i)
      
      #Prep the data for ANIr calculation
      LR[[2]][, anir_contribution := V3*(V2-V1+1)]
      LR[[2]]$contig <- LR[[1]]$V1[findInterval(LR[[2]]$V1, LR[[1]]$V3, left.open = T)+1]
      
      #Different data will have different cutoffs for in-group. This script is designed to just be run once and produce an entire scaffold for a user to work from.
      ANIr_in_group_levels <- seq(99.5, 90, -0.5)
      
      #Gets ANIr for whole genome, at each in-group level
      whole_genome_stats <- data.table(ANI_in_group_minimum = ANIr_in_group_levels, ANIr = unlist(lapply(ANIr_in_group_levels, function(x){
        
        sum(LR[[2]][V3>=x,anir_contribution])/sum(LR[[2]][V3>=x,V2-V1+1])
        
      })))
      
      whole_genome_stats <-  cbind(data.table("Whole_Genome"), data.table(do.call(rbind, list(whole_genome_stats$ANIr))))
      
      names(whole_genome_stats) = c("Contig", paste0("ANIr_min_", ANIr_in_group_levels))
      
      setkey(LR[[2]], "contig")
      
      #Spotting drops in ANIr in a given contig may be useful for spotting gene loss, etc.
      per_contig_stats <- data.table(contig = unique(LR[[2]]$contig) , do.call(cbind, lapply(ANIr_in_group_levels, function(x){
        
        (LR[[2]][V3>=x, sum(anir_contribution), by = key(LR[[2]])]$V1)/(LR[[2]][V3>=x, sum(V2-V1+1), by = key(LR[[2]])]$V1)
        
      })))
      
      LR[[2]][,anir_contribution:=NULL]
      LR[[2]][,contig:=NULL]
      
      names(per_contig_stats) = c("Contig", paste0("ANIr_min_", ANIr_in_group_levels))
      
      per_contig_stats <- rbind(whole_genome_stats, per_contig_stats)
      
      fwrite(per_contig_stats, paste0(i, "_ANIr_values.tsv"), sep = "\t")
      
      rm(whole_genome_stats, per_contig_stats)
      
      bin_cap <- max(LR[[1]]$V3)
      
      update_CD <- function(start, end){
        coverage_and_depth[start:end] <<- coverage_and_depth[start:end] + 1L
        return(NA)
      }
      
      ANIr_in_group_levels <- rev(ANIr_in_group_levels)
      LR[[2]] <- LR[[2]][V3>=ANIr_in_group_levels[1], .SD]
      
      
      LR[[2]][, id_bin := findInterval(V3, ANIr_in_group_levels)]
      setkey(LR[[2]], "id_bin")
      
      LR[[2]][, id := .I]
      
      bin_id <- unique(LR[[2]]$id_bin)
      
      LR[[2]] <- LR[[2]][, list(list(.SD)), by = "id_bin"]$V1
      
      coverage_and_depth <- rep(0L, max(LR[[1]]$V3))
      
      tad_indices <- seq(0, .25, by = 0.025)
      tad_indices <- round(length(coverage_and_depth) * tad_indices)
      
      
      
      TADs <- lapply(rev(LR[[2]]), function(j){
        
        j[, update_CD(V1, V2), by = id]
        
        tmp <- sort(coverage_and_depth)
        
        tad_values <- unlist(lapply(tad_indices, function(y){
          
          mean(tmp[y:(bin_cap-y)])
          
        }))
        
        return(tad_values)
        
      })
      
      rm(LR)
      
      tad_table <- data.table(cbind(rev(ANIr_in_group_levels), do.call(rbind, TADs)))
      
      names(tad_table) = c("Minimum_Percent_Identity", paste0("TAD_", (1-seq(0, .25, by = 0.025))*100))
      
      fwrite(tad_table, paste0(i, "_TAD_Stats.tsv"), sep = "\t")
      
      return(NA)
      
    })
    
  }else{
  
  library(doParallel, lib.loc = library_location)
    
  cl <- makeCluster(min(detectCores(), threads, length(prefixes)))
  
  clusterExport(cl = cl, varlist = "get_lim_rec", envir = environment())
  
  print(paste("Working on:", min(detectCores(), threads, length(prefixes)), "threads."))
  
  registerDoParallel(cl)
  
  foreach(i = prefixes, .packages = "data.table") %dopar% {
    
    LR <- get_lim_rec(i)
    
    #Prep the data for ANIr calculation
    LR[[2]][, anir_contribution := V3*(V2-V1+1)]
    LR[[2]]$contig <- LR[[1]]$V1[findInterval(LR[[2]]$V1, LR[[1]]$V3, left.open = T)+1]
    
    #Different data will have different cutoffs for in-group. This script is designed to just be run once and produce an entire scaffold for a user to work from.
    ANIr_in_group_levels <- seq(99.5, 90, -0.5)
    
    #Gets ANIr for whole genome, at each in-group level
    whole_genome_stats <- data.table(ANI_in_group_minimum = ANIr_in_group_levels, ANIr = unlist(lapply(ANIr_in_group_levels, function(x){
      
      sum(LR[[2]][V3>=x,anir_contribution])/sum(LR[[2]][V3>=x,V2-V1+1])
      
    })))
    
    whole_genome_stats <-  cbind(data.table("Whole_Genome"), data.table(do.call(rbind, list(whole_genome_stats$ANIr))))
    
    names(whole_genome_stats) = c("Contig", paste0("ANIr_min_", ANIr_in_group_levels))
    
    setkey(LR[[2]], "contig")
    
    #Spotting drops in ANIr in a given contig may be useful for spotting gene loss, etc.
    per_contig_stats <- data.table(contig = unique(LR[[2]]$contig) , do.call(cbind, lapply(ANIr_in_group_levels, function(x){
      
      (LR[[2]][V3>=x, sum(anir_contribution), by = key(LR[[2]])]$V1)/(LR[[2]][V3>=x, sum(V2-V1+1), by = key(LR[[2]])]$V1)
      
    })))
    
    LR[[2]][,anir_contribution:=NULL]
    LR[[2]][,contig:=NULL]
    
    names(per_contig_stats) = c("Contig", paste0("ANIr_min_", ANIr_in_group_levels))
    
    per_contig_stats <- rbind(whole_genome_stats, per_contig_stats)
    
    fwrite(per_contig_stats, paste0(i, "_ANIr_values.tsv"), sep = "\t")
    
    rm(whole_genome_stats, per_contig_stats)
    
    bin_cap <- max(LR[[1]]$V3)
    
    update_CD <- function(start, end){
      coverage_and_depth[start:end] <<- coverage_and_depth[start:end] + 1L
      return(NA)
    }
    
    ANIr_in_group_levels <- rev(ANIr_in_group_levels)
    LR[[2]] <- LR[[2]][V3>=ANIr_in_group_levels[1], .SD]
    
    
    LR[[2]][, id_bin := findInterval(V3, ANIr_in_group_levels)]
    setkey(LR[[2]], "id_bin")
    
    LR[[2]][, id := .I]
    
    bin_id <- unique(LR[[2]]$id_bin)
      
    LR[[2]] <- LR[[2]][, list(list(.SD)), by = "id_bin"]$V1
    
    coverage_and_depth <- rep(0L, max(LR[[1]]$V3))
    
    tad_indices <- seq(0, .25, by = 0.025)
    tad_indices <- round(length(coverage_and_depth) * tad_indices)
    
    
    
    TADs <- lapply(rev(LR[[2]]), function(j){
      
      j[, update_CD(V1, V2), by = id]
      
      tmp <- sort(coverage_and_depth)
      
      tad_values <- unlist(lapply(tad_indices, function(y){
        
        mean(tmp[y:(bin_cap-y)])
        
      }))
      
      return(tad_values)
      
    })

    rm(LR)
    
    tad_table <- data.table(cbind(rev(ANIr_in_group_levels), do.call(rbind, TADs)))
    
    names(tad_table) = c("Minimum_Percent_Identity", paste0("TAD_", (1-seq(0, .25, by = 0.025))*100))
    
    fwrite(tad_table, paste0(i, "_TAD_Stats.tsv"), sep = "\t")
    
    return(NA)
  }
  
  stopCluster(cl)
  
  }
  
}

prefixes <- list.files(pattern = ".lim")
prefixes <- substr(prefixes, 1, nchar(prefixes)-4)

recplot_4_statistics(prefixes, threads = threads)









