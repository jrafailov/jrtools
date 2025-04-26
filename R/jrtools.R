#' @name hunt
#' @title grab your bows and arrows cuz it's hunting season
#'
#' @description
#' Hunt function to assess local memory usage (mskilab server)
#'
#' @param user user you want to "hunt down"
#' @return A data table of user and memory usage in GB.
#' @import data.table
#' @import ps
#' @import tidyverse
#' @import data.tree
#' @export
#'
#' @examples
#' hunt()
hunt = function(huntdown = NULL){
  #if(tree & is.null(huntdown)) {stop("user must be specified for PID tree")}

  ps = ps() %>%
    data.table() %>%
    transmute(pid,
           ppid,
           task = name,
           username,
           status,
           resident_mem_gb = rss,
           virtual_mem_gb = vms
           )

  ps.tree = ps %>%
    FromDataFrameNetwork()

  print(ps.tree, "task", "username", "status", "resident_mem_gb", "virtual_mem_gb")

  ps_all <- ps[,.(total_memory_GB = round(sum(resident_mem_gb / (1024 * 1024), na.rm = TRUE) / 1024, 4)), by = username][order(-total_memory_GB)]

   if(!is.null(huntdown)){
    ps_user = ps[grep(huntdown, username),.(username, pid, ppid, status, resident_mem_gb = round(resident_mem_gb / (1024^3), 2), virtual_mem_gb = round(virtual_mem_gb / (1024^3), 2))][order(-resident_mem_gb)]
    ps_user_tasks = ps[pid %in% ps_user$ppid | pid %in% ps_user$pid]
  }
  #if(tree){
   # ps_ppid = unique(ps_user$ppid)
  #}
  if(is.null(huntdown)){
    return(ps_all)
  } else{
    return(ps_user)
  }
}


#' @name concat_file_paths
#' @title let's parallelize concat some data.tables
#'
#' @description
#' function to parallelize and concatenate a list of vector file paths that contain rds that hold only data.tables
#'
#' @param filepaths vector of rds filepaths
#' @param cores number of cores to mclapply over lol
#'
#' @return data.table
#' @export
concat_file_paths = function(filepaths, cores = 1){
  filepaths = filepaths[!is.na(filepaths)]
  mclapply(filepaths, function(path){
    #print(path)
    data = readRDS(path)
    if(!is.data.table(data)){
      data <- data %>% as.data.table()
    }
    return(data)
  }, mc.cores = cores) %>% rbindlist() -> to_return
  return(to_return)
}

#' @name concat_file_paths
#' @title let's parallelize concat some data.tables
#'
#' @description
#' function to parallelize and concatenate a list of vector file paths that contain rds that hold only data.tables
#'
#' @param filepaths vector of rds filepaths
#' @param cores number of cores to mclapply over lol
#'
#' @return data.table
#' @export
concat_file_paths = function(filepaths, cores = 1){
  filepaths = filepaths[!is.na(filepaths)]
  mclapply(filepaths, function(path){
    #print(path)
    data = readRDS(path)
    if(!is.data.table(data)){
      data <- data %>% as.data.table()
    }
    return(data)
  }, mc.cores = cores) %>% rbindlist() -> to_return
  return(to_return)
}

calculate_MD_tag <- function(seq, cigar, reference) {
  # Initialize variables
  md_tag <- ""
  seq_pos <- 1  # Position in the aligned sequence
  ref_pos <- 1  # Position in the reference sequence
  num_matches <- 0  # Number of matching bases

  # Split the CIGAR string into operations and lengths
  cigar_ops <- str_split(gsub("\\d", "", cigar), "") %>% unlist
  cigar_lengths <- as.integer(unlist(str_extract_all(cigar, "\\d+")))


  # Process each CIGAR operation
  for (i in seq_along(cigar_ops)) {
    op <- cigar_ops[i]
    length <- cigar_lengths[i]

    if (op %in% c("M", "X", "=", "D", "N", "P")) {
      # Match, mismatch, deletion, skip, padding
      aligned_seq <- substr(seq, seq_pos, seq_pos + length - 1)
      ref_seq <- substr(reference, ref_pos, ref_pos + length - 1)

      for (j in 1:length) {
        seq_base <- substr(aligned_seq, j, j)
        ref_base <- substr(ref_seq, j, j)

        if (seq_base != ref_base) {
          # Mismatch  # Reset match counter
          md_tag <- paste0(md_tag, num_matches, ref_base)
          num_matches <- 0
        } else {
          # Match
          num_matches <- num_matches + 1
        }

        seq_pos <- seq_pos + 1
        ref_pos <- ref_pos + 1
      }
    } else if (op %in% c("I", "S")) {
      # Insertion, soft clipping
      seq_pos <- seq_pos + length
    } else if (op %in% c("H")) {
      # Hard clipping
    }
  }

  return(paste("MD:Z:", md_tag))
}

md <- function(seq, cigar, reference){
  md_tag = ""
  seq_pos <- 1  # Position in the aligned sequence
  ref_pos <- 1  # Position in the reference sequence
  num_matches <- 0  # Number of matching bases

  # Split the CIGAR string into operations and lengths
  cigar_ops <- str_split(gsub("\\d", "", cigar), "") %>% unlist
  cigar_lengths <- as.integer(unlist(str_extract_all(cigar, "\\d+")))

  for(i in seq_along(cigar_ops)){
    op =  cigar_ops[i]
    length = cigar_lengths[i]

    if(op %in% c("I", "S")){
      seq_pos = seq_pos + length
    } else if (op %in% c("M", "X", "=", "D", "N", "P")) {
      aligned_seq = substr(seq, seq_pos, seq_pos + length - 1)
      ref_seq = substr(reference, ref_pos, ref_pos + length - 1)

      for(j in 1:length){
        seq_base <- substr(aligned_seq, j, j)
        ref_base <- substr(ref_seq, j, j)

        if(seq_base == ref_base && j == length){
          md_tag <- paste0(md_tag, num_matches, ref_base)
        } else if (seq_base == ref_base) {
          num_matches = num_matches + 1
        } else {

        }
      }

    }

  }

}
