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

md <- function(seq, cigar, reference) {
  md_tag = ""
  seq_pos <- 1 # Position in the aligned sequence
  ref_pos <- 1 # Position in the reference sequence
  num_matches <- 0 # Number of matching bases

  # Split the CIGAR string into operations and lengths
  cigar_ops <- str_split(gsub("\\d", "", cigar), "") %>% unlist()
  cigar_lengths <- as.integer(unlist(str_extract_all(cigar, "\\d+")))

  for (i in seq_along(cigar_ops)) {
    op = cigar_ops[i]
    length = cigar_lengths[i]

    if (op %in% c("I", "S")) {
      seq_pos = seq_pos + length
    } else if (op %in% c("M", "X", "=", "D", "N", "P")) {
      aligned_seq = substr(seq, seq_pos, seq_pos + length - 1)
      ref_seq = substr(reference, ref_pos, ref_pos + length - 1)

      for (j in 1:length) {
        seq_base <- substr(aligned_seq, j, j)
        ref_base <- substr(ref_seq, j, j)

        if (seq_base == ref_base && j == length) {
          md_tag <- paste0(md_tag, num_matches, ref_base)
        } else if (seq_base == ref_base) {
          num_matches = num_matches + 1
        } else {

        }
      }
    }
  }
}

parsesnpeff = function (
  vcf,
  snpeff_path = system.file("extdata", "snpeff_scripts", package = "multiplicity"),
  tumor_id = NULL,
  normal_id = NULL,
  filterpass = TRUE,
  coding_alt_only = TRUE, 
  geno = NULL,
  gr = NULL,
  keepfile = FALSE,
  altpipe = FALSE, 
  debug = FALSE,
  verbose = FALSE,
  bcftools = "bcftools", # FIXME: hardcoded for now.
  canonical_only = FALSE
  ) {
  if (debug)
    browser()
  tmp.path = tempfile(pattern = "tmp_", fileext = ".vcf.gz")
  if (!keepfile)
    on.exit(unlink(tmp.path))
  try2({
    # catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
    catcmd = paste(bcftools, "view -Ov")
    ## Redundant logic is for backwards compatibility with workflows
    ## in which a snpeff path is already provided.
    ## Fallback option encoded here in case it doesn't.
    onepline_path1 = paste0(snpeff_path, "/vcfEffOnePerLine.pl")
    onepline_path2 = paste0(snpeff_path, "/scripts/vcfEffOnePerLine.pl")
    onepline_path3 = system.file("extdata", "snpeff_scripts", "vcfEffOnePerLine.pl", package = "multiplicity")
    onepline = onepline_path3
    is_path1_same_as_path3 = identical(onepline_path1, onepline_path3)
    if (!is_path1_same_as_path3 && file.exists(onepline_path1)) {
      onepline = onepline_path1
    } else if (!is_path1_same_as_path3 && file.exists(onepline_path2)) {
      onepline = onepline_path2
    }
    # onepline = paste0(snpeff_path, "/vcfEffOnePerLine.pl")
    snpsift_path1 = paste0(snpeff_path, "/SnpSift.jar")
    snpsift_path2 = paste0(snpeff_path, "/scripts/SnpSift.jar")
    snpsift_path3 = system.file("extdata", "snpeff_scripts", "SnpSift.jar", package = "multiplicity")
    is_path1_same_as_path3 = identical(snpsift_path1, snpsift_path3)
    snpsift = snpsift_path3
    if (!is_path1_same_as_path3 && file.exists(snpsift_path1)) {
      snpsift = snpsift_path1
    } else if (!is_path1_same_as_path3 && file.exists(snpsift_path2)) {
      snpsift = snpsift_path2
    }
    if (verbose)(message(paste0("applying SnpSift to VCF: ", vcf)))
    if (coding_alt_only) {
      if(verbose)(message("Coding alterations only."))
      filt = paste0("java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar ",
        # snpeff_path, "/SnpSift.jar ",
        snpsift, " ",
        "filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\"")
      if (filterpass) {
        if(verbose) (message("Coding alterations only and FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, filt, bcftools, tmp.path)  
      } else {
        if(verbose) (message("Coding alterations only."))
          cmd = sprintf("%s %s | %s | %s | bgzip -c > %s", catcmd, vcf, onepline, filt, tmp.path)
      }
    } else {
      filt = ""
      if (filterpass){
        if(verbose)(message("FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, bcftools, tmp.path)
      } else {
        if(verbose)(message("All alterations included."))
        cmd = sprintf(paste(catcmd, "%s | %s | bgzip -c > %s"), vcf, onepline, tmp.path)
      }
    }
    if(verbose)(message("Performing command."))
    system(cmd)
    if(verbose)(message("SnpSift successfully applied!"))
  })
  if (!altpipe)
    out = grok_vcf(tmp.path, long = TRUE, geno = geno, gr = gr)
  else {
    if (verbose)(message(paste0("reading in SnpSift VCF.")))
    vcf = VariantAnnotation::readVcf(tmp.path)
    rr = MatrixGenerics::rowRanges(vcf)
    rr$REF = as.character(rr$REF)
    vcf_geno_lst = VariantAnnotation::geno(vcf)
    vcf_info = VariantAnnotation::info(vcf)
    #browser()
    ann = as.data.table(data.table::tstrsplit(unlist(vcf_info$ANN),
      "\\|"))[, 1:15, with = FALSE, drop = FALSE]
    fn = c("allele", "annotation", "impact", "gene", "gene_id",
      "feature_type", "feature_id", "transcript_type", 
      "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos", 
      "protein_pos", "distance")
    data.table::setnames(ann, fn)
    
    # Extract additional INFO fields if available
    additional_fields = c("CpG", "METH_PROB", "V_HAT", "MUT_ALT_EST", "METH_ALT_EST", "STATUS")
    available_fields = intersect(additional_fields, names(vcf_info))
    
    if (length(available_fields) > 0) {
      if (verbose)(message(paste0("Found additional INFO fields: ", paste(available_fields, collapse = ", "))))
      additional_info = vcf_info[available_fields]
    } else {
      additional_info = NULL
    }
    
    if ("AD" %in% names(vcf_geno_lst)) {
      if (verbose)(message(paste0("parsing AD field in VCF.")))
      vcf.ncol <- ncol(vcf_geno_lst$AD)
      if (verbose)(message(paste0(vcf.ncol, " columns found in the VCF.")))
      vcf.names <- colnames(vcf_geno_lst$AD)
      if(vcf.ncol > 1) {
        if (verbose)(message(paste0("parsing tumor and normal alt/ref counts.")))
        ## grab the last item in the array if ids not specified... presumably tumor
        tumor_col_id = vcf.ncol ## assume last column bydefault
        normal_col_id = NA_integer_
        is_id_matchable = function(id) {
          is_null_or_na = is.null(id) || is.na(id) || identical(id, "NA")
          return(!is_null_or_na)
        }
        if (is_id_matchable(tumor_id) && any(tumor_id %in% vcf.names)) {
          tumor_col_id = base::match(tumor_id, vcf.names)
        }
        if (is_id_matchable(normal_id) && any(normal_id %in% vcf.names)) {
          normal_col_id = base::match(normal_id, vcf.names)
        }
        adep = data.table::transpose(vcf_geno_lst$AD[,tumor_col_id])
        adep = as.data.table(adep)
        adep = base::subset(
          adep,
          select = 1:2
        )
        data.table::setnames(adep, c("ref", "alt"))
        adep$normal.ref = NA_integer_
        adep$normal.alt = NA_integer_
        if (!is.na(normal_col_id)) {
          adep.n = data.table::transpose(vcf_geno_lst$AD[,normal_col_id])
          adep$normal.ref = adep.n[[1]]
          adep$normal.alt = adep.n[[2]]
        }
        gt = vcf_geno_lst$GT
      } else {
        if (verbose)(message(paste0("parsing only tumor alt/ref counts.")))
        adep = vcf_geno_lst$AD[, vcf.ncol]
        adep = data.table::transpose(adep)
        adep = setnames(as.data.table(adep), c("ref", "alt"))
        gt = vcf_geno_lst$GT
      }
    } else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in%
                     c(names(vcf_geno_lst)))) {
      if(verbose)(message("parsing AU/GU/CU/TAR/TIR fields in VCF."))
      this.col = dim(vcf_geno_lst[["AU"]])[2]
      d.a = vcf_geno_lst[["AU"]][, , 1, drop = F][, this.col, 1]
      d.g = vcf_geno_lst[["GU"]][, , 1, drop = F][, this.col, 1]
      d.t = vcf_geno_lst[["TU"]][, , 1, drop = F][, this.col, 1]
      d.c = vcf_geno_lst[["CU"]][, , 1, drop = F][, this.col, 1]
      mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
      refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
      refid = ifelse(!isSNV(vcf), NA_integer_, refid)
      ## Note this assumes that multiallelic splitting via bcftools norm -m -any is run upstream (this breaks some fields if run after vcfOneLine.pl)
        ## otherwise you can't just unlist the ALT field if multiallelics are present
      altid = match(as.character(unlist(VariantAnnotation::fixed(vcf)$ALT)), colnames(mat))
      altid = ifelse(!isSNV(vcf), NA_integer_, altid)
      refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
      altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
      this.icol = dim(vcf_geno_lst[["TAR"]])[2]
      refindel = d.tar = vcf_geno_lst[["TAR"]][, , 1, drop = F][, this.icol, 1]
      altindel = d.tir = vcf_geno_lst[["TIR"]][, , 1, drop = F][, this.icol, 1]
      try2({
        n.d.a = vcf_geno_lst[["AU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.g = vcf_geno_lst[["GU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.t = vcf_geno_lst[["TU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.c = vcf_geno_lst[["CU"]][, , 1, drop = F][, this.col - 1, 1]
        n.mat = cbind(A = n.d.a, G = n.d.g, T = n.d.t, C = n.d.c)
        n.refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(n.mat))
        n.refid = ifelse(!isSNV(vcf), NA_integer_, refid)
        ## Note this assumes that multiallelic splitting via bcftools norm -m -any is run upstream (this breaks some fields if run after vcfOneLine.pl)
        ## otherwise you can't just unlist the ALT field if multiallelics are present
        n.altid = match(as.character(unlist(VariantAnnotation::fixed(vcf)$ALT)), colnames(n.mat))
        n.altid = ifelse(!isSNV(vcf), NA_integer_, altid)
        n.refsnv = n.mat[cbind(seq_len(nrow(n.mat)), refid)]
        n.altsnv = n.mat[cbind(seq_len(nrow(n.mat)), altid)]
        n.refindel = n.d.tar = vcf_geno_lst[["TAR"]][, , 1, drop = F][, this.icol - 1, 1]
        n.altindel = n.d.tir = vcf_geno_lst[["TIR"]][, , 1, drop = F][, this.icol - 1, 1]
      })
      adep = data.table(ref = coalesce(refsnv, refindel),
        alt = coalesce(altsnv, altindel))
      try2({
        adep.n = data.table(normal.ref = coalesce(n.refsnv, n.refindel),
          normal.alt = coalesce(n.altsnv, n.altindel))
        adep =  adep %>% cbind(adep.n)
      })
      gt = data.frame(GT = rep_len("", NROW(vcf)))
    } else {
      message("ref and alt count columns not recognized")
      adep = NULL
      gt = NULL
    }
    ## cbinding S4Vector DataFrame objects works much faster with dataframes
    ## Need to figure out problems with data.table conversions as some point
    mcols_to_cbind = list(mcols(rr), data.table::setDF(ann), data.table::setDF(adep))
    
    if (!is.null(additional_info)) {
      mcols_to_cbind = append(mcols_to_cbind, list(data.table::setDF(as.data.table(additional_info))))
    }
    
    mcols_to_cbind = append(mcols_to_cbind, list(gt = gt[, 1]))
    
    mcols(rr) = do.call(cbind, mcols_to_cbind)
    
    rr = S4Vectors::expand(rr, "ALT")
    rr$ALT = as.character(rr$ALT)
    
    # Add canonical transcript filtering
    if (canonical_only) {
      if (verbose)(message("Filtering to canonical transcripts only."))
      rr = dt2gr(gr2dt(rr)[, .SD[1], by = .(seqnames, start, end, REF, ALT, gene)])
    }
    
    out = rr
  }
  this.env = environment()
  return(this.env$out)
} 



#' @name rand.string
#' @title make a random string
#'
#' @return random string
#' @author Someone from Stackoverflow
rand.string <- function(n = 1, length = 12) {
  randomString <- c(1:n) # initialize vector
  for (i in 1:n) {
    randomString[i] <- paste(
      sample(c(0:9, letters, LETTERS),
        length,
        replace = TRUE
      ),
      collapse = ""
    )
  }
  return(randomString)
}

#' @name try2
#' @title wrapper around tryCatch - robust to parallel:: functions; A slightly more robust version of try that works within the parallel:: set of functions that pre-deploy a cluster.
try2 <- function(expr, ..., finally) {
  tryCatch(expr,
    error = function(e) {
      msg <- structure(paste(conditionMessage(e), conditionCall(e), sep = "\n"), class = "err")
      cat("Error: ", msg, "\n\n")
      return(msg)
    },
    finally = finally,
    ... = ...
  )
}

#' @name normalize_path
#' @title easy function that returns NULL if file path set to /dev/null
normalize_path <- function(path) {
  is_character = is.character(path)
  is_length_one = NROW(path) == 1
  is_a_potential_path = is_character && is_length_one
  is_invalid_input = is_character && !is_length_one
  out = path
  if (is_a_potential_path && !file.exists(path)) stop("Provided path does not exist")
  if (is_invalid_input) stop("Provided invalid paths (not length == 1)")
  if (!(is_a_potential_path)) return(out)
  if (path == "/dev/null") return(NULL)
  ## After all this.. this will return variable "path" (could be single length character or NULL)
  return(out)
}

