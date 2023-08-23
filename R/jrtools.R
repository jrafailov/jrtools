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
#' @export
#'
#' @examples
#' hunt()
hunt = function(huntdown = NULL){
  ps = ps() %>% data.table()
  ps_all <- ps[,.(total_memory_GB = round(sum(rss / (1024 * 1024), na.rm = TRUE) / 1024, 4)), by = username][order(-total_memory_GB)]
  if(!is.null(huntdown)){
    ps_user = ps[grep(huntdown, username),.(username, pid, ppid, name, status, resident_mem_gb = round(rss / (1024^3), 2), virtual_mem_gb = round(rss / (1024^3), 2), created)][order(-resident_mem_gb)]
  }
  if(is.null(huntdown)){
    return(ps_all)
  } else{
    return(ps_user)
  }
}
