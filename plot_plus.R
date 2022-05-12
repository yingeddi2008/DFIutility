pdf_date <- function(filename, height, width, onefile){
  pdf(paste0(filename,"_", format(Sys.Date(), format = "%Y%m%d"),".pdf"),
      height = height,
      width = width,
      onefile = onefile)
}
