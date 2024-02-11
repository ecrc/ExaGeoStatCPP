# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CreateSymbolTable.R
# @version 1.1.0
# @author David Helmy
# @date 2024-02-09

read_symbols_from_object_file <- function(f){
    nm <- Sys.getenv("UserNM")
    if(!nzchar(nm)) {
        ## reasonable to assume nm is on the path
        nm <- Sys.which("nm")
        if(nzchar(nm)) nm <- shQuote(nm)
    }
    if(!nzchar(nm)) {
        warning("this requires 'nm' to be on the PATH")
        return()
    }
    if(!(file.size(f))) return()
    s <- strsplit(system(sprintf("%s -Pg %s", nm, shQuote(f)),
                         intern = TRUE),
                  " +")
    ## Cannot simply rbind() this because elements may have 2-4 entries.
    n <- length(s)
    tab <- matrix("", nrow = n, ncol = 4L)
    colnames(tab) <- c("name", "type", "value", "size")
    ## Compute desired i and j positions in tab.
    i <- rep.int(seq_len(n), lengths(s))
    j <- unlist(lapply(s, seq_along))
    tab[n * (j - 1L) + i] <- unlist(s)
    tab
}
list_files_with_extension <- function(directory_path, file_extension) {
  # Get a list of all files in the directory and its subdirectories
  all_files <- list.files(directory_path, recursive = TRUE, full.names = TRUE)

  # Filter files based on the provided file extension
  selected_files <- all_files[grep(paste0("\\.", file_extension, "$"), all_files, ignore.case = TRUE)]

  return(selected_files)
}

# Define the arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("\n\n\n\n")
  stop("Please provide correct arguments")
}

directory_path <- (args[1])
file_path <- (args[2])

file_extension <- "o"


objects <- list_files_with_extension(directory_path, file_extension)
tables <- lapply(objects, read_symbols_from_object_file)
names(tables) <- objects
saveRDS(tables, file = file_path, version = 2)