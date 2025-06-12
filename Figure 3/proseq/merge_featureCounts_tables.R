library(data.table)

# List all featureCounts output files
files <- list.files(pattern = "*.featureCounts.txt")

# Function to process each file
process_file <- function(file) {
  # Read the file to determine the number of comment lines
  lines <- readLines(file)
  comment_lines <- sum(grepl("^#", lines))
  
  # Read the file, skipping the comment lines to get the header and data
  dt <- fread(file, skip = comment_lines, header = TRUE)
  
  # Ensure "Geneid" column is present
  if (!("Geneid" %in% names(dt))) {
    stop("Geneid column not found in file: ", file)
  }
  
  # Select "Geneid" and the specific count columns
  dt_selected <- dt[, c("Geneid", grep("sorted.bam$", names(dt), value = TRUE)), with = FALSE]
  
  # Rename the count columns to retain only the relevant part of their names
  new_count_col_names <- gsub("\\.sorted\\.bam$", "", names(dt_selected)[-1])  # Remove '.sorted.bam' suffix
  setnames(dt_selected, old = names(dt_selected)[-1], new = new_count_col_names)
  
  return(dt_selected)
}

# Apply the function to all files and merge
merged_data <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), lapply(files, process_file))

# Write the merged matrix to a file
fwrite(merged_data, "merged_featureCounts_matrix.txt", sep = "\t")

