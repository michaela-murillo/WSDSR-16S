#  ___        __  __ _____ ____    _____           _           _   
# |__ \      |  \/  |_   _|  _ \  |  __ \         (_)         | |  
#    ) ______| \  / | | | | |_) | | |__) _ __ ___  _  ___  ___| |_ 
#   / |______| |\/| | | | |  _ <  |  ___| '__/ _ \| |/ _ \/ __| __|
#  / /_      | |  | |_| |_| |_) | | |   | | | (_) | |  __| (__| |_ 
# |____|     |_|  |_|_____|____/  |_|   |_|  \___/| |\___|\___|\__|
#                                                _/ |              
#                                                |__/               


# Comprehensive code for automated reservoir water analysis
# Designed to analyse genus-level kraken2/bracken data by batch

# Output list
#   [] QC report of kraken2/bracken data
#   [] QC report graphs
#   [x] Graph of top 10 genera found in each sample
#   [x] Excel file of cyanobacteria
#   [x] Relative abundance graphs of 2-MIB producing cyanobacteria
#   [] Excel file of actinomycetes
#   [] Relative abundance graphs of actinomycetes
#   [] Krona plots of each sample

# Improvements awaiting implementation
#   [] Designated color for every 2-MIB producer cyanobacteria for consistency
#   [] Make Top 10 plots comparative tables instead of bar charts - cleaner
#   [] Create NA result checks for taxizedb query
#   [] For bacteria excel files, save each location in separate sheets



# ==============================================================================
# Update local taxizedb db (must do every new day of analysis!)
# ==============================================================================

# library(taxizedb)
# print("Updating local taxizedb database...")
# db_download_ncbi(overwrite = TRUE)
# print("Database update successful.")



# ==============================================================================
# Specify input and analysis parameters
# ==============================================================================

# Path to dir with raw kraken2-bracken .txt files + where output is placed       note: must replace \ with /
filepath <- "C:/Users/Mikki/Documents/WSD/ONS_week1-4/serviceRun_batch1-4/raw_batch4"

# Path to dir for storing cyanobacteria list                                     note: must replace \ with \\
cyanolist_path <- "C:\\Users\\Mikki\\Documents\\WSD\\ONS_week1-4\\serviceRun_batch1-4\\raw_batch4"

# Batch name
batchName <- "Batch004"

# Sample collection date
samplingdate <- "20231123"


# --------------------------------------------------------
# Choose reservoir (PC = Plover Cove, TLC = Tai Lam Chung)
# --------------------------------------------------------

# reservoir <- "PC"
reservoir <- "TLC"


# -----------------------
# Barcode to sample names                                                        note: must place in ascending order of barcode no.
# -----------------------

# barcodeNames <- c(                                                              # for PC
#   "Intake Tower Surface",
#   "Intake Tower -3 MPD",
#   "Intake Tower -8 MPD",
#   "Portal AJ")

barcodeNames <- c(                                                               # for TLC
  "Portal L",
  "Draw Off Tower Surface",
  "Draw Off Tower 12M",
  "Draw Off Tower Bottom")


# ------------------------------------------------------------------------------

# Set relative abundance cutoff
cutoffValue_relab = 0.001

# Known cyanobacterial 2-MIB producer list
mibProducers <- c(
  "Anabaena",
  "Cyanobium",
  "Cylindrospermopsis",
  "Dolichospermum",
  "Leptolyngbya",
  "Microcystis",
  "Oscillatoria",
  "Phormidium",
  "Planktothrix",
  "Pseudanabaena",
  "Synechococcus")



# ==============================================================================
# Prerequisites
# ==============================================================================

# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(taxizedb)
library(writexl)

# Set working directory
setwd(filepath)

# Read and rename raw kraken2/bracken output
inputlist = list.files(pattern = "*.txt$")
for (i in 1:length(inputlist))
  assign(barcodeNames[i], 
         read.delim(inputlist[i], header = TRUE))

# 'Barcode' each df by location (i.e. df name)
object_names <- ls()
for (name in object_names) {
  if (is.data.frame(get(name))) {
    df <- get(name)
    df$Location <- name
    assign(name, df)
  }}

rm(df)


# ==============================================================================
# kraken2/bracken file cleanup
# ==============================================================================

# Apply relative abundance cutoff
df_names <- ls()
for (df_name in df_names) {
  if (is.data.frame(get(df_name))) {
    original_df <- get(df_name)
    subset_df <- original_df[original_df$fraction_total_reads >= cutoffValue_relab, ]
    new_df_name <- paste0("relabCut", " ", df_name)
    assign(new_df_name, subset_df, envir = .GlobalEnv)
  }}

# Clean out redundant variables
rm(original_df)
rm(subset_df)



# ==============================================================================
# Graph of top 10 genera found in each sample
# ==============================================================================

# Create new dfs containing only top 10 genera in each location
df_list <- ls(pattern = "^relabCut")
for (df_name in df_list) {
  sorted_df <- get(df_name)[order(-get(df_name)$fraction_total_reads), ]
  top10_df <- sorted_df[1:10, ]
  new_df_name <- paste0("top10 ", df_name)
  assign(new_df_name, top10_df)
}

# Pool into one df
df_list <- ls(pattern = "^top10 relabCut")
top10_comp <- data.frame()
for (df_name in df_list) {
  df <- get(df_name)
  top10_comp <- rbind(top10_comp, df)
}

# Clean out redundant variables
rm(sorted_df)
rm(top10_df)
rm(df)


# ----------
# Make graph
# ----------

# Changeable elements
plotTitle <- paste0(batchName, ": ", reservoir, " ", samplingdate,
                    "\n Top 10 Bacterial Genera By Location")
plotCaption <- paste("Relative abundance cutoff = ", cutoffValue_relab)


# Define bar labels
top10_comp$barlabel <- paste(top10_comp$name, 
                             top10_comp$fraction_total_reads, 
                             sep = ": ")

# Plot
plot_top10 <-
  ggplot(data = top10_comp, 
         aes(x = Location, 
             y = fraction_total_reads,
             fill = name,
             label = fraction_total_reads)) +
  geom_bar(stat = "identity",
           color = "grey") +
  geom_label_repel(mapping = aes(label = barlabel),
                   size = 3,
                   position = position_stack(vjust = 0.5),
                   color = "white",
                   segment.alpha = 0.5,
                   segment.color = "black",
                   max.overlaps = Inf) +
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  ggtitle(plotTitle) +
  scale_x_discrete("Source",
                   limits = barcodeNames) +
  scale_y_continuous("Relative Abundance") +
  labs(fill = "Genus", caption = plotCaption) +
  theme(plot.title = element_text(hjust = 0.5))

plot_top10


# Save as image
ggsave(paste0(as.character(samplingdate), "_",
              batchName, "_",
              reservoir, 
              "_Top10.jpg"),
       plot_top10,
       path = filepath,
       width = 6.188*2,
       height = 4.375*2)



# ==============================================================================
# taxizedb query: phylum and class queries of every genus
# ==============================================================================

# Filter out the correct dfs for query
objects <- ls()
filtered_objects <- objects[grep("^relabCut", objects)]

# Convert into list of dfs
dfs_list <- lapply(filtered_objects, get)

# rbind into a single df
monsterquery <- do.call(rbind, dfs_list)

# Put taxids into vector
query_list <- c(monsterquery$taxonomy_id)


# Phylum query -----------------------------------------------------------------

# run phylum query on locally downloaded NCBI db
print("Running taxizedb query at for phyla...")
query_output <- taxa_at(query_list,
                        rank = "phylum",
                        db = "ncbi",
                        missing = "lower",
                        verbose = TRUE,
                        warn = TRUE)
print("Query complete.")

# convert output from list into df
query_output_df <- do.call(rbind, query_output)

# add column of G-specific taxid to query output df
query_output_df$taxonomy_id <- rownames(query_output_df)

# add phylum information to original df
phylum_output <- merge(x = monsterquery, 
                       y = query_output_df, 
                       by = "taxonomy_id", 
                       all.x = TRUE)



# ==============================================================================
# Filter for cyanobacteria
# ==============================================================================

# subset for cyanobacteria only
output_cyano <- subset(phylum_output, name.y == "Cyanobacteriota")

# rename/delete some variables for clarity
output_cyano <- rename(output_cyano, !!"genus" := name.x)
output_cyano <- rename(output_cyano, !!"phylum" := name.y)
output_cyano <- select(output_cyano, -taxonomy_lvl, -rank, -id)
output_cyano$source <- gsub("relabCut ", "", output_cyano$source)

# export as excel
cyanolist_path_full <- paste0(cyanolist_path, 
                              "\\", as.character(samplingdate), 
                              "_", batchName,
                              "_", reservoir, 
                              "_cyanoList.xlsx")

write_xlsx(output_cyano, 
           path = cyanolist_path_full)



# ==============================================================================
# Filter for known 2-MIB producing cyanobacteria
# ==============================================================================

mibproducing_cyano <- output_cyano[output_cyano$genus %in% mibProducers, ]       # note: in cyanoList output, should rearrange so there's a column denoting which are known 2-MIB producers



# ==============================================================================
# Plot generation
# ==============================================================================

# Changeable elements
plotTitle <- paste0(batchName, ": ", reservoir, " ", samplingdate)
plotCaption <- paste("Relative abundance cutoff = ", cutoffValue_relab)

# Plot parameters
finalplot <-
ggplot(data = mibproducing_cyano, 
       aes(x = source, 
           y = fraction_total_reads,
           fill = genus,
           label = fraction_total_reads)) +
  geom_bar(stat = "identity",
           color = "grey") +
  geom_label_repel(mapping = aes(label = fraction_total_reads),
                   size = 3,
                   position = position_stack(vjust = 0.5),
                   color = "white",
                   segment.alpha = 0.5,
                   segment.color = "black",
                   max.overlaps = Inf) +
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  ggtitle(plotTitle) +
  scale_x_discrete("Source",
                   limits = barcodeNames) +
  scale_y_continuous("Relative Abundance") +
  labs(fill = "Cyanobacteria \n(Known 2-MIB Producers)",
       caption = plotCaption) +
  theme(plot.title = element_text(hjust = 0.5))

finalplot

# Save as image
ggsave(paste0(as.character(samplingdate), "_", reservoir, "_mibProducers.jpg"),
       finalplot,
       path = filepath,
       width = 6.188*1.75,
       height = 4.375*1.75)
