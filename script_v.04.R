# Comprehensive code for automated metagenomic analysis
# Designed to analyse genus-level kraken2/bracken data by batch

# Output list
#   [] Analysis statistics
#   [] Excel of all bacteria in samples post relAb-cutoff
#   [x] Graph of top 10 genera found in each sample
#   [x] Excel file of cyanobacteria
#   [x] Relative abundance graphs of cyanobacteria
#   [x] Relative abundance graphs of 2-MIB producing cyanobacteria
#   [x] Excel file of actinomycetes
#   [x] Relative abundance graphs of actinomycetes

# Improvements awaiting implementation (sorted by priority)
#   [] Stop analysis if new cyanobacterial genus that hasn't been checked for 2-MIB production is detected
#   [] Create NA result checks for taxizedb query
#   [~] Designated color for every 2-MIB producer cyanobacteria for consistency
#   [] For bacteria excel files, save each location in separate sheets
#   [] Add column to 'cyanobacteria' excel (Boolean variable called 'mib_producer')
#   [] Make Top 10 plots comparative tables in addition to bar charts



# Update local taxizedb db (must do every new day of analysis!) ================

# library(taxizedb)
# print("Updating local taxizedb database...")
# db_download_ncbi(overwrite = TRUE)
# print("Database update successful.")



# Specify input and analysis parameters ========================================


# Path to dir for input and output files
# Note: must replace \ with /
filepath <- "C:/Users/Michaela/Documents/Rprojects/WSD_16S_analysis/raw_batch58"

# Path to dir for storing excel files
# Note: must replace \ with \\
excelOutput_path <- "C:\\Users\\Michaela\\Documents\\Rprojects\\WSD_16S_analysis\\raw_batch58"

# Batch name
batchName <- "Batch058"

# Sample collection date
samplingdate <- "20240410"



# Choose reservoir (PC = Plover Cove, TLC = Tai Lam Chung) ---------------------

# reservoir <- "PC"
reservoir <- "TLC"
# reservoir <- "SP"



# Barcode to sample names ------------------------------------------------------
# note: must place in ascending order of barcode no!

# barcodeNames <- c(                                                             # for PC
#   "Intake Tower Surface",
#   "Intake Tower -3 MPD",
#   "Intake Tower -8 MPD",
#   "Portal AJ"
#   )

barcodeNames <- c(                                                             # for TLC
  "Portal L",
  "Draw Off Tower Surface",
  "Draw Off Tower 12M",
  "Draw Off Tower Bottom")

# barcodeNames <- c(                                                               # for SP
#   "Draw Off Tower Surface",
#   "Draw Off Tower Middle",
#   "Draw Off Tower Bottom")

# barcodeNames <- c(
#   "Bottom",
#   "Portal L",
#   "Surface",
#   "12M")



# Set relative abundance cutoff ------------------------------------------------
cutoffValue_relab = 0.001

# Known cyanobacterial 2-MIB producer list
mib_db <- read.delim("C:/Users/Michaela/Documents/Rprojects/WSD_16S_analysis/mib_producer_db.txt")

# temporary workaround: custom color scheme for mib-producing cyanobacteria
mib_db_customPalette <- c("Anabaena" = "#f8766d",
                          "Cyanobium" = "#e88526",
                          "Cylindrospermopsis" = "#d39200",
                          "Dolichospermum" = "#b79f00",
                          "Hapalosiphon" = "#93aa00",
                          "Hyella" = "#5eb300",
                          "Leptolyngbya" = "#00ba38",
                          "Lyngbya" = "#00bf74",
                          "Microcoleus" = "#00c19f",
                          "Microcoleus pseudatumnalis" = "#00bfc4",
                          "Microcystis" = "#00b9e3",
                          "Odorella" = "#00adfa",
                          "Oscillatoria" = "#619cff",
                          "Phormidium" = "#ae87ff",
                          "Planktothricoides" = "#db72fb",
                          "Planktothrix" = "#f564e3",
                          "Pseudanabaena" = "#ff61c3",
                          "Synechococcus" = "#ff699c")



# Prerequisites ================================================================

# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(rstudioapi)
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
    assign(name, df)}}

# clean redundant variables
rm(df)



# kraken2/bracken file cleanup =================================================

# Apply relative abundance cutoff
df_names <- barcodeNames
for (df_name in df_names) {
  if (is.data.frame(get(df_name))) {
    original_df <- get(df_name)
    subset_df <- original_df[original_df$fraction_total_reads >= cutoffValue_relab, ]
    new_df_name <- paste0("relabCut", " ", df_name)
    assign(new_df_name, subset_df, envir = .GlobalEnv)}}


# Change relAb values from decimal to percentage
df_list <- ls(pattern = "^relabCut")
for (df_name in df_list) {
  df <- get(df_name)
  df$fraction_total_reads <- df$fraction_total_reads * 100
  new_df_name <- paste0("trans ", df_name)
  assign(new_df_name, df, envir = .GlobalEnv)}

cutoffValue_relab <- cutoffValue_relab * 100


# Clean out redundant variables
rm(original_df)
rm(subset_df)



# Graph of top 10 genera found in each sample ==================================

# Create new dfs containing only top 10 genera in each location
df_list <- ls(pattern = "^trans relabCut")
for (df_name in df_list) {
  sorted_df <- get(df_name)[order(-get(df_name)$fraction_total_reads), ]
  top10_df <- sorted_df[1:10, ]
  new_df_name <- paste0("top10 ", df_name)
  assign(new_df_name, top10_df)}


# remove placeholder variable to avoid bugged graph
rm(top10_df)


# Pool into one df
df_list <- ls(pattern = "^top10")
top10_comp <- data.frame()
for (df_name in df_list) {
  df <- get(df_name)
  top10_comp <- rbind(top10_comp, df)}


# Clean out redundant variables
rm(sorted_df)
rm(df)



# Make graph -------------------------------------------------------------------

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
  scale_x_discrete("Location",
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



# taxizedb query: phylum and class queries =====================================

# Filter out the correct dfs for query
objects <- ls()
filtered_objects <- objects[grep("^trans relabCut", objects)]

# Pool into list of dfs
dfs_list <- lapply(filtered_objects, get)

# rbind into a single df
monsterquery <- do.call(rbind, dfs_list)

# Put variable 'taxids' into vector
query_list <- c(monsterquery$taxonomy_id)


# PHASE 1: Phylum-level query --------------------------------------------------

# run phylum query on locally downloaded NCBI db
print("Running taxizedb query for phyla...")
query_output <- taxa_at(query_list,
                        rank = "phylum",
                        db = "ncbi",
                        missing = "lower",
                        verbose = TRUE,
                        warn = TRUE)
print("Phylum query complete.")

# convert output from list into df
query_output_df <- do.call(rbind, query_output)

# add column of G-specific taxid to query output df
query_output_df$taxonomy_id <- rownames(query_output_df)

# add phylum information to original df
phylum_output <- merge(x = monsterquery, 
                       y = query_output_df, 
                       by = "taxonomy_id", 
                       all.x = TRUE)


# PHASE 2: Class-level query ---------------------------------------------------

# Clean up phylum output
phylum_output <- select(phylum_output, -taxonomy_lvl, -rank, -id)

# run class query on locally downloaded NCBI db
print("Running taxizedb query for class...")
query_output <- taxa_at(query_list,
                        rank = "class",
                        db = "ncbi",
                        missing = "lower",
                        verbose = TRUE,
                        warn = TRUE)
print("Class query complete.")

# convert output from list into df
query_output_df <- do.call(rbind, query_output)

# add column of G-specific taxid to query output df
query_output_df$taxonomy_id <- rownames(query_output_df)

# add phylum information to original df
class_output <- merge(x = phylum_output, 
                       y = query_output_df, 
                       by = "taxonomy_id", 
                       all.x = TRUE)



# Clean up final output --------------------------------------------------------
taxizedb_output <- rename(class_output, !!"genus" := name.x)
taxizedb_output <- rename(taxizedb_output, !!"phylum" := name.y)
taxizedb_output <- rename(taxizedb_output, !!"class" := name)
taxizedb_output <- select(taxizedb_output, -rank, -id)


# remove redundant variables
rm(class_output)
rm(dfs_list)
rm(monsterquery)
rm(phylum_output)
rm(query_output)
rm(query_output_df)



# Subset for cyanoabacteriota (phylum) =========================================

# subset for cyanobacteria only
subset_cyano <- subset(taxizedb_output, phylum == "Cyanobacteriota")


# create path for saving excel
cyanoList_path_full <- paste0(excelOutput_path, 
                              "\\", as.character(samplingdate), 
                              "_", batchName,
                              "_", reservoir, 
                              "_cyanobacteria.xlsx")

# save excel of cyanobacteria
write_xlsx(subset_cyano, 
           path = cyanoList_path_full)



# Subset for actinomycetes (class) =============================================

# subset for actinomycetes only
subset_actinomycetes <- subset(taxizedb_output, class == "Actinomycetes")


# create path for saving excel
actinomyceteList_path_full <- paste0(excelOutput_path, 
                              "\\", as.character(samplingdate), 
                              "_", batchName,
                              "_", reservoir, 
                              "_actinomycetes.xlsx")

# save excel of actinomycetes
write_xlsx(subset_actinomycetes, 
           path = actinomyceteList_path_full)



# Plot: relative abundance of cyanobacteria ====================================

# Changeable elements
plotTitle <- paste0(batchName, ": ", reservoir, " ", samplingdate,
                    "\n Relative Abundance of Cyanobacteria")
plotCaption <- paste("Relative abundance cutoff = ", cutoffValue_relab)

# Plot parameters
plot_cyano <-
  ggplot(data = subset_cyano, 
         aes(x = Location, 
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
  scale_x_discrete("Location",
                   limits = barcodeNames) +
  scale_y_continuous("Relative Abundance (%)") +
  labs(fill = "Genus",
       caption = plotCaption) +
  theme(plot.title = element_text(hjust = 0.5))

plot_cyano

# Save as image
ggsave(paste0(as.character(samplingdate), "_", 
              batchName, "_",
              reservoir, "_cyanobacteria.jpg"),
       plot_cyano,
       path = filepath,
       width = 6.188*2,
       height = 4.375*2)



# Plot: relative abundance of 2-MIB producing cyanobacteria ====================

# Subset known mib-producing cyanobacteria into new df
mib_db_true <- subset(mib_db, mib_production == TRUE)
mibproducing_cyano <- subset_cyano[subset_cyano$genus %in% mib_db_true$genus, ]
# mibproducing_cyano <- subset_cyano[subset_cyano$genus %in% mibProducers, ]


# Changeable elements
plotTitle <- paste0(batchName, ": ", reservoir, " ", samplingdate,
                    "\n Relative Abundance of Known 2-MIB Producing Cyanobacteria")
plotCaption <- paste("Relative abundance cutoff = ", cutoffValue_relab)


# (Optional) custom order for x-axis
customOrder <- c("Portal L", "Surface", "12M", "Bottom")
customRename <- c("Portal L,", "Draw Off Tower Surface",
                  "Draw Off Tower 12M", "Draw Off Tower Bottom")


# Plot parameters
plot_cyano_mibProducers <-
  ggplot(data = mibproducing_cyano, 
       aes(x = Location, 
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
  scale_x_discrete("Location",
                   limits = barcodeNames) +
  scale_y_continuous("Relative Abundance (%)") +
  labs(fill = "Genus",
       caption = plotCaption) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = mib_db_customPalette)

plot_cyano_mibProducers

# Save as image
ggsave(paste0(as.character(samplingdate), "_", 
              batchName, "_",
              reservoir, "_cyano_mibProducers.jpg"),
       plot_cyano_mibProducers,
       path = filepath,
       width = 6.188*1.75,
       height = 4.375*1.75)



# Plot: relative abundance of actinomycetes ====================================

# Changeable elements
plotTitle <- paste0(batchName, ": ", reservoir, " ", samplingdate,
                    "\n Relative Abundance of Actinomycetes")
plotCaption <- paste("Relative abundance cutoff = ", cutoffValue_relab)

# Plot parameters
plot_actinomycetes <-
  ggplot(data = subset_actinomycetes, 
         aes(x = Location, 
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
  scale_x_discrete("Location",
                   limits = barcodeNames) +
  scale_y_continuous("Relative Abundance (%)") +
  labs(fill = "Genus",
       caption = plotCaption) +
  theme(plot.title = element_text(hjust = 0.5))

plot_actinomycetes

# Save as image
ggsave(paste0(as.character(samplingdate), "_", 
              batchName, "_",
              reservoir, "_actinomycetes.jpg"),
       plot_actinomycetes,
       path = filepath,
       width = 6.188*1.75,
       height = 4.375*1.75)



# Analysis statistics and report ===============================================


# total reads per raw kraken file
# % of total reads that are cyanobacteria
# % of total reads taht are cyano mib producers
# % of total reads that are actinomycetes



# Excel for internal report ====================================================

# subset streptomyces from actinomycetes df
streptomyces <- subset(subset_actinomycetes, taxonomy_id == 1883)
  
# merge strep + mib-producing cyano data into 1 df
internal_report_merged <- rbind(mibproducing_cyano, streptomyces)

# create path for saving excel
internal_report_path <- paste0(excelOutput_path, 
                             "\\", as.character(samplingdate), 
                             "_", batchName,
                             "_", reservoir, 
                             "_internalReport.xlsx")

# save excel for internal report
write_xlsx(internal_report_merged, 
           path = internal_report_path)




# Closing actions ==============================================================

# Save current R script to working dir
system_date <- format(Sys.Date(), "%Y%m%d")

file_name <- paste0("Rscript_16S_",
                    system_date, "_",
                    batchName,".R")

destination_file_path <- file.path(filepath, file_name)

unsaved_script <- getSourceEditorContext()$contents

writeLines(unsaved_script, destination_file_path)


# Successful analysis message
print("Analysis complete ~")

