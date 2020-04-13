
# Compare fraction spliced in each of the experiments
# and methods

setwd("../")

library(diem)

#=========================================
#=========================================

infn <- "data/raw/splice_frctn/all.splice_fraction.txt"
sf <- read.table(infn, header = TRUE, row.names = 1)

infn <- "data/raw/splice_frctn/midpoints.txt"
types <- read.table(infn, header = FALSE, row.names = 1)
types[,1] <- 100 * types[,1]

# Labels
at_labs <- c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
samples <- c("adpcyte", "mouse_nuclei_2k", at_labs)
expr <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))
names(samples) <- expr
labels <- c("DiffPA", "Mouse Brain", at_labs)
names(labels) <- samples
mthds <- c("DIEM", "EmptyDrops", "Quantile")

data_feat <- c("total_counts", 
               "n_genes", 
               "score.debris", 
               "pct.mt", 
               "MALAT1", 
               "SpliceFrctn")

sf_calls <- list()

for (i in 1:length(expr)){
    s <- samples[i]
    e <- names(samples)[i]

    # Read in results
    diem_fn <- paste0("data/processed/", e, "/diem/", s, ".diem_sce.rds")
    e_diem <- readRDS(diem_fn)
    ed_fn <- paste0("data/processed/", e, "/emptydrops/", s, ".emptydrops_counts.rds")
    e_ed <- readRDS(ed_fn)
    quant_fn <- paste0("data/processed/", e, "/quantile/", s, ".quantile_counts.rds")
    e_quant <- readRDS(quant_fn)

    test_data <- e_diem@test_data

    # Add droplet ID prefixes  
    rownames(e_diem@test_data) <- paste(samples[i], rownames(e_diem@test_data), sep="_")
    diem_si <- intersect(rownames(e_diem@test_data), rownames(sf))
    e_diem@test_data <- e_diem@test_data[diem_si,]

    colnames(e_ed$counts) <- paste(samples[i], colnames(e_ed$counts), sep="_")
    ed_si <- intersect(colnames(e_ed$counts), rownames(sf))
    e_ed$counts <- e_ed$counts[,ed_si]

    colnames(e_quant$counts) <- paste(samples[i], colnames(e_quant$counts), sep="_")
    q_si <- intersect(colnames(e_quant$counts), rownames(sf))
    e_quant$counts <- e_quant$counts[,q_si]

    test_drop <- rownames(e_diem@test_data)

    diem_call <- rep("Debris", length(test_drop)); names(diem_call) <- test_drop
    diem_call <- e_diem@test_data$Call
    ed_call <- rep("Debris", length(test_drop)); names(ed_call) <- test_drop
    ed_call[colnames(e_ed$counts)] <- "Clean"
    quant_call <- rep("Debris", length(test_drop)); names(quant_call) <- test_drop
    quant_call[colnames(e_quant$counts)] <- "Clean"

    calls <- list(diem_call, ed_call, quant_call)

    all_dfl <- lapply(1:3, function(j){
                      datf1 <- data.frame("Droplet" = test_drop, 
                                          "Method" = rep(mthds[j], length(test_drop)), 
                                          "Exprm" = rep(expr[i], length(test_drop)),
                                          "Sample" = rep(samples[i], length(test_drop)), 
                                          "Label" = rep(labels[i], length(test_drop)),
                                          "Call" = calls[[j]])
                      datf1 <- cbind(datf1, e_diem@test_data[,data_feat])
                      datf1[,"Truth"] <- rep("Nuclear", nrow(datf1))
                      datf1[datf1[,"SpliceFrctn"] >= types[samples[i],1], "Truth"] <- "Background"
                      return(datf1) })
    all_df <- do.call(rbind, all_dfl)
    rownames(all_df) <- NULL
                      
    sf_calls[[s]] <- all_df
}

all_sf <- do.call(rbind, sf_calls)
rownames(all_sf) <- NULL

write.table(all_sf, 
            "data/processed/meta_data.txt", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE, 
            sep = "\t")

