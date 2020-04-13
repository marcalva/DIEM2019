
# Process spliced fraction reads

setwd("../")



old_ids <- c("DiffPA18", "mouse_nuclei_2k", "SAT3", "SAT5", "SAT1", "SAT4", "SAT6", "SAT2")
lab_ids <- c("adpcyte", "mouse_nuclei_2k", "AT1", "AT2", "AT3", "AT4", "AT5", "AT6")

datfl <- lapply(1:length(old_ids), function(a){
                i <- old_ids[a]
                ifn <- paste0("data/raw/splice_frctn/", i, ".splice_fraction.txt")
                tb <- read.table(ifn, header = TRUE, row.names = 1, sep = ",", stringsAsFactors=FALSE)
                rownames(tb) <- gsub(old_ids[a], lab_ids[a], rownames(tb))
                return(tb)
})

datf <- do.call(rbind, datfl)

rownames(datf) <- sub("-", "_", rownames(datf))


write.table(datf, 
            paste0("data/raw/splice_frctn/all.splice_fraction.txt"), 
            row.names = TRUE, 
            col.names = NA, 
            sep = "\t", 
            quote = FALSE)

