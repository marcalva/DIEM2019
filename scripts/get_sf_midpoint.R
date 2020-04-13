
# Get midpoint from GMM of SF

setwd("../")

library(mixtools)

old_ids <- c("DiffPA18", "mouse_nuclei_2k", "SAT3", "SAT5", "SAT1", "SAT4", "SAT6", "SAT2")
lab_ids <- c("adpcyte", "mouse_nuclei_2k", "AT1", "AT2", "AT3", "AT4", "AT5", "AT6")
dir_ids <- c("adpcyte", "mouse_nuclei_2k", rep("atsn", 6))

sf <- read.table("data/raw/splice_frctn/all.splice_fraction.txt", 
                 header = TRUE, 
                 row.names = 1, 
                 sep = "\t", 
                 stringsAsFactors = FALSE)

set.seed(1)

midpoints <- list()

for (i in 1:8){
    s <- lab_ids[i]
    d <- dir_ids[i]
    scefn <- paste0("data/processed/", d, "/diem/", s, ".diem_sce.rds")
    sce <- readRDS(scefn)
    datf <- sf

    test_data <- sce@test_data
    #keep <- rownames(test_data)[test_data$n_genes >= 100]
    #test_data <- test_data[keep,]
    rownames(test_data) <- paste0(s, "_", rownames(test_data))
    
    si <- intersect(rownames(datf), rownames(test_data))
    test_data <- test_data[si,]
    datf <- datf[si,]

    ret <- normalmixEM(datf, k=2, mu = c(.3, .7))
    #ret <- normalmixEM(datf, k=2)

    if (ret$mu[1] < ret$mu[2]){
        mu1 <- ret$mu[1]
        mu2 <- ret$mu[2]
        v1 <- ret$sigma[1]^2
        v2 <- ret$sigma[2]^2
    } else {
        mu1 <- ret$mu[2]
        mu2 <- ret$mu[1]
        v1 <- ret$sigma[2]^2
        v2 <- ret$sigma[1]^2
    }


    a <- (v1 - v2)/(v1*v2)
    b <- 2*(v2*mu1 - v1*mu2)/(v1*v2)
    c <- (v1*mu2^2 - v2*mu1^2)/(v1*v2) - log(v1/v2)
    s1 <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
    s2 <- (-b - sqrt(b^2 - 4*a*c))/(2*a)

    print(s)
    if (s1 < mu2 & s1 > mu1){
        print(s1)
        midpoints[[s]] <- s1
    } else{
        print(s2)
        midpoints[[s]] <- s2
    }
}

midpoints <- do.call(rbind, midpoints)

write.table(midpoints, "data/raw/splice_frctn/midpoints.txt", row.names = T, col.names = F, sep = "\t")


