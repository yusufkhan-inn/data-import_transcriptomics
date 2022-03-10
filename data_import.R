## RNA-seq pipeline #############
# Import and quality-check data #


# Source functions and libraries 
source("./R/lib_fun.R")



# Notes:

# Read counts and transcripts per million are stored in txt/tsv files 
# the following script imports them and store in suitable formats (as RDA-files (ignored by git))

# Methods/counts

# RSEM: readcount (RC) & transcript per million (TPM)
# Kallisto: RC & TPM
# Salmon: RC & TPM
# HISAT2: RC
# STAR: RC



###### RSEM  ############

rsem_rc <- read_delim("./seqDATA/rsem_geneid_expectedcount_matrix.txt", delim = "\t")
rsem_tpm <- read_delim("./seqDATA/rsem_geneid_tpm_matrix.txt", delim = "\t")



# Load leg-condition data
legs_condition <- read.csv("./data/oneThreeSetLeg.csv", sep = ";") %>%
  pivot_longer(names_to = "sets", values_to = "leg", multiple:single) %>%
  print()

# Extract sample information from colnames in the count matrix
design.df <- data.frame(sample = colnames(rsem_tpm)[-1]) %>%
  rowwise() %>%
  mutate(sample = as.character(sample), 
         subject = str_split_fixed(sample, "w", 2)[1], 
         leg = substrRight(sample, 1), 
         time = gsub(subject, "", sample), 
         time = gsub(leg, "", time)) %>%
  inner_join(legs_condition) %>%
  dplyr::select(sample, subject, time, sets) %>%
  print()

#### Kallisto ###########
  
  
kallisto_rc <- read_delim("./seqDATA/kallisto_geneid_readcount_all_samples.tsv", delim = "\t")
kallisto_tpm <- read_delim("./seqDATA/kallisto_geneid_tpm_all_samples.tsv", delim = "\t")
  
#### Salmon ###########


salmon_rc <- read_delim("./seqDATA/salmon_geneid_readcount_all_samples.tsv", delim = "\t")
salmon_tpm <- read_delim("./seqDATA/salmon_geneid_tpm_all_samples.tsv", delim = "\t")

#### Hisat2 ########

hisat2_rc <- read_delim("./seqDATA/hisat2-htseq_geneid_expectedcount_matrix1.txt", delim = "\t")



#### Star ############


star_rc <- read_delim("./seqDATA/star-htseq_geneid_expectedcount_matrix1.txt", delim = "\t")



########### Create DEGlist of all count methods #################

###### Filter and normalization ###########

# Removes genes with zero count in all samples
nonzero_rsem     <- rowSums(rsem_rc[,-1]) > 0
nonzero_kallisto <- rowSums(kallisto_rc[,-1]) > 0
nonzero_salmon   <- rowSums(salmon_rc[,-1]) > 0
nonzero_star     <- rowSums(star_rc[,-1]) > 0
nonzero_hisat2   <- rowSums(hisat2_rc[,-1]) > 0

# Percentage of genes with non-zero count in each method


nonzero <- data.frame(method = c("rsem", "kallisto", "salmon", "star", "hisat2"),
             nonzero = c(sum(nonzero_rsem, na.rm = TRUE) / length(nonzero_rsem), 
                       sum(nonzero_kallisto, na.rm = TRUE) / length(nonzero_kallisto),
                       sum(nonzero_salmon, na.rm = TRUE) / length(nonzero_salmon),
                       sum(nonzero_star, na.rm = TRUE) / length(nonzero_star),
                       sum(nonzero_hisat2, na.rm = TRUE) / length(nonzero_hisat2))) %>%
  print()
  

### Save data for preliminary report ...
saveRDS(nonzero, file = "./data/derivedData/prel_analysis/nonzero.RDS")




# Remove zero counts
ct_rsem     <- rsem_rc[nonzero_rsem,]
ct_kallisto <- kallisto_rc[nonzero_kallisto,]
ct_salmon   <- salmon_rc[nonzero_salmon,]
ct_star     <- star_rc[nonzero_star,]
ct_hisat2   <- hisat2_rc[nonzero_hisat2,]

# Create matrix 
y_rsem     <- as.matrix(data.frame(ct_rsem)[,-1])
y_kallisto <- as.matrix(data.frame(ct_kallisto)[,-1])
y_salmon   <- as.matrix(data.frame(ct_salmon)[,-1])
y_star     <- as.matrix(data.frame(ct_star)[,-1])
y_hisat2   <- as.matrix(data.frame(ct_hisat2)[,-1])




# Set rownames (gene ID)  #


rownames(y_rsem)      <- data.frame(ct_rsem)[,1]
rownames(y_kallisto)  <- data.frame(ct_kallisto)[,1]
rownames(y_salmon)    <- data.frame(ct_salmon)[,1]
rownames(y_star)      <- data.frame(ct_star)[,1]
rownames(y_hisat2)    <- data.frame(ct_hisat2)[,1]


## Store DGELists in a list of lists (!)

dge_lists <- list(rsem = DGEList(y_rsem),    
                  kallisto = DGEList(y_kallisto), 
                  salmon = DGEList(y_salmon),
                  star = DGEList(y_star),    
                  hisat = DGEList(y_hisat2))



### Set sample information ###

for(i in 1:length(dge_lists)) {
  
  dge_lists[[i]]$samples <- dge_lists[[i]]$samples %>%
    mutate(sample = rownames(.)) %>%
    rowwise() %>%
    mutate(sample = as.character(sample), 
           subject = str_split_fixed(sample, "w", 2)[1], 
           leg = substrRight(sample, 1), 
           time = gsub(subject, "", sample), 
           time = gsub(leg, "", time)) %>%
    inner_join(legs_condition) %>%
    dplyr::select(group, lib.size, norm.factors, subject, sex, time, sets, sample) %>%
    ungroup() %>%
    mutate(time = factor(time, levels = c("w0", "w2pre", "w2post", "w12")), 
           sets = factor(sets, levels = c("single", "multiple"))) %>%
    `rownames<-` (.$sample) %>%
    dplyr::select(-sample) %>%
    data.frame()
  
  
  
}



#### Filter away low expression genes #############

# Create logCPM-plot from all methods #



# Calculate L and M values used for cutoff in each method #

LM.fun <- function(lib){
  
  L <- mean(lib) * 1e-6
  M <- median(lib) * 1e-6
  
  return(data.frame(L = L, M = M))
  
}


# These df combines information about average CPM from each method together 
# with a threshold for low expression genes (cutoff = 10/M + 2/L) where
# L is the mean library size and M is the median.

aCPM.cutoff <- list()


for(i in 1:length(dge_lists)) {
  
  aCPM.cutoff[[i]] <- data.frame(method = names(dge_lists[i]), 
                                 LM.fun(dge_lists[[i]]$samples$lib.size))
  
}

aCPM.cutoff <- bind_rows(aCPM.cutoff) %>%
  mutate(cutoff = log2(10/M + 2/L)) %>% 
  print()



aCPM <- list()

for(i in 1:length(dge_lists)) {
  
  aCPM[[i]] <- data.frame(logCPM = aveLogCPM(dge_lists[[i]]), 
                          method = names(dge_lists[i]))
}

 aCPM <- bind_rows(aCPM)



# Plot of distribution of average log2-CPM. The distribution should be bimodal with a low-abundance peak and
# high abundance peak representing expressed genes. The threshold should separate the two peaks

# Density plot prior to filtering  
denPlot1 <- aCPM %>%
  ggplot(aes(x = logCPM, fill = method, group = method)) + 
  geom_segment(data = aCPM.cutoff, aes(x = cutoff, y = 0, xend = cutoff, yend = 1), 
               color="red", linetype="dashed") +
  geom_density(alpha = 0.2, stat = "density") + scale_x_continuous(limits = c(-5, 15), expand = c(0,0)) +
  facet_wrap(~ method, scales = "free")


###### Filtering by expression edgeR ################
# Note: this also adds the design matrix to the list.
# The design matrix must be explicitly retrieved in downstream functions,
# e.g. design = dge_lists[[i]]$design


for(i in 1:length(dge_lists)) {
  
  
  dge_lists[[i]]$design <- model.matrix(~ time + time:sets, data = dge_lists[[i]]$samples)
  
  
  keep.exprs <- filterByExpr(dge_lists[[i]], design = dge_lists[[i]]$design)
  
  dge_lists[[i]] <- dge_lists[[i]][keep.exprs, ,keep.lib.sizes = FALSE]
  
}



aCPM.postfilter <- list()

for(i in 1:length(dge_lists)) {
  
  aCPM.postfilter[[i]] <- data.frame(logCPM = aveLogCPM(dge_lists[[i]]), 
                          method = names(dge_lists[i]))
}

aCPM.postfilter <- bind_rows(aCPM.postfilter)


logCPM_data <- rbind(data.frame(aCPM, filtering = "pre"), 
                     data.frame(aCPM.postfilter, filtering = "post"))


# Save data for preliminary analysis 

saveRDS(logCPM_data, file = "./data/derivedData/prel_analysis/logCPM_data.RDS")





######### Calculate normalization factor ###############


for(i in 1:length(dge_lists)) {
  
  dge_lists[[i]] <- calcNormFactors(dge_lists[[i]], method = "TMM")

}


### Distrubutions of different normalization factors per method

nf <- list() 

for(i in 1:length(dge_lists)) {
  
  nf[[i]] <- data.frame(nf = dge_lists[[i]]$samples$norm.factors, method = names(dge_lists[i]), 
                        subject = dge_lists[[i]]$samples$subject, 
                        time = dge_lists[[i]]$samples$time, 
                        sets = dge_lists[[i]]$samples$sets)
  
}



nf <- bind_rows(nf)

saveRDS(nf, file = "./data/derivedData/prel_analysis/norm_factors.RDS")



#### Exploratory plots #####
i <- 1

lcpm <- cpm(dge_lists[[i]], log = TRUE)

## From glimma package 

glMDSPlot(lcpm, labels=paste(dge_lists[[i]]$samples$sex), 
          groups=paste(dge_lists[[i]]$samples$time, sep = "_"),
          launch=TRUE)


# save DGELists for other methods #######################

saveRDS(dge_lists, file = "./data/derivedData/dge_lists/dge_list.RDS")






