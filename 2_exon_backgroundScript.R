library(readr)
library(here)
library(dplyr)
library(ggplot2)
directorio <- here()


# load data ---------------------------------------------------------------
#exones intrones, inferidos desde trinotate supertranscripts
gtf_trinity_exons <- read_delim(paste0(directorio,"/data/gtf_trinity_exons.csv"), 
                                delim = "\t", escape_double = FALSE, 
                                col_types = cols(exon_start = col_integer(), 
                                exon_end = col_integer()), trim_ws = TRUE)


Dinv_CNS_match_Dlaeve <- read_delim(paste0(directorio,"/data/Dinv_CNS_match_Dlaeve.txt"), 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, col_types = cols(X2 = col_skip(), 
                                                                        X3 = col_integer(), X4 = col_integer(), 
                                                                        X5 = col_integer(), X6 = col_skip(), 
                                                                        X7 = col_skip(), X8 = col_skip(), 
                                                                        X9 = col_skip()), trim_ws = TRUE)

colnames(Dinv_CNS_match_Dlaeve) <- c("qseqid", "length", "qlen", "slen", "evalue")

Dret_CNS_match_Dlaeve <- read_delim(paste0(directorio,"/data/Dret_CNS_match_Dlaeve.txt"), 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, col_types = cols(X2 = col_skip(), 
                                                                        X3 = col_integer(), X4 = col_integer(), 
                                                                        X5 = col_integer(), X6 = col_skip(), 
                                                                        X7 = col_skip(), X8 = col_skip(), 
                                                                        X9 = col_skip()), trim_ws = TRUE)



colnames(Dret_CNS_match_Dlaeve) <- c("qseqid", "length", "qlen", "slen", "evalue")

StrandedAssGCcontent <- read_table(paste0(directorio,"/data/StrandedAssGCcontent"), 
                                   col_types = cols(content = col_skip()))

colnames(StrandedAssGCcontent) <- c("transcript_id", "gc_content")

### rfam results for non coding rna families
rfam <- read_table("data/dlaevis.tbl", col_names = FALSE, 
                   col_types = cols(X2 = col_skip(), X5 = col_skip(), 
                                    X6 = col_skip(), X7 = col_skip(), 
                                    X8 = col_skip(), X9 = col_skip(), 
                                    X12 = col_skip(), X14 = col_skip(), 
                                    X15 = col_skip(), X17 = col_skip(), 
                                    X18 = col_skip(), X19 = col_skip()), 
                   skip = 2)

rfam <- rfam %>% filter(str_detect(X1, "TRINITY"))

colnames(rfam) <- c("transcript_id", "nonCodingtype", "rfamily", "rfam_strand", "rfam_truncated", "gc", "rfam_evalue")
rfam <- rfam %>% group_by(transcript_id) %>% filter(row_number(desc((rfam_evalue)))== 1, rfam_evalue< 1e-05)


#get the number of exons per transcript
#
gtf_trinity_exons <- gtf_trinity_exons %>% group_by(transcript_id) %>% mutate(number_of_exons = n(), exon_length = (1 + exon_end - exon_start), mean_exon_length = mean(exon_length)) %>%
  ungroup() %>% group_by(gene_id) %>% mutate(n_isoforms = n_distinct(transcript_id), gene_length = max(exon_end)) 

gtf_trinity_transcripts <- gtf_trinity_exons %>% ungroup() %>% group_by(transcript_id) %>% 
  summarise(gene_id = unique(gene_id),
            transcript_id = unique(transcript_id),
            n_isoforms = unique(n_isoforms),
            transcript_length = sum(exon_length),
            gene_length = unique(gene_length),
            mean_exon_length = unique(mean_exon_length),
            number_of_exons = unique(number_of_exons)
            ) %>% ungroup() %>% group_by(gene_id) %>% mutate(all_isoforms = sum(transcript_length), 
                                                             coverage = all_isoforms/(n_isoforms*gene_length))

#join with transcript-level classification: non-coding, swissprot, etc..
gtf_trinity_transcripts <- inner_join(transcript_classification, gtf_trinity_transcripts) %>% 
  select(transcript_id, gene_id, rna_type, rna_type_2, rna_type_joint, total_length, transcript_length, gene_length, mean_exon_length, number_of_exons, coverage)

#join with RNA-seq data from amputation experiments
gtf_trinity_transcripts <- select(results_full_amp_degs, transcript_id, baseMean, log2FoldChange, padj) %>% right_join(gtf_trinity_transcripts)

#join with strandedness info
gtf_trinity_transcripts <- left_join(gtf_trinity_transcripts,ss_analysisAMP)



# Join Transcripts with additional sources of data ------------------------

consolidated_transcript_df <- inner_join(gtf_trinity_transcripts, StrandedAssGCcontent)

#join with the blast data for conservation compared to invadens and reticulatum transcriptomes
consolidated_transcript <- Dinv_CNS_match_Dlaeve %>% transmute(transcript_id = qseqid, cons = 100 * length /qlen) %>%
  group_by(transcript_id) %>% summarise(Invadens_conservation = max(cons)) %>% 
  right_join(consolidated_transcript_df)


consolidated_transcript <- Dret_CNS_match_Dlaeve %>% transmute(transcript_id = qseqid, cons = 100 * length /qlen) %>%
  group_by(transcript_id) %>% summarise(Reticulatum_conservation = max(cons)) %>% 
  right_join(consolidated_transcript)

#NAs can be interpreted as zero conservation
consolidated_transcript <- consolidated_transcript %>% replace_na(list(Reticulatum_conservation = 0, Invadens_conservation = 0))

consolidated_transcript <- left_join(consolidated_transcript, rfam)

consolidated_transcript <- consolidated_transcript %>% mutate(rna_type_2 = coalesce(rna_type_2, nonCodingtype))

#Load the maker annotation data
exonByTxMaker <- read_csv("data/exonByTx.csv", 
                          col_types = cols(n_exons = col_integer(), 
                                           fivep = col_integer(), threep = col_integer()))
exonByTxMaker <- separate(exonByTxMaker, col = Trinity, sep = "to ", into = c("trash", "transcript_id")) 
exonByTxMaker$trash <- NULL

exonByTxMaker <- exonByTxMaker %>% select(transcript_id, exons, fivep, threep, n_exons, Swissprot) %>% left_join(consolidated_transcript)
# Dimensionality reduction and clustering ---------------------------------
############# Dimensionality reduction experiments
matKmeans <- consolidated_transcript %>% 
  transmute(Reticulatum_conservation, Invadens_conservation, log1p(baseMean), log1p(total_length), log1p(mean_exon_length), number_of_exons, coverage, gc_content) %>% 
  scale()


kn <- 1:20
wo <- sapply(kn, function(k){wa <-kmeans(matKmeans, centers = k, iter.max = 10)
return(wa$tot.withinss)
  })

clustering <- kmeans(matKmeans, centers = 10, iter.max = 10)

pca <- prcomp(matKmeans, scale = F, center = F)
pca_df <- as.data.frame(pca$x)
pca_df$cluster <- clustering$cluster
pca_df$type <- consolidated_transcript$rna_type

consolidated_transcript$cluster <- clustering$cluster
library(GGally)
parallel_df <- as.data.frame(clustering$centers)
parallel_df$cluster <- rownames(parallel_df)
parallel_df <- clean_names(parallel_df)

# Save dataframe ----------------------------------------------------------
write_csv(consolidated_transcript, "data/transcriptome_consolidated.csv")
write_csv(exonByTxMaker, "data/maker_exons.csv")

