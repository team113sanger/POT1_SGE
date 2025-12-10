#DEPENDENCIES - START
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(DEGreport)
library(pheatmap)
library(GGally)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(cowplot)
library(ggpubr)
library(ggridges)
library(mclust)
#DEPENDENCIES - END


###############----------------------------------------------------###########################################################
###############-----------COUNTS_START-----------------------------###########################################################
###############----------------------------------------------------###########################################################

#DIRECTORY ORGANISATION AND READING IN THE DATA

#Data preparation
setwd("/Users/aw28/Documents/POT1_analysis")
#######SGE READ PREPARATION 
#makes lists of file locations to unzip, then run the output list of commands using termincal
gen.count.locations <- function(exon,sg) {
  paste0("cd /Users/aw28/Documents/POT1_analysis/POT1_SGE_organised/pot1_exon_",exon,"_",sg,"/results/pycroquet
gunzip *.query_to_library_counts.tsv.gz")
}
#Function to apply the above function to listed exons and designated sgRNA_A
a_locs<-lapply (c("2","3_2","4","6","7","8","10","11","12"), function(y) {gen.count.locations(exon=y, sg="a")
})
#Function to apply the above function to listed exons and designated sgRNA_B
b_locs<-lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(y) {gen.count.locations(exon=y, sg="b")
})
#make a file of all the files directories to gunzip
lapply(a_locs, write, "./a_locs.txt", append=TRUE)
lapply(b_locs, write, "./b_locs.txt", append=TRUE)

#######PLASMID READ PREPARATION 
#in terminal
'cd /Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised'
'gunzip *.query_to_library_counts.tsv.gz'






#reads need to be merged for plasmid libraries for 3_2a/b, 4a/b and 6a/b - import combine then export

one_3_2_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/3_2a_lib__44030_1_101.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_3_2_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/3_2a_lib__44030_2_101.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_3_2_a<-merge(one_3_2_a,two_3_2_a,by="X.id")
plasmid_3_2_a<-plasmid_3_2_a %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")

one_3_2_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/3_2b_lib__44030_1_102.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_3_2_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/3_2b_lib__44030_2_102.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_3_2_b<-merge(one_3_2_b,two_3_2_b,by="X.id")
plasmid_3_2_b<-plasmid_3_2_b %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")

one_4_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/4a_lib__44030_1_69.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_4_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/4a_lib__44030_2_69.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_4_a<-merge(one_4_a,two_4_a,by="X.id")
plasmid_4_a<-plasmid_4_a %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")

one_4_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/4b_lib__44030_1_70.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_4_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/4b_lib__44030_2_70.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_4_b<-merge(one_4_b,two_4_b,by="X.id")
plasmid_4_b<-plasmid_4_b %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")

one_6_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/6a_lib__44030_1_31.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_6_a<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/6a_lib__44030_2_31.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_6_a<-merge(one_6_a,two_6_a,by="X.id")
plasmid_6_a<-plasmid_6_a %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")

one_6_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/6b_lib__44030_1_32.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)
two_6_b<-read.table("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/6b_lib__44030_2_32.query_to_library_counts.tsv",skip = 2, header =T, comment.char = "", fill=TRUE)

plasmid_6_b<-merge(one_6_b,two_6_b,by="X.id")
plasmid_6_b<-plasmid_6_b %>% select(c(1,5,9)) %>% mutate(PLASMID = rowSums(.[2:3])) %>% select(c(1,4)) #%>% rename(id = "X.id")


#Function to read in separate count.csv files, $PATH will need to be changed. 
read.counts <- function(exon,sg) {
  #directories 
  dir=paste0("/Users/aw28/Documents/POT1_analysis/POT1_SGE_organised/pot1_exon_",exon,"_",sg,"/results/pycroquet/")
  dir_plasmid=paste0("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/")
  out=paste0("/Users/aw28/Documents/POT1_analysis/", "E", exon, "_SG", sg, "_count_frame.csv")
  #makes a list of count file names held in the dir 
  temp = list.files(dir, pattern="*rep_1.query_to_library_counts.tsv|rep_2.query_to_library_counts.tsv|rep_3.query_to_library_counts.tsv")
  temp_plasmid = list.files(dir_plasmid, pattern=paste0("Exon","_",exon,"_",sg,"_merged.query_to_library_counts.tsv"))
  #makes a list of full paths to each file
  file_all_count=paste0(dir,temp)
  file_all_count_plasmid=paste0(dir_plasmid,temp_plasmid)
  #reads the data from these paths into a list of dataframes
  myfiles = lapply(file_all_count, read.table, skip = 2, header =T, comment.char = "")
  myfiles_plasmid = lapply(file_all_count_plasmid, read.table, skip = 2, header =T, comment.char = "", fill=TRUE)
  myfiles_plasmid = lapply(myfiles_plasmid,"[",c(1,5))
  #makes a large dataframe with all counts - n.b. be careful of the ordering of the new column names make sure agrees to the order in the lists above - this version is ascending numerical
  df <- myfiles %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("rep"))
  df_plasmid <- myfiles_plasmid %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("merged"))
  df<- df %>% left_join(df_plasmid, by=c("X.id"))
  df <- df %>% `colnames<-` (c("id", "D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3","D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","PLASMID")) %>% replace(is.na(.), 0)
  #saves dataframe to output location and re-imports it as a named df with names taken from definitions above
  write.csv(df, file=out, row.names = FALSE)
  x=paste0("count_table_", "E", exon, "_SG", sg)
  assign(x,read.csv(out), envir = .GlobalEnv)
}

###############-----------RUN FOR ALL LIBRARIES_STRART-------------###########################################################
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("2","8","10","11","12"), function(y) {read.counts(exon=y, sg="a")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("1","5","9_1","9_2","13"), function(y) {read.counts(exon=y, sg="b")
})
###############-----------RUN FOR ALL LIBRARIES_END----------------###########################################################


#run 3_2a/b, 4a/b and 6a/b separately as plasmids (and runs) from merged lanes 
#run 3_2a/b, 4a/b and 6a/b separately as plasmids (and runs) from merged lanes
#run 3_2a/b, 4a/b and 6a/b separately as plasmids (and runs) from merged lanes

#Function to read in separate count.csv files, $PATH will need to be changed - WITHOUT PLASMID COMBINATION
read.counts <- function(exon,sg) {
  #directories 
  dir=paste0("/Users/aw28/Documents/POT1_analysis/POT1_SGE_organised/pot1_exon_",exon,"_",sg,"/results/pycroquet/")
  #dir_plasmid=paste0("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/")
  out=paste0("/Users/aw28/Documents/POT1_analysis/", "E", exon, "_SG", sg, "_count_frame.csv")
  #makes a list of count file names held in the dir 
  temp = list.files(dir, pattern="*query_to_library_counts.tsv|query_to_library_counts.tsv|query_to_library_counts.tsv")
  #temp_plasmid = list.files(dir_plasmid, pattern=paste0("Exon","_",exon,"_",sg,"_merged.query_to_library_counts.tsv"))
  #makes a list of full paths to each file
  file_all_count=paste0(dir,temp)
  #file_all_count_plasmid=paste0(dir_plasmid,temp_plasmid)
  #reads the data from these paths into a list of dataframes
  myfiles = lapply(file_all_count, read.table, skip = 2, header =T, comment.char = "")
  #myfiles_plasmid = lapply(file_all_count_plasmid, read.table, skip = 2, header =T, comment.char = "", fill=TRUE)
  #myfiles_plasmid = lapply(myfiles_plasmid,"[",c(1,5))
  #makes a large dataframe with all counts - n.b. be careful of the ordering of the new column names make sure agrees to the order in the lists above - this version is ascending numerical
  df <- myfiles %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("rep"))
  df<-df %>%  mutate(D10R1 = rowSums(.[2:3])) %>%
    mutate(D10R2 = rowSums(.[4:5])) %>%
    mutate(D10R3 = rowSums(.[6:7])) %>%
    mutate(D14R1 = rowSums(.[8:9])) %>%
    mutate(D14R2 = rowSums(.[10:11])) %>%
    mutate(D14R3 = rowSums(.[12:13])) %>%
    mutate(D21R1 = rowSums(.[14:15])) %>%
    mutate(D21R2 = rowSums(.[16:17])) %>%
    mutate(D21R3 = rowSums(.[18:19])) %>%
    mutate(D4R1 = rowSums(.[20:21])) %>%
    mutate(D4R2 = rowSums(.[22:23])) %>%
    mutate(D4R3 = rowSums(.[24:25])) %>%
    mutate(D7R1 = rowSums(.[26:27])) %>%
    mutate(D7R2 = rowSums(.[28:29])) %>%
    mutate(D7R3 = rowSums(.[30:31])) %>% select(c(1,32:46))
  #df_plasmid <- myfiles_plasmid %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("merged"))
  #df<- df %>% left_join(df_plasmid, by=c("X.id"))
  plasmid_file<-paste0("plasmid_",exon,"_",sg)
  plasmid_file_actual <- get(plasmid_file)
  df<- df %>% left_join(plasmid_file_actual, by=c("X.id"))
  df <- df %>% `colnames<-` (c("id", "D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3","D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","PLASMID")) %>% replace(is.na(.), 0)
  #saves dataframe to output location and re-imports it as a named df with names taken from definitions above
  write.csv(df, file=out, row.names = FALSE)
  x=paste0("count_table_", "E", exon, "_SG", sg)
  assign(x,read.csv(out), envir = .GlobalEnv)
}

###############-----------RUN FOR ALL LIBRARIES_STRART-------------###########################################################
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("3_2","4","6"), function(y) {read.counts(exon=y, sg="a")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("3_2","4","6"), function(y) {read.counts(exon=y, sg="b")
})
###############-----------RUN FOR ALL LIBRARIES_END----------------###########################################################

#3_1_a, 14_b have no plasmid available, 7, 0 targetons need to be added by Sofia 
#3_1_a, 14_b have no plasmid available, 7, 0 targetons need to be added by Sofia 
#3_1_a, 14_b have no plasmid available, 7, 0 targetons need to be added by Sofia 

#Function to read in separate count.csv files, $PATH will need to be changed. 
read.counts <- function(exon,sg) {
  #directories 
  dir=paste0("/Users/aw28/Documents/POT1_analysis/POT1_SGE_organised/pot1_exon_",exon,"_",sg,"/results/pycroquet/")
  #dir_plasmid=paste0("/Users/aw28/Documents/POT1_analysis/POT1_PLASMID_organised/")
  out=paste0("/Users/aw28/Documents/POT1_analysis/", "E", exon, "_SG", sg, "_count_frame.csv")
  #makes a list of count file names held in the dir 
  temp = list.files(dir, pattern="*rep_1.query_to_library_counts.tsv|rep_2.query_to_library_counts.tsv|rep_3.query_to_library_counts.tsv")
  #temp_plasmid = list.files(dir_plasmid, pattern=paste0("Exon","_",exon,"_",sg,"_merged.query_to_library_counts.tsv"))
  #makes a list of full paths to each file
  file_all_count=paste0(dir,temp)
  #file_all_count_plasmid=paste0(dir_plasmid,temp_plasmid)
  #reads the data from these paths into a list of dataframes
  myfiles = lapply(file_all_count, read.table, skip = 2, header =T, comment.char = "")
  #myfiles_plasmid = lapply(file_all_count_plasmid, read.table, skip = 2, header =T, comment.char = "", fill=TRUE)
  #myfiles_plasmid = lapply(myfiles_plasmid,"[",c(1,5))
  #makes a large dataframe with all counts - n.b. be careful of the ordering of the new column names make sure agrees to the order in the lists above - this version is ascending numerical
  df <- myfiles %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("rep"))
  #df_plasmid <- myfiles_plasmid %>% purrr::reduce(full_join, by=c("X.id")) %>% dplyr::select(contains("X.id")|contains("merged"))
  #df<- df %>% left_join(df_plasmid, by=c("X.id"))
  df <- df %>% `colnames<-` (c("id", "D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3","D4R1","D4R2","D4R3","D7R1","D7R2","D7R3")) %>% replace(is.na(.), 0)
  #saves dataframe to output location and re-imports it as a named df with names taken from definitions above
  write.csv(df, file=out, row.names = FALSE)
  x=paste0("count_table_", "E", exon, "_SG", sg)
  assign(x,read.csv(out), envir = .GlobalEnv)
}

###############-----------RUN FOR ALL LIBRARIES_STRART-------------###########################################################
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("3_1"), function(y) {read.counts(exon=y, sg="a")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("14"), function(y) {read.counts(exon=y, sg="b")
})


###############----------------------------------------------------###########################################################
###############-----------COUNTS_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_START---------###########################################################
###############----------------------------------------------------###########################################################

#Define directories where VaLiAnT VCFs are stored
dir_vcf=paste0("/Users/aw28/Documents/POT1_analysis/february_2024/repeat/repeats_exclusions/valiant_repeat_april_24/vcfs/")
#vcf_out=paste0("/Users/aw28/Documents/POT1_analysis/vep_input_vcfs/")
vcf_out=paste0("/Users/aw28/Documents/POT1_analysis/february_2024/repeat/repeats_exclusions/valiant_repeat_april_24/vep_input_vcfs/")

#makes a list of count file names held in the dir 
temp_vcf = list.files(dir_vcf, pattern="*.vcf")
#makes a list of full paths to each file
file_all_vcf=paste0(dir_vcf,temp_vcf)
#reads the data from these paths into a list of dataframes
my_vcf_files = sapply(file_all_vcf, read.table, skip = 9, header =F, comment.char = "#", simplify=FALSE)
#names columns
colnames_vcf <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
my_vcf_files = lapply(my_vcf_files, setNames, colnames_vcf)
#gets list of directory names so one knows which targetons are which dataframe
vcf_names <-names(my_vcf_files)
#gets file names from full paths
vcf_names_base <-basename(vcf_names)
rm(vcf_names)

#########################################################################
#FUNCTIONS
#########################################################################

#Function to extract the oligo ID from the INFO field and place in the ID column
wrangle.vcf<-function(number) {
  tmp_vcf_data<-my_vcf_files[[number]]
  tmp_vcf_data[c('INFO1', 'INFO2', 'INFO3')] <- str_split_fixed(tmp_vcf_data$INFO, ';', 3)
  tmp_vcf_data$INFO2<-str_remove(tmp_vcf_data$INFO2, "SGE_OLIGO=")
  tmp_vcf_data[ ,c('INFO1', 'INFO3', 'ID')] <- list(NULL)
  tmp_vcf_data$ID = tmp_vcf_data$INFO2
  tmp_vcf_data$INFO2 = NULL
  tmp_vcf_data <- tmp_vcf_data %>% relocate(ID, .before = REF)
  colnames_vcf_out <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  setNames(tmp_vcf_data, colnames_vcf_out)
}

#run the above function over all VCFs
my_vcf_files = sapply (c(1:34), function(z) {wrangle.vcf(number=z)
}, simplify = FALSE, USE.NAMES=TRUE)
#name the dataframes based on the original VCF file 
names(my_vcf_files) <- vcf_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_vcf<-my_vcf_files[[object_name]]
  suffix<-names(my_vcf_files)[object_name]
  write.table(table_output_vcf, file=paste0(vcf_out,"VEP_INPUT_", suffix), sep = "\t", row.names = FALSE, quote = FALSE)
}
#run the above function over all VCFs to output appropriately named VEP INPUTS 
lapply (c(1:34), function(zy) {make.tables(object_name=zy)
})

###############----------------------------------------------------###########################################################
###############-----------VCF_CLEANING_FOR_VEP_INPUT_END-----------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------COUNT_FILTERING_INPUT_START--------------###########################################################
###############----------------------------------------------------###########################################################

#directoy of metafiles - ONLY A at the moment (expand 13 to 26 if want to do all)
dir_meta=paste0("/Users/aw28/Documents/POT1_analysis/february_2024/repeat/repeats_exclusions/valiant_repeat_april_24/meta/")
#makes a list of count file names held in the dir 
temp_meta = list.files(dir_meta, pattern="*.csv")
#makes a list of full paths to each file
file_all_meta=paste0(dir_meta,temp_meta)
#reads the data from these paths into a list of dataframes
my_meta_files = sapply(file_all_meta, read.csv, header =T, comment.char = "#", simplify=FALSE)
#meta names, full paths
meta_names <-names(my_meta_files)
#gets file names from full paths
meta_names_base <-basename(meta_names)
rm(meta_names)
#name the dataframes based on the original VCF file 
names(my_meta_files) <- meta_names_base
#function to output each VCF as a separate dataframe
make.tables<-function(object_name) {
  table_output_meta<-my_meta_files[[object_name]]
  identifier<-names(my_meta_files)[object_name]
  names(table_output_meta)[names(table_output_meta) == 'oligo_name'] <- 'id'
  assign(identifier, table_output_meta, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named metas into the global environment
lapply (c(1:34), function(zy) {make.tables(object_name=zy)
})

#merge meta and counts_SGA
E2_SGA_master<-merge(count_table_E2_SGa, chr7_124870880_124871127_minus_sgRNA_w2_a_meta.csv, by="id", all.x=FALSE)
E3_1_SGA_master<-merge(count_table_E3_1_SGa, chr7_124863272_124863520_minus_sgRNA_w3.1_a_meta.csv, by="id", all.x=FALSE)
E3_2_SGA_master<-merge(count_table_E3_2_SGa, chr7_124863468_124863717_minus_sgRNA_w3.2_a_meta.csv, by="id", all.x=FALSE)
E4_SGA_master<-merge(count_table_E4_SGa,chr7_124858890_124859138_minus_sgRNA_w4_a_meta.csv, by="id", all.x=FALSE)
E6_SGA_master<-merge(count_table_E6_SGa,chr7_124851781_124852030_minus_sgRNA_w6_a_meta.csv, by="id", all.x=FALSE)
E8_SGA_master<-merge(count_table_E8_SGa,chr7_124842757_124843006_minus_sgRNA_w8_a_meta.csv, by="id", all.x=FALSE)
E10_SGA_master<-merge(count_table_E10_SGa,chr7_124835214_124835463_minus_sgRNA_w10_a_meta.csv, by="id", all.x=FALSE)
E11_SGA_master<-merge(count_table_E11_SGa,chr7_124829216_124829465_minus_sgRNA_w11_a_meta.csv, by="id", all.x=FALSE)
E12_SGA_master<-merge(count_table_E12_SGa,chr7_124827112_124827359_minus_sgRNA_w12_a_meta.csv, by="id", all.x=FALSE)


#merge meta and counts_SGB

E1_SGB_master<-merge(count_table_E1_SGb,chr7_124892164_124892411_minus_sgRNA_w1_b_meta.csv, by="id", all.x=FALSE)
E3_2_SGB_master<-merge(count_table_E3_2_SGb,chr7_124863468_124863717_minus_sgRNA_w3.2_b_meta.csv, by="id", all.x=FALSE)
E4_SGB_master<-merge(count_table_E4_SGb,chr7_124858890_124859138_minus_sgRNA_w4_b_meta.csv, by="id", all.x=FALSE)
E5_SGB_master<-merge(count_table_E5_SGb,chr7_124852932_124853177_minus_sgRNA_w5_b_meta.csv, by="id", all.x=FALSE)
E6_SGB_master<-merge(count_table_E6_SGb,chr7_124851781_124852030_minus_sgRNA_w6_b_meta.csv, by="id", all.x=FALSE)
E9_1_SGB_master<-merge(count_table_E9_1_SGb,chr7_124840935_124841182_minus_sgRNA_w9.1_b_meta.csv, by="id", all.x=FALSE)
E9_2_SGB_master<-merge(count_table_E9_2_SGb,chr7_124840984_124841232_minus_sgRNA_w9.2_b_meta.csv, by="id", all.x=FALSE)
E13_SGB_master<-merge(count_table_E13_SGb,chr7_124825226_124825474_minus_sgRNA_w13_b_meta.csv, by="id", all.x=FALSE)
E14_SGB_master<-merge(count_table_E14_SGb,chr7_124823919_124824168_minus_sgRNA_w14_b_meta.csv, by="id", all.x=FALSE)


#E3_1_A and 14B do not have plasmid data available so add a column named plasmid with no counts
E14_SGB_master$PLASMID<-NA
E14_SGB_master <- E14_SGB_master %>% relocate(PLASMID, .before = species)
E3_1_SGA_master$PLASMID<-NA
E3_1_SGA_master <- E3_1_SGA_master %>% relocate(PLASMID, .before = species)

#remove rows that have fewer than 10 counts total across all 9 samples - OPTIONAL REMOVE VARIANTS IN CONSTANT REGTIONS
filter.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master")
  filtered_counts <- get(master_input)
  filtered_counts$rowsum <- rowSums(filtered_counts[, c("D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3")])
  filtered_counts<- filtered_counts %>% filter(rowsum>=10)
  #optional_remove high counts
  #filtered_counts<- filtered_counts %>% filter(rowsum<=10000)
  #OPTIONAL - REMOVE VARIANTS PAM MUTATION CODONS
  #filtered_counts<- filtered_counts %>% filter(pam_mut_sgrna_id %in% "" | is.na(pam_mut_sgrna_id))
  #OPTIONAL - REMOVE VARIANTS IN CONSTANT REGIONS - CAN REMOVE THIS STEP IF ADAPTERS STRIPPED SAME AS SEQUENCING PRIMERS
  filtered_counts<- filtered_counts %>% filter(vcf_var_in_const==0 | is.na(vcf_var_in_const))
  tablename<-paste0(master_input,"_filtered_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}


#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(y) {filter.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(y) {filter.counts(exon=y, sg="B")
})


###############----------------------------------------------------###########################################################
###############-----------COUNT_FILTERING_INPUT_END----------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------VEP_START--------------------------------###########################################################
###############----------------------------------------------------###########################################################

#mergeing vep outputs to create vep + meta
#Define directories where VEP OUTPUTS are stored
dir_vep=paste0("/Users/aw28/Documents/POT1_analysis/february_2024/repeat/repeats_exclusions/valiant_repeat_april_24/vep_outputs/")
#makes a list of count file names held in the dir 
temp_vep = list.files(dir_vep, pattern="*.tsv")
#makes a list of full paths to each file
file_all_vep=paste0(dir_vep,temp_vep)
#reads the data from these paths into a list of dataframes
my_vep_files = sapply(file_all_vep, read.table, skip = 92, header =T, comment.char = "", simplify=FALSE)
#meta names, full paths
vep_names <-names(my_vep_files)
#gets file names from full paths
vep_names_base <-basename(vep_names)
rm(vep_names)
#name the dataframes based on the original VEP file 
names(my_vep_files) <- vep_names_base
#function to output each VCF as a separate dataframe
make.vep.tables<-function(object_name) {
  table_output_vep<-my_vep_files[[object_name]]
  identifier<-names(my_vep_files)[object_name]
  names(table_output_vep)[names(table_output_vep) == 'X.Uploaded_variation'] <- 'id'
  table_output_vep <- table_output_vep %>% filter(Feature %in% "ENST00000357628.8")
  assign(paste0("VEP_OUTPUT_",identifier), table_output_vep, envir = .GlobalEnv)
}
#run the above function over all metas to output appropriately named metas into the global environment
lapply (c(1:34), function(zy) {make.vep.tables(object_name=zy)
})

#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGA 
E2_SGA_master_annotated<-merge(E2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124870880_124871127_minus_sgRNA_w2_a_pam.tsv, by="id", all.x=TRUE)
E3_1_SGA_master_annotated<-merge(E3_1_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863272_124863520_minus_sgRNA_w3.1_a_pam.tsv , by="id", all.x=TRUE)
E3_2_SGA_master_annotated<-merge(E3_2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863468_124863717_minus_sgRNA_w3.2_a_pam.tsv, by="id", all.x=TRUE)
E4_SGA_master_annotated<-merge(E4_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124858890_124859138_minus_sgRNA_w4_a_pam.tsv, by="id", all.x=TRUE)
E6_SGA_master_annotated<-merge(E6_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124851781_124852030_minus_sgRNA_w6_a_pam.tsv, by="id", all.x=TRUE)
E8_SGA_master_annotated<-merge(E8_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124842757_124843006_minus_sgRNA_w8_a_pam.tsv, by="id", all.x=TRUE)
E10_SGA_master_annotated<-merge(E10_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124835214_124835463_minus_sgRNA_w10_a_pam.tsv, by="id", all.x=TRUE)
E11_SGA_master_annotated<-merge(E11_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124829216_124829465_minus_sgRNA_w11_a_pam.tsv, by="id", all.x=TRUE)
E12_SGA_master_annotated<-merge(E12_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124827112_124827359_minus_sgRNA_w12_a_pam.tsv, by="id", all.x=TRUE)
#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGB
E1_SGB_master_annotated<-merge(E1_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124892164_124892411_minus_sgRNA_w1_b_pam.tsv, by="id", all.x=TRUE)
E3_2_SGB_master_annotated<-merge(E3_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863468_124863717_minus_sgRNA_w3.2_b_pam.tsv, by="id", all.x=TRUE)
E4_SGB_master_annotated<-merge(E4_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124858890_124859138_minus_sgRNA_w4_b_pam.tsv, by="id", all.x=TRUE)
E5_SGB_master_annotated<-merge(E5_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124852932_124853177_minus_sgRNA_w5_b_pam.tsv, by="id", all.x=TRUE)
E6_SGB_master_annotated<-merge(E6_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124851781_124852030_minus_sgRNA_w6_b_pam.tsv, by="id", all.x=TRUE)
E9_1_SGB_master_annotated<-merge(E9_1_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124840935_124841182_minus_sgRNA_w9.1_b_pam.tsv, by="id", all.x=TRUE)
E9_2_SGB_master_annotated<-merge(E9_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124840984_124841232_minus_sgRNA_w9.2_b_pam.tsv, by="id", all.x=TRUE)
E13_SGB_master_annotated<-merge(E13_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124825226_124825474_minus_sgRNA_w13_b_pam.tsv, by="id", all.x=TRUE)
E14_SGB_master_annotated<-merge(E14_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124823919_124824168_minus_sgRNA_w14_b_pam.tsv, by="id", all.x=TRUE)


###############----------------------------------------------------###########################################################
###############-----------VEP_END----------------------------------###########################################################
###############----------------------------------------------------###########################################################

###############----------------------------------------------------###########################################################
###############-----------DESEQ2_START-----------------------------###########################################################
###############----------------------------------------------------###########################################################

#filter the above tables to get normalization tables
normalization.counts<-function(exon, sg){
  master_input<-paste0("E",exon,"_SG",sg,"_master_annotated")
  filtered_counts <- get(master_input)
  filtered_counts <-filtered_counts[!duplicated(filtered_counts$mseq), ]
  filtered_counts<- filtered_counts %>% filter(Consequence %in% "synonymous_variant" | Consequence %in% "intron_variant")
  tablename<-paste0("E",exon,"_SG",sg,"_normalization_counts")
  assign(tablename, filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(y) {normalization.counts(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(y) {normalization.counts(exon=y, sg="B")
})



#make normalization matricies
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(y) {deseq.norm.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(y) {deseq.norm.input(exon=y, sg="B")
})


#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D7R3","D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}

#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(y) {deseq.count.input(exon=y, sg="A")
})
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(y) {deseq.count.input(exon=y, sg="B")
})



#read in experimental designs SGA
annot_a<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_continuous.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)



#########################################################################
#FUNCTIONS_DESEQ_START
#########################################################################

#########################################################################
#OUTPUT FUNCTIONS START
#########################################################################

#PCA SCREE PLOT
plotPCA.hk <- function (object, intgroup = "condition", ntop = 500, pc_1 = 1, pc_2 = 2, returnData = FALSE, scree = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc_1], PC2 = pca$x[, pc_2], group = group,
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc_1:pc_2]
    return(d)
  }
  if (scree) {
    xx<- barplot(round(percentVar, digits = 2)*100, names.arg=c(1:length(percentVar)),xlab="PC",ylab="% Variance",ylim=c(0,100), main="Scree Plot")
    text(x = xx, y = round(percentVar, digits = 4)*100, label = round(percentVar, digits = 4)*100, pos = 3, cex = 0.8, col = "black")
  }
  else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color="condition", shape="type", label = "sample")) + scale_shape_manual(values=seq(0,127)) + geom_point(size = 3) + xlab(paste0("PC",pc_1,": ", round(percentVar[pc_1] * 100, digits = 2), "% variance")) + ylab(paste0("PC",pc_2,": ", round(percentVar[pc_2] * 100, digits=2), "% variance")) + coord_fixed() + geom_text_repel(size=3)
  }
}

#needed for scatter pairs plots
ggpairs_ext <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) +
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt)
}

# Function to estimate size factors (typically run for control oligos only)
estimate_control_size_factors <- function ( countData = countData, colData = colData, design = design, minRowSum = 10, ref = NULL ) {
  print( "Estimating control size factors...")
  dds <- DESeqDataSetFromMatrix( countData = countData, 
                                 colData = colData, 
                                 design = as.formula( design ) )
  
  dds <- dds[ rowSums( counts( dds ) ) > minRowSum, ]
  
  if ( ! is.null( ref ) ) {
    dds$condition <- relevel( dds$condition, ref = ref )
  }
  
  control_size_factors <- sizeFactors( estimateSizeFactors( dds ) )
  
  return( control_size_factors )
}

#function to add text to column names
appendDataFrameColumns<-function(df, prefix='', suffix='', sep='')
{
  colnames(df) <- paste(prefix, colnames(df), suffix, sep=sep)
  
  return(df)
}

#########################################################################
#OUTPUT FUNCTIONS END
#########################################################################

#the main deseq funtion - setwd() to where you want all the files to be outputted
run.deseq<-function(exon_assay, sg_assay){
  x_normalization_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_normalization_counts_MATRIX")
  x_filtered_counts = paste0("E",exon_assay,"_","SG",sg_assay,"_master_annotated_filtered_counts_MATRIX")
  x_normalization_counts=get(x_normalization_counts)
  x_filtered_counts =get(x_filtered_counts )
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a,
                                        ref = "D4",
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a, design = ~condition)
  dds$condition <- factor(dds$condition, levels=c("D4", "D7", "D10","D14","D21"))
  dds$condition <- relevel(dds$condition, ref = "D4")
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq
  dds <- DESeq(dds)
  rld <- rlog(dds)
  res<-results(dds)
  res<-as.data.frame(res)
  z_score <- assay(rld) %>% as.matrix() %>% t() %>% scale() %>% t() %>% as.data.frame()
  colnames(z_score) <- paste0(colnames(z_score), "_z_score")
  z_score$mseq <- row.names(z_score)  
  row.names(z_score) <- NULL
  z_score <-  z_score %>% select(mseq, everything())
  #shirnkage type set to normal, apeglm would be better for RNASeq - but have selected no shrinkage is table
  table_wald <- degComps(dds, combs = "condition", contrast = list("condition_D7_vs_D4","condition_D10_vs_D4" , "condition_D14_vs_D4", "condition_D21_vs_D4"), alpha = 0.05, skip = FALSE, type = "normal", pairs = FALSE, fdr = "default")
  
  #Summary Table
  #production of summary table - RAW - change to shrunken if want LFC shrinkage - apeglm is most appropriate for RNASeq
  summary <- purrr::reduce(c(deg(table_wald[[1]], "raw") %>% appendDataFrameColumns(suffix="_D4_D7") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[2]], "raw") %>% appendDataFrameColumns(suffix="_D4_D10") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[3]], "raw") %>% appendDataFrameColumns(suffix="_D4_D14") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list(),
                             deg(table_wald[[4]], "raw") %>% appendDataFrameColumns(suffix="_D4_D21") %>% as.data.frame() %>% rownames_to_column(var="mseq") %>% list()), left_join, by="mseq") %>% select(-c("baseMean_D4_D10","baseMean_D4_D14","baseMean_D4_D21")) %>% dplyr::rename(baseMean = "baseMean_D4_D7") %>% data.frame()
  #bind the z_scores to the summary table
  summary<-merge(summary, z_score, by="mseq", all.x=TRUE)
  
  #OUTPUT PLOTS AND TABLES 
  ident<-paste0("E",exon_assay,"_","SG",sg_assay)
  #SCREE
  pdf(paste0(ident,"_scree.pdf"))
  plotPCA.hk(rld,intgroup=c("condition", "type"), returnData=FALSE,pc_1=1, pc_2=2, scree=TRUE)
  dev.off()
  #DISPERSION
  pdf(paste0(ident,"_dispersion.pdf"))
  plotDispEsts(dds, ylim =c(1e-4,2e1))
  dev.off()
  #HEATMAP
  sampleDistMatrix <- as.matrix( dist( t( assay(rld) ) ) )
  pdf(paste0(ident,"_heatmap.pdf"))
  pheatmap(sampleDistMatrix, trace="none", col=colorRampPalette(rev(brewer.pal(9, "Blues")) )(255), adjRow = c(1,1))
  dev.off()
  remove(sampleDistMatrix)
  ###Scatter plot, checking replicate consistency
  pdf(paste0(ident,"_scatter_matrix.pdf"), width=8, height=8)
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D4")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D7")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D10")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count" ) + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D14")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  print(ggpairs(as.data.frame(assay(rld)) %>% select(starts_with("D21")), lower = list(continuous = wrap(ggpairs_ext, pts=list(size=0.4,colour="black"), smt=list(method="lm", se=F, size=0.2, colour="blue"))), title = "Regularized-log Transformed Read Count") + theme(plot.title = element_text(hjust = 0.5, size=13, face="bold")))
  dev.off()
  
  #PCA PLOTS
  pcaData <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pca<-ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()+
    theme_classic()
  ggsave(paste0(ident,"_PCA.pdf"), pca, height=6, width=8)
  
  #Export rld
  write.table(as.data.frame(assay(rld)), sep="\t",file=paste0(ident,"_rld.txt"), col.names=NA)
  #Export normalized read count
  write.table(counts(dds,normalized=TRUE), sep="\t",file=paste0(ident,"_norm_count.txt"), col.names=NA)
  #Export size factor
  write.table(dds@assays@data@listData %>% as.data.frame(),sep="\t",file=paste0(ident,"_normalization_table.txt"), col.names=NA)
  #Export full table
  write.table(dds@rowRanges@elementMetadata@listData %>% as.data.frame() ,sep="\t",file=paste0(ident,"_disper_table.txt"), col.names=NA)
  
  #continuous DESEQ to add SGE RATE 
  # Get control size factors
  SF <- estimate_control_size_factors(  countData = x_normalization_counts,
                                        colData = annot_a_continuous,
                                        design = "~ condition")
  #Prep DESEQ
  dds <- DESeqDataSetFromMatrix(countData = x_filtered_counts, colData = annot_a_continuous, design = ~condition)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds) <- SF
  #Run Deseq - SHRUNKEN LFC - CHANGE TO table_wald[[2]] if want shrunken LFC, type="apgelm" is appropriate for RNASeq perhaps not SGE
  dds <- DESeq(dds)
  table_wald <- degComps(dds, combs = "condition", alpha = 0.05, skip = FALSE, type = "normal", pairs = FALSE, fdr = "default")
  rate <- as.data.frame(table_wald[[1]])%>%rownames_to_column(var="mseq")%>%appendDataFrameColumns(suffix="_continuous")
  colnames(`rate`)[colnames(`rate`) == "mseq_continuous"] <- "mseq"
  summary<-merge(summary, rate, by="mseq", all.x=TRUE)
  
  
  #Substract the LFC median of the synonymous_variant and intron_variant
  adj<- as.data.frame(x_normalization_counts)
  adj<-adj %>% mutate(mseq=rownames(adj)) %>% select(mseq, everything()) %>% remove_rownames() %>% select(1)
  adj <- adj %>% left_join(summary %>%select("mseq","log2FoldChange_D4_D7","log2FoldChange_D4_D10","log2FoldChange_D4_D14","log2FoldChange_D4_D21", "log2FoldChange_continuous")) %>% summarise(median_D4_D7=median(log2FoldChange_D4_D7),median_D4_D10=median(log2FoldChange_D4_D10),median_D4_D14=median(log2FoldChange_D4_D14),median_D4_D21=median(log2FoldChange_D4_D21), median_continuous=median(log2FoldChange_continuous))
  median_scaled<- summary %>% mutate(median_D4_D7=adj$median_D4_D7) %>% mutate(median_D4_D10=adj$median_D4_D10)%>% mutate(median_D4_D14=adj$median_D4_D14)%>% mutate(median_D4_D21=adj$median_D4_D21) %>% mutate(median_continuous=adj$median_continuous) %>% mutate(adj_lfc_D4_D7=log2FoldChange_D4_D7-median_D4_D7) %>% mutate(adj_lfc_D4_D10=log2FoldChange_D4_D10-median_D4_D10)%>% mutate(adj_lfc_D4_D14=log2FoldChange_D4_D14-median_D4_D14) %>% mutate(adj_lfc_D4_D21=log2FoldChange_D4_D21-median_D4_D21) %>% mutate(adj_lfc_continuous=log2FoldChange_continuous-median_continuous) %>% mutate(adj_score_D4_D7=adj_lfc_D4_D7/lfcSE_D4_D7) %>% mutate(adj_score_D4_D10=adj_lfc_D4_D10/lfcSE_D4_D10)%>% mutate(adj_score_D4_D14=adj_lfc_D4_D14/lfcSE_D4_D14) %>% mutate(adj_score_D4_D21=adj_lfc_D4_D21/lfcSE_D4_D21) %>% mutate(adj_score_continuous=adj_lfc_continuous/lfcSE_continuous)                              
  median_scaled<-median_scaled %>% select("mseq","median_D4_D7","median_D4_D10","median_D4_D14", "median_D4_D21", "median_continuous","adj_lfc_D4_D7", "adj_lfc_D4_D10", "adj_lfc_D4_D14", "adj_lfc_D4_D21", "adj_lfc_continuous", "adj_score_D4_D7", "adj_score_D4_D10", "adj_score_D4_D14", "adj_score_D4_D21", "adj_score_continuous")
  
  
  #recalculate p values based on z score produced from median scaling, standard error does not need to be recalculated after the scaling as this is a simple translation (ie moving the y axis up and down)
  median_scaled<-median_scaled %>% mutate(uncombined_two_tailed_p_D4_D7= pnorm(abs(adj_score_D4_D7),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D7 = p.adjust(uncombined_two_tailed_p_D4_D7, method = "BH")) %>% 
    mutate(uncombined_two_tailed_p_D4_D10= pnorm(abs(adj_score_D4_D10),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D10 = p.adjust(uncombined_two_tailed_p_D4_D10, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_D4_D14= pnorm(abs(adj_score_D4_D14),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D14 = p.adjust(uncombined_two_tailed_p_D4_D14, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_D4_D21= pnorm(abs(adj_score_D4_D21),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_D4_D21 = p.adjust(uncombined_two_tailed_p_D4_D21, method = "BH")) %>%
    mutate(uncombined_two_tailed_p_continuous= pnorm(abs(adj_score_continuous),lower.tail = FALSE) *2) %>% mutate(uncombined_BH_FDR_continuous = p.adjust(uncombined_two_tailed_p_continuous, method = "BH"))
  summary<-merge(summary, median_scaled, by="mseq", all.x=TRUE)
  #merge summary to annotation file - this is the final output file
  annotation_file <-paste0(ident,"_master_annotated")
  annotation_file <-get(annotation_file)
  OUT<-summary
  OUT$SG<-sg_assay
  OUT$IDENT<-ident
  OUT$EXON_GROUP<-exon_assay
  OUT<-merge(OUT, annotation_file[ , c("id","mseq","pam_mut_sgrna_id")],  by="mseq", all.x=TRUE)
  OUT$combined_name<-paste0(OUT$id,"_",OUT$EXON_GROUP)
  OUT<-OUT %>% relocate(SG, .after = mseq) %>% relocate(IDENT, .after = SG) %>% relocate(EXON_GROUP, .after = IDENT) %>% relocate(combined_name, .after = mseq) %>% relocate(pam_mut_sgrna_id, .after = EXON_GROUP) %>% relocate(id, .after = mseq)
  #make multiple outputs - one full like before and one input for the combination
  OUT_full_anot<-merge(OUT, annotation_file,  by="id", all.x=TRUE)
  names(OUT_full_anot)[names(OUT_full_anot) == 'pam_mut_sgrna_id.x'] <- 'pam_mut_sgrna_id'
  names(OUT_full_anot)[names(OUT_full_anot) == 'mseq.x'] <- 'mseq'
  OUT_full_anot$pam_mut_sgrna_id.y<-NULL
  OUT_full_anot$mseq.y<-NULL
  write.csv(OUT, file=(paste0("./",ident,"_OUT.csv")), row.names = FALSE)
  write.csv(OUT_full_anot, file=(paste0("./",ident,"_OUT_full_anot.csv")), row.names = FALSE)
  #name and output files to global environment
  #deseqname<-deparse(substitute(x_filtered_counts))
  deseqname<-paste0(ident,"_OUT")
  deseqname_full<-paste0(ident,"_OUT_full_anot")
  assign(deseqname, OUT, envir = .GlobalEnv)
  assign(deseqname_full, OUT_full_anot, envir = .GlobalEnv)
}

#########################################################################
#FUNCTIONS_DESEQ_END
#########################################################################

setwd("/Users/aw28/Documents/POT1_analysis/february_2024/final_run_april_2024")

#RUN ALL OF THE ANALYSIS AND QC OUTPUT !!!!!!!!!!!!!!!!! ######START############################
#Function to apply to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(xx) {run.deseq(exon_assay=xx,sg_assay="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(xx) {run.deseq(exon_assay=xx,sg_assay="B")
})

###############----------------------------------------------------###########################################################
###############-----------DESEQ2_END-------------------------------###########################################################
###############----------------------------------------------------###########################################################


#########################################################################
#RUNS THAT REQUIRE REPLICATE EXLCUSION- START 
#########################################################################


setwd("/Users/aw28/Documents/POT1_analysis/february_2024/final_run_april_2024/repeats_exclusions")


#E9_1_SGB_R3
#E9_1_SGB_R3
#E9_1_SGB_R3

#read in experimental designs SGA
annot_a<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_9_1.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_continuous_9_1.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)


#MAKE NEW NORMALIZATION MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make normalization matricies - EXCLUDED
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D7R1", "D7R2", "D10R1", "D10R2",  "D14R1", "D14R2",  "D21R1", "D21R2")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("9_1"), function(y) {deseq.norm.input(exon=y, sg="B")
})

#MAKE NEW COUNT MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1", "D4R2", "D7R1", "D7R2", "D10R1", "D10R2",  "D14R1", "D14R2",  "D21R1", "D21R2")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("9_1"), function(y) {deseq.count.input(exon=y, sg="B")
})

#RUN DESEQ AGAIN WITH THE EXLCUDED COUNT AND NORMALIZATION MATRICIES AND NEW EXPERIMENTAL DESIGN
lapply (c("9_1"), function(xx) {run.deseq(exon_assay=xx,sg_assay="B")
})



#E14_SGB_R1
#E14_SGB_R1
#E14_SGB_R1



#read in experimental designs SGA
annot_a<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_14.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_continuous_14.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)


#MAKE NEW NORMALIZATION MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make normalization matricies - EXCLUDED
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R2", "D4R3", "D7R2", "D7R3", "D10R2", "D10R3",  "D14R2", "D14R3",  "D21R2", "D21R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("14"), function(y) {deseq.norm.input(exon=y, sg="B")
})

#MAKE NEW COUNT MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R2", "D4R3", "D7R2", "D7R3", "D10R2", "D10R3",  "D14R2", "D14R3",  "D21R2", "D21R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("14"), function(y) {deseq.count.input(exon=y, sg="B")
})

#RUN DESEQ AGAIN WITH THE EXLCUDED COUNT AND NORMALIZATION MATRICIES AND NEW EXPERIMENTAL DESIGN
lapply (c("14"), function(xx) {run.deseq(exon_assay=xx,sg_assay="B")
})


#E10_SGA_D7R3_REMOVED
#E10_SGA_D7R3_REMOVED
#E10_SGA_D7R3_REMOVED


#read in experimental designs SGA
annot_a<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_10.csv")
annot_a<-as.matrix.data.frame(annot_a)
annot_a_continuous<-read.csv("/Users/aw28/Documents/POT1_analysis/POT1_november_23/annot_a_continuous_10.csv")
annot_a_continuous$condition<-as.numeric(annot_a_continuous$condition)


#MAKE NEW NORMALIZATION MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make normalization matricies - EXCLUDED
deseq.norm.input <- function(exon, sg){
  x_normalization_counts <- paste0("E",exon,"_SG",sg,"_normalization_counts")
  x_normalization_counts <-get(x_normalization_counts)
  deseq_input_noramlization_counts <- x_normalization_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3")
  deseq_input_noramlization_counts <- deseq_input_noramlization_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_noramlization_counts <- as.matrix.data.frame(deseq_input_noramlization_counts)
  normdfname<-paste0("E",exon,"_SG",sg,"_normalization_counts_MATRIX")
  assign(normdfname, deseq_input_noramlization_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_A
lapply (c("10"), function(y) {deseq.norm.input(exon=y, sg="A")
})

#MAKE NEW COUNT MATRIX THAT HAS UNDESIRED REPLICATES EXLUDED
#make count matricies
deseq.count.input <- function(exon,sg){
  x_filtered_counts <- paste0("E",exon,"_SG",sg,"_master_annotated")
  x_filtered_counts <- get(x_filtered_counts)
  #line below removes duplicates within the targeton file based on mseq - some mutators produce the same mseq but oligo library contains one instance
  x_filtered_counts <-x_filtered_counts[!duplicated(x_filtered_counts$mseq), ]
  deseq_input_filtered_counts <- x_filtered_counts %>% dplyr::select("mseq", "D4R1","D4R2","D4R3","D7R1","D7R2","D10R1","D10R2","D10R3","D14R1","D14R2","D14R3","D21R1","D21R2","D21R3")
  deseq_input_filtered_counts <- deseq_input_filtered_counts %>% remove_rownames %>% column_to_rownames(var="mseq")
  deseq_input_filtered_counts <- as.matrix.data.frame(deseq_input_filtered_counts)
  countdfname<-paste0("E",exon,"_SG",sg,"_master_annotated_filtered_counts_MATRIX")
  assign(countdfname, deseq_input_filtered_counts, envir = .GlobalEnv)
}
#Function to apply the above function to listed exons and designated sgRNA_B
lapply (c("10"), function(y) {deseq.count.input(exon=y, sg="A")
})

#RUN DESEQ AGAIN WITH THE EXLCUDED COUNT AND NORMALIZATION MATRICIES AND NEW EXPERIMENTAL DESIGN
lapply (c("10"), function(xx) {run.deseq(exon_assay=xx,sg_assay="A")
})


#########################################################################
#RUN REPEATS THAT REQUIRE REPLICATE EXCLUSION - END 
#########################################################################



#########################################################################
#EXON 7 - TESTING
#########################################################################







###############----------------------------------------------------###########################################################
###############------QC plots for editing-START--------------------###########################################################
###############----------------------------------------------------###########################################################

#positional effect QC generation 
position.plas.plot<-function(exon,sg){
  locus_input <- paste0("E",exon,"_SG",sg,"_OUT_full_anot")
  locus_input <- get(locus_input)
  locus_input_plasmid<-locus_input
  locus_input_plasmid$D4_MEAN <-(locus_input_plasmid$D4R1 + locus_input_plasmid$D4R2 + locus_input_plasmid$D4R3)/3
  locus_input_plasmid$D4_DIV_PLAS<-locus_input_plasmid$D4_MEAN/locus_input_plasmid$PLASMID
  
  locus_input_plasmid$D4R1_DIV_PLAS<-locus_input_plasmid$D4R1/locus_input_plasmid$PLASMID
  locus_input_plasmid$D4R2_DIV_PLAS<-locus_input_plasmid$D4R2/locus_input_plasmid$PLASMID
  locus_input_plasmid$D4R3_DIV_PLAS<-locus_input_plasmid$D4R3/locus_input_plasmid$PLASMID
  
  pam_df<-locus_input_plasmid %>% filter(!is.na(pam_mut_sgrna_id))
  pam_df<-pam_df %>% filter(!pam_mut_sgrna_id %in% "")
  locus_input_plasmid$mut_type<-as.factor(locus_input_plasmid$mut_type)
  levels(locus_input_plasmid$mut_type)[levels(locus_input_plasmid$mut_type)==''] <- 'other'
  
  pos_plot_mean<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate Mean")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("Average Day_4 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_mean)
  
  pos_plot_D4R1<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R1_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R1_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 1")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R1 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R1)
  
  pos_plot_D4R2<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R2_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R2_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 2")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R2 SGE Counts / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R2)
  
  pos_plot_D4R3<-ggplot(locus_input_plasmid, aes(x=mut_position, y=D4R3_DIV_PLAS, color = factor(mut_type, levels=c("non","syn","mis","other")))) +
    scale_y_continuous(trans='log2')+
    #ylim(-10, 50)+
    theme_classic()+
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          legend.title=element_blank(),
          axis.ticks.x = element_blank())+
    geom_point(aes(), size = 1, alpha=0.5) +
    #geom_point(aes(), shape = 1,size = 2,colour = "black",stroke = .25)+
    scale_color_brewer(palette = "Set1")+
    geom_point(data=pam_df, 
               aes(x=mut_position,y=D4R3_DIV_PLAS), 
               colour='black',
               size=1, 
               shape=3, 
               alpha=0.5)+
    ggtitle("Replicate 3")+
    xlab("GRCh38 Genomic Coordinate") +
    ylab("D4R3 / Plasmid Counts")+
    geom_smooth(data = locus_input_plasmid, method = "loess", se = FALSE, span=0.5, color='blue')
  print(pos_plot_D4R3)
  
  p<-plot_grid(pos_plot_mean, pos_plot_D4R1, pos_plot_D4R2, pos_plot_D4R3, nrow = 1, ncol = 4, label_size = 12)
  
  pdf(paste0("E",exon,"_SG",sg,"_positional_effect_D4_DIV_PLASMID.pdf"), width=16, height=4, pointsize = 1)
  
  print(p)
  
  dev.off() 
}

#Function to apply to listed exons and designated sgRNA_A
lapply (c("2","3_2","4","6","8","10","11","12"), function(xxx) {position.plas.plot(exon=xxx,sg="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(xxx) {position.plas.plot(exon=xxx,sg="B")
})

###############----------------------------------------------------###########################################################
###############------QC plots for editing-END----------------------###########################################################
###############----------------------------------------------------###########################################################





###############----------------------------------------------------###########################################################
###############-----------library tiling calculations-START -------###########################################################
###############----------------------------------------------------###########################################################

#combine and count 

#the above DESeq2 scripts use unique mseq per targeton - mseqs will be duplicated if not made unique as. This is because different 'id' are not unique for sequence 'mseq'.
#mutator function used by valiant is part of the 'id' nomenclature, so the same variant at the 'mseq' level is produced by varied means and will produce a duplicated mseq at the annotation stage
#the mseqs are expanded again at the final annotation stage in the deseq process and all ids will again be present in the below bindings

#tiled libraries
total_single_pot1 <- do.call("rbind.fill", list(E1_SGB_OUT,
                                                E2_SGA_OUT,
                                                E4_SGA_OUT,
                                                E5_SGB_OUT,
                                                E6_SGB_OUT,
                                                E8_SGA_OUT,
                                                E10_SGA_OUT,
                                                E11_SGA_OUT,
                                                E12_SGA_OUT,
                                                E13_SGB_OUT,
                                                E14_SGB_OUT))
total_single_pot1$process<-"non_tiled"
#single libraries
total_tiled_pot1 <- do.call("rbind.fill", list(E3_1_SGA_OUT, 
                                               E3_2_SGA_OUT,
                                               E9_1_SGB_OUT,
                                               E9_2_SGB_OUT))

total_tiled_pot1$process<-"tiled"                                                    
#can bind together libraries - the 'id' field produced by valiant in these dataframes are guide library agnostic (ie. the duplicate variant will will have the same id)
total_pot1 <- do.call("rbind.fill", list(total_single_pot1,total_tiled_pot1))

#this produces a data frame with all variants processed - not duplicated for 'id' but will be duplicated for mseq
total_unique_id_count<-total_pot1[!duplicated(total_pot1$id), ]
#this produces a data frame with all variants processed - not duplicated for 'mseq' - what the experimental flasks would have contained as collated unique libraries
total_unique_id_count<-total_pot1[!duplicated(total_pot1$mseq), ]



total_single_pot1_pam_full<-total_pot1 %>%
  group_by(mseq) %>%
  dplyr::summarize(dup_mseq= n(), 
                   PAM_status_full_arrayed = paste(pam_mut_sgrna_id, collapse = ','))%>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_b")~ "B", TRUE~ "-")) %>%
  mutate(PAM_status_simplified=case_when(str_detect(PAM_status_full_arrayed,"_a")~ "A", TRUE~ PAM_status_simplified)) %>%
  select(c(1,2,4))

total_single_pot1_pam_corrected<-merge(total_pot1,total_single_pot1_pam_full,by="mseq",all.x=TRUE)

total_single_pot1_pam_corrected<-total_single_pot1_pam_corrected  %>% arrange(EXON_GROUP,id,SG)




#select the necessary columns from the uncombined data
uncombined_sub <- total_single_pot1_pam_corrected
#add a column and populate with description of this data as being from an uncombined source, that is one guide was used to get this data
uncombined_sub$Variant_duplication<-NA
uncombined_sub$dataset_process<-"not_combined"
#change the names of the following columns, so that they can be bound to the combined dataframe
uncombined_sub<-uncombined_sub %>% dplyr::rename(Variant_Sources = "SG")
uncombined_sub<-uncombined_sub %>% dplyr::rename(PAM_status = "pam_mut_sgrna_id")
uncombined_sub<-uncombined_sub %>% dplyr::rename(PAM_status_simplified_collated = "PAM_status_simplified")
uncombined_sub<-uncombined_sub %>% dplyr::rename(dup_mseq_collated = "dup_mseq")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D7 = "uncombined_two_tailed_p_D4_D7")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D10 = "uncombined_two_tailed_p_D4_D10")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D14 = "uncombined_two_tailed_p_D4_D14")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_D4_D21 = "uncombined_two_tailed_p_D4_D21")
uncombined_sub<-uncombined_sub %>% dplyr::rename(two_tailed_p_continuous = "uncombined_two_tailed_p_continuous")

#turn mseq in uncombined sub into mutseq_a and mutseq_b in order to bind tables
uncombined_sub <- uncombined_sub %>% mutate(mseq_b=case_when(Variant_Sources == "B" ~ mseq, TRUE ~ "NA")) %>%
  mutate(mseq_a=case_when(Variant_Sources == "A" ~ mseq, TRUE ~ "NA"))
#remove the mseq field from the combined as it is now redundant
uncombined_sub$mseq<-NULL  



targeton_input<-uncombined_sub %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"1")~"2",TRUE~"-")) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"2")~"3",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"3")~"4",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"4")~"5",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"5")~"6",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"6")~"7",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"8")~"9",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"9")~"10",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"10")~"11",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"11")~"12",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"12")~"13",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"13")~"14",TRUE~REGION)) %>%
  mutate(REGION=case_when(str_detect(EXON_GROUP,"14")~"15",TRUE~REGION)) 

targeton_collapse <- targeton_input %>%
  group_by(id,REGION) %>%
  dplyr::summarize(variant_n_instances_between_targetons= n(), 
                   #exon_specific_names = paste(combined_name, collapse = '|'),
                   source_dataset = paste(dataset_process, collapse ="|"),
                   source_targetons = paste(EXON_GROUP, collapse = '|'),
                   source_libraries = paste(Variant_Sources, collapse = '|'),
                   variant_n_instances_between_libraries = paste(Variant_duplication, collapse = '|'),
                   PAM_status_simplified_collated = paste(PAM_status_simplified_collated, collapse = '|'), 
                   PAM_status = paste(PAM_status, collapse='|'),
                   dup_mseq_collated = paste(dup_mseq_collated, collapse='|'),
                   arrayed_mseq_a = paste(mseq_a, collapse = '|'), 
                   arrayed_mseq_b = paste(mseq_b, collapse = '|'), 
                   
                   lfcSE_D4_D7 = paste(lfcSE_D4_D7, collapse = ','), 
                   lfcSE_D4_D10 = paste(lfcSE_D4_D10, collapse = ','), 
                   lfcSE_D4_D14 = paste(lfcSE_D4_D14, collapse = ','), 
                   lfcSE_D4_D21 = paste(lfcSE_D4_D21, collapse = ','), 
                   lfcSE_continuous = paste(lfcSE_continuous, collapse = ','), 
                   
                   adj_lfc_D4_D7 = paste(adj_lfc_D4_D7, collapse = ','),
                   adj_lfc_D4_D10 = paste(adj_lfc_D4_D10, collapse = ','),
                   adj_lfc_D4_D14 = paste(adj_lfc_D4_D14, collapse = ','),
                   adj_lfc_D4_D21 = paste(adj_lfc_D4_D21, collapse = ','),
                   adj_lfc_continuous = paste(adj_lfc_continuous, collapse = ',')) %>%
  
  ungroup() %>%
  #ungroup the key data into primary and secondary targetons fields
  separate(adj_lfc_D4_D7, sep=",", c("primary_tg_adj_lfc_D4_D7","secondary_tg_adj_lfc_D4_D7")) %>%
  separate(adj_lfc_D4_D10, sep=",", c("primary_tg_adj_lfc_D4_D10","secondary_tg_adj_lfc_D4_D10")) %>%
  separate(adj_lfc_D4_D14, sep=",", c("primary_tg_adj_lfc_D4_D14","secondary_tg_adj_lfc_D4_D14")) %>%
  separate(adj_lfc_D4_D21, sep=",", c("primary_tg_adj_lfc_D4_D21","secondary_tg_adj_lfc_D4_D21")) %>%
  separate(adj_lfc_continuous, sep=",", c("primary_tg_adj_lfc_continuous","secondary_tg_adj_lfc_continuous")) %>%
  
  separate(lfcSE_D4_D7, sep=",", c("primary_tg_lfcSE_D4_D7","secondary_tg_lfcSE_D4_D7")) %>%
  separate(lfcSE_D4_D10, sep=",", c("primary_tg_lfcSE_D4_D10","secondary_tg_lfcSE_D4_D10"))%>%
  separate(lfcSE_D4_D14, sep=",", c("primary_tg_lfcSE_D4_D14","secondary_tg_lfcSE_D4_D14"))%>%
  separate(lfcSE_D4_D21, sep=",", c("primary_tg_lfcSE_D4_D21","secondary_tg_lfcSE_D4_D21")) %>%
  separate(lfcSE_continuous, sep=",", c("primary_tg_lfcSE_continuous","secondary_tg_lfcSE_continuous"))

#make the newly created fields numeric in order to do weighted calculations
targeton_collapse[, c(13:32)] <- sapply(targeton_collapse[, c(13:32)], as.numeric)
#make the number of targeton source field character based to get subsequent string detection to work
targeton_collapse$variant_n_instances_between_targetons<-as.character(targeton_collapse$variant_n_instances_between_targetons)

#PERFORM THE WEIGHTED CALCULATIONS
targeton_processed <-targeton_collapse %>%
  mutate(weight_tg1_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D7)^2)) %>%
  mutate(weight_tg1_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D10)^2)) %>%
  mutate(weight_tg1_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D14)^2)) %>%
  mutate(weight_tg1_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_D4_D21)^2)) %>%
  mutate(weight_tg1_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(primary_tg_lfcSE_continuous)^2)) %>%
  
  mutate(weight_tg2_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D7)^2)) %>%
  mutate(weight_tg2_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D10)^2)) %>%
  mutate(weight_tg2_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D14)^2)) %>%
  mutate(weight_tg2_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_D4_D21)^2)) %>%
  mutate(weight_tg2_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~as.numeric(NA), TRUE~ 1/(secondary_tg_lfcSE_continuous)^2)) %>%
  
  mutate(sum_of_weight_D4_D7=rowSums(cbind(weight_tg1_D4_D7,weight_tg2_D4_D7), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D10=rowSums(cbind(weight_tg1_D4_D10,weight_tg2_D4_D10), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D14=rowSums(cbind(weight_tg1_D4_D14,weight_tg2_D4_D14), na.rm=FALSE))%>%
  mutate(sum_of_weight_D4_D21=rowSums(cbind(weight_tg1_D4_D21,weight_tg2_D4_D21), na.rm=FALSE))%>%
  mutate(sum_of_weight_continuous=rowSums(cbind(weight_tg1_continuous,weight_tg2_continuous), na.rm=FALSE))%>%
  
  mutate(SE_bind_D4_D7 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D7, TRUE~ (sum_of_weight_D4_D7)^(-0.5))) %>%
  mutate(SE_bind_D4_D10 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D10, TRUE~ (sum_of_weight_D4_D10)^(-0.5))) %>%
  mutate(SE_bind_D4_D14 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D14, TRUE~ (sum_of_weight_D4_D14)^(-0.5))) %>%
  mutate(SE_bind_D4_D21 = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_D4_D21, TRUE~ (sum_of_weight_D4_D21)^(-0.5))) %>%
  mutate(SE_bind_continuous = case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_lfcSE_continuous, TRUE~ (sum_of_weight_continuous)^(-0.5))) %>%
  
  mutate(weighted_primary_tg_LFC_D4_D7= weight_tg1_D4_D7*primary_tg_adj_lfc_D4_D7) %>% mutate(weighted_secondary_tg_LFC_D4_D7= weight_tg2_D4_D7*secondary_tg_adj_lfc_D4_D7) %>%
  mutate(weighted_primary_tg_LFC_D4_D10= weight_tg1_D4_D10*primary_tg_adj_lfc_D4_D10) %>% mutate(weighted_secondary_tg_LFC_D4_D10= weight_tg2_D4_D10*secondary_tg_adj_lfc_D4_D10) %>%
  mutate(weighted_primary_tg_LFC_D4_D14= weight_tg1_D4_D14*primary_tg_adj_lfc_D4_D14) %>% mutate(weighted_secondary_tg_LFC_D4_D14= weight_tg2_D4_D14*secondary_tg_adj_lfc_D4_D14) %>%
  mutate(weighted_primary_tg_LFC_D4_D21= weight_tg1_D4_D21*primary_tg_adj_lfc_D4_D21) %>% mutate(weighted_secondary_tg_LFC_D4_D21= weight_tg2_D4_D21*secondary_tg_adj_lfc_D4_D21) %>%
  mutate(weighted_primary_tg_LFC_continuous= weight_tg1_continuous*primary_tg_adj_lfc_continuous) %>% mutate(weighted_secondary_tg_LFC_continuous= weight_tg2_continuous*secondary_tg_adj_lfc_continuous) %>%
  
  mutate(sum_of_weighted_LFC_D4_D7=rowSums(cbind(weighted_primary_tg_LFC_D4_D7,weighted_secondary_tg_LFC_D4_D7), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D7=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D7, TRUE~ sum_of_weighted_LFC_D4_D7/sum_of_weight_D4_D7)) %>%
  mutate(sum_of_weighted_LFC_D4_D10=rowSums(cbind(weighted_primary_tg_LFC_D4_D10,weighted_secondary_tg_LFC_D4_D10), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D10=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D10, TRUE~ sum_of_weighted_LFC_D4_D10/sum_of_weight_D4_D10)) %>%
  mutate(sum_of_weighted_LFC_D4_D14=rowSums(cbind(weighted_primary_tg_LFC_D4_D14,weighted_secondary_tg_LFC_D4_D14), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D14=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D14, TRUE~ sum_of_weighted_LFC_D4_D14/sum_of_weight_D4_D14)) %>%
  mutate(sum_of_weighted_LFC_D4_D21=rowSums(cbind(weighted_primary_tg_LFC_D4_D21,weighted_secondary_tg_LFC_D4_D21), na.rm=FALSE)) %>%
  mutate(processed_LFC_D4_D21=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_D4_D21, TRUE~ sum_of_weighted_LFC_D4_D21/sum_of_weight_D4_D21)) %>%
  mutate(sum_of_weighted_LFC_continuous=rowSums(cbind(weighted_primary_tg_LFC_continuous,weighted_secondary_tg_LFC_continuous), na.rm=FALSE)) %>%
  mutate(processed_LFC_continuous=case_when(str_detect(variant_n_instances_between_targetons,"1")~primary_tg_adj_lfc_continuous, TRUE~ sum_of_weighted_LFC_continuous/sum_of_weight_continuous)) #%>%

#check that no duplicate ids remain
duplicated <- targeton_processed  %>% 
  group_by(id) %>%
  filter(n() > 1)

#now make a new mseq identifier for the appropriate mseq if either mseq is ok then can use a or b annotation

#we now need to derive the correct annotation source for each variant
#we can use the mseq as an identifier for the variant that can be made unique so that BH FDR calculation is correct for the number of unique observations
#we also need to annotate correctly, taking into account whether the variant scores were derived from library a or library b in the case of the PPE codons

#first = is the variant derived from library a or library b alone? if yes then the correct annotation source for each variant is simply the library used:
targeton_processed_req_anno<-targeton_processed %>% mutate(required_annotation_source=case_when(source_dataset=="not_combined"~source_libraries,TRUE~"-")) %>%
  relocate(required_annotation_source,.after=id) %>%
  #second = if the variant was observed in library a and library b did the variant fall in a PPE codon? if yes then the variant calculations were based on the alternative library alone, so get that annotation
  mutate(required_annotation_source=case_when(source_dataset=="combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #third = is the variant observed in both libraries and also in an overlapping region of a tiled targeton and a PPE codon? if yes the annotation will again be from the alternative libary to the PPE containing library 
  mutate(required_annotation_source=case_when(source_dataset=="combined|combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined|combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #fourth = is the variant observed in one targeton where library a and b werd combined, and in one targeton where the only one library was used, and does the variant fall in a PPE? if so then take the non-PPE annotaiton source  
  mutate(required_annotation_source=case_when(source_dataset=="combined|not_combined" & str_detect(PAM_status_simplified_collated,"B")~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(source_dataset=="combined|not_combined" & str_detect(PAM_status_simplified_collated,"A")~"B",TRUE~required_annotation_source)) %>%
  #fifth = is the variant not so far allocated an annotation source? therefore it is not a PPE codon. Is the variant only observed in library a or library b, if so then that library is the correct source  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="A"~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="B"~"B",TRUE~required_annotation_source)) %>%
  #sixth = is the variant observed in only one library a or b? and is the variant in an overlapping region between targetons? if yes, then again the annotaiton source is that library (a or b)  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="A|A"~"A",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(required_annotation_source=="-" & source_libraries=="B|B"~"B",TRUE~required_annotation_source)) %>%
  #finally = if the variant has so far not ben allocated an annotation source, then either library a or library b mseq can be used as an annotation as they have identical VEP and VaLiAnT annotations  
  mutate(required_annotation_source=case_when(required_annotation_source=="-" ~ "either",TRUE~required_annotation_source))



#updated approach to annotation for POT1

targeton_processed_req_anno<-targeton_processed %>% 
  mutate(required_annotation_source="-") %>%
  mutate(required_annotation_source=case_when(str_detect(source_libraries,"B")~"B",TRUE~required_annotation_source)) %>%
  mutate(required_annotation_source=case_when(str_detect(source_libraries,"A")~"A",TRUE~required_annotation_source)) %>%
  relocate(required_annotation_source,.after=id)



#seperate out the arrayed mseq field (that were combined into one string at the weighted overlapping tile step, and separated by "|" with targeton 1 on the left, targeton 2 on the right)
targeton_processed_req_anno_2<-targeton_processed_req_anno %>% separate(arrayed_mseq_a, sep='\\|',c("mseq_a_tg1","mseq_a_tg2")) %>%
  separate(arrayed_mseq_b, sep='\\|',c("mseq_b_tg1","mseq_b_tg2")) 

#now create a field with the appropriate mseq for the variant - most variants are not in overlapping region so take targeton 1 (either a or b) for the overlapping take also take targeton 1 as the unqiue identifer
targeton_processed_req_anno_mseq<-targeton_processed_req_anno_2 %>% mutate(mseq_for_filter=case_when(required_annotation_source=="A" ~ mseq_a_tg1,TRUE~"-")) %>%
  mutate(mseq_for_filter=case_when(required_annotation_source=="B" ~ mseq_b_tg1,TRUE~mseq_for_filter)) %>%
  mutate(mseq_for_filter=case_when(required_annotation_source=="either" ~ mseq_a_tg1,TRUE~mseq_for_filter)) %>%
  relocate(mseq_for_filter,.after=id)

#make unique for the mseq
targeton_processed_req_anno_mseq_unique<-targeton_processed_req_anno_mseq[!duplicated(targeton_processed_req_anno_mseq$mseq_for_filter), ]

#check that no duplicate mseqs remain
duplicated_mseq <- targeton_processed_req_anno_mseq_unique  %>% 
  group_by(mseq_for_filter) %>%
  filter(n() > 1)

#produce the final statistical model now that the variant rows are truly unique 
final_stats_frame<-targeton_processed_req_anno_mseq_unique %>%
  #PERFORM THE CALCULATION OF Z SCORES > P VALUES > BH_FDR - this is the final statistics that will be used 
  mutate(processed_Z_D4_D7=processed_LFC_D4_D7/SE_bind_D4_D7) %>% mutate(two_tailed_p_D4_D7= pnorm(abs(processed_Z_D4_D7),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D7 = p.adjust(two_tailed_p_D4_D7, method = "BH")) %>%
  mutate(processed_Z_D4_D10=processed_LFC_D4_D10/SE_bind_D4_D10) %>% mutate(two_tailed_p_D4_D10= pnorm(abs(processed_Z_D4_D10),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D10 = p.adjust(two_tailed_p_D4_D10, method = "BH")) %>%
  mutate(processed_Z_D4_D14=processed_LFC_D4_D14/SE_bind_D4_D14) %>% mutate(two_tailed_p_D4_D14= pnorm(abs(processed_Z_D4_D14),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D14 = p.adjust(two_tailed_p_D4_D14, method = "BH")) %>%
  mutate(processed_Z_D4_D21=processed_LFC_D4_D21/SE_bind_D4_D21) %>% mutate(two_tailed_p_D4_D21= pnorm(abs(processed_Z_D4_D21),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_D4_D21 = p.adjust(two_tailed_p_D4_D21, method = "BH")) %>%
  mutate(processed_Z_continuous=processed_LFC_continuous/SE_bind_continuous) %>% mutate(two_tailed_p_continuous= pnorm(abs(processed_Z_continuous),lower.tail = FALSE) *2) %>% mutate(processed_BH_FDR_continuous = p.adjust(two_tailed_p_continuous, method = "BH"))

#metrics calculation
final_stats_frame$functional_classification<-NA
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous<0.01 & final_stats_frame$processed_LFC_continuous<0] <- 'depleted'
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous<0.01 & final_stats_frame$processed_LFC_continuous>0] <- 'enriched'
final_stats_frame$functional_classification[final_stats_frame$processed_BH_FDR_continuous>=0.01] <- 'unchanged'

#check the categories
table(final_stats_frame$functional_classification)
#check there are no NA rows in new classification column
is_na_class<-final_stats_frame %>% filter(functional_classification %in% NA)


###############----------------------------------------------------###########################################################
###############------QC plots for editing-START--------------------###########################################################
###############----------------------------------------------------###########################################################



###############----------------------------------------------------###########################################################
###############---Expansion of data frame to full annotation-START-###########################################################
###############----------------------------------------------------###########################################################



#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGA 
E2_SGA_master_annotated_META_VEP<-merge(E2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124870880_124871127_minus_sgRNA_w2_a_pam.tsv, by="id", all.x=TRUE)
E3_1_SGA_master_annotated_META_VEP<-merge(E3_1_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863272_124863520_minus_sgRNA_w3.1_a_pam.tsv , by="id", all.x=TRUE)
E3_2_SGA_master_annotated_META_VEP<-merge(E3_2_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863468_124863717_minus_sgRNA_w3.2_a_pam.tsv, by="id", all.x=TRUE)
E4_SGA_master_annotated_META_VEP<-merge(E4_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124858890_124859138_minus_sgRNA_w4_a_pam.tsv, by="id", all.x=TRUE)
E6_SGA_master_annotated_META_VEP<-merge(E6_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124851781_124852030_minus_sgRNA_w6_a_pam.tsv, by="id", all.x=TRUE)
E8_SGA_master_annotated_META_VEP<-merge(E8_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124842757_124843006_minus_sgRNA_w8_a_pam.tsv, by="id", all.x=TRUE)
E10_SGA_master_annotated_META_VEP<-merge(E10_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124835214_124835463_minus_sgRNA_w10_a_pam.tsv, by="id", all.x=TRUE)
E11_SGA_master_annotated_META_VEP<-merge(E11_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124829216_124829465_minus_sgRNA_w11_a_pam.tsv, by="id", all.x=TRUE)
E12_SGA_master_annotated_META_VEP<-merge(E12_SGA_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124827112_124827359_minus_sgRNA_w12_a_pam.tsv, by="id", all.x=TRUE)
#merge filtered master tables with VEP annotation tables and syn_filter
#merge filtered counts + meta with vep output to give filtered annotated dataframes SGB
E1_SGB_master_annotated_META_VEP<-merge(E1_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124892164_124892411_minus_sgRNA_w1_b_pam.tsv, by="id", all.x=TRUE)
E3_2_SGB_master_annotated_META_VEP<-merge(E3_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124863468_124863717_minus_sgRNA_w3.2_b_pam.tsv, by="id", all.x=TRUE)
E4_SGB_master_annotated_META_VEP<-merge(E4_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124858890_124859138_minus_sgRNA_w4_b_pam.tsv, by="id", all.x=TRUE)
E5_SGB_master_annotated_META_VEP<-merge(E5_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124852932_124853177_minus_sgRNA_w5_b_pam.tsv, by="id", all.x=TRUE)
E6_SGB_master_annotated_META_VEP<-merge(E6_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124851781_124852030_minus_sgRNA_w6_b_pam.tsv, by="id", all.x=TRUE)
E9_1_SGB_master_annotated_META_VEP<-merge(E9_1_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124840935_124841182_minus_sgRNA_w9.1_b_pam.tsv, by="id", all.x=TRUE)
E9_2_SGB_master_annotated_META_VEP<-merge(E9_2_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124840984_124841232_minus_sgRNA_w9.2_b_pam.tsv, by="id", all.x=TRUE)
E13_SGB_master_annotated_META_VEP<-merge(E13_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124825226_124825474_minus_sgRNA_w13_b_pam.tsv, by="id", all.x=TRUE)
E14_SGB_master_annotated_META_VEP<-merge(E14_SGB_master_filtered_counts, VEP_OUTPUT_VEP_OUTPUT_chr7_124823919_124824168_minus_sgRNA_w14_b_pam.tsv, by="id", all.x=TRUE)



#FUNCTION to add an additional field which is ID and SG combined - to label a variant annotation as either from SGA or SGB source - which is needed subsequently 
#label annotations with unique id
label <- function(exon, sg){
  input<-paste0("E",exon,"_SG",sg,"_master_annotated_META_VEP")
  input<-get(input)
  input$SG<-sg
  input$exon_group<-exon
  input$required_annotation_id<-paste0(input$id,"_",input$SG)
  countdfname<-paste0("E",exon,"_SG",sg,"_META_VEP")
  assign(countdfname, input, envir = .GlobalEnv)
}

#Function to apply to listed exons and designated sgRNA_A
lapply (c("2","3_1","3_2","4","6","8","10","11","12"), function(xx) {label(exon=xx,sg="A")
})
#Function to apply to listed exons and designated sgRNA_B
lapply (c("1","3_2","4","5","6","9_1","9_2","13","14"), function(xx) {label(exon=xx,sg="B")
})


total_META_VEP <- do.call("rbind.fill", list(E2_SGA_META_VEP,
                                             E3_1_SGA_META_VEP,
                                             E3_2_SGA_META_VEP,
                                             E4_SGA_META_VEP,
                                             # E6_SGA_META_VEP,
                                             E8_SGA_META_VEP,
                                             E10_SGA_META_VEP,
                                             E11_SGA_META_VEP,
                                             E12_SGA_META_VEP,
                                             
                                             E1_SGB_META_VEP,
                                             #E3_2_SGB_META_VEP,
                                             #E4_SGB_META_VEP,
                                             E5_SGB_META_VEP,
                                             E6_SGB_META_VEP,
                                             E9_1_SGB_META_VEP,
                                             E9_2_SGB_META_VEP,
                                             E13_SGB_META_VEP,
                                             E14_SGB_META_VEP))



final_stats_frame_to_merge<-final_stats_frame %>% dplyr::rename(mseq= "mseq_for_filter")

expanded_test<-merge(final_stats_frame_to_merge, total_META_VEP, by="mseq")


duplicated_id <- expanded_test  %>% 
  group_by(id.y) %>%
  filter(n() > 1)

write.csv(duplicated_id,"./duplicated_id.csv", row.names = FALSE)


#is_na_expanded_test<-expanded_test %>% filter(is.na(HGVSc))
#write.csv(is_na_expanded_test,"./is_na_expanded_test.csv", row.names = FALSE)



#There are some mseqs associated with an ID that we don't want because they are an meseq that is a duplicate of a mseq that is annotated as having a pam edit, but isn't itself annotated as a pam containing mseq
#this happens with 1del deletions, and in some cases 2del if the pam site is in a 2 del location 
#some ids are represented by multiple mseqs 
#these mseqs will have remained distinct in the annotation process
#for library a and library b the required annotation id will have pulled the correct mseq(s)
#some required annotation ids will have pulled mseqs that are duplications

#specific dupliciations to remove - these are in the 3_1, 3_2 overlap region - all are unchanged

expanded_filtered<-expanded_test %>% filter(!id.y %in% "ENST00000357628.8.ENSG00000128513.16_chr7:124863497_124863499_TAC>GAC_snvre_rc") %>%
  filter(!id.y %in% "ENST00000357628.8.ENSG00000128513.16_chr7:124863497_124863499_TAC>CAC_snvre_rc") %>% 
  filter(!id.y %in% "ENST00000357628.8.ENSG00000128513.16_chr7:124863497_124863499_TAC>AAC_snvre_rc") #%>% 

#filter(!id.y %in%"ENST00000357628.8.ENSG00000128513.16_chr7:124863497_T>G_snv_rc") %>% 
#filter(!id.y %in%"ENST00000357628.8.ENSG00000128513.16_chr7:124863497_T>C_snv_rc") %>% 
# filter(!id.y %in%"ENST00000357628.8.ENSG00000128513.16_chr7:124863497_T>A_snv_rc")






#stats_unique<-final_stats_frame_to_merge[!duplicated(final_stats_frame_to_merge$mseq), ]

#expanded_unique<-expanded_filtered[!duplicated(expanded_filtered$mseq), ]


clean_annotated_dataset<-expanded_filtered[, c(
  "mseq",
  "HGVSc",
  "HGVSp",
  "id.y",
  #"required_annotation_source",
  "REGION",
  #"variant_n_instances_between_targetons",
  #"source_dataset",
  "source_targetons",
  "source_libraries",
  #"variant_n_instances_between_libraries",
  "PAM_status_simplified_collated",
  "PAM_status",
  #"dup_mseq_collated",
  "SE_bind_D4_D7",
  "SE_bind_D4_D10",
  "SE_bind_D4_D14",
  "SE_bind_D4_D21",
  "SE_bind_continuous",
  "processed_LFC_D4_D7",
  "processed_LFC_D4_D10",
  "processed_LFC_D4_D14",
  "processed_LFC_D4_D21",
  "processed_LFC_continuous",
  "processed_Z_D4_D7",
  "two_tailed_p_D4_D7",
  "processed_BH_FDR_D4_D7",
  "processed_Z_D4_D10",
  "two_tailed_p_D4_D10",
  "processed_BH_FDR_D4_D10",
  "processed_Z_D4_D14",
  "two_tailed_p_D4_D14",
  "processed_BH_FDR_D4_D14",
  "processed_Z_D4_D21",
  "two_tailed_p_D4_D21",
  "processed_BH_FDR_D4_D21",
  "processed_Z_continuous",
  "two_tailed_p_continuous",
  "processed_BH_FDR_continuous",
  "functional_classification",
  "D4R1",
  "D4R2",
  "D4R3",
  "D7R1",
  "D7R2",
  "D7R3",
  "D10R1",
  "D10R2",
  "D10R3",
  "D14R1",
  "D14R2",
  "D14R3",
  "D21R1",
  "D21R2",
  "D21R3",
  "PLASMID",
  "gene_id",
  "transcript_id",
  "ref_chr",
  "ref_strand",
  "ref_start",
  "ref_end",
  "mut_position",
  "ref",
  "new",
  "ref_aa",
  "alt_aa",
  "mut_type",
  "mutator",
  "oligo_length",
  "pam_mut_sgrna_id",
  "mave_nt",
  "vcf_var_in_const",
  "Consequence",
  "cDNA_position",
  "CDS_position",
  "Protein_position",
  "Amino_acids",
  "Codons",
  "Existing_variation",
  "SIFT",
  "PolyPhen",
  "EXON",
  "INTRON",
  "DOMAINS",
  "HGVS_OFFSET",
  "AF",
  "gnomAD_AF",
  "gnomAD_AFR_AF",
  "gnomAD_AMR_AF",
  "gnomAD_ASJ_AF",
  "gnomAD_EAS_AF",
  "gnomAD_FIN_AF",
  "gnomAD_NFE_AF",
  "gnomAD_OTH_AF",
  "gnomAD_SAS_AF",
  "CLIN_SIG",
  "gnomAD_FLAG",
  "gnomAD_AF.1",
  "ClinVar",
  "ClinVar_CLNSIG",
  "ClinVar_CLNREVSTAT",
  "dbSNP")]


#Consequence simplification
clean_pot1_v1<- clean_annotated_dataset  %>% mutate(slim_consequence=case_when(str_detect(Consequence,"stop_gained") ~ "stop_gained")) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"UTR_variant") ~ "UTR", TRUE ~ slim_consequence)) %>%
  mutate(slim_consequence=case_when(str_detect(Consequence,"frameshift_variant") ~ "frameshift", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(mutator,"inframe") ~ "codon_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"intron_variant") ~ "intron", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"missense_variant") ~ "missense", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"synonymous_variant") ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_acceptor_variant") ~ "splice_acceptor", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"splice_donor_variant") ~ "splice_donor", TRUE ~ slim_consequence)) %>%
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"start_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"start_lost") & mutator== "snv" |
                                      str_detect(Consequence,"start_lost") & mutator== "custom" |
                                      str_detect(Consequence,"start_lost") & mutator== "ala"~ "start_lost", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_lost") & mutator== "snvre" |
                                      str_detect(Consequence,"stop_lost") & mutator== "snv" |
                                      str_detect(Consequence,"stop_lost") & mutator== "custom" |
                                      str_detect(Consequence,"stop_lost") & mutator== "ala"~ "stop_lost", TRUE ~ slim_consequence)) %>% 
  
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"stop_retained_variant") & mutator=="snvre" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="snv" |
                                      str_detect(Consequence,"stop_retained_variant") & mutator=="custom" ~ "synonymous", TRUE ~ slim_consequence)) %>% 
  
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_deletion") & mutator=="custom" ~ "clinical_inframe_deletion", TRUE ~ slim_consequence)) %>% 
  mutate(slim_consequence=case_when(str_detect(Consequence,"inframe_insertion") & mutator=="custom" ~ "clinical_inframe_insertion", TRUE ~ slim_consequence)) 
#relocate new column
clean_pot1_v1 <- clean_pot1_v1 %>% relocate(slim_consequence, .after = Consequence)

#add LFC_continuous as a new column called functional socre
clean_pot1_v1$functional_score <- clean_pot1_v1$processed_LFC_continuous
#relocate new column
clean_pot1_v1 <- clean_pot1_v1 %>% relocate(functional_score, .after = functional_classification)


#FEBRUARY 2024 BIND 

#clinvar updated binding

env_pot1_clinvar_all_120324<-pot1_clinvar_all_120324


#bind clinvar release 1 star only - 'one_star_040923' is a dataframe downloaded on 4th September 2023 which contains >=1* review status variants
env_pot1_clinvar_all_120324$CLIN_SIG<-gsub("\\s*\\([^\\)]+\\)","",as.character(env_pot1_clinvar_all_120324$Germline.classification))
env_pot1_clinvar_all_120324$HGVSc<-gsub("\\s*\\([^\\)]+\\)","",as.character(env_pot1_clinvar_all_120324$Name))
env_pot1_clinvar_all_120324$HGVSc<-str_replace(env_pot1_clinvar_all_120324$HGVSc, "NM_015450.3", "ENST00000357628.8")
env_pot1_clinvar_all_120324<-env_pot1_clinvar_all_120324 %>% select(c(10,26,27))

env_pot1_clinvar_all_120324<-env_pot1_clinvar_all_120324 %>% filter(str_detect(HGVSc,"ENST00000357628.8"))



env_pot1_clinvar_onestar_120324<-pot1_clinvar_onestar_120324
env_pot1_clinvar_onestar_120324$clinvar_onestar<-"Y"

env_pot1_clinvar_onestar_120324$CLIN_SIG<-gsub("\\s*\\([^\\)]+\\)","",as.character(env_pot1_clinvar_onestar_120324$Germline.classification))
env_pot1_clinvar_onestar_120324$HGVSc<-gsub("\\s*\\([^\\)]+\\)","",as.character(env_pot1_clinvar_onestar_120324$Name))
env_pot1_clinvar_onestar_120324$HGVSc<-str_replace(env_pot1_clinvar_onestar_120324$HGVSc, "NM_015450.3", "ENST00000357628.8")
filtered_env_pot1_clinvar_onestar_120324<-env_pot1_clinvar_onestar_120324 %>% select(c(26,28))

filtered_env_pot1_clinvar_onestar_120324<-filtered_env_pot1_clinvar_onestar_120324 %>% filter(str_detect(HGVSc,"ENST00000357628.8"))


merged_clinvar<-merge(env_pot1_clinvar_all_120324,filtered_env_pot1_clinvar_onestar_120324, by="HGVSc", all.x=TRUE)

merged_clinvar<-merged_clinvar %>% mutate(clinvar_onestar=case_when(is.na(clinvar_onestar) ~ "N", TRUE~ "Y"))



march_pot1_24_clinvar<-merged_clinvar 
clean_pot1_v2<-merge(clean_pot1_v1,march_pot1_24_clinvar, by="HGVSc",all.x=TRUE)





#dataframe for the plotting on database information (clinvar and gnomad)
databases<-clean_pot1_v2 %>% mutate(is_in_clinvar=case_when(is.na(VariationID) ~ "N", TRUE ~ "Y")) %>%
  mutate(is_in_gnomad=case_when(gnomAD_AF == "-" ~ "N", TRUE ~ "Y")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" ~ "Unobserved", TRUE ~ "ClinVar only")) %>%
  mutate(present_accession=case_when(is_in_clinvar == "N" & is_in_gnomad =="Y" ~ "gnomAD only", TRUE ~present_accession)) %>%
  mutate(present_accession=case_when(is_in_clinvar == "Y" & is_in_gnomad =="Y" ~ "ClinVar and gnomAD", TRUE ~present_accession)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Uncertain significance" ~ "Uncertain significance")) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(is.na(CLIN_SIG.y) ~ "Unobserved", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Pathogenic/Likely pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Likely pathogenic" ~ "Pathogenic/Likely pathogenic", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Likely benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Benign/Likely benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Benign" ~ "Benign/Likely benign", TRUE ~ClinVar_CLNSIG_slim)) %>%
  #mutate(ClinVar_CLNSIG_slim=case_when(ClinVar_CLNSIG=="risk_factor" ~ "Risk Factor", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(ClinVar_CLNSIG_slim=case_when(CLIN_SIG.y=="Conflicting classifications of pathogenicity" ~ "Conflicting interpretation", TRUE ~ClinVar_CLNSIG_slim)) %>%
  mutate(present_accession_filter=paste0(present_accession,"_",mseq))
#reorder the dataframe by present_accession factor
database_levels<-c("ClinVar and gnomAD","ClinVar only","gnomAD only","Unobserved")
#re-order by present_accession to set levels
databases<-databases %>% mutate(present_accession =  factor(present_accession, levels = database_levels)) %>%
  arrange(present_accession)
#save full dataframe
databases_expanded <- databases


#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
databases_expanded<-databases_expanded %>% arrange(-mut_position) %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)

databases_expanded_cleaned<-databases_expanded %>% filter(!HGVSc %in% "-") %>% dplyr::rename(id= "id.y")


unique<-databases_expanded_cleaned[!duplicated(databases_expanded_cleaned$mseq), ]



######CLUSTERING and CLASSIFICATION
#clustering
cluster.frame<-unique 

#stratify for variants that have FDR>0.01
unchanged<-cluster.frame %>% filter(functional_classification %in% "unchanged")
enriched<-cluster.frame  %>% filter(functional_classification %in% "enriched") 
depleted<-cluster.frame  %>% filter(functional_classification %in% "depleted")


library(mclust)
set.seed(0)
EEV <- Mclust(depleted%>% dplyr::select(c("processed_LFC_D4_D21","processed_LFC_D4_D7")), G=2, modelNames = "EEV")
EEV_plot<-plot(EEV, what="classification", xlim=c(-5,2), ylim=c(-5,2), xlab="D4D21", ylab="D4D7", symbols=c(3,0,4,2))
dev.off()

all_depleted <- cbind(depleted, EEV$z %>% as.data.frame() %>% `colnames<-` (c("C1_z","C2_z")),EEV$classification %>% as.data.frame() %>% `colnames<-` (c("Cluster")), EEV$uncertainty %>% as.data.frame() %>% `colnames<-` (c("uncertainty")))
enriched$Cluster<-3
unchanged$Cluster<-4

rates<- do.call("rbind.fill", list(all_depleted, enriched, unchanged))
rates$Cluster<-as.factor(rates$Cluster)

levels(rates$Cluster)[levels(rates$Cluster)=='1'] <- 'strong depleted'
levels(rates$Cluster)[levels(rates$Cluster)=='2'] <- 'weak depleted'
levels(rates$Cluster)[levels(rates$Cluster)=='3'] <- 'enriched'
levels(rates$Cluster)[levels(rates$Cluster)=='4'] <- 'unchanged'

rates$Cluster <- relevel(rates$Cluster, "unchanged", "strong depleted", "weak depleted", "enriched")


rates<-rates %>% select(2,109:112)


databases_expanded_cleaned_rate<-merge(databases_expanded_cleaned,rates,by="mseq")


#re-order by mutator according to set levels
clean_one<-databases_expanded_cleaned_rate %>% arrange(-mut_position) %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)



#FLAGGING OF primer bindig site proximity PPE codons

flagged_version<-clean_one %>% mutate(primer_flag=case_when(mut_position %in% 124892388:124892384~ "Y", TRUE~ "N")) %>%
  mutate(primer_flag=case_when(mut_position %in% 124870909:124870905~ "Y", TRUE~ primer_flag)) %>%
  mutate(primer_flag=case_when(mut_position %in% 124859117:124859113~ "Y", TRUE~ primer_flag)) %>%  
  mutate(primer_flag=case_when(mut_position %in% 124825254:124825250~ "Y", TRUE~ primer_flag))


#simplify the pam status field and create a flag field for variants that only come from one dataset and are pam variant codon
flagged_version<- flagged_version %>% 
  mutate(pam_codon=case_when(str_detect(PAM_status,"sgRNA") ~ "Y", TRUE ~ "N"))


#column clean up 

flagged_version<-flagged_version %>% relocate(35,36,.after=id) 

flagged_version<-flagged_version %>% relocate(111,.after=5) 

flagged_version<-flagged_version %>% relocate(113,114,.after=7) 

flagged_version_cleaned <- flagged_version[ -c(67,69,70,72,97,100,101,102) ]


library(openxlsx)
write.xlsx(flagged_version_cleaned, 'pot1_expanded.xlsx')

















#4th april - wait for prashant to re-run


#prashant added splice ai values and fixed chorm pos ref alt, one variant needed to be removed as could not be 
#mapped to GRCH38 due to indel and containing PPE mutation this variant stripped out: "ENST00000357628.8.ENSG00000128513.16_chr7:124892332_124892347_clinvar_20201107_rc"



#new file with splice ai and correct mappings is:

pot1_expanded_fix_updated<-pot1_expanded_updated_fix



#some counts were in the annotation file which is present bound to the final dataset that were excluded at the analysis stage
#make these values NA as they were not included in the calculations

flagged_version_cleaned_NAS<-pot1_expanded_fix_updated %>% 
  mutate(D4R3=ifelse(str_detect(source_targetons,"9_1"),NA,D4R3)) %>% 
  mutate(D7R3=ifelse(str_detect(source_targetons,"9_1"),NA,D7R3)) %>% 
  mutate(D10R3=ifelse(str_detect(source_targetons,"9_1"),NA,D10R3)) %>% 
  mutate(D14R3=ifelse(str_detect(source_targetons,"9_1"),NA,D14R3)) %>% 
  mutate(D21R3=ifelse(str_detect(source_targetons,"9_1"),NA,D21R3)) %>% 
  
  mutate(D4R1=ifelse(str_detect(source_targetons,"14"),NA,D4R1)) %>% 
  mutate(D7R1=ifelse(str_detect(source_targetons,"14"),NA,D7R1)) %>% 
  mutate(D10R1=ifelse(str_detect(source_targetons,"14"),NA,D10R1)) %>% 
  mutate(D14R1=ifelse(str_detect(source_targetons,"14"),NA,D14R1)) %>% 
  mutate(D21R1=ifelse(str_detect(source_targetons,"14"),NA,D21R1)) %>% 
  
  mutate(D7R3=ifelse(str_detect(source_targetons,"10"),NA,D7R3))



#merge EVE and ddG scores




#EVOLUTIONARY BOXPLOT

#GET THE EVE SCORES FROM VEP using HGVSc to pull 
eve_hgvsc<-clean_one_unique %>% select("HGVSc")
write.table(eve_hgvsc, file = "eve_hgvsc.txt", sep = "\t", row.names=FALSE, col.names = FALSE, quote=FALSE)

#import the VEP output
eve_pot1<-ZkZpM15Zt43ZfEkZ

#filter the data
evo_EVE_2<-eve_pot1 %>% filter (Feature  %in% "ENST00000357628.8") %>% select('X.Uploaded_variation',"EVE_CLASS","EVE_SCORE","CADD_PHRED","CADD_RAW") %>% dplyr::rename(HGVSc = "X.Uploaded_variation")
eve_EVE_2<-merge(clean_one_unique,evo_EVE_2,by="HGVSc",all.x=TRUE)
eve_EVE_2<- eve_EVE_2[!duplicated(eve_EVE_2$mseq), ]
eve_EVE_2<- eve_EVE_2[!duplicated(eve_EVE_2$HGVSc), ]
eve_EVE_2$EVE_SCORE<-as.numeric(eve_EVE_2$EVE_SCORE)

#only snvs

#eve_EVE_3<-eve_EVE_2 %>% filter(!PAM_status_simplified_collated %in% "A" | !PAM_status_simplified_collated %in% "B")

EVE_TO_MERGE<-eve_EVE_2 %>% select(HGVSc, EVE_SCORE, EVE_CLASS,CADD_PHRED,CADD_RAW) %>% filter(!is.na(EVE_SCORE)) %>% filter(!EVE_SCORE %in% "-")

#DELTA_TO_MERGE<-delta_merged %>% select(HGVSc,ddG)

#scores_to_merge<-merge(EVE_TO_MERGE,DELTA_TO_MERGE,by="HGVSc",all.x = TRUE)

complete_scores<-merge(flagged_version_cleaned_NAS,EVE_TO_MERGE,by="HGVSc",all.x = TRUE)

#re-order by mutator according to set levels
complete_scores<-complete_scores %>% arrange(X) 

#final column organisation

coulmn_names<-as.data.frame(colnames(complete_scores))

write.csv(coulmn_names,"./coulmn_names.csv", row.names = FALSE)





# FINAL re-organisation


complete_final<-complete_scores[, c(
  "mseq",
  "HGVSc",
  "HGVSp",
  "id",
  "fixed_pos",
  "fixed_ref",
  "fixed_alt",
  "Cluster",
  "functional_score",
  "EXON",
  "INTRON",
  "slim_consequence",
  "cDNA_position",
  "CDS_position",
  "Protein_position",
  "Amino_acids",
  "ref_aa",
  "alt_aa",
  "Codons",
  "gene_id",
  "transcript_id",
  "REGION",
  "source_targetons",
  "primer_flag",
  "pam_codon",
  "mutator",
  "mave_nt",
  "Consequence",
  "DOMAINS",
  "HGVS_OFFSET",
  "processed_LFC_D4_D7",
  "processed_LFC_D4_D10",
  "processed_LFC_D4_D14",
  "processed_LFC_D4_D21",
  "processed_LFC_continuous",
  "processed_BH_FDR_D4_D7",
  "processed_BH_FDR_D4_D10",
  "processed_BH_FDR_D4_D14",
  "processed_BH_FDR_D4_D21",
  "processed_BH_FDR_continuous",
  "SE_bind_D4_D7",
  "SE_bind_D4_D10",
  "SE_bind_D4_D14",
  "SE_bind_D4_D21",
  "SE_bind_continuous",
  "processed_Z_D4_D7",
  "processed_Z_D4_D10",
  "processed_Z_D4_D14",
  "processed_Z_D4_D21",
  "processed_Z_continuous",
  "two_tailed_p_D4_D7",
  "two_tailed_p_D4_D10",
  "two_tailed_p_D4_D14",
  "two_tailed_p_D4_D21",
  "two_tailed_p_continuous",
  "C1_z",
  "C2_z",
  "uncertainty",
  "D4R1",
  "D4R2",
  "D4R3",
  "D7R1",
  "D7R2",
  "D7R3",
  "D10R1",
  "D10R2",
  "D10R3",
  "D14R1",
  "D14R2",
  "D14R3",
  "D21R1",
  "D21R2",
  "D21R3",
  "PLASMID",
  "SIFT",
  "PolyPhen",
  "EVE_SCORE",
  "EVE_CLASS",
  "CADD_PHRED",
  "CADD_RAW",
  #"ddG",
  "DS_AG",
  "DS_AL",
  "DS_DG",
  "DS_DL",
  "DP_AG",
  "DP_AL",
  "DP_DG",
  "DP_DL",
  "REF_DS_AG",
  "REF_DS_AL",
  "REF_DS_DG",
  "REF_DS_DL",
  "ALT_DS_AG",
  "ALT_DS_AL",
  "ALT_DS_DG",
  "ALT_DS_DL",
  "AF",
  "gnomAD_AF",
  "gnomAD_AFR_AF",
  "gnomAD_AMR_AF",
  "gnomAD_ASJ_AF",
  "gnomAD_EAS_AF",
  "gnomAD_FIN_AF",
  "gnomAD_NFE_AF",
  "gnomAD_OTH_AF",
  "gnomAD_SAS_AF",
  "gnomAD_FLAG",
  "dbSNP",
  "VariationID",
  "is_in_clinvar",
  "clinvar_onestar",
  "CLIN_SIG.y",
  "ClinVar_CLNSIG_slim",
  "is_in_gnomad",
  "present_accession",
  "present_accession_filter")]



#rename columns 
complete_final<-complete_final %>% dplyr::rename(alt="fixed_alt")
complete_final<-complete_final %>% dplyr::rename(BH_FDR_continuous="processed_BH_FDR_continuous")
complete_final<-complete_final %>% dplyr::rename(BH_FDR_D4_D10="processed_BH_FDR_D4_D10")
complete_final<-complete_final %>% dplyr::rename(BH_FDR_D4_D14="processed_BH_FDR_D4_D14")
complete_final<-complete_final %>% dplyr::rename(BH_FDR_D4_D21="processed_BH_FDR_D4_D21")
complete_final<-complete_final %>% dplyr::rename(BH_FDR_D4_D7="processed_BH_FDR_D4_D7")
complete_final<-complete_final %>% dplyr::rename(clinical_significance="CLIN_SIG.y")
complete_final<-complete_final %>% dplyr::rename(clinical_significance_slim="ClinVar_CLNSIG_slim")
complete_final<-complete_final %>% dplyr::rename(exon="EXON")
complete_final<-complete_final %>% dplyr::rename(functional_classificaiton="Cluster")
complete_final<-complete_final %>% dplyr::rename(functional_score="functional_score")
complete_final<-complete_final %>% dplyr::rename(intron="INTRON")
complete_final<-complete_final %>% dplyr::rename(LFC_continuous="processed_LFC_continuous")
complete_final<-complete_final %>% dplyr::rename(LFC_D4_D10="processed_LFC_D4_D10")
complete_final<-complete_final %>% dplyr::rename(LFC_D4_D14="processed_LFC_D4_D14")
complete_final<-complete_final %>% dplyr::rename(LFC_D4_D21="processed_LFC_D4_D21")
complete_final<-complete_final %>% dplyr::rename(LFC_D4_D7="processed_LFC_D4_D7")
complete_final<-complete_final %>% dplyr::rename(pos="fixed_pos")
complete_final<-complete_final %>% dplyr::rename(ref="fixed_ref")
complete_final<-complete_final %>% dplyr::rename(SE_continuous="SE_bind_continuous")
complete_final<-complete_final %>% dplyr::rename(SE_D4_D10="SE_bind_D4_D10")
complete_final<-complete_final %>% dplyr::rename(SE_D4_D14="SE_bind_D4_D14")
complete_final<-complete_final %>% dplyr::rename(SE_D4_D21="SE_bind_D4_D21")
complete_final<-complete_final %>% dplyr::rename(SE_D4_D7="SE_bind_D4_D7")
complete_final<-complete_final %>% dplyr::rename(target_region="REGION")
complete_final<-complete_final %>% dplyr::rename(target_region_data_source="source_targetons")
complete_final<-complete_final %>% dplyr::rename(variation_id="VariationID")
complete_final<-complete_final %>% dplyr::rename(vep_consequence="Consequence")
complete_final<-complete_final %>% dplyr::rename(vep_consequence_slim="slim_consequence")
complete_final<-complete_final %>% dplyr::rename(z_continuous="processed_Z_continuous")
complete_final<-complete_final %>% dplyr::rename(z_D4_D10="processed_Z_D4_D10")
complete_final<-complete_final %>% dplyr::rename(z_D4_D14="processed_Z_D4_D14")
complete_final<-complete_final %>% dplyr::rename(z_D4_D21="processed_Z_D4_D21")
complete_final<-complete_final %>% dplyr::rename(z_D4_D7="processed_Z_D4_D7")

complete_final<-complete_final %>% dplyr::rename(Ref_AA="ref_aa")
complete_final<-complete_final %>% dplyr::rename(Alt_AA="alt_aa")



#make chrom pos ref alt
complete_final$chrom_pos_ref_alt<-paste0("7_",complete_final$pos,"_",complete_final$ref,"_",complete_final$alt)


# rearrange
complete_final<-complete_final %>% relocate(chrom_pos_ref_alt,.after=id)

#complete_final$pos<-NULL
#complete_final$ref<-NULL
#complete_final$alt<-NULL

complete_final<-complete_final %>% relocate(mseq,.after=transcript_id)
complete_final<-complete_final %>% relocate(id,.after=mseq)

# correct targeton naming system to match exon numbers

complete_final_fix<-complete_final %>% 
  mutate(target_region_rep="-") %>%
  mutate(target_region_rep=ifelse(target_region %in% "2","6",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "3","7",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "4","8",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "5","9",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "6","10",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "7","11",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "9","13",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "10","14",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "11","15",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "12","16",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "13","17",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "14","18",target_region_rep)) %>% 
  mutate(target_region_rep=ifelse(target_region %in% "15","19",target_region_rep))

complete_final_fix<-complete_final_fix %>% relocate(target_region_rep,.after=target_region)

complete_final_fix_again<-complete_final_fix %>% 
  mutate(target_region_data_source_rep="-") %>%
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "1","6",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "2","7",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "3_1","8_1",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "3_1|3_2","8-1|8_2",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "3_2","8_2",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "4","9",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "5","10",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "6","11",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "8","13",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "9_1","14_1",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "9_1|9_2","14_1|14_2",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "9_2","14_2",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "10","15",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "11","16",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "12","17",target_region_data_source_rep)) %>% 
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "13","18",target_region_data_source_rep)) %>%
  mutate(target_region_data_source_rep=ifelse(target_region_data_source %in% "14","19",target_region_data_source_rep))


complete_final_fix_again<-complete_final_fix_again %>% relocate(target_region_data_source_rep,.after=target_region_data_source)


complete_final_fix_again$target_region<-NULL
complete_final_fix_again$target_region_data_source<-NULL


complete_final_fix_again<-complete_final_fix_again %>% dplyr::rename(target_region="target_region_rep")
complete_final_fix_again<-complete_final_fix_again %>% dplyr::rename(target_region_data_source="target_region_data_source_rep")



#THRESHOLDS
complete_final_fix_again<-complete_final_fix_again %>% mutate(functional_score_below_threshold = case_when(functional_score < -0.03359 ~ "T", TRUE ~ "F"))
complete_final_fix_again<-complete_final_fix_again %>% mutate(functional_score_below_threshold = case_when(functional_classificaiton %in% "unchanged" ~ NA, TRUE ~ functional_score_below_threshold))
complete_final_fix_again<-complete_final_fix_again %>% mutate(functional_score_below_threshold = case_when(functional_classificaiton %in% "enriched" ~ NA, TRUE ~ functional_score_below_threshold))


complete_final_fix_again<-complete_final_fix_again %>% relocate(functional_score_below_threshold,.before=primer_flag)



##need to address N/A in amino acid column - could change in excel (asn>alanine single letter is N/A which R and excel think are NA)

complete_final_fix_again<-complete_final_fix_again %>% 
  mutate(Amino_acids=case_when(Codons %in% "AAT/GCC" ~ "N/A", TRUE ~ Amino_acids)) %>%
  mutate(Amino_acids=case_when(Codons %in% "AAC/GCC" ~ "N/A", TRUE ~ Amino_acids))




library(openxlsx)
#output final expanded dataframe
write.xlsx(complete_final_fix_again, 'pot1_sge_expanded_annotations.xlsx')



#output unique data frame
complete_final_unique<-complete_final_fix_again[!duplicated(complete_final_fix_again$mseq), ]

complete_final_unique<-complete_final_unique[!duplicated(complete_final_unique$HGVSc), ]

write.xlsx(complete_final_unique, 'pot1_sge.xlsx')



####### ENDS 


## ADDITION OF MISSING PLASMID SEQUENCING DATA - FINAL DATAFRAME - START

#add in missing plasmid data for exon 14b (now 19) and exon 3_1a (now 8_1)

pre_plasmid_complete<-complete_final_fix_again

plasmid_3_1_a<-pot1_exon_3_1_a_lib.counts

plasmid_3_1_a<-plasmid_3_1_a %>% select(NAME,reads_pot1_exon_3_1_a_lib)
plasmid_3_1_a<-plasmid_3_1_a %>% dplyr::rename(id="NAME") %>% dplyr::rename(PLASMID_3="reads_pot1_exon_3_1_a_lib")

ADDED_3<-merge(pre_plasmid_complete,plasmid_3_1_a,by="id",all.x=TRUE)


plasmid_14_b<-POT1_14_B_PLASMID_LIB.lib_counts

plasmid_14_b<-plasmid_14_b %>% select(NAME,COUNT)
plasmid_14_b<-plasmid_14_b %>% dplyr::rename(id="NAME") %>% dplyr::rename(PLASMID_14="COUNT")

ADDED_3_14<-merge(ADDED_3,plasmid_14_b,by="id",all.x=TRUE)

#final frame
FINAL_expanded_plus_plasmid<-ADDED_3_14 %>% mutate(PLASMID=case_when(is.na(PLASMID) & str_detect(target_region_data_source,"8_1") ~ PLASMID_3, TRUE ~ PLASMID)) %>%
  mutate(PLASMID=case_when(is.na(PLASMID) & str_detect(target_region_data_source,"8-1") ~ PLASMID_3, TRUE ~ PLASMID)) %>%
  mutate(PLASMID=case_when(is.na(PLASMID) & str_detect(target_region_data_source,"19") ~ PLASMID_14, TRUE ~ PLASMID))


FINAL_expanded_plus_plasmid$PLASMID_3<-NULL
FINAL_expanded_plus_plasmid$PLASMID_14<-NULL

FINAL_expanded_plus_plasmid<-FINAL_expanded_plus_plasmid %>% relocate(id,.after=mseq)


FINAL_expanded_plus_plasmid$mapping<-paste0(FINAL_expanded_plus_plasmid$Ref_AA,"_",FINAL_expanded_plus_plasmid$Protein_position,"_",FINAL_expanded_plus_plasmid$Alt_AA)


#add delta delta G scores

delta<-test

delta$mapping<-paste0(delta$aa_wt,"_",delta$uniprot_pos,"_",delta$aa_mt)

delta_merged<-merge(FINAL_expanded_plus_plasmid,delta,by="mapping",all.x=TRUE)

delta_merged$mapping<-NULL


delta_merged<-delta_merged %>% select(!119:128)
delta_merged<-delta_merged %>% select(!120:121)


FINAL_expanded_plus_plasmid<-delta_merged %>% relocate(ddG,.after=CADD_RAW)




#reorder the dataframe by mutator factor
x_levels<-c("snv","1del","custom","2del0","2del1","inframe","stop","ala","snvre")

#re-order by mutator according to set levels
FINAL_expanded_plus_plasmid<-FINAL_expanded_plus_plasmid %>% arrange(-pos) %>% mutate(mutator =  factor(mutator, levels = x_levels)) %>%
  arrange(mutator)




#UPDATE EVE AND DDG VALUES - MERGE BY HGVSC 







#output final expanded dataframe
write.xlsx(FINAL_expanded_plus_plasmid, 'pot1_sge_expanded_annotations_updated.xlsx')



#output unique data frame
complete_final_unique<-FINAL_expanded_plus_plasmid[!duplicated(FINAL_expanded_plus_plasmid$mseq), ]

complete_final_unique<-complete_final_unique[!duplicated(complete_final_unique$HGVSc), ]

write.xlsx(complete_final_unique, 'pot1_sge_updated.xlsx')


## ADDITION OF MISSING PLASMID SEQUENCING DATA - FINAL DATAFRAME - END




########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
########################## END #########################################################
















