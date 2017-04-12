library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
synapseLogin()

inchi<-read.table(synGet("syn8682658")@filePath, header = TRUE)

## pull drug data and filter for human targets, and eliminate drugs with
## 0 quantitative effects measured
drugdat <- synTableQuery("SELECT * FROM syn7341038")
drugdat <- as.data.frame(drugdat@values)
drugdat.filtered <- filter(drugdat, Organism == "Homo sapiens")

## obtain data to map Uniprot id with Hugo Genes
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "dec2016.archive.ensembl.org")
bm <-
  getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"),
        mart = mart)
bm <- filter(bm, uniprot_swissprot != "")
colnames(bm) <- c("Hugo_Gene", "Uniprot_accession_numbers")

## pull evotec data to map compound names to structure IDs
compound.data <- synTableQuery("SELECT * FROM syn8118065")
compound.data <- as.data.frame(compound.data@values)
compound.data <-
  compound.data[!duplicated(compound.data$Structure_ID),]

compounds <- full_join(drugdat.filtered, compound.data)

inchi$smiles <- gsub("%23", "#", inchi$smiles)
colnames(inchi) <- c("inchikey", "Original_molecule_SMILES")

evotec<-full_join(compounds,inchi)
evotec<-full_join(evotec,bm)
evotec<-dplyr::select(evotec, Structure_ID, Original_molecule_SMILES, inchikey, Hugo_Gene, Uniprot_accession_numbers, Protein_names, Organism, Min_activity_operator,
      MinActivity_nM, MeanActivity_nM, MedianActivity_nM, Geom_meanActivity_nM, Max_activity_operator, MaxActivity_nM, N_quantitative, N_qualitative, N_inactive,
      Supplier, Supplier_ID, Supplier_Molname, Supplier_Data_1, Supplier_Data_2, Supplier_Data_3)

write.table(evotec, "evotec_with_inchikey.txt", sep = "\t")
synStore(File("evotec_with_inchikey.txt", parentId = "syn8682571"), used = c("syn8118065","syn7341038","syn8682658"), 
         executed = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/map_inchi_to_evotec_by_smiles.R")


