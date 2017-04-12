library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
library(webchem)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/inchi_evotec_map.R"

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

## some uniprot ids do not successfully map to hugo genes, note them
## here, these will not show up in analysis
unannotated <-
  filter(drugdat.filtered, is.na(drugdat.filtered$Hugo_Gene))

## pull evotec data to map compound names to structure IDs
compound.data <- synTableQuery("SELECT * FROM syn8118065")
compound.data <- as.data.frame(compound.data@values)
compound.data <-
  compound.data[!duplicated(compound.data$Structure_ID),]

compounds <- full_join(drugdat.filtered, compound.data)
compounds$Original_molecule_SMILES <- gsub("#", "%23", compounds$Original_molecule_SMILES)
smiles<-unique(compounds$Original_molecule_SMILES)

data<-cir_query(unique(smiles[1:2000]), representation = "stdinchikey", first = TRUE, verbose = TRUE)
data1<-cir_query(unique(smiles[2001:4000]), representation = "stdinchikey", first = TRUE, verbose = TRUE)
data2<-cir_query(unique(smiles[4001:6000]), representation = "stdinchikey", first = TRUE, verbose = TRUE)
data3<-cir_query(unique(smiles[6001:8000]), representation = "stdinchikey", first = TRUE, verbose = TRUE)
data4<-cir_query(unique(smiles[8001:8145]), representation = "stdinchikey", first = TRUE, verbose = TRUE)
#one structure has a messed up smiles!
data5<-cir_query(unique(smiles[8147:8482]), representation = "stdinchikey", first = TRUE, verbose = TRUE)

data<-data.frame(data)
data1<-data.frame(data1)
data2<-data.frame(data2)
data3<-data.frame(data3)
data4<-data.frame(data4)
data5<-data.frame(data5)

colnames(data) <- "inchi"
colnames(data1) <- "inchi"
colnames(data2) <- "inchi"
colnames(data3) <- "inchi"
colnames(data4) <- "inchi"
colnames(data5) <- "inchi"

alldat<-rbind(data,data1,data2,data3,data4,data5)
alldat$smiles <- rownames(alldat)

write.table(alldat, "inchikey_smiles_map.txt", sep = "\t", row.names = FALSE)
synStore(File("inchikey_smiles_map.txt", parentId = "syn8682571"), executed = this.file, used = c("syn8118065","syn7341038"))
