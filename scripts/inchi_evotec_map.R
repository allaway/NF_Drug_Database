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
alldat$inchikey <- gsub("InChIKey=", "", alldat$inchikey)

##59 smiles didnt map to inchikeys, these were mapped with pubchem ID Exchange service 
unmapped<-dplyr::select(filter(alldat, is.na(inchikey)),smiles)
write.table(unmapped, "unmapped.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
nowmapped<-read.table("unmappedNowMapped.txt", comment.char = "", sep = "\t")
colnames(nowmapped)<-c("smiles", "inchikey")
nowmapped<-dplyr::select(nowmapped, inchikey, smiles)
alldat<-bind_rows(filter(alldat, !is.na(inchikey)), nowmapped)

#nine remain unmapped (really noncanonical smiles?) - map manually with pubchem structure search
alldat$inchikey[alldat$smiles=="CC(C)[C@@H]1NC(=O)[C@H](CCCCN)NC(=O)[C@@H](CC2=CNC3=C2C=CC=C3)NC(=O)[C@H](CC2=CC=C(O)C=C2)NC(=O)[C@H](C)N(C)C(=O)[C@H](CC2=CC=CC=C2)NC1=O"] <- "NPJIOCBFOAHEDO-AVWFULIKSA-N"
alldat$inchikey[alldat$smiles=="BrCCCBr.CN(C)CCCCCCN(C)C"] <- "KZKAYEGOIJEWQB-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="[H]N1C(=O)C(C%23N)=C(NC2=CC(C)=CC=C2)SC11CCCC1"] <- "PZKAUBRBLGYPFT-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="C1=CC2=CC3=C(C=CC=C3)C=C2C=C1"] <- "MWPLVEDNUUSJAV-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="CC1=NC(C(=O)NCC(O)=O)=C(O)C2=CC=C(OC3=CC=CC=C3)C=C12"] <- "YOZBGTLTNGAVFU-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="[Na+].[Na+].[Na+].CC1=C(O)C(C=O)=C(COP([O-])([O-])=O)C(N=NC2=CC=C(C=C2)C([O-])=O)=N1"] <- "PBGKQYOUZYMMOR-UHFFFAOYSA-K"
alldat$inchikey[alldat$smiles=="CC1=NC2=C(C=C(C=C2N1CC1=CC=CC(=C1C)C(F)(F)F)N1CCOCC1)C(O)=O"] <- "XTKLTGBKIDQGQL-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="CC1=NN(C(C)=C1C=NN1CCN(CC2=CC=CC=C2)CC1)C1=CC=CC=C1"] <- "FOORCIAZMIWALX-UHFFFAOYSA-N"
alldat$inchikey[alldat$smiles=="ClC1=C(N=C2NN=C(NC(=O)C3CC3)C2=C1)C1=CC=CC=C1"] <- "DOKPTPBHIRPUBR-UHFFFAOYSA-N"


write.table(alldat, "inchikey_smiles_map.txt", sep = "\t", row.names = FALSE)
synStore(File("inchikey_smiles_map.txt", parentId = "syn8682571"), executed = this.file, used = c("syn8118065","syn7341038", "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi", "https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds&query_type=structure&query_subtype=identity"))


