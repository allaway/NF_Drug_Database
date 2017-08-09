library(rcdk)
library(fingerprint)
library(synapseClient)
library(plyr)
library(dplyr)
library(parallel)
synapseLogin()

evo<-synTableQuery("SELECT * FROM syn7341038")
evo <- evo@values
evo$Uniprot_accession_numbers <- gsub(",.+$", "", evo$Uniprot_accession_numbers)
evo <- filter(evo, Organism == "Homo sapiens")

## obtain data to map Uniprot id with Hugo Genes
library(UniProt.ws)
up <- UniProt.ws(taxId=9606)
keys <- unique(evo$Uniprot_accession_numbers)
columns <- c("HGNC")
kt <- "UNIPROTKB"
res <- select(up, keys, columns, kt)
colnames(res) <- c("Uniprot_accession_numbers", "HGNC")
evo <- left_join(evo, res)

##get smiles files
evotec.smiles <- synTableQuery("SELECT * FROM syn8118065")@values
evotec.smiles <- evotec.smiles %>% distinct() %>% group_by(Structure_ID) %>% sample_n(1)
ncats<-read.table("ncatssmiles.txt",header = F, sep = "\t", comment.char = "", quote = "") %>% filter(V2 != "") %>% distinct()

ncats.mol <- parse.smiles(as.character(unique(ncats$V2)))
evotec.mol <- parse.smiles(as.character(evotec.smiles$Original_molecule_SMILES))

fp.nca <- lapply(ncats.mol, get.fingerprint, type='extended')
fp.evo <- lapply(evotec.mol, get.fingerprint, type = 'extended')

#to find bad smiles strings if get this err: Error in FUN(X[[i]], ...) : 
#must supply an IAtomContainer or something coercable to it

for(i in 1:8457){
   if(is.null(attr(evotec.mol[[i]], "jclass"))) print(i)
}
sims <- mclapply(fp.nca, function(i){
  sim <- sapply(fp.evo, function(j){
    distance(i, j)
  })
  bar <- as.data.frame(sim)
  bar$match <- rownames(bar)
  bar <- bar %>% top_n(1, sim) %>% sample_n(1)
}, mc.cores=detectCores())


##calculating inter-evotec sims for other purposes
names(fp.evo) <- evotec.smiles$Structure_ID
sims2 <- as.data.frame(fp.sim.matrix(fp.evo, method = "tanimoto"))
colnames(sims2) <- names(fp.evo)
sims2$Structure_ID <- as.integer(names(fp.evo))
map <- dplyr::select(evotec.smiles, Original_molecule_SMILES, Structure_ID)
sims2 <- left_join(sims2, map)

write.table(sims2, "evotecsims.txt")
synStore(File("evotecsims.txt", parentId = "syn7287882"), used = c("syn8118065"), executed = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/molecular_mapping.R")

temp <- ldply(sims)
temp$original <- temp$.id
temp$evo.match<- temp$match
temp <- select(temp, original, evo.match, sim)

colnames(temp)[1] <- "ncats_mipe"
colnames(ncats)[2] <- "ncats_mipe"

map <- left_join(ncats, temp) %>% distinct() %>% group_by(V1) %>% top_n(1, sim) %>% sample_n(1) %>% ungroup()

colnames(evotec.smiles) <- c("Original_molecule_SMILES", "Structure_ID")
colnames(map)[2] <- "Original_molecule_SMILES"
map <- left_join(map, evotec.smiles)

map <- left_join(map, evo) %>% filter(N_quantitative > 0 | N_qualitative > 0) %>% filter(sim >= 0.9)



outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

outersect(map$V1, ncats$V1)
