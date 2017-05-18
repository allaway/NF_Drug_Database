library(synapseClient)
library(dplyr)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/map_pNF_to_evotec.R"

##take these to pubchem id exchange, uploaded to synapse, read in
ncats_inchikey<-read.table(synGet("syn8690496")@filePath, header = FALSE, sep = "\t", comment.char = "")
colnames(ncats_inchikey) <- c("NCGC.SID", "inchikey")
evotec<-read.table(synGet("syn8685703")@filePath, header = TRUE, sep = "\t", comment.char = "") %>% dplyr::select(inchikey, Hugo_Gene, Structure_ID)

evotec <- evotec %>% semi_join(evotec %>% count(Structure_ID) %>% filter(n < 250))

nf2.ncats<-read.table(synGet("syn8314523")@filePath, header = TRUE, sep = "\t", comment.char = "") %>% select(-Protocol.Name)

##schwannoma

hs01<-nf2.ncats %>% filter(Cell.line == "HS01") %>% select(-Cell.line)
hs11<-nf2.ncats %>% filter(Cell.line == "HS12") %>% select(-Cell.line)

colnames(hs01) <- c("NCGC.SID", "Cell.Type", "NF2.Status.HS01", "AC50.HS01", "CRC.HS01", 
                    "Max.Resp.HS01", "Log.AC50.uM.HS01", "Sample.Name", "AUC.HS01", "AUC.Fit.HS01", "Gene.Symbol")
colnames(hs11) <- c("NCGC.SID", "Cell.Type", "NF2.Status.HS12", "AC50.HS12", "CRC.HS12", 
                    "Max.Resp.HS12", "Log.AC50.uM.HS12", "Sample.Name", "AUC.HS12", "AUC.Fit.HS12", "Gene.Symbol")

schwannoma <- full_join(hs01, hs11)

schwannoma$AUC.Ratio <- (schwannoma$AUC.HS01/schwannoma$AUC.HS12)

schwannoma.up <- full_join(schwannoma, ncats_inchikey) %>% full_join(evotec) %>% filter(AUC.Ratio>=1.25) %>% 
  select(NCGC.SID, Hugo_Gene, AUC.Ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA")

nodes <- data.frame(c(unique(as.character(schwannoma.up$Hugo_Gene)), unique(schwannoma.up$NCGC.SID)))
write.table(nodes,"schwannoma_up_nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

edges <- schwannoma.up
write.table(edges,"schwannoma_up_edges.txt", sep= "\t", row.names = FALSE, quote = FALSE)


schwannoma.down <- full_join(schwannoma, ncats_inchikey) %>% full_join(evotec) %>% filter(AUC.Ratio<=0.75) %>% 
  select(NCGC.SID, Hugo_Gene, AUC.Ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA")

nodes <- data.frame(c(unique(as.character(schwannoma.down$Hugo_Gene)), unique(schwannoma.down$NCGC.SID)))
write.table(nodes,"schwannoma_down_nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

edges <- schwannoma.down
write.table(edges,"schwannoma_down_edges.txt", sep= "\t", row.names = FALSE, quote = FALSE)

##meningioma

syn5<-nf2.ncats %>% filter(Cell.line == "Syn5") %>% select(-Cell.line)
syn1<-nf2.ncats %>% filter(Cell.line == "Syn1") %>% select(-Cell.line)

colnames(syn5) <- c("NCGC.SID", "Cell.Type", "NF2.Status.SYN5", "AC50.SYN5", "CRC.SYN5", 
                    "Max.Resp.SYN5", "Log.AC50.uM.SYN5", "Sample.Name", "AUC.SYN5", "AUC.Fit.SYN5", "Gene.Symbol")
colnames(syn1) <- c("NCGC.SID", "Cell.Type", "NF2.Status.SYN1", "AC50.SYN1", "CRC.SYN1", 
                    "Max.Resp.SYN1", "Log.AC50.uM.SYN1", "Sample.Name", "AUC.SYN1", "AUC.Fit.SYN1", "Gene.Symbol")

meningioma <- full_join(syn5, syn1)

meningioma$AUC.Ratio <- (meningioma$AUC.SYN5/meningioma$AUC.SYN1)

meningioma.up <- full_join(meningioma, ncats_inchikey) %>% full_join(evotec) %>% filter(AUC.Ratio>=1.25) %>% 
  select(NCGC.SID, Hugo_Gene, AUC.Ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA")

nodes <- data.frame(c(unique(as.character(meningioma.up$Hugo_Gene)), unique(meningioma.up$NCGC.SID)))
write.table(nodes,"meningioma_up_nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

edges <- meningioma.up
write.table(edges,"meningioma_up_edges.txt", sep= "\t", row.names = FALSE, quote = FALSE)

meningioma.down <- full_join(meningioma, ncats_inchikey) %>% full_join(evotec) %>% filter(AUC.Ratio<=0.75) %>% 
  select(NCGC.SID, Hugo_Gene, AUC.Ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA")

nodes <- data.frame(c(unique(as.character(meningioma.down$Hugo_Gene)), unique(meningioma.down$NCGC.SID)))
write.table(nodes,"meningioma_down_nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

edges <- meningioma.down
write.table(edges,"meningioma_down_edges.txt", sep= "\t", row.names = FALSE, quote = FALSE)

