library(synapseClient)
library(dplyr)
library(igraph)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/map_pNF_to_evotec.R"

pNF.ncats<-read.table(synGet("syn5522642")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"")
write.table(pNF.ncats %>% dplyr::dplyr::select(NCGC.SID), "ncats-pNF-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

##take these to pubchem id exchange, uploaded to synapse, read in
ncats_inchikey<-read.table(synGet("syn8690496")@filePath, header = FALSE, sep = "\t", comment.char = "")
colnames(ncats_inchikey) <- c("NCGC.SID", "inchikey")

evotec<-read.table(synGet("syn8685703")@filePath, header = TRUE, sep = "\t", comment.char = "") %>% dplyr::select(inchikey, Hugo_Gene, Structure_ID)

pNF.ncats1<-read.table(synGet("syn5522642")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF02.3_2l_LAC50` = LAC50, `ipNF02.3_2l_MAXR` = MAXR, `ipNF02.3_2l_FAUC` = FAUC)

pNF.ncats2<-read.table(synGet("syn5522643")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF02.8_LAC50` = LAC50, `ipNF02.8_MAXR` = MAXR, `ipNF02.8_FAUC` = FAUC)

pNF.ncats3<-read.table(synGet("syn5522644")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF05.5_MC_LAC50` = LAC50, `ipNF05.5_MC_MAXR` = MAXR, `ipNF05.5_MC_FAUC` = FAUC)

pNF.ncats4<-read.table(synGet("syn5522645")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF05.5_SC_LAC50` = LAC50, `ipNF05.5_SC_MAXR` = MAXR, `ipNF05.5_SC_FAUC` = FAUC)

pNF.ncats5<-read.table(synGet("syn5522646")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF06.2A_LAC50` = LAC50, `ipNF06.2A_MAXR` = MAXR, `ipNF06.2A_FAUC` = FAUC)

pNF.ncats6<-read.table(synGet("syn5522648")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF95.11C_LAC50` = LAC50, `ipNF95.11C_MAXR` = MAXR, `ipNF95.11C_FAUC` = FAUC)

pNF.ncats7<-read.table(synGet("syn5522647")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(SID, LAC50,MAXR,FAUC) %>% 
  rename(`ipNF95.11b_C_T_LAC50` = LAC50, `ipNF95.11b_C_T_MAXR` = MAXR, `ipNF95.11b_C_T_FAUC` = FAUC, NCGC.SID = SID)

pNF.ncats8<-read.table(synGet("syn5522649")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC, name, target) %>% 
  rename(`ipNF95.6_LAC50` = LAC50, `ipNF95.6_MAXR` = MAXR, `ipNF95.6_FAUC` = FAUC)

#pNF.ncats9<-read.table(synGet("syn8556314")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"") %>% 
#  dplyr::select(NCGC.SID, LAC50,MAXR,FAUC) %>% 
#  rename(`ipNF95.11bC_LAC50` = LAC50, `ipNF95.11bC_MAXR` = MAXR, `ipNF95.11bC_FAUC` = FAUC)

pNF.all <- full_join(pNF.ncats1, pNF.ncats2) %>% full_join(pNF.ncats3) %>% full_join(pNF.ncats4) %>% full_join(pNF.ncats5) %>% 
  full_join(pNF.ncats6) %>% full_join(pNF.ncats7) %>% full_join(pNF.ncats8) 

pNF.LAC50 <- dplyr::select(pNF.all, NCGC.SID, target, ipNF02.3_2l_LAC50, ipNF02.8_LAC50, ipNF05.5_MC_LAC50, ipNF05.5_SC_LAC50,
                           ipNF06.2A_LAC50, ipNF95.11C_LAC50, ipNF95.11b_C_T_LAC50, ipNF95.6_LAC50)

pNF.MAXR <- dplyr::select(pNF.all, NCGC.SID, target, ipNF02.3_2l_MAXR, ipNF02.8_MAXR, ipNF05.5_MC_MAXR, ipNF05.5_SC_MAXR,
                          ipNF06.2A_MAXR, ipNF95.11C_MAXR, ipNF95.11b_C_T_MAXR, ipNF95.6_MAXR)

pNF.FAUC <- dplyr::select(pNF.all, NCGC.SID, target, ipNF02.3_2l_FAUC, ipNF02.8_FAUC, ipNF05.5_MC_FAUC, ipNF05.5_SC_FAUC,
                          ipNF06.2A_FAUC, ipNF95.11C_FAUC, ipNF95.11b_C_T_FAUC, ipNF95.6_FAUC)


pNF.FAUC.1 <- dplyr::select(pNF.FAUC, target, NCGC.SID, ipNF05.5_SC_FAUC, ipNF02.3_2l_FAUC)
pNF.FAUC.1$ratio <- (pNF.FAUC.1$ipNF02.3_2l_FAUC/pNF.FAUC.1$ipNF05.5_SC_FAUC)
pNF.FAUC.1 <- full_join(pNF.FAUC.1, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio>=1.5) %>% 
  select(NCGC.SID, Hugo_Gene, ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA")

nodes <- data.frame(c(unique(as.character(pNF.FAUC.1$Hugo_Gene)), unique(pNF.FAUC.1$NCGC.SID)))
color<- c(rep("0", 412), rep("1", 19))
nodes <- cbind(nodes,color)
write.table(nodes,"nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

edges <- pNF.FAUC.1
write.table(edges,"edges.txt", sep= "\t", row.names = FALSE, quote = FALSE)

write.table(ncats_join, "NTAP_ipNF95.11bC_MIPE_qHTS_evotec_targets.txt", sep = "\t")
synStore(File("NTAP_ipNF95.11bC_MIPE_qHTS_evotec_targets.txt", 
              parentId="syn8690652"), used = c("syn5522642","syn8690496", "syn8685703", "syn8556314"), executed = this.file)

pNF.LAC50.1 <- dplyr::select(pNF.LAC50, target, NCGC.SID, ipNF05.5_SC_LAC50, ipNF02.3_2l_LAC50)
pNF.LAC50.1$ratio <- (pNF.LAC50.1$ipNF02.3_2l_LAC50/pNF.LAC50.1$ipNF05.5_SC_LAC50)
pNF.LAC50.1 <- full_join(pNF.LAC50.1, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio>=1.5) %>% 
  select(NCGC.SID, Hugo_Gene, ratio) %>% filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00347946-01"
                                                & NCGC.SID != "NCGC00347951-01" & NCGC.SID != "NCGC00347957-01")

nodes <- data.frame(c(unique(as.character(pNF.LAC50.1$Hugo_Gene)), unique(pNF.LAC50.1$NCGC.SID)))
color<- c(rep("0", 412), rep("1", 19))
nodes <- cbind(nodes,color)
write.table(nodes,"nodes.txt", sep= "\t", row.names = FALSE, quote = FALSE)

net <- pNF.LAC50.1
write.table(net,"net.txt", sep= "\t", row.names = FALSE, quote = FALSE)

