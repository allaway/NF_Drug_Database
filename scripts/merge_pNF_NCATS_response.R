library(synapseClient)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/map_pNF_to_evotec.R"

pNF.ncats<-read.table(synGet("syn5522642")@filePath, header = TRUE, sep = ",", comment.char = "", quote = "\"")
write.table(pNF.ncats %>% dplyr::select(NCGC.SID), "ncats-pNF-ID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
up <- full_join(pNF.FAUC.1, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio>2) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00346937-01")
net.up <- up %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.up, "netup.txt", row.names = FALSE, quote = FALSE)

down <- full_join(pNF.FAUC.1, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio<0.5) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00346937-01")
net.down <- down %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.down, "netdown.txt", row.names = FALSE, quote = FALSE)

pNF.FAUC.2 <- dplyr::select(pNF.FAUC, target, NCGC.SID, ipNF02.3_2l_FAUC, ipNF95.11b_C_T_FAUC)
pNF.FAUC.2$ratio <- (pNF.FAUC.2$ipNF02.3_2l_FAUC/pNF.FAUC.2$ipNF95.11b_C_T_FAUC)
up <- full_join(pNF.FAUC.2, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio>1.2) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00347952-01" &
           NCGC.SID != "NCGC00346964-01")
net.up <- up %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.up, "netup.txt", row.names = FALSE, quote = FALSE)

down <- full_join(pNF.FAUC.2, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio<0.6) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00346937-01")
net.down <- down %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.down, "netdown.txt", row.names = FALSE, quote = FALSE)

pNF.FAUC.3 <- dplyr::select(pNF.FAUC, target, NCGC.SID, ipNF02.8_FAUC, ipNF95.6_FAUC)
pNF.FAUC.3$ratio <- (pNF.FAUC.3$ipNF02.8_FAUC/pNF.FAUC.3$ipNF95.6_FAUC)
up <- full_join(pNF.FAUC.3, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio>1.2) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00346971-01" &
           NCGC.SID != "NCGC00346964-01" & NCGC.SID != "NCGC00347958-01" & NCGC.SID != "NCGC00347303-01" & NCGC.SID != "NCGC00346942-01")
net.up <- up %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.up, "netup.txt", row.names = FALSE, quote = FALSE)

down <- full_join(pNF.FAUC.3, ncats_inchikey) %>% full_join(evotec) %>% filter(ratio<0.6) %>% 
  filter(Hugo_Gene != "NA" & Hugo_Gene != "" & NCGC.SID != "NA" & NCGC.SID != "NCGC00347952-01" & NCGC.SID != "NCGC00346937-01" &
           NCGC.SID != "NCGC00347958-01" & NCGC.SID != "NCGC00347278-01")
net.down <- down %>% select(NCGC.SID, Hugo_Gene, ratio)
write.table(net.down, "netdown.txt", row.names = FALSE, quote = FALSE)

