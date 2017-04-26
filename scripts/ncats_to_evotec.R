library(synapseClient)
library(dplyr)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF_Drug_Database/master/scripts/ncats_to_evotec.R"

ncats<-read.table(synGet("syn8314523")@filePath, header = TRUE, sep = "\t", comment.char = "")
evotec<-read.table(synGet("syn8685703")@filePath, header = TRUE, sep = "\t", comment.char = "")

write.table(ncats %>% dplyr::select(Sample.ID), "ncatsID.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

##take these to pubchem id exchange, uploaded to synapse, read in
ncats_inchikey<-read.table(synGet("syn8685814")@filePath, header = FALSE, sep = "\t", comment.char = "")
colnames(ncats_inchikey) <- c("Sample.ID", "inchikey")

ncats_join<-left_join(ncats,distinct(ncats_inchikey))
ncats_join<-left_join(ncats_join,distinct(dplyr::select(evotec, inchikey, Hugo_Gene)))
ncats_join<-distinct(filter(ncats_join, Hugo_Gene != ""))

write.table(ncats_join, "NF2_ncats_with_evotec_targets.txt", sep = "\t")
synStore(File("NF2_ncats_with_evotec_targets.txt", parentId = "syn8682571"), executed = this.file, used = c("syn8685814", "syn8314523", "syn8685703", "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"))
