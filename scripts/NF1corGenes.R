x<-read.table(synGet("syn9772360")@filePath, header = TRUE)
up <- as.data.frame(filter(x, cors > 0.4)$gene_id)
down <- as.data.frame(filter(x, cors < -0.4)$gene_id)

write.table(up, "up.txt", col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
write.table(down, "down.txt", col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t")
