gff <- read.table("/home/jberry/Danforth/Bart_Lab/Shiny-cassava_atlas/data/tme7_200703_falcon_phase0_renamed.gff",sep = "\t",stringsAsFactors = F)
gff <- gff[gff$V3 == "gene",]
gff[1:5,1:9]
gff_df <- data.frame(
    "ID"=unlist(lapply(strsplit(gff$V9,"ID="),function(i) strsplit(i[2],";")[[1]][1])),
    "Size" = abs(gff$V5-gff$V4),
    "annot"=unlist(lapply(strsplit(gff$V9,"Note="),function(i) strsplit(i[2],";")[[1]][1])),
    "GO_maker"=unlist(lapply(strsplit(gff$V9,"Ontology_term="),function(i) strsplit(i[2],";")[[1]][1])),
    stringsAsFactors = F
)
head(gff_df)

# There are only 549 genes in the at list that match the gff, not worth it
#atGOs <- data.frame(data.table::fread("tme7_falcon_200703__arabidopsis_go_terms.csv",stringsAsFactors = F),stringsAsFactors = F)
#colnames(atGOs) <- c("ID","GO_at")
#atGOs$ID_g <- as.character(sapply(atGOs$ID,function(i) strsplit(i,"[.]")[[1]][1]))
#atGOs[atGOs$ID_g == "tme7p0_S01300",]
#test <- do.call(rbind,lapply(split(atGOs,atGOs$ID_g),function(i) aggregate(data = i,GO_at~ID_g,FUN=function(i) paste(i,collapse = ","))))
#rownames(test) <- NULL
#colnames(test)[1] <- "ID"
#gff_df <- plyr::join(gff_df,test, by="ID")
#head(gff_df)
#gff_df$GO_at
#gff_df$ID %in% test$ID
#sum(!is.na(gff_df$GO_at))

gcm <- read.csv("/home/jberry/Danforth/Bart_Lab/Shiny-cassava_atlas/data/gene_count_matrix.csv",stringsAsFactors = F)
gcm <- gcm[!is.na(unlist(lapply(strsplit(gcm$gene_id,"[|]"),function(i) i[2]))),]
gcm <- gcm[!duplicated(as.character(unlist(lapply(strsplit(gcm$gene_id,"[|]"),function(i) i[2])))),]
rownames(gcm) <- as.character(unlist(lapply(strsplit(gcm$gene_id,"[|]"),function(i) i[2])))
gcm <- gcm[,-1]
gcm <- gcm[rownames(gcm) %in% gff_df$ID,]

gcm[1:5,1:5]
dim(gcm)
dim(gff_df)
mean(gff_df$ID %in% rownames(gcm))
mean(rownames(gcm) %in% gff_df$ID)
rownames(gcm)[!(rownames(gcm) %in% gff_df$ID)]

"tme7p0_08G00495" %in% gff_df$ID

####################################################################################################
# Create DEseq2 Object and get FPKMs
####################################################################################################
metadata <- data.frame("sample_name"=colnames(gcm),
                       "tissue"=sapply(colnames(gcm),function(i) strsplit(i,"[.]")[[1]][2]),
                       stringsAsFactors = F,row.names = NULL)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = gcm,
                              colData = metadata,
                              design = ~ tissue)
mcols(dds)$basepairs <- gff_df$Size
rfpkm <- data.frame(fpkm(dds))
rfpkm$ID <- rownames(rfpkm)
rownames(rfpkm) <- NULL

rfpkm <- plyr::join(rfpkm,gff_df,by="ID")
rfpkm$annot <- gsub("\n"," ",rfpkm$annot)
rfpkm$annot <- as.character(unlist(lapply(strsplit(rfpkm$annot,"_Phase0"),function(i) i[1])))
rfpkm$GO_maker <- gsub("\n"," ",rfpkm$GO_maker)
rfpkm$GO_maker <- as.character(unlist(lapply(strsplit(rfpkm$GO_maker," "),function(i) i[1])))
rfpkm$GO_maker[is.na(rfpkm$GO_maker)] <- ""
colnames(rfpkm) <- c(as.character(sapply(c("FEC","Fibrous.root","Lateral.bud","Leaf","Mid.Vein","OES","Petiole","RAM","SAM","Stem"),function(i) paste0(i,c(0,1,2)))),"Storage.root0","Storage.root1","gene_id","size","annot","go")
write.table(rfpkm,"gf_data_p0_withreps.txt",row.names = F,quote = F,sep = "\t")

rfpkm_avg <- data.frame(cbind(sapply(c("FEC","Fibrous.root","Lateral.bud","Leaf","Mid.Vein","OES","Petiole","RAM","SAM","Stem","Storage.root"),function(i) rowMeans(rfpkm[,stringr::str_detect(colnames(rfpkm),i)])),rfpkm[,c("gene_id","size","annot","go")]),stringsAsFactors = F)
write.table(rfpkm_avg,"gf_data_p0.txt",row.names = F,quote = F,sep = "\t")
head(rfpkm_avg)
