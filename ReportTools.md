# ReportingTools 

ReportingTools can be used to generated nice HTML output of DESeq2 Results with boxplot. The following code demonstrates how to add human gene symbols in the HTML output when working with real data

```R
library(tximport)
library(dplyr)
library("Biostrings")
library(readr)
library(org.Hs.eg.db)
library(DESeq2)


txi <- tximport(files, type="salmon", tx2gene=tx2gene, txOut=FALSE, varReduce = TRUE)

# Deseq2 on two way analysis

condition <- factor(c(rep("control",8),rep("knockout",7)))
ddsTwoWay <- DESeqDataSetFromTximport(txi, as.data.frame(condition), ~condition)
colData(ddsTwoWay)$condition <- factor(colData(ddsTwoWay)$condition, levels=c("control","knockout"))
ddsTwoWay <- DESeq(ddsTwoWay)
TwoWayResults <- results(ddsTwoWay)

# Add Gene Symbols and chromosome location information to the DESeq2 Results
TwoWayResults$symbol <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(TwoWayResults))
                               , "SYMBOL", "ENSEMBL", multiVals="first")
TwoWayResults$CHR <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(TwoWayResults))
                            , "CHR", "ENSEMBL", multiVals="first")
TwoWayResults$CHRSTART <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(TwoWayResults))
                                 , "CHRLOC", "ENSEMBL", multiVals="first")
TwoWayResults$CHREND <- mapIds(org.Hs.eg.db, sub("\\..*", "", rownames(TwoWayResults))
                               , "CHRLOCEND", "ENSEMBL", multiVals="first")

# For generate HTML Report using ReportingTool

library(ReportingTools)
DummyTwoWay <- ddsTwoWay

# no version numbers for ENSEMBL
rownames(DummyTwoWay) <- sub("\\..*", "", rownames(DummyTwoWay)) 

# Replace ENSEMBL with ENTREZID with org.Hs.eg.db
rownames(DummyTwoWay) <- mapIds(org.Hs.eg.db, rownames(DummyTwoWay),"ENTREZID","ENSEMBL")

# Remove data with rowname=NA since ReportTools cannot generate while having NA values
DummyTwoWayNOTNA <- DummyTwoWay[!is.na(rownames(DummyTwoWay)),] 

# Run DESeq results on only non-missing DESeq2 dataset
res2NotNA <- results(DummyTwoWayNOTNA)

# specify titles and output directory of the HTML report
des2Report <- HTMLReport(shortName ='TwoWay w Existing GeneName', title ='Differential Expression 
                         between Control versus Knockout with Existing GeneName', reportDirectory = "./reports")

# generate report using publish(DEseq2 results, reportformat, DESeq2 Dataset, ...)
# annotation.db gives us the Gene Symbol and ID of the given gene, it does not always work well
# if performing a three-way analysis with contrast, specify the contrast results in the results option, 
# and specify the three-way DESeq2 Dataset in the DESeq2 Dataset option
# when using contrast, I haven't figured out a way to include annotation.db within the report
publish(res2NotNA, des2Report,DummyTwoWayNOTNA, pvalueCutoff=0.05, make.plot=TRUE,
        annotation.db="org.Hs.eg.db", factor=colData(DummyTwoWayNOTNA)$condition, reportDir="./reports")
finish(des2Report)

# Same procedure can be used for rowname=NA observations but we cannot use annotation.db in publish().
```

