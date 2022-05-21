
dir.create("data")
# 5-way OIS meta-analysis: identifies meta-OIS DEGs from 5 meta analysis methods
load(file.path(getwd(), "../data/OIS_dataset_processed_20220405.rda"))
.metaGene(
  OIS_dataset_processed_20220405,
  .05,
  nperm = 1000,
  bedArgs = list(label.sig = c(0.01, 2)),
  seed = 76543217
)
dir(file.path(getwd()))[grep("_metaVennRes.RData", dir(file.path(getwd())))]
ois.tar = vennresults.temp$FiveMetaMethodsOverlap$a31

# 5-way quiescence meta-analysis: identifies meta-quiescence DEGs from 5 meta analysis methods
load(file.path(
  getwd(),
  "../data/OIS_di_quiescence_dataset_processed_20220405.rda"
))
.metaGene(
  OIS_di_quiescence_dataset_processed_20220405,
  .05,
  nperm = 1000,
  bedArgs = list(label.sig = c(0.01, 2)),
  seed = 76543217
)
qui.tar = vennresults.temp$FiveMetaMethodsOverlap$a31

# Meta-analytic unsupervised learning through fuzzy clustering on times series and bypass datasets.
load(file.path(
  getwd(),
  "../data/OIS_dataset_processed_timeSeries_20220405.rda"
))
load(file.path(
  getwd(),
  "../data/OIS_dataset_processed_bypass_20220405.rda"
))
load(file.path(
  getwd(),
  "../data/OIS_dataset_processed_quiescence_20220405.rda"
))
load(file.path(getwd(), "../data/HPA_2022_Jan9.rda"))
print(unique(unlist(strsplit(HPA[, c("Protein.class")], ', '))))
entrez.tar <-
  mapg(HPA[grep("Predicted secreted proteins", HPA[, 8]), "Gene"]) 
inx = na.omit(entrez.tar)

d.merged.ois.timeSeries <- OIS_dataset_processed_timeSeries_20220405
d.merged.ois.bypass <- OIS_dataset_processed_bypass_20220405
tmp = d.merged.ois.timeSeries$GSE59522
time.labels = c("0 hr",
                "2 hrs",
                "8 hrs",
                "24 hrs",
                "2 days",
                "4 days",
                "6 days",
                "8 days")
tmp = d.merged.ois.timeSeries$GSE144397oisRas3oxygen
time.labels = c("0 hr", "24 hrs", "48 hrs", "72 hrs", "94 hrs", "144 hrs")
tmp = d.merged.ois.timeSeries$GSE144397oisRaf3oxygen
time.labels = c("0 hr", "12 hrs", "24 hrs", "48 hrs", "72 hrs", "96 hrs")
tmp = d.merged.ois.bypass$GSE2487
time.labels = c("0 hr", "OIS", "OIS-bypass")
x = tmp$x
y = data.frame(tmp$y)
rownames(y) = colnames(x)
inx2 = intersect(rownames(x), inx)

# Meta-analytic unsupervised learning through fuzzy clustering on times series and bypass datasets. 
softClusteringCenterNum = 4 # The optimal parameter K needs to be determined by observing the graph and was manually set to the following : softClusteringCenterNum = 4 
tpg.res <-
  tpg(
    x[inx2,],
    y[, 1, drop = F],
    "soft clustering",
    softClusteringCenterNum = 4,
    seed = 1234,
    m = 2,
    time.labels = time.labels,
    labels.angle = 45,
    PCAOverlapPlot = F
  )

fuzzy.GSE59522 <- tpg.res # # The above function must be repeated for each dataset
#fuzzy.GSE144397oisRas3oxygen <- tpg.res
#fuzzy.GSE144397oisRaf3oxygen <- tpg.res
#fuzzy.GSE2487 <- tpg.res
# fuzzy.res = list(fuzzy.GSE59522=fuzzy.GSE59522, fuzzy.GSE144397oisRas3oxygen=fuzzy.GSE144397oisRas3oxygen, fuzzy.GSE144397oisRaf3oxygen=fuzzy.GSE144397oisRaf3oxygen, fuzzy.GSE2487=fuzzy.GSE2487)

load(file.path(getwd(), "../data/OIS_fuzzyRes_20220405.rda"))
fuzzy.res <- OIS_fuzzyRes_20220405
vec.list = list(up.vec.list = list(4, 2, 2, 4),
                down.vec.list = list(3, 3, 4, 1))
# vec.list=list(up.vec.list = list(4,2,2, c(4,6)),  down.vec.list = list(3,3,4, c(1,2)))

fuzzy.GSE59522 <- fuzzy.res$fuzzy.GSE59522
fuzzy.GSE144397oisRas3oxygen <-
  fuzzy.res$fuzzy.GSE144397oisRas3oxygen
fuzzy.GSE144397oisRaf3oxygen <-
  fuzzy.res$fuzzy.GSE144397oisRaf3oxygen
fuzzy.GSE2487 <- fuzzy.res$fuzzy.GSE2487
inx.updown <- list()
for (i in 1:length(vec.list)) {
  vec = vec.list[[i]]
  inx.updown[[i]] <-
    Reduce(
      intersect,
      list(
        fuzzy.GSE59522$normalized$genes[which(fuzzy.GSE59522$df$module %in% vec[[1]])] ,
        fuzzy.GSE144397oisRas3oxygen$normalized$genes[which(fuzzy.GSE144397oisRas3oxygen$df$module %in% vec[[2]])],
        fuzzy.GSE144397oisRaf3oxygen$normalized$genes[which(fuzzy.GSE144397oisRaf3oxygen$df$module %in% vec[[3]])],
        fuzzy.GSE2487$normalized$genes[which(fuzzy.GSE2487$df$module %in% vec[[4]])]
      )
    )
}
hugoUp = mapg(inx.updown[[1]])
hugoDown = mapg(inx.updown[[2]])

inxQuiescence = qui.tar
intersectQInxUp = intersect(inxQuiescence, c(hugoUp))
intersectQInxUp
intersectQInxDown = intersect(inxQuiescence, c(hugoDown))
intersectQInxDown
if (length(intersectQInxUp) != 0)
  hugoUp = setdiff(hugoUp, intersectQInxUp)
if (length(intersectQInxDown) != 0)
  hugoDown = setdiff(hugoDown, intersectQInxDown)
inxFinal.up = intersect(c(hugoUp),
                        OISres$ois.vennresults$FiveMetaMethodsOverlap$a31)
inxFinal.up
inxFinal.down = intersect(c(hugoDown),
                          OISres$ois.vennresults$FiveMetaMethodsOverlap$a31)
inxFinal.down
overlap4 = list(inxFinal.up = inxFinal.up, inxFinal.down = inxFinal.down)


# Identification of the core features through supervised learning from OISP
load(file.path(getwd(), "../data/OIS_ML_dataset_processed_20220405.rda"))
authors = c(
  "GSE46801 (Pawlikowski JS et al., 2013)",
  "GSE33710 (Benhamed M et al ., 2012)",
  "GSE75207 (Tordella L et al., 2016)",
  "GSE60652 (Takebayashi S et al., 2015)",
  "GSE54402 (Nelson DM et al., 2014)",
  "GSE19864 (Chicas A et al., 2010)",
  "GSE59522 (Young AR et al., 2009)",
  "GSE2487 (Collado M et al., 2005)",
  "GSE144397 (Martínez-Zamudio RI., al., 2020)",
  "GSE144397 (Martínez-Zamudio RI., al., 2020)"
)
devCohorts = c(0, rep(1, 6), 0, 1, 0)
res <-
  biPDSmlv2(
    OIS_ML_dataset_processed_20220405,
    CVal = "LOSOCV",
    GlobalOp = "epsgo",
    devCohorts = devCohorts,
    verbosePlotTable = F,
    alphas.seq = alphas.seq,
    authors = authors,
    path = "cv"
  )
res$coefDf


# OISP validation in TCGA cancer patient samples with high frequency of Ras/Raf/MEK mutations
load('~/DB/Firehose/pancanTCGA/TCGAentrezInsClass202111.RData')
TCGA.pancan@pheno$cancerType <-
  toTCGAAcronym(TCGA.pancan@pheno$X_primary_disease)

load(file.path(
  getwd(),
  "../data/OIS_TCGA_processed_RasRafMEK_20220405.rda"
))
RasRafMEK <- OIS_TCGA_processed_RasRafMEK_20220405
# Exclude silent mutation
tmpp <- RasRafMEK
tmpp <- lapply(tmpp, function(x) {
  x <-
    lapply(x, function(x) {
      if (length(inx <-
                 which(x$Variant_Classification == "Silent")) != 0)  {
        x = x[-inx,]
        x
      } else {
        x = x
        x
      }
    })
})
RasRafMEK <- tmpp
RasRafMEK.all.maf <-
  do.call("rbind", lapply(RasRafMEK, function(x) {
    if (dim(x$maf)[2] == 3)  {
      x$maf$AAChange = NA
      colnames(x$maf) =  c("Hugo_Symbol",
                           "Variant_Classification",
                           "Tumor_Sample_Barcode",
                           "AAChange")
      x$maf
    } else {
      colnames(x$maf) =  c("Hugo_Symbol",
                           "Variant_Classification",
                           "Tumor_Sample_Barcode",
                           "AAChange")
      x$maf
    }
  }))
RasRafMEK.all.maf$Tumor_Sample_Barcode <-
  substr(RasRafMEK.all.maf$Tumor_Sample_Barcode, 1, 15)
RasRafMEK.all.maf.collapse <- column.merge2(RasRafMEK.all.maf, 3)


show.line.paired = T
minInx = 12
onlypaired = T
excludeRASRAFMEK = T
tarGene.specified.vec <-
  c("1029", "4319", "1469", "2069", "6387", "2192")

indGplotList = indGplotList2 = list()

for (k in tarGene.specified.vec) {
  tarGene.specified <- k
  for (tar in names(RasRafMEK)) {
    inx <-
      grep(paste0(substr(
        unique(RasRafMEK[[tar]]$maf$Tumor_Sample_Barcode), 1, 12
      ), collapse = "|"), colnames(TCGA.pancan@x))
    inx2 = colnames(TCGA.pancan@x)[inx]
    
    if (excludeRASRAFMEK) {
      excludeInx = rownames(TCGA.pancan@pheno)[which(TCGA.pancan@pheno$cancerType %in% tar)]
      inx2 = setdiff(excludeInx, inx2)
    }
    
    if (onlypaired) {
      inx3 = inx2[duplicated(substr(unique(inx2), 1, 12))]
      grepInx = paste0(substr(unique(inx3), 1, 12), collapse = "|")
      if (grepInx != "") {
        inx2 = inx2[grep(grepInx, inx2)]
      } else {
        inx2 = NULL
      }
    }
    
    if (is.null(inx2))
      next
    expr <- TCGA.pancan@x[, inx2]
    
    if (dim(expr)[2] < minInx)
      next
    anno <-
      data.frame(
        SampleID = colnames(expr),
        ExperimentID = "",
        SampleType = TCGA.pancan@pheno[inx2, "sample_type"]
      )
    if (length(inx4 <-
               which(
                 !anno$SampleType %in% c("Solid Tissue Normal", "Primary Tumor")
               )) != 0) {
      anno <- anno[-inx4,]
      expr = expr[, -inx4]
    }
    
    tmp <-
      list(expr.list = list(tmp.expr = expr),
           anno.list = list(tmp.anno = anno))
    tmp$anno.list$tmp.anno$SampleType <-
      gsub("Solid Tissue Normal",
           "  Solid Tissue Normal",
           tmp$anno.list$tmp.anno$SampleType)
    tmp$anno.list$tmp.anno$SampleType <-
      gsub("Primary Tumor",
           " Primary Tumor",
           tmp$anno.list$tmp.anno$SampleType)
    
    indGplotList[tar] <-
      indGplot(
        tmp,
        tarGene.specified,
        entrezHugo = T,
        x.axis.label.angle = 90,
        ggarrange.ncol = 1,
        ggarrange.nrow = 1,
        combinedPlot = F,
        hide.ns = T,
        line = F,
        eb.batch.correct = F,
        legend.position = "none",
        show.line.paired = show.line.paired,
        subtitle = tar,
        col = .col("c2")(8)[c(5, 3, 7, 7)]
      )
    
    tar.meltSec.pre$Variant_Classification <-
      ifelse(
        tar.meltSec.pre$y %in% c(" Primary Tumor", "Metastatic"),
        RasRafMEK.all.maf.collapse[paste0(tar.meltSec.pre$pairdID, "-01"), c("Variant_Classification")],
        NA
      )
    tar.meltSec.pre$Oncogene <-
      ifelse(
        tar.meltSec.pre$y  %in% c(" Primary Tumor", "Metastatic"),
        RasRafMEK.all.maf.collapse[paste0(tar.meltSec.pre$pairdID, "-01"), "Hugo_Symbol"],
        NA
      )
    tar.meltSec.pre$Mutation <-
      paste0(tar.meltSec.pre$Oncogene,
             " (",
             tar.meltSec.pre$Variant_Classification,
             ")")
    tar.meltSec.pre$Mutation <-
      gsub("NA \\(NA\\)", NA, tar.meltSec.pre$Mutation)
    tar.meltSec.pre$Mutation <- factor(tar.meltSec.pre$Mutation)
    
    kk = mapg(k)
    ggg <-
      gboxplot(
        value ~ y,
        tar.meltSec.pre,
        1,
        3,
        add = "none",
        method = "wilcox",
        ref.group = "first",
        hide.ns = T,
        col = .col("c2")(8)[c(5, 3, 4, 4)],
        angle = 90,
        p.size = 4
      ) + labs(
        colg = '',
        title = tar,
        subtitle = "",
        x = "",
        y = substitute(paste(italic(kk)), list(kk = kk))
      )
    require(ggnewscale)
    indGplotList[[tar]] <- ggg + new_scale_color()
  }
  indGplotList2[[k]] <- indGplotList
}

indGplotList2[[1]]
do.call(ggpubr::ggarrange, c(
  do.call("c", indGplotList2),
  list(
    ncol = 5,
    nrow = 2,
    legend = "none"
  )
))  
