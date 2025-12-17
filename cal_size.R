################################
### Calibration size
################################
library(SingleCellExperiment)
library(scater)
library(scran)
library(VGAM)
library(igraph)
library(ontoProc)
library(MerfishData)
library(ExperimentHub)
library(kableExtra)
library(scConform)
library(BiocParallel)
library(tidyr)
library(ggplot2)

`%notin%` = Negate(`%in%`)

cl <- getOnto("cellOnto", "2023")

# Load data

spe_baysor <- MouseIleumPetukhov2021(
  segmentation = "baysor",
  use.images = FALSE, use.polygons = FALSE
)
# Load ontology
cl <- getOnto("cellOnto", "2023")

tags <- c(
  "CL:0009022", # Stromal
  "CL:0000236", # B cell
  "CL:0009080", # Tuft
  "CL:1000411", # Endothelial
  "CL:1000335", # Enterocyte
  "CL:1000326", # Goblet
  "CL:0002088", # ICC
  "CL:0009007", # Macrophage + DC
  "CL:1000343", # Paneth
  "CL:0000669", # Pericyte
  "CL:1000278", # Smooth Muscle
  "CL:0009017", # Stem + TA
  "CL:0000492", # T (CD4+)
  "CL:0000625", # T (CD8+)
  "CL:0017004" # Telocyte
)
opi <- graph_from_graphnel(onto_plot2(cl, tags))

## Delete instances from other ontologies
sel_ver <- V(opi)$name[c(grep("CARO", V(opi)$name), grep("BFO", V(opi)$name))]
opi1 <- opi - sel_ver

## Rename vertex to match annotations
V(opi1)$name[grep("CL:0000236", V(opi1)$name)] <- "B cell"
V(opi1)$name[grep("CL:1000411", V(opi1)$name)] <- "Endothelial"
V(opi1)$name[grep("CL:1000335", V(opi1)$name)] <- "Enterocyte"
V(opi1)$name[grep("CL:1000326", V(opi1)$name)] <- "Goblet"
V(opi1)$name[grep("CL:0002088", V(opi1)$name)] <- "ICC"
V(opi1)$name[grep("CL:0009007", V(opi1)$name)] <- "Macrophage + DC"
V(opi1)$name[grep("CL:1000343", V(opi1)$name)] <- "Paneth"
V(opi1)$name[grep("CL:0000669", V(opi1)$name)] <- "Pericyte"
V(opi1)$name[grep("CL:1000278", V(opi1)$name)] <- "Smooth Muscle"
V(opi1)$name[grep("CL:0009017", V(opi1)$name)] <- "Stem + TA"
V(opi1)$name[grep("CL:0009022", V(opi1)$name)] <- "Stromal"
V(opi1)$name[grep("CL:0000492", V(opi1)$name)] <- "T (CD4+)"
V(opi1)$name[grep("CL:0000625", V(opi1)$name)] <- "T (CD8+)"
V(opi1)$name[grep("CL:0017004", V(opi1)$name)] <- "Telocyte"
V(opi1)$name[grep("CL:0009080", V(opi1)$name)] <- "Tuft"

## Add the edge from connective tissue cell and telocyte and delete redundant
## nodes
opi1 <- add_edges(opi1, c("connective\ntissue cell\nCL:0002320", "Telocyte"))
gr <- as_graphnel(opi1)
opi2 <- opi1 - c(
  "somatic\ncell\nCL:0002371", "contractile\ncell\nCL:0000183",
  "native\ncell\nCL:0000003"
)

V(opi2)$name <- trimws(gsub("CL:.*|\\n", " ", V(opi2)$name))

gr1 <- as_graphnel(opi2)

## Plot the final ontology
attrs <- list(node = list(shape = "box", fontsize = 15, fixedsize = FALSE))

spe_baysor$cell_type <- spe_baysor$leiden_final
spe_baysor$cell_type[spe_baysor$cell_type %in% c(
  "B (Follicular, Circulating)",
  "B (Plasma)"
)] <- "B cell"
spe_baysor$cell_type[grep("Enterocyte", spe_baysor$cell_type)] <- "Enterocyte"
spe_baysor <- spe_baysor[, spe_baysor$cell_type %notin% c(
  "Removed",
  "Myenteric Plexus"
)]

spe_baysor <- addPerCellQCMetrics(spe_baysor)

# Randomly select 500 cells
set.seed(1703)
train <- sample(1:ncol(spe_baysor), 500)

# Training data
spe_train <- spe_baysor[, train]
spe_other <- spe_baysor[, -train]

# get HVGs
spe_train <- logNormCounts(spe_train)
v <- modelGeneVar(spe_train)
hvg <- getTopHVGs(v, n = 50)

rm(spe_baysor)

# Extract counts and convert data into a data.frame format
df_train <- as.data.frame(t(as.matrix(counts(spe_train[hvg, ]))))
df_other <- as.data.frame(t(as.matrix(counts(spe_other[hvg, ]))))
df_train$Y <- spe_train$cell_type
df_other$Y <- spe_other$cell_type

# Fit multinomial model
fit <- vglm(Y ~ .,
            family = multinomial(refLevel = "B cell"),
            data = df_train
)

# split data
set.seed(1237)
cal <- sample(1:nrow(df_other), 1000)
df_cal <- df_other[cal,]
df_test <- df_other[-cal,]

# See performance on the test data
p_test <- predict(fit, newdata=df_test, type="response")
pr_class <- apply(p_test, 1, function(row) colnames(p_test)[which.max(row)])

conf.pred.n <- function(B, n){
  cvg.crc <- cvg.conf <- rep(NA, B)
  for(i in 1:B){
    cal <- sample(1:nrow(df_other), n)
    df.cal <- df_other[cal,]
    df.test <- df_other[-cal,]
    p.cal <- predict(fit, newdata=df.cal, type="response")
    p.test <- predict(fit, newdata=df.test, type="response")
    sets_cf <- getPredictionSets(p.test, p.cal, df.cal$Y, follow_ontology = FALSE,
                                 resample = FALSE, return_sc = FALSE, 
                                 labels = unique(df_cal$Y))
    cvg.conf[i] <- mean(mapply(function(vec, lst) vec %in% lst, df.test$Y, sets_cf))
    
    sets_graphs <- getPredictionSets(p.test, p.cal, df.cal$Y, follow_ontology = TRUE,
                                     resample = FALSE, return_sc = FALSE, onto = opi2,
                                     labels = unique(df_cal$Y),
                                     BPPARAM = MulticoreParam(workers = 6))
    cvg.crc[i] <- mean(mapply(function(vec, lst) vec %in% lst, df.test$Y, sets_graphs))
    if(i%%10==0) cat(i)
  }
  return(list(cvg.conf=cvg.conf, cvg.crc=cvg.crc))
}

set.seed(1221)
n100 <- conf.pred.n(B=100, n=100)
save(n100, file="n100.RData")
n500 <- conf.pred.n(B=100, n=500)
save(n500, file="n500.RData")
n1000 <- conf.pred.n(B=100, n=1000)
save(n1000, file="n1000.RData")
