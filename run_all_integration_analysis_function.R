input <- commandArgs(trailingOnly = TRUE)[1]
curdir <- "~/Labolatory/JST_Mouse_Human/Tabula/Integration"

setwd(curdir)
library(here)
set_here(curdir)
renv::restore()

library(tidyverse)
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)
source(paste(curdir, "integration_analysis_function.R", sep = "/"))

read_tsv("~/Labolatory/JST_Mouse_Human/Tabula/Integration/Hs_Mm_one2one_genelist.tsv",
         show_col_types = FALSE) -> Hs_Mm
input <- "Liver"
outdir <- paste(curdir, input, sep = "/")

# ヒト
readRDS(file=paste(outdir, "/TS_", input, ".rds", sep = "")) -> Hs_seurat
Hs_cell <- celltype_pick(Hs_seurat)
Hs_seurat %>%
  extract_one2one_cell(Hs_cell, "Gene_name") -> Hs_seurat.one2one.table

# マウス
readRDS(file=paste(outdir, "/TM_", input, ".rds", sep = "")) -> Mm_seurat
Mm_cell <- celltype_pick(Mm_seurat)
Mm_seurat %>%
  extract_one2one_cell(Mm_cell, "Mouse_gene_name") -> Mm_seurat.one2one.table

# 遺伝子選択
Hs_seurat.one2one.table %>%
  inner_join(Mm_seurat.one2one.table, by = "Gene_name") %>%
  dplyr::select("Gene_name") -> gene_list

# Seuratオブジェクトの作成
Hs_meta_col <- c("donor", "method", "cell_ontology_class", "free_annotation", "gender")
Mm_meta_col <- c("mouse.id", "FACS.selection", "cell_ontology_class", "subtissue", "mouse.sex")
create_one2one_object(Hs_seurat.one2one.table, Hs_seurat,
                      "Hs", Hs_cell, Hs_meta_col) -> Hs_seurat.one2one
create_one2one_object(Mm_seurat.one2one.table, Mm_seurat,
                      "Mm", Mm_cell, Mm_meta_col) -> Mm_seurat.one2one

# Variable featureの選択
genes.use = Reduce(intersect, list(VariableFeatures(Hs_seurat.one2one),
                                   VariableFeatures(Mm_seurat.one2one),
                                   rownames(Hs_seurat.one2one),
                                   rownames(Mm_seurat.one2one)))

# Seurat integration
Hs_Mm.list <- list(Hs = Hs_seurat.one2one, Mm = Mm_seurat.one2one)
features <- SelectIntegrationFeatures(object.list = Hs_Mm.list)
Hs_Mm.anchors <- FindIntegrationAnchors(object.list = Hs_Mm.list,
                                        anchor.features = features)

# this command creates an 'integrated' data assay
Hs_Mm.seurat.combined <- IntegrateData(anchorset = Hs_Mm.anchors)
DefaultAssay(Hs_Mm.seurat.combined) <- "integrated"
Hs_Mm.seurat.combined <- ScaleData(Hs_Mm.seurat.combined)
## Run PCA and UMAP
Hs_Mm.seurat.combined <- RunPCA(Hs_Mm.seurat.combined, do.print=FALSE)
ElbowPlot(object = Hs_Mm.seurat.combined, ndims = 50)

# UMAP＋可視化
npcs <- 30
Hs_Mm.seurat.combined <- FindNeighbors(Hs_Mm.seurat.combined,
                                       reduction = "pca", dims = 1:npcs)
Hs_Mm.seurat.combined <- FindClusters(Hs_Mm.seurat.combined,
                                      resolution = 0.5)
Hs_Mm.seurat.combined <- RunUMAP(Hs_Mm.seurat.combined, dims = 1:npcs)

DimPlot(Hs_Mm.seurat.combined, group.by = "cell_ontology_class",
        split.by = "SP", repel = T, label = T)  + NoLegend() -> integrated.cell_ontology.Hs_Mm
DimPlot(Hs_Mm.seurat.combined,
        split.by = "SP", group.by = "ident",
        repel = T, label = T)  + NoLegend() -> integrated.ident
DimPlot(Hs_Mm.seurat.combined, group.by = "SP",
        repel = T, label = T, label.box = T) + NoLegend() -> integrated.Hs_Mm
DimPlot(Hs_Mm.seurat.combined, group.by = "cell_ontology_class",
        repel = T, label = T) + NoLegend() -> integrated.cell_ontology
integrated.cell_ontology.Hs_Mm +
  integrated.ident +
  {integrated.Hs_Mm + integrated.cell_ontology} -> integrated.out
outfile.integrated <- paste(outdir, "/Figure/Hs_Mm_", input, ".seurat_integrate.celltype_pick.ident.png", sep = "")
ggsave(outfile.integrated, integrated.out, width = 12, height = 20)

# 統合データの保存
paste(outdir, "/RDS/", "Hs_Mm_", input, ".seurat.combined.rds", sep = "") %>%
  saveRDS(Hs_Mm.seurat.combined, file = .)

# 統合データの細胞クラスター情報を元のデータにくっつける＆統合データの細胞クラスターを利用した、マーカー遺伝子の抽出
# マウス
Mm_seurat.picked <- picked_gene_list(Mm_seurat, Mm_cell, "Mm")
Mm_seurat.picked.markers <- find_merker(Mm_seurat.picked, "Mm")
paste(outdir, "/RDS/Mm_", input, "picked.rds", sep = "") %>%
  saveRDS(Mm_seurat.picked, file = .)
# ヒト
Hs_seurat.picked <- picked_gene_list(Hs_seurat, Hs_cell, "Hs")
Hs_seurat.picked.markers <- find_merker(Hs_seurat.picked, "Hs")
paste(outdir, "/RDS/Hs_", input, "picked.rds", sep = "") %>%
  saveRDS(Hs_seurat.picked, file = .)

# ヒト・マウスの種特異的な遺伝子のリストを用意
# 系統特異的な遺伝子リストを用意
Hs_only <- read_tsv(paste(curdir, "Gene_name/Hs.Gene_name.species_specific.v2.txt", sep = "/"),
                    show_col_types = FALSE)
Mm_only <- read_tsv(paste(curdir, "/Gene_name/Mm.Gene_name.species_specific.v2.txt", sep = "/"),
                    show_col_types = FALSE)

# 対応関係表に含まれている遺伝子リスト
Hs_all.gene <- read_tsv(paste(curdir, "Gene_name/Hs.Gene_name.homology_type.txt", sep = "/"),
                        show_col_types = FALSE) %>%
                  select("Gene_name") %>%
                  distinct
Mm_all.gene <- read_tsv(paste(curdir, "/Gene_name/Mm.Gene_name.homology_type.txt", sep = "/"),
                        show_col_types = FALSE) %>%
                  select("Gene_name") %>%
                  distinct

# ヒトとマウスについて、マーカー遺伝子のうち種特異的なものを抽出
# マウス特異的マーカー遺伝子
paste(outdir, "/Marker_gene/Mm_", input, ".seurat_integrate.celltype_pick.Mm_only.markergene.tsv",
      sep = "") %>%
  ss_markers(Mm_seurat.picked.markers, Mm_only, 0.05, .)
# ヒトオーソログをもつマウスマーカー遺伝子
paste(outdir, "/Marker_gene/Mm_", input, ".seurat_integrate.celltype_pick.Mm_ortholog.markergene.tsv",
      sep = "") %>%
  ortholog_markers(Mm_seurat.picked.markers, Mm_only, Mm_all.gene, 0.05, .)

# ヒト特異的マーカー遺伝子
paste(outdir, "/Marker_gene/Hs_", input, ".seurat_integrate.celltype_pick.Hs_only.markergene.tsv",
      sep = "") %>%
  ss_markers(Hs_seurat.picked.markers, Hs_only, 0.05, .)
# マウスオーソログをもつヒトマーカー遺伝子
paste(outdir, "/Marker_gene/Hs_", input, ".seurat_integrate.celltype_pick.Hs_ortholog.markergene.tsv",
      sep = "") %>%
  ortholog_markers(Hs_seurat.picked.markers, Hs_only, Hs_all.gene, 0.05, .)

# 種特異的な遺伝子の発現量をプロット
# 統合データのメタ情報を抽出
Hs_Mm.seurat.combined@meta.data %>%
  rownames_to_column(var = "cell") ->
  Hs_Mm.seurat.meta.data

# マウス特異的な遺伝子の発現量合計を計算
sum_target_gene_data(Mm_seurat, Mm_only) -> Mm_seurat.Mm_only.sum
# ヒト特異的な遺伝子の発現量合計を計算
sum_target_gene_data(Hs_seurat, Hs_only) -> Hs_seurat.Hs_only.sum

# 種特異的な遺伝子の発現量合計を統合データのメタ情報と結合
merge_Hs_Mm_analyzed(Mm_seurat.Mm_only.sum, Hs_seurat.Hs_only.sum,
                     Hs_Mm.seurat.meta.data) -> seurat.only.sum.meta

# 系統特異的な遺伝子の合計をプロット
seurat.only.sum.meta %>%
    ggplot(aes(x = seurat_clusters, y = SUM, fill = seurat_clusters)) +
    geom_violin() +
    facet_wrap( ~ SP) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 5) +
    theme_classic() +
    NoLegend() -> integrated.sum

# マーカー遺伝子のうち、種特異的な遺伝子が占める割合をプロットする
# マウス
marker_gene_ratio(Mm_seurat.picked.markers, Mm_only) %>%
    mutate(SP = "Mm") -> Mm_seurat.picked.counts
# ヒト
marker_gene_ratio(Hs_seurat.picked.markers, Hs_only) %>%
    mutate(SP = "Hs") -> Hs_seurat.picked.counts

# Cluster数
Hs_seurat.picked.counts %>%
  bind_rows(Mm_seurat.picked.counts) %>%
  dplyr::select("cluster") %>%
  distinct %>%
  dim -> num.cluster

# ヒトとマウスの統合
Hs_seurat.picked.counts %>%
  bind_rows(Mm_seurat.picked.counts) %>%
  transform(cluster = factor(cluster, levels = 0:(num.cluster-1))) ->
  seurat.picked.count

# ヒトとマウスの種特異的なマーカー遺伝子数をプロット
seurat.picked.count %>%
  ggplot(aes(x=cluster, y = n.y, fill = cluster)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~ SP) +
    theme_classic() +
    NoLegend() -> integrated.counts
integrated.cell_ontology.Hs_Mm + integrated.ident +
  integrated.sum + integrated.counts -> integrated_figure.sum
paste(outdir,  "/Figure/Hs_Mm_", input, ".seurat_integrate.celltype_pick.only_gene_sum.png",
      sep = "") %>%
  ggsave(integrated_figure.sum, width = 12, height = 20)

# Ratioで計算する
# 全対応表に含まれている遺伝子の合計を細胞ごとに計算
# マウス
sum_target_gene_data(Mm_seurat, Mm_all.gene) -> Mm_seurat.all_sum
# ヒト
sum_target_gene_data(Hs_seurat, Hs_all.gene) -> Hs_seurat.all_sum

merge_Hs_Mm_analyzed(Mm_seurat.all_sum, Hs_seurat.all_sum,
                     Hs_Mm.seurat.meta.data) -> seurat.all.sum.meta

# 系統特異的な遺伝子の合計をプロット
seurat.only.sum.meta %>%
  inner_join(seurat.all.sum.meta, by = c("cell", "seurat_clusters", "SP")) %>%
    mutate(ratio = SUM.x / SUM.y) %>%
    ggplot(aes(x = seurat_clusters, y = ratio, fill = seurat_clusters)) +
      geom_violin() +
      facet_wrap( ~ SP) +
      geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.0005) +
      theme_classic() +
      NoLegend() -> integrated.gene.ratio

# ヒトとマウスのマーカー遺伝子にしめる種特異的な遺伝子の割合をプロット
seurat.picked.count %>%
  mutate(ratio = n.y / n.x) %>%
  ggplot(aes(x=cluster, y = ratio, fill = cluster)) +
    geom_bar(stat = "identity") +
    facet_wrap( ~ SP) +
    theme_classic() +
    NoLegend() -> integrated.count.ratio

integrated.cell_ontology.Hs_Mm + integrated.ident +
  integrated.gene.ratio + integrated.count.ratio -> integrated_figure.ratio

paste(outdir,  "/Figure/Hs_Mm_", input, ".seurat_integrate.celltype_pick.only_gene_ratio.png",
      sep = "") %>%
ggsave(integrated_figure.ratio, width = 12, height = 20)

