# celltype pick
celltype_pick <- function(seurat_object){
  result <- c()
  for(cell_type in unique(seurat_object@meta.data$cell_ontology_class)){
    seurat_object@meta.data %>%
      dplyr::filter(cell_ontology_class == cell_type) %>%
      nrow -> n.cell_type

    if(n.cell_type > 200){
      seurat_object@meta.data %>%
        dplyr::filter(cell_ontology_class == cell_type) %>%
        slice_sample(n=200) %>%
        rownames %>%
        c(result, .) -> result
    }else{
      seurat_object@meta.data %>%
        dplyr::filter(cell_ontology_class == cell_type) %>%
        rownames %>%
        c(result, .) -> result
    }
  }
  return(result)
}

extract_one2one_cell <- function(x, cells.list, target){
  x@assays$RNA@counts %>%
    as.data.frame() %>%
    dplyr::select(all_of(cells.list)) %>%
    rownames_to_column(var= target) %>%
    inner_join(Hs_Mm, by = target) %>%
    dplyr::select(c(-Mouse_gene_name))
}
create_one2one_object <- function(seurat.one2one.table, seurat_object,
                                  sp, cells.list, meta_col){
  sp.ident <- rep(sp, length(cells.list))
  # メタデータの用意
  seurat_object@meta.data[cells.list,] %>%
    dplyr::select(meta_col) %>%
    rownames_to_column(var = "cell") %>%
    mutate(SP = sp.ident) ->
    sp.meta.data

  # Seuratオブジェクトに変換
  seurat.one2one.table %>%
    inner_join(gene_list, by = "Gene_name") %>%
    column_to_rownames(var = "Gene_name") %>%
    dplyr::select(all_of(cells.list)) %>%
    as.matrix %>%
    CreateSeuratObject() -> seurat.one2one

  # Seuratオブジェクトにメタ情報を追加
  seurat.one2one@meta.data %>%
    rownames_to_column(var = "cell") %>%
    inner_join(sp.meta.data, by = "cell") %>%
    column_to_rownames(var = "cell") -> seurat.one2one@meta.data

  # 正規化
  seurat.one2one <- SCTransform(seurat.one2one, variable.features.n = 2000)
  return(seurat.one2one)
}


# 統合データの細胞クラスター情報を元のデータにくっつける
picked_gene_list <- function(seurat_object, cells.list, sp){
  seurat_object@assays$RNA@counts[, cells.list] %>%
    CreateSeuratObject() -> seurat.picked
  Hs_Mm.seurat.combined@meta.data %>%
    filter(SP == sp) %>%
    select(seurat_clusters) ->
    seurat.picked[["integrated_seurat_cluster"]]

  seurat.picked <- SCTransform(seurat.picked,
                               variable.features.n = 2000)
  return(seurat.picked)
}

# 統合データの細胞クラスターを利用した、マーカー遺伝子の抽出
find_merker <- function(seurat.picked, sp){
  # Ident情報を細胞型情報に変換
  Hs_Mm.seurat.combined@meta.data %>%
    filter(SP == sp) %>%
    select(seurat_clusters) %>%
    t -> cell.ident
  cell.ident %>%
    as.character() -> cell_idents
  Hs_Mm.seurat.combined@meta.data %>%
    filter(SP == sp) %>%
    rownames_to_column(var = "cell") %>%
    select(cell) %>%
    unlist %>%
    as.vector -> cell.list
  names(cell_idents) <- cell.list
  cell_idents %>%
    unique %>%
    as.numeric %>%
    sort %>%
    as.character -> cell_levels

  cell.ident <- factor(cell_idents, levels = cell_levels)
  seurat.picked@active.ident <- cell.ident
  FindAllMarkers(seurat.picked, only.pos = TRUE, min.pct = 0.25,
                 logfc.threshold = 0.25) -> seurat.picked.markers
  outable_marker <- paste(outdir, "/Marker_gene/", sp, "_", input,
                            ".seurat_integrate.celltype_pick_markergene.tsv",
                            sep = "")

  seurat.picked.markers %>%
    filter(p_val_adj < 0.05) %>%
    write.table(file = outable_marker, sep = "\t",
                row.names = F, quote = F)
  return(seurat.picked.markers)
}

# 種特異的マーカー遺伝子
ss_markers <- function(picked.markers, only.gene_list,
                       threshold = 0.05, outfile){
  picked.markers %>%
    inner_join(only.gene_list, by = c("gene" = "Gene_name")) %>%
    filter(p_val_adj < threshold) %>%
    write.table(file = outfile, sep = "\t", row.names = F, quote = F)
}

# オーソログをもつマーカー遺伝子
ortholog_markers <- function(picked.markers, only.gene_list, all.gene,
                             threshold = 0.05, outfile){
  picked.markers %>%
    anti_join(only.gene_list, by = c("gene" = "Gene_name")) %>%
    inner_join(all.gene, by = c("gene" = "Gene_name")) %>%
    filter(p_val_adj < threshold) %>%
    write.table(file = outfile, sep = "\t", row.names = F, quote = F)
}

# 系統特異的な遺伝子の合計を細胞ごとに計算
sum_target_gene_data <- function(x, target.gene_list, analysis = sum)
  x@assays$RNA@data %>%
    as.data.frame %>%
    rownames_to_column("Gene") %>%
    inner_join(target.gene_list, by = c("Gene" = "Gene_name")) %>%
    column_to_rownames(var = "Gene") %>%
    apply(2, analysis)

# 遺伝子の発現量解析データを統合データのメタ情報と結合
merge_Hs_Mm_analyzed <- function(Mm, Hs, integrated.meta.data, key = "cell"){
  analyzed_expression <- c(Mm, Hs)
  data.frame(cell = names(analyzed_expression),
             SUM = analyzed_expression) %>%
    as_tibble %>%
    inner_join(integrated.meta.data, by = key)
}

# マーカー遺伝子のうち、種特異的な遺伝子が占める割合をプロットする
marker_gene_ratio <- function(picked.markers, target_gene.list,
                              threshold = 0.05){
  picked.markers %>%
    inner_join(target_gene.list, by = c("gene" = "Gene_name")) %>%
    filter(p_val_adj < threshold) %>%
    dplyr::count(cluster) -> picked.only_markers.counts
  picked.markers %>%
    filter(p_val_adj < threshold) %>%
    dplyr::count(cluster) %>%
    inner_join(picked.only_markers.counts,
               by = "cluster")
}

