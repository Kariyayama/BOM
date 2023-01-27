extract_Hs_marker_gene <- function(input, clstr,
                                   Hs_all.marker.ident, Mm_all.marker.ident){
  # 3.4 表を結合
  homology %>%
      inner_join(Hs_all.marker.ident,
                 by = c("Gene_name" = "gene")) %>%
      inner_join(Mm_all.marker.ident,
                 by = c("Mouse_gene_name" = "gene")) %>%
      select(c(1:13)) %>%
      mutate(share = "shared") -> Hs_Mm.conserved_marker.ident

  homology %>%
      inner_join(Hs_all.marker.ident,
                 by = c("Gene_name" = "gene")) %>%
      anti_join(Hs_Mm.conserved_marker.ident,
                by = c("Gene_stable_ID", "Gene_name")) %>%
      mutate(share = "Hs_only") -> Hs_all_uniq_marker.ident

  Hs_all.marker.ident %>%
      anti_join(Hs_Mm.conserved_marker.ident,
                by = c("gene" = "Gene_name")) %>%
      anti_join(Hs_all_uniq_marker.ident,
                by = c("gene" = "Gene_name")) %>%
      mutate(share = "No_homology_info") -> Hs_all_noorthology_marker.ident

  # 共通遺伝子をデータフレームに保存
  bind_rows(Hs_Mm.conserved_marker.ident,
            Hs_all_uniq_marker.ident,
            Hs_all_noorthology_marker.ident) ->  Hs_all.marker.ident.homology
  paste(input, "/Marker_gene/Hs_Mm_", input,
        ".seurat_integrate.cluster_", clstr, ".homology.tsv", sep = "") %>%
    write.table(Hs_all.marker.ident.homology, file = ., quote = F, row.names = F, sep = "\t")
  return(Hs_all.marker.ident.homology)
}

# マーカー遺伝子を対応関係ごとに抽出
annotate_Hs_marker_to_homology <- function(input){
  human_table <- tibble()
  # 表の読み込み
  paste(input, "/Marker_gene/Hs_", input ,
        ".seurat_integrate.celltype_pick_markergene.tsv", sep = "") %>%
    read_tsv(show_col_types = FALSE) -> Hs_all.marker
  paste(input, "/Marker_gene/Mm_", input,
        ".seurat_integrate.celltype_pick_markergene.tsv", sep = "") %>%
    read_tsv(show_col_types = FALSE) -> Mm_all.marker
  Hs_all.marker %>%
    select(cluster) %>%
    distinct() %>%
    unlist -> Hs_idents

  # 1. クラスター選択
  for(clstr in Hs_idents){
    # 2. マーカー遺伝子抽出
    Hs_all.marker %>%
      filter(cluster == clstr) -> Hs_all.marker.ident
    Mm_all.marker %>%
      filter(cluster == clstr) %>%
      select("gene") -> Mm_all.marker.ident

    extract_Hs_marker_gene(input, clstr,
                           Hs_all.marker.ident,
                           Mm_all.marker.ident) -> Hs_all.marker.ident.homology
    Hs_all.marker.ident.homology %>%
        select("Gene_name", "gene", "Homology_type", "share") %>%
        distinct%>%
        group_by(share) %>%
        dplyr::count(Homology_type) %>%
        mutate(type = paste(share, Homology_type, sep = "_")) %>%
        mutate(cluster = clstr) %>%
        bind_rows(human_table, .) -> human_table
  }
  return(human_table)
}

# クラスターのラベルをSeuratオブジェクトから取り出す
extract_cluster_label <- function(input, sp){
  paste(input, "/RDS/", sp, "_", input, "picked.rds", sep = "") %>%
    readRDS() -> seurat.picked
  seurat.picked@meta.data %>%
    rownames_to_column(var = "cell")  %>%
    select(integrated_seurat_cluster, customclassif) %>%
    distinct %>%
    arrange(integrated_seurat_cluster)
}
