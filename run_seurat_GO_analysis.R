# 解析用の関数
filter_cluster <- function(x, ident){
  x %>%
    filter(cluster == ident) %>%
    select(gene) %>% unlist %>% as.character()
}

# GO解析の結果の出力ファイル名を作成
createGOanalysis_outfile <- function(sp, ident,
                                     testtype = "FISHER", target){
  paste(outdir, "/PANTHER/", sp, "/", sp, "_PANTHER_ident_",
        ident, "_", test_type, "_BP_", target, "_gene.tsv", sep="")
}

# GO解析
run_GOanalysis_for_a_ident <- function(marker, sp, label, ident, taxon.id){
  marker %>%
    filter_cluster(., ident) -> marker.genes
  # Outfile
  outfile <- createGOanalysis_outfile(sp, ident,
                                       test_type, label)
  # GO analysis
  if(length(marker.genes) > 0){
    rba_panther_enrich(marker.genes, taxon.id, BP, test_type = test_type,
                       cutoff = cutoff) -> marker.result
    write.table(marker.result$result, file = outfile,
                quote = F, row.names = F, sep="\t")
  }else{
    "number_in_list\tfold_enrichment\tfdr\texpected \tnumber_in_reference\tpValue\tplus_minus\tterm.id\tterm.label" %>%
      write(file = outfile)
  }
}

# マーカー遺伝子情報を読み込み
read_marker_gene_list <- function(sp, target = "_"){
    paste(outdir, "/Marker_gene/", sp, "_", input,
        ".seurat_integrate.celltype_pick",
        target, "markergene.tsv", sep = "") %>%
    read_tsv(show_col_types = FALSE)
}

# GO解析を特定の種全体に行う
run_GOanalysis_for_all_ident <- function(sp, taxon.id, cluster_list){
  # マーカー検出結果の取得
  paste(".", sp, "_only.", sep = "") %>%
      read_marker_gene_list(sp, .) -> only.marker
  paste(".", sp, "_ortholog.", sep = "") %>%
      read_marker_gene_list(sp, .)-> ortholog.marker

  for(cluster.id in cluster_list){
    print(cluster.id)
    # SS genes
    run_GOanalysis_for_a_ident(only.marker,
                               sp, "SS", cluster.id, taxon.id)
    # Ortholog genes
    run_GOanalysis_for_a_ident(ortholog.marker,
                               sp, "ortholog", cluster.id, taxon.id)
  }
}

make_GO_short <- function(acc, short){
    result <- c()
    for(i in 1:length(acc)){
          if(nchar(unlist(short[i])) > 20){
              result <- c(result,
                          paste(substr(short[i], 1, 17),
                                "... (", acc[i], ")",
                                sep=""))
          }else{
              result <- c(result,
                          paste(short[i], " (", acc[i], ")",
                                sep=""))
          }
    }
    return(result)
}

read_data <- function(sp, cluster.id, target = "all"){
  paste(outdir, "/PANTHER/", sp, "/", sp, "_PANTHER_ident_",
        cluster.id, "_", test_type, "_BP_", target, "_gene.tsv", sep="") %>%
  read_tsv(show_col_types = FALSE) %>%
    filter(plus_minus == "+")
}

# データの加工、NAを1に変換する
na_to_1 <- function(x) as.numeric(replace(x, is.na(x), 1))

# GO解析の結果をplot
plot_gocompare <- function(out_table){
    out_table %>%
      mutate(GO_short_acc = make_GO_short(term.id, term.label)) %>%
      ggplot() +
        theme_classic() +
        theme(axis.text.x = plot_para)
}

# 種間比較
merge_different_species <- function(hs, mm){
    hs %>%
      full_join(mm, by = c("term.id", "term.label")) %>%
      dplyr::rename(human = fdr.x, mouse = fdr.y) %>%
      mutate(human = na_to_1(human), mouse = na_to_1(mouse)) %>%
      dplyr::select(c("term.id", "term.label", "human", "mouse")) %>%
      pivot_longer(cols=-c(term.id, term.label),
                   names_to = "SP", values_to = "FDR") %>%
      mutate(fdr_log = -log10(FDR))
}

# homology_typeごと
merge_different_homology_type <- function(SS, ortholog){
    SS %>%
        full_join(ortholog, by = c("term.id", "term.label")) %>%
        dplyr::rename(SS = fdr.x, ortholog = fdr.y) %>%
        mutate(SS = na_to_1(SS), ortholog = na_to_1(ortholog)) %>%
        dplyr::select(c("term.id", "term.label", "SS", "ortholog")) %>%
        pivot_longer(cols=-c(term.id, term.label),
                     names_to = "ortholog_type", values_to = "FDR") %>%
        mutate(fdr_log = -log10(FDR))
}

plot_SS_ortholog_compare <- function(sp){
  result <- data.frame()
  read_marker_gene_list(sp, "_") %>%
    select(cluster) %>%
    distinct %>% unlist %>%
    as.character -> cluster_list

  for(ident in cluster_list){
    print(ident)
    read_data(sp, ident, "SS") -> all.ident
    read_data(sp, ident, "ortholog") -> ortholog.ident
    outfile <- paste(outdir, "/PANTHER/Comparison/All/", sp,
                     "compare_ident_", ident, "_SS-ortholog.png", sep= "")
    merge_different_homology_type(all.ident, ortholog.ident) -> out_table
    result <- out_table %>%
      mutate(cluster.id = ident) %>%
      bind_rows(result, .)

    print(out_table)
    if(dim(out_table)[1] > 0){
      plot_gocompare(out_table) +
        geom_bar(aes(x=GO_short_acc, y = fdr_log, fill = ortholog_type),
                 stat = "identity", position = "dodge") +
        geom_hline(aes(yintercept = threshold), linetype = "dashed") -> p1

      ggsave(outfile, p1, width = width(dim(out_table)[1]), height = height,
             limitsize = FALSE)
    }
  }
  return(result)
}
