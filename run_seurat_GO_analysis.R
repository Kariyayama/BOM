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

classify_GO_term <- function(ident, sp = "Hs"){
  read_data(sp, ident, "SS") %>%
      mutate(type = "SS") -> only.ident
  read_data(sp, ident, "ortholog") %>%
      mutate(type = "ortholog") -> ortholog.ident

  # SSとOrthologどっちにもhitしたID
  only.ident %>%
      inner_join(ortholog.ident, by = "term.id") %>%
      dplyr::select(term.id)-> both_GO.id

  # Both, SS, Orthologを抽出する
  rbind(only.ident, ortholog.ident) %>%
      inner_join(both_GO.id, by = "term.id") %>%
      mutate(group = "Both") -> both_GO
  anti_join(only.ident, both_GO.id, by = "term.id") %>%
      mutate(group = "SS") -> SS_only_GO
  anti_join(ortholog.ident, both_GO.id, by = "term.id") %>%
      mutate(group = "ortholog") -> ortholog_only_GO

  rbind(both_GO, SS_only_GO, ortholog_only_GO) -> integrated
  if(dim(integrated)[1] > 0){
      integrated %>%
          mutate(fdr_log = -log10(fdr)) %>%
          mutate(GO_short = cut_GO_short(term.id, term.label, 50)) %>%
          mutate(SP = sp)
  }else{
      data.frame()
  }
}

plot_GO_separate <- function(x){
    x %>%
        ggplot() +
        theme_classic() +
        theme(axis.text.x = plot_para,
              axis.title.y = element_blank(),
              legend.position = "none") +
        facet_grid(group ~ ., scales = "free", space = "free") +
        geom_bar(aes(x=fdr_log, y = reorder(GO_short, fdr_log), fill = type),
                 stat = "identity", position = "dodge") +
        ggtitle(unique(x$SP))
}

create_plot <- function(x){
    if(dim(x)[1] > 0){
        x %>%
            plot_GO_separate
    }else{
        ggplot(x)
  }
}

extract_term.num <- function(x){
    if(dim(x)[1] > 0){
        x %>%
            dplyr::select("term.id") %>%
            distinct %>%
            dplyr::count() %>%
            unlist %>%
            return
    }else{
        return(0)
    }
}

decide_height <- function(x, y){
    extract_term.num(x) -> height.x
    extract_term.num(y) -> height.y
    if(height.x >= height.y) height <- 5 + height.x * 0.1
    if(height.y > height.x) height <- 5 + height.y * 0.1
    if(height > 49){
        return(49)
    }else{
        return(height)
    }
}

# GO解析の結果をplot
plot_gocompare <- function(out_table){
    out_table %>%
        mutate(GO_short_acc = make_GO_short(term.id, term.label)) %>%
        ggplot() +
            theme_classic() +
            theme(axis.text.x = plot_para)
}

cut_GO_short <- function(acc, short, length = 40){
    result <- c()
    for(i in 1:length(acc)){
          if(nchar(unlist(short[i])) > length){
              result <- c(result,
                          paste(substr(short[i], 1, (length - 3)),
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
