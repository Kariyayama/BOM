go_analysis <- function(marker.genes, taxon.id){
 if(length(marker.genes) > 0){
    rba_panther_enrich(marker.genes, taxon.id, BP, test_type = test_type,
                       cutoff = cutoff) -> marker.result
    write.table(marker.result$result, file = outfile,
                quote = F, row.names = F, sep="\t")
  }else{
    header %>%
      write(file = outfile)
  }
}

default_table <- function(type){
  tibble(number_in_list = 0,
        fold_enrichment = 0,
        fdr = 1,
        expected = 0,
        number_in_reference = 0,
        plus_minus = "",
        term.id = "",
        term.label = "",
        GO_short = "",
        group = "",
        type = type,
        fdr_log = 0)}

read_GOanalysis_result <- function(input, clstr, type){
  paste(outdir, "/PANTHER/Hs_shared-only/",
        "Hs_Mm_", input, ".seurat_integrate.cluster_", clstr, "_", type, ".tsv",
        sep = "") -> filename
  if(file.exists(filename)){
    filename %>%
      read_tsv(show_col_types = FALSE) %>%
      mutate(type = type) -> outtable
    if(dim(outtable)[1] > 0){
      outtable
    }else{
      default_table(type)
    }
  }else{
    default_table(type)
  }
}

integrate_species_sharing <- function(group1, group2, label, lngth=50){
  # SSとOrthologどっちにもhitしたID
  group1 %>%
      inner_join(group2, by = "term.id") %>%
      dplyr::select(term.id) -> both_GO.id

  # Both, SS, Orthologを抽出する
  group1 %>%
      inner_join(both_GO.id, by = "term.id") %>%
      mutate(group = "Both") -> both_GO
  anti_join(group1, both_GO.id, by = "term.id") %>%
      mutate(group = label) -> group1_only_GO

  rbind(both_GO, group1_only_GO) -> integrated
  if(dim(integrated)[1] > 0){
      integrated %>%
          mutate(fdr_log = -log10(fdr)) %>%
          mutate(GO_short = cut_GO_short(term.id, term.label, lngth))
  }else{
      default_table
  }
}

merged_shared_only <- function(input, clstr, n = 20){
    shared.ident <- tibble()
    read_GOanalysis_result(input, clstr, "shared") -> shared.ident
    only.ident <- tibble()
    read_GOanalysis_result(input, clstr, "Hs_only") -> only.ident

    integrated.shared <- tibble()
    if(dim(shared.ident)[1] > 0)
      integrate_species_sharing(shared.ident, only.ident, "Shared", 50) %>%
        group_by(group) %>%
        slice_head(n=n) %>%
        transform(GO_short = reorder(GO_short, fdr_log)) %>%
        as_tibble -> integrated.shared

    integrated.only <- tibble()
    if(dim(only.ident)[1] > 0)
      integrate_species_sharing(only.ident, shared.ident, "Hs_only", 50) %>%
        group_by(group) %>%
        slice_head(n=n) %>%
        transform(GO_short = reorder(GO_short, fdr_log)) %>%
        as_tibble -> integrated.only

    bind_rows(integrated.shared, integrated.only)
}
