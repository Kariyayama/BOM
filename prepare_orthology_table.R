# 遺伝子関係表を作成
read_tsv("221110_Hs_Mm_gene_name_table.txt", show_col_types = FALSE) %>%
  dplyr::select(c("Gene_name", "Mouse_gene_name", "Mouse_homology_type")) %>%
  filter(Mouse_homology_type == "ortholog_one2one") %>%
  dplyr::select(c("Gene_name", "Mouse_gene_name")) %>%
  distinct -> Hs_Mm

Hs_Mm %>%
  dplyr::count(Gene_name) %>%
  dplyr::filter(n < 2) -> Hs_list

Hs_Mm %>%
  dplyr::count(Mouse_gene_name) %>%
  dplyr::filter(n < 2) -> Mm_list

Hs_Mm %>%
  inner_join(Hs_list, by = "Gene_name") %>%
  inner_join(Mm_list, by = "Mouse_gene_name") -> Hs_Mm

Hs_Mm %>%
  write.table("~/Labolatory/JST_Mouse_Human/Tabula/Integration/Hs_Mm_one2one_genelist.tsv",
              sep = "\t", row.names = F, quote = F)


