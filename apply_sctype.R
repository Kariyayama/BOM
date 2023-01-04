extract_cellmarker <- function(cl, meta.data){
    meta.data %>%
        filter(integrated_seurat_cluster == cl) %>%
        rownames -> es.max.list
    es.max[, es.max.list] %>%
        rowSums %>%
        sort(decreasing = !0) -> es.max.cl

    data.frame(cluster = cl,
               type = names(es.max.cl),
               scores = es.max.cl,
               ncells = sum(clusters==cl)
               ) %>%
        head(10)
}
