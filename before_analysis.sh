dir=`pwd`
mkdir $1
cd $1
mkdir Figure
mkdir Marker_gene
mkdir -p PANTHER/Hs
mkdir -p PANTHER/Mm
mkdir -p PANTHER/Comparison/All
mkdir -p PANTHER/Comparison/Marker_gene
mkdir RDS
ln -s ~/Labolatory/JST_Mouse_Human/Tabula/Tabula_sapiens/TS_${1}.rds ./
ln -s ~/Labolatory/JST_Mouse_Human/Tabula/Tabula_muris/tabula-muris/TM_${1}.rds ./
cd $dir

# Rscript run_all_integration_analysis_function.R $1
# Rscript run_all_seurat_GO_analysis.R $1
