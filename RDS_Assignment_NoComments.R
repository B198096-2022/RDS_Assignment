

#Acquire data, check for zeros, log transform 

data_all <- read.table(file= "data_all.csv",header=T,sep=",")

gene_annotation <- read.table(file= "gene_annotation.csv",header=T,sep=",")

sample_annotation <- read.table(file= "sample_annotation.csv",header=T,sep=",")

data_all_df <- as.data.frame(data_all[1:200, 2:13])

zero_test <- function(df, find=0) {
  all_data_matrix <- as.matrix(df)
  count <- 0
  zero_positions <- c()
  for (i in all_data_matrix){
    count = count + 1
    if (i == find){
      print("Zero value found at position:")
      print(count)
      zero_positions <- c(zero_positions, count)
    }
  }
  return(zero_positions)
}


mimimum_test <- function(df, min=1) {
  all_data_matrix <- as.matrix(df)
  count <- 0
  min_positions <- c()
  for (i in all_data_matrix){
    count = count + 1
    if (i < min){
      print("Filter value:")
      print(i)
      print("Found at position:")
      print(count)
      min_positions <- c(min_positions, count)
    }
  }
  return(min_positions)
}



change_values <- function(df, find=0, change=1) {
  all_data_matrix <- as.matrix(df)
  find_matrix <- all_data_matrix <= find
  all_data_matrix[find_matrix] <- change
  all_data_df <- as.data.frame(all_data_matrix)
  return(all_data_df)
}


zeros <- zero_test(data_all_df)
filter_response = 'n'

if (is.null(zeros) == FALSE){
  print("Zero values found in data")
  filter_response <- readline(prompt="Filter Data to change Zeros to Ones? (y/n): ")
}


if (filter_response == "y"|| filter_response == "yes") {
  data_all_df <- change_values(all_data_df)
}


data_all_log <- log(data_all_df)


#Annotate data and format for pheatmap 
gene_list <- read.delim("genelist_7.txt")

gene_vec <- unlist(gene_list)

uniq_gene_vec <- sort(unique(gene_vec))

selected_genes <- data_all_log[uniq_gene_vec,1:12]

gene_names <- as.data.frame(gene_annotation[uniq_gene_vec,3:4])

labeled_data <- cbind(gene_names,selected_genes)

labeled_data_matrix <- as.matrix(labeled_data[,3:14])

rownames(labeled_data_matrix) = gene_names[,2]


#Generate pheatmaps 
library(pheatmap)

genetype_annotate <- data.frame("Type" = labeled_data$Type)

rownames(genetype_annotate) = rownames(labeled_data_matrix)

colnames(genetype_annotate) = "Gene Type"

treatment_annotate <- data.frame(Treatment = sample_annotation$TreatmentGroup)

rownames(treatment_annotate) = colnames(labeled_data_matrix)

colnames(treatment_annotate) = 'Treatment Group'


pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, 
         annotation_row=genetype_annotate, 
         annotation_col=treatment_annotate, 
         main = "Log(expression) of selected genes", 
         treeheight_col = 10, treeheight_row = 30)


pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, 
         annotation_row=genetype_annotate, 
         annotation_col=treatment_annotate, 
         main = "Log(expression) of selected genes", 
         cluster_cols = FALSE, treeheight_row = 30)



library(grid)


setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=1.0, name="vp", just=c("right","top"))), action="prepend")
pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, annotation_row=genetype_annotate, annotation_col=treatment_annotate, main = "Log(expression) of selected genes", treeheight_col = 5, treeheight_row = 30)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample Name", x=0.32, y=0.1, gp=gpar(fontsize=16))
grid.text("Gene Name", x=-0.05, rot=90, gp=gpar(fontsize=16))


setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=1, name="vp", just=c("right","top"))), action="prepend")
pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, annotation_row=genetype_annotate, annotation_col=treatment_annotate, main = "Log(expression) of selected genes", cluster_cols = FALSE, treeheight_row = 30)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample Name", x=0.32, y=0.1, gp=gpar(fontsize=16))
grid.text("Gene Name", x=-0.05, rot=90, gp=gpar(fontsize=16))





