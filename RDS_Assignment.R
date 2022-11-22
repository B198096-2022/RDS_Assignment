#######################################################################
#                                                                     #
#                   Using R for Data Science Assignment               #
#                                                                     #
#                 B198096  November 2022 Submitted version.           #
#                                                                     #
#######################################################################

#This pipeline is designed to generate a pretty heat map of a subset of 
#Gene expression data from a data file of expression values, using associated 
#Annotation files to label the data for useful interpretation. 

#Note that before the pipeline begins processing the data I need to import 
#The data into the working directory. I performed this in the terminal 
#command line at the beginning of the project. The command for this is below
#I commented ou the command, since it no longer needs to be run to 
#Perform the function of this pipeline, and the necessary files will 
#Be included in the zipped folder with this script

#Retrieving the files 
#cp /shared_files/RDS_assignment_files/s2249132_files.zip .
#unzip s2249132_files.zip



##########################################################
#                                                        #
#. Import data, check for Zeros, and log transform it    #
#                                                        #
##########################################################


#First I need to read the data files into objects for use in the R session 

#all_data is read in as a data table from the file data_all.csv that contains
#All of the results from the experiment 
#The variable containing the table is named to match the file name 
#The table includes 13 columns, X column is indexed to match the 
#indexed gene annotation table, the other 12 columns label the A through L to 
#indicate whicih sample the data is from 
data_all <- read.table(file= "data_all.csv",header=T,sep=",")

#gene_annotation is read in as a data table from the file gene_annotation.csv 
#that contains the gene type (XA, XB, XC) and 
#the long name for every gene in the experiment 
#The variable containing the table is named to match the file name 
gene_annotation <- read.table(file= "gene_annotation.csv",header=T,sep=",")


#sample_annotation is read in as a data table from the file sample_annotation.csv 
#that contains the treatment group assignments for the 12 samples a through L
#The variable containing the table is named to match the file name 
sample_annotation <- read.table(file= "sample_annotation.csv",header=T,sep=",")


#I am now converting the data table into a data frame for easier manipulation 
#I am extracting all rows (1:200) and all columns except for the indexing column
#Because the data frame is already indexed by row number. So the data frame 
#Is made from row 1:200 and columns 2:13 of the data table 
data_all_df <- as.data.frame(data_all[1:200, 2:13])




#These are a set of tests that the can be implemented to find 0s and decimals 
#in the data, since zeros cause problems for log transformations 

#Test to see if there are any zeros
#Make a matrix to scan through 
#The counter helps identify the position in the matrix
#The for loop will iterate through every value in the matrix
#If a given value is zero the function will print the value and position 
#The function then returns a list of the zero positions 
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


#This test is identical to the above function but scans for values that are 
#Less than one, but can be adjusted to set a minimum threshold if the user
#Wishes to specify a different minimum read threshold 
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


#Changing zero values to one 
#One solution to the zero value problem is changing all zeros to one 
#The below function will scan the data for zero and change them to a one 
#But the function can be specified to scan for a minimum threshold and can 
#Can be specified to change found values to a given value

#The arguments are the function are desired data frame (df), 
#The desired filtering threshold (find), default set to zero
#The desired value that the excluded values should be changed to (change), 
#Default set to 1 
#The function transforms the data frame into a matrix, then generates a T/F 
#Matrix with the positions of the values that fit filtering criteria   
#Then this T/F matrix is used to change every position that meets filter criteria
#To the desired change value
#Then convert matrix back to data frame and return it  
change_values <- function(df, find=0, change=1) {
  all_data_matrix <- as.matrix(df)
  find_matrix <- all_data_matrix <= find
  all_data_matrix[find_matrix] <- change
  all_data_df <- as.data.frame(all_data_matrix)
  return(all_data_df)
}
 

#The variable zero is populated with the matrix position of
#Every zero value in the data 

zeros <- zero_test(data_all_df)
filter_response = 'n'

#IF there are no zero values then zero is NULL
#So this if statement will evauate whether there are any zero values, 
#And if there are it will inform the user and then ask if they want to filter 
#Them out

if (is.null(zeros) == FALSE){
  print("Zero values found in data")
  filter_response <- readline(prompt="Filter Data to change Zeros to Ones? (y/n): ")
}

#If the user says yes then the change_values function is called and 
#Changes all zeros to ones 
if (filter_response == "y"|| filter_response == "yes") {
  data_all_df <- change_values(all_data_df)
}


#There are no zeros in the data used for this pipeline so these functions 
#And if statements are not used

#Now I am converting all of the values in the data frame into a log scale 
data_all_log <- log(data_all_df)


##########################################################
#                                                        #
#Time to organize the data and create a labeled data set #
#                                                        #
##########################################################

#Start by reading the text file of the gene list into a list variable 
#This is the personalized set of genes that I am supposed to be analyzing
gene_list <- read.delim("genelist_7.txt")

#Now using the unlist() function to convert the list into a vector  
gene_vec <- unlist(gene_list)

#I am first using the unique() function to 
#make sure that the list is non-redundant 
#Then I am using the sort() function to put the genes in order 
uniq_gene_vec <- sort(unique(gene_vec))

#Now I am filtering the data frame, pulling only the rows that are 
#specified by the uniq_gene_vec (and pulling all columns)
#selected_genes is teh resulting data frame with  log transformed data 
#For all of the genes in my personalized gene set 
selected_genes <- data_all_log[uniq_gene_vec,1:12]
selected_genes

#Now I am making a data frame with all of the gene names 
#This data frame possesses the gene type and long name for all genes in 
#My personalized gene set 
gene_names <- as.data.frame(gene_annotation[uniq_gene_vec,3:4])

#Now I am just binding the gene names to the date 
labeled_data <- cbind(gene_names,selected_genes)
labeled_data

#Lastly, I am making a matrix from the data frame
#This is because the pheatmap function uses matrices as input
labeled_data_matrix <- as.matrix(labeled_data[,3:14])

#I am then using the LongNames of the genes to label the rows of the matrix
rownames(labeled_data_matrix) = gene_names[,2]




##########################################################
#                                                        #
#                 Generating Heat Maps                   #
#                                                        #
##########################################################

#Need to download the pheatmap library 
library(pheatmap)

#GENE ANNOTATION
#This is generating a data frame with a single column pulled from the
#"Type" column of the labeled_data data frame 
genetype_annotate <- data.frame("Type" = labeled_data$Type)

#I am assigning the row names of the genetype_annotate as the same row names from 
#The labeled_data_matrix
rownames(genetype_annotate) = rownames(labeled_data_matrix)
#And then assigning the column name of this data frane as Gene Type
colnames(genetype_annotate) = "Gene Type"


#SAMPLE ANNOTATION 

#This is generating a data frame with a single column pulled 
#from the sample_annotation data frame TreatmentGroup column 
treatment_annotate <- data.frame(Treatment = sample_annotation$TreatmentGroup)

#Importantly, the row names of this data frame need to match the column names 
#of the matrix used to generate the heatmap
rownames(treatment_annotate) = colnames(labeled_data_matrix)
colnames(treatment_annotate) = 'Treatment Group'


#HEATMAP 

#This is now making the heat map from the labeled_data_matrix
#I am scaling it to the rows, setting font size at 6 so that the gene names
#Are not overlapping AND can fit within the graphical window 
#Specifying the row and column annotations made above 
#And labeling the heat map 
#lastly, I specified the height of the clusteing trees, trimming them for 
#Aesthetics 
clustered <- pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, 
                      annotation_row=genetype_annotate, 
                      annotation_col=treatment_annotate, 
                      main = "Log(expression) of selected genes", 
                      treeheight_col = 10, treeheight_row = 30)

#This is the same heat map except it is not clustering the samples 
#This change is made with  cluster_cols = FALSE
nocluster <- pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, 
                      annotation_row=genetype_annotate, 
                      annotation_col=treatment_annotate, 
                      main = "Log(expression) of selected genes", 
                      cluster_cols = FALSE, treeheight_row = 30)

#Now I am adding labels to the X and Y axis (columns and rows, Samples and Genes)
#this requires the grid library to paste the plot onto a grid and add features 
#on top of that grid 
#I am simply putting the heat map onto a grid, 
#Then pasting text on top of it (Sample Name and Gene Name)
#Then specifying the find size and grid location of these text labels 
#The pheatmap function is identical to the clustered heatmap above 
library(grid)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=1.0, name="vp", just=c("right","top"))), action="prepend")
pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, annotation_row=genetype_annotate, annotation_col=treatment_annotate, main = "Log(expression) of selected genes", treeheight_col = 5, treeheight_row = 30)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample Name", x=0.32, y=0.1, gp=gpar(fontsize=16))
grid.text("Gene Name", x=-0.05, rot=90, gp=gpar(fontsize=16))


#Then I am doing the same for the non-sampel clustering heat map 
#The pheatmap function is identical to the nocluster heatmap above 
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=1, name="vp", just=c("right","top"))), action="prepend")
pheatmap(labeled_data_matrix, scale='row', fontsize_row = 6, annotation_row=genetype_annotate, annotation_col=treatment_annotate, main = "Log(expression) of selected genes", cluster_cols = FALSE, treeheight_row = 30)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample Name", x=0.32, y=0.1, gp=gpar(fontsize=16))
grid.text("Gene Name", x=-0.05, rot=90, gp=gpar(fontsize=16))



##########################################################
#                                                        #
#                       The End! :)                      #
#                                                        #
##########################################################

