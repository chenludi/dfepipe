#2022/02/12
# I added an allernative column name for length_fsfile and length_inputvcf because I used a different
# column name in the fsmut files for conserved sites data.
#This is a script used to concate a list of FS file into one by adding frequency count in each bin.
# The first argument is the output directory, 2: output filename 3. input directory 4.input files split by "\n"
# Usage example:
# input_files=$(ls fsmut.smer.has_qui.s.nc_promoters_states_hg38.chr*)
# Rscript $script . test.tab . $input_files
args = commandArgs(trailingOnly=TRUE)
# Read output file from arguments
dir_output=args[1]
output=args[2]

# Read input file from arguments
dir_input=args[3]
print(paste0("dir_input:",dir_input))
setwd(dir_input)
print(paste0("inputfile1:",args[4]))
table_sum=NULL
i=NULL
theta_all=0
for (i in c(4:length(args))){
  subfsfile=basename(args[i])
  print(subfsfile)
  table_subfsfile=read.table(subfsfile,header = TRUE, sep = "\t")
  # print(table_subfsfile[1,1])  
  theta_all=theta_all+table_subfsfile[1,'avepersite_diploi_Roulette_mutrate']*table_subfsfile[1,'length_adj_inputbed']
  table_sum=rbind(table_sum,table_subfsfile)
}

# check if input files have the same projection number
if (length(unique(table_sum[,'projection_num'])) !=1){
  print("Error: The projection number is not the same in the subfiles.")
  q()
}

#
colnames_table_sum=colnames(table_sum)

fsfile=paste0(table_sum[1,'fsfile'],".etal")
# Because I renamed the column, I need to add an option to try another column name

if('length_fsfile' %in% colnames_table_sum) {
  length_fsfile=sum(table_sum[,'length_fsfile'])

} else if('sites_fsfile' %in% colnames_table_sum)
{
  length_fsfile=sum(table_sum[,'sites_fsfile'])

}else{
  print("Error: The column length_fsfile or sites_fsfile is not in the subfiles.")
}

if('length_inputvcf' %in% colnames_table_sum) {
  length_inputvcf=sum(table_sum[,'length_inputvcf'])
} else if ('sites_vcf' %in% colnames_table_sum) {
  length_inputvcf=sum(table_sum[,'sites_vcf'])
}else {
     print("Error: The column length_inputvcf or sites_inputvcf is not in the subfiles.")
  }


ratio_keptsites=length_fsfile/length_inputvcf
length_inputbed=sum(table_sum[,'length_inputbed'])
length_adj_inputbed=sum(table_sum[,'length_adj_inputbed'])
projection_num=as.character(unique(table_sum[,'projection_num']))
mutbed=paste0(table_sum[1,'mutbed'],".etal")
length_withmutinfo=sum(table_sum[,'length_withmutinfo'])
inputbed=paste0(table_sum[1,'inputbed'],".etal")
ratio_withmutinfo=length_withmutinfo/length_inputbed
avepersite_diploi_Roulette_mutrate=theta_all/length_adj_inputbed

out_colnames=colnames(table_subfsfile)

out_line=c(fsfile,length_fsfile,length_inputvcf,
           ratio_keptsites,length_inputbed,length_adj_inputbed,
           projection_num,mutbed,length_withmutinfo,
           inputbed,length_inputbed,ratio_withmutinfo,
           avepersite_diploi_Roulette_mutrate)
table_out=rbind(out_line,table_sum)

# write a file out
setwd(dir_output)
write.table(table_out,output,sep="\t",quote = F,row.names = F)


