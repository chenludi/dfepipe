# 2024/02/09
#!/usr/bin/env Rscript
# This is a script used to concate a list of FS file into one by adding frequency count in each bin.
# The first argument is the output directory, 2: output filename 3. input directory 4.input files split by "\n"
# Usage example:
# input_files=$(ls smer.has_qui.s.nc_promoters_states_hg38.mask.bisnp.YRI_2504_108.1kGP.chr*)
# Rscript $script . test.fs . $input_files
args = commandArgs(trailingOnly=TRUE)

read_input_sfs_headtail=function(input_file)  {
  ## Reads input SFS in Dadi Format
  this_file = file(input_file)
  on.exit(close(this_file))
  headtext = readLines(this_file)[1]
  tailtext= readLines(this_file)[3]
  out=c(headtext,tailtext)
  return(out)
} 
read_input_sfs_l2 = function(input_file)  {
  ## Reads input SFS in Dadi Format
  this_file = file(input_file)
  on.exit(close(this_file))
  sfs_string = readLines(this_file)[2]
  output_sfs = as.numeric(unlist(strsplit(sfs_string, ' ')))
  return(output_sfs)
}

# Read output file from arguments
dir_output=args[1]
output=args[2]
dir_input=args[3]


#print(dir_output)
#print(output)
print(dir_input)


setwd(dir_input)

# Concate FS files
fs_l2s=NULL
i=NULL
for (i in c(4:length(args))){
    subfsfile=basename(args[i])
    print(subfsfile)
    if (i==4){
    fs_headtail=read_input_sfs_headtail(subfsfile)
   # print(fs_headtail)
    fs_head=fs_headtail[1]
   # print(fs_head)
    fs_tail=fs_headtail[2]
   # print(fs_tail)
    }
    fs_headtail_i=read_input_sfs_headtail(subfsfile)
    # Check if the projection number and population matches in FS files to concate.
    if(identical(fs_headtail_i,fs_headtail)){
      fs_l2_i=read_input_sfs_l2(subfsfile)
        fs_l2s=rbind(fs_l2s,fs_l2_i)
    }else{print("Error: head or tail text not equal")}
}
i=NULL
out_fs_l2=NULL
for (i in c(1:ncol(fs_l2s))){
  out_fs_l2_coli=sum(fs_l2s[1:nrow(fs_l2s),i])
  out_fs_l2=c(out_fs_l2,out_fs_l2_coli)
}

# Write the concated FS file out.
Lfs_l2= paste(out_fs_l2, collapse = " ")
#write
#print(paste0("Writing :",output))
output_cotent=c(fs_head,Lfs_l2,fs_tail)
setwd(dir_output)
writeLines(output_cotent, output)





