#!/usr/bin/env python 
from common_pipeline_functions import comment_or_blank
from common_pipeline_functions import get_fastq_files
from common_pipeline_functions import create_slurm_header
from common_pipeline_functions import find_matching_reads
from common_pipeline_functions import write_command_lists_to_file
from common_pipeline_functions import read_key_word_input_file
#from common_pipeline_functions import construct_hisat_command
from common_pipeline_functions import construct_hisat_array_command
from common_pipeline_functions import convert_list_to_bash_array
from common_pipeline_functions import dereference_array_string
from common_pipeline_functions import construct_cufflinks_array_command
from common_pipeline_functions import construct_feature_counts_command
from common_pipeline_functions import construct_cuffnorm_command
import os 

def construct_samtools_commands(list_of_sam_files,deref_sam_array,threads=3):
	"""Puts the arguments needed to create a sorted bam file from the
	list of sam files. 
	Returns a tuple where the first item is a list of the commands
	and the second is a list of the bam files"""
	samtools_program="samtools"
	task_id="$SLURM_ARRAY_TASK_ID"
	commands=[]
	bam_files=[]
	bam_files_array_name="sorted_bam_files_array"
	command=[samtools_program, "sort", "-@", str(threads),"-o"]
	for sam in list_of_sam_files:
		
		name_stem,file_type=sam.split(".")
		sorted_out_bam=name_stem+"_sorted.bam"
		bam_files.append(sorted_out_bam)
	bam_array=convert_list_to_bash_array(bam_files,bam_files_array_name)
	command.append(dereference_array_string(bam_files_array_name,task_id))
	command.append(deref_sam_array)
	command_str=" ".join(command)
	return(bam_array,command_str,bam_files)


def convert_list_of_list_to_bash_array(list_of_lists,array_name):
	array_element_strings=[]
	array_prefix="declare -a "+array_name+"="+"("
	for l in list_of_lists:
		array_element='"'+" ".join(l)+'"'
		array_element_strings.append(array_element)

	array_content_string=" ".join(array_element_strings)
	array_string=array_prefix+array_content_string+")"
	return array_string
	
def convert_list_to_bash_array(input_list,array_name):
	array_element_strings=[]
	array_prefix="declare -a "+array_name+"="+"("
	for i in input_list:
		array_element='"'+i+'"'
		array_element_strings.append(array_element)
	array_content_string=" ".join(array_element_strings)
	array_string=array_prefix+array_content_string+")"
	return array_string

def write_out_hisat_array_slurm_file(arguments_file,outfile="hisat_array_commands.sh"):
	"""Uses data from input file to write out the hisat commands as an array job"""
	
	arguments_dic=read_key_word_input_file(arguments_file,write_output_dir=True)
	#hisat_array_name="hisat_array"
	#sam_array="sam_array" I don't think I need this for right now
	samtools_sort_array_name="samtools_sort_array"
	#bam_array="bam_array" I don't think I need this for anything right now
	task_id="$SLURM_ARRAY_TASK_ID"
	workdir=arguments_dic["main_output_dir"]
	print workdir
	#Create the directory for the alignment files
	alignment_dir="Alignments"
	os.mkdir(workdir+"/"+alignment_dir)
	#Dictionary pointing from sample names to reads
	read_dict=find_matching_reads(get_fastq_files(arguments_dic["sample_dir"]))
	#List of slurm options to be written on their own lines
        slurm_commands=create_slurm_header(workdir=arguments_dic["main_output_dir"],
                                        time=":".join([arguments_dic["hours"],arguments_dic["minutes"],"00"]),
                                        ntasks=arguments_dic["threads"],
                                        mem=arguments_dic["mem"],
                                        mail_user=arguments_dic["email_to_send_output_to"],submit_to_rra=arguments_dic["submit_to_rra"])



	hisat_command_input,sam_files=construct_hisat_array_command(index_files_path_and_basename=arguments_dic["genome_index"],
                                                        read_dict=read_dict,
                                                        number_of_processors=arguments_dic["threads"],
                                                        input_dict=arguments_dic,alignment_dir=alignment_dir)
	


	read1_array,read2_array,hisat_output_array,hisat_command,deref_hisat_outfiles=hisat_command_input
	#print hisat_command_input
	bam_array,samtools_sort_command,bam_files=construct_samtools_commands(sam_files,deref_hisat_outfiles,threads=arguments_dic["threads"])
	arguments_dic["alignment_files"]=bam_files
	job_size,cufflinks_input_array_string,cufflinks_outputdirs_array_string,cufflinks_command=construct_cufflinks_array_command(arguments_dic)
	gff_file=arguments_dic["gff"]
	if "gff_annotation" in arguments_dic:
		id_value=arguments_dic["gff_annotation"]
	else:
		 id_value="gene_id"
	featurecounts_command_list=construct_feature_counts_command(gff_file,bam_files,id_value=id_value)
	featurecounts_command=" ".join(featurecounts_command_list)
	normalization_groups=arguments_dic["normalization_groups"]
	normalization_method=arguments_dic["normalization_method"]
	cuffnorm_command_list=construct_cuffnorm_command(gff_file,bam_files,alignment_dir+"/",groupings=normalization_groups,normalization_method=normalization_method)
	cuffnorm_command=" ".join(cuffnorm_command_list)
	narray=len(sam_files)
	##May want to include a warning here if only one job is allowed to run at a time
	narray_arg="#SBATCH --array=0-"+str(narray-1)+"%"+str(narray)
	slurm_commands.append(narray_arg)
	slurm_commands.append("#SBATCH --output=rna-seq_pipeline.out")
	slurm_commands.append("module purge")
	slurm_commands.append("module load apps/hisat2/2.1.0")
	slurm_commands.append("module load apps/samtools/1.3.1")
	slurm_commands.append("module load apps/cufflinks/2.2.1")
	slurm_commands.append("module load apps/subread/1.6.3")
	
	list_of_commands=[[read1_array],["\n"],[read2_array],["\n"],[hisat_output_array],["\n"],[bam_array],["\n"],[cufflinks_outputdirs_array_string],
["\n"],[hisat_command],["\n"],[samtools_sort_command],["\n"],[cufflinks_command],["\n"],[featurecounts_command],["\n"],[cuffnorm_command]]
	write_command_lists_to_file(slurm_commands=slurm_commands,list_of_command_lists=list_of_commands,outfile=outfile)
	return list_of_commands

#arguments_file="/shares/biocomputing_ii/Bitbucket/Hisat_Pipeline/Sample_Data/sample_hisat_array_input.txt"


if __name__=="__main__":
	import sys
	import subprocess
	args=sys.argv
	program,infile,outfile=args
	write_out_hisat_array_slurm_file(infile,outfile)
	subprocess.call(["sbatch",outfile])
