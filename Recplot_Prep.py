#!/usr/bin/env python3

from optparse import OptionParser
import sys
import re

#Produces the rec file for blast output.
def handle_blast(adjuster, reads, prefix):
	print("Treating reads as BLAST tabular format... ", end="", flush=True)
	
	blast = open(reads)
	rec = open(prefix+".rec", "w")
	
	for line in blast:
		segment = line.split("\t")
		
		ref = segment[1]

		pct_id = segment[2]
		
		pos1 = int(segment[8])
		pos2 = int(segment[9])
		
		start = min(pos1, pos2)+adjuster[ref]
		end = start+(max(pos1, pos2)-min(pos1, pos2))
		name = segment[0]
		
		rec.write(str(start)+"\t"+str(end)+"\t"+pct_id+"\t"+name+"\n")
		
	blast.close()
	rec.close()
	
	print("done!")
	return(0)

#Produces the rec file for magic blast output. Functionally identical to the regular blast function.
def handle_magic_blast(adjuster, reads, prefix):
	print("Treating reads as Magic-BLAST format... ", end="", flush=True)
	magic = open(reads)
	rec = open(prefix+".rec", "w")
	
	
	
	for line in magic:
		if line.startswith("#"):
			pass
		
		segment = line.split("\t")
		
		ref = segment[1]

		pct_id = segment[2]
		
		pos1 = int(segment[8])
		pos2 = int(segment[9])
		
		start = min(pos1, pos2)+adjuster[ref]
		end = start+(max(pos1, pos2)-min(pos1, pos2))
		name = segment[0]
		
		rec.write(str(start)+"\t"+str(end)+"\t"+pct_id+"\t"+name+"\n")
	
	magic.close()
	rec.close()
	print("done!")
	return(0)

#Produces the rec file for a SAM format set of reads. Requires the MD:Z: flag to be filled in the sam file, which is default behavior for bowtie2.
def handle_sam(adjuster, reads, prefix):
	print("Treating reads as SAM format... ", end="", flush=True)
	
	sam = open(reads)
	rec = open(prefix+".rec", "w")
	
	for line in sam:
		if "MD:Z:" not in line:
			continue
		else :
			#This flag is neither guaranteed to be present or in the same column. Has to be found.
			mdz_seg = re.search('MD:Z:(.+)\t', line).group(1).split("\t")[0]
			match_count = re.findall('[0-9]+', mdz_seg)
			
			sum=0
			for num in match_count:
				sum+=int(num)
				
			mismatch_count = len(''.join([i for i in mdz_seg if not i.isdigit()]))
			
			pct_id = (sum/(sum+mismatch_count))*100

			segment = line.split("\t")
			ref = segment[2]
			
			start = int(segment[3])+adjuster[ref]
			
			end = start+len(segment[9])
			name = segment[0]
			
			rec.write(str(start)+"\t"+str(end)+"\t"+str(pct_id)+"\t"+name+"\n")
	
	sam.close()
	rec.close()
	
	print("done!")
	return(0)
	
#Gets the dimensions of a counts matrix for a recplot. Possibly useful in a later script to handle very big datasets. Not ready.
def prep_matrix(id_step, bin_width, starts, ends):
	print("Preparing counts matrix... ", end="", flush=True)
	
	#Y axis of the recplot matrix
	identity_seq = [100]
	id = 100-id_step
	while id > 70 :
		identity_seq.append(id)
		id -= id_step
	identity_seq.append(70)
	
	bin_seq = []
	
	#X axis of the recplot matrix. Divides contig runs into chunks of roughly bin_width bp; as close as possible without making bins smaller than this size.
	for i in range(0, len(starts)):
		num_chunks = (ends[i]-starts[i]+1)//bin_width
		chunk_size = (ends[i]-starts[i]+1)/num_chunks
		#print(num_chunks, chunk_size)
		
		bin_seq.append(starts[i])
		
		current = starts[i]+chunk_size
		for j in range(1, num_chunks-1):
			bin_seq.append(current)
			current += chunk_size
			
		bin_seq.append(ends[i])
		
	for i in range(0, len(bin_seq)):
		bin_seq[i] = round(bin_seq[i])
	
	
	
	print("done!")
	return(0)
	
#Reads a FASTA file and gets the lengths of its sequences. Sorts the contigs by length, descending, and then provides start/stop locations assuming a concatenation of these contigs in order from long to short. 
#Produces the .LIM file using the supplied prefix, and returns the contigs and their starts and ends for producing the .REC file
def read_contigs(contig_file_name, out_prefix):
	print("Reading contigs... ", end="", flush=True)
	
	current_contig = ""
	output = {}
	
	fh = open(contig_file_name)
	
	for line in fh:
		if line[0] == ">":
			current_contig = line[1:].strip()
			output[current_contig] = 0
		else :
			output[current_contig] += len(line.strip())
	
	fh.close()
	
	contig_names = []
	contig_length = []
	
	for key in sorted(output, key = output.get, reverse = True):
			contig_names.append(key)
			contig_length.append(output[key])
			
			contig_start = [1]
			contig_ends = [contig_length[0]]
			
	for i in range(1, len(contig_length)):
		contig_start.append(contig_ends[i-1]+1)
		contig_ends.append(contig_start[i]+contig_length[i])
		
	fh = open(out_prefix+".lim", "w")
	
	for i in range(0, len(contig_names)):
		fh.write(contig_names[i]+"\t"+str(contig_start[i])+"\t"+str(contig_ends[i])+"\n")
	
	fh.close()
		
	print("done!")
	
	return(contig_names, contig_start, contig_ends)

#A function for reading CLI
def opts():
	parser = OptionParser()
	parser.add_option("-c", "--contigs", dest="contigs", help = "This should be a FASTA file containing all and only the contigs that you would like to be part of your recruitment plot.")
	parser.add_option("-r", "--reads", dest="reads", help = "This should be a file with reads aligned to your contigs in any of the following formats: tabular BLAST(outfmt 6), SAM, or Magic-BLAST")
	parser.add_option("-f", "--format", dest="format", default="blast", help="The format of the reads file (write 'blast', 'sam', or 'magic'). Defaults to tabular BLAST.")
	#parser.add_option("-g", "--genes", dest="genes", default = "", help = "Optional GFF3 file containing gene starts and stops to be use in the recruitment plot.")
	#parser.add_option("-i", "--ID-step", dest="id_step", default = 0.5, help = "Percent identity bin width. Default 0.5.")
	#parser.add_option("-w", "--bin-width", dest="bin_width", default = 1000, help = "Approximate genome bin width in bp. Default 1000.")
	parser.add_option("-o", "--output", dest="out_file_name", default = "recruitment_plot", help = "Prefix of results to be output. Default: 'recruitment_plot'")
	return(parser.parse_args())
	
#Program starts here. Acquires file names, format of reads, and prefix for output. Creates the .lim file with read_contigs, then calls the appropriate read handler to produce the .rec file.
def main():
	options, args = opts()
	
	contigs = options.contigs
	reads = options.reads
	format = options.format
	#genes = options.genes
	#step = options.id_step
	#width = options.bin_width
	prefix = options.out_file_name
	
	names, starts, ends = read_contigs(contigs, prefix)
	#prep_matrix(step, width, starts, ends)
	adjust = {}
	for i in range(0, len(names)):
		adjust[names[i]] = starts[i]-1
	
	if format == "blast":
		handle_blast(adjust, reads, prefix)
		
	if format == "magic":
		handle_magic_blast(adjust, reads, prefix)
		
	if format == "sam":
		handle_sam(adjust, reads, prefix)

#Just runs main.
if __name__ == "__main__":main()