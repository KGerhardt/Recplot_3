#!/usr/bin/env python3

import sys
from optparse import OptionParser
import re
import bisect
import numpy as np

#reads a fasta file, gets the lengths of each sequence, and returns a dictionary each with a single sequence name, 
#its appropriate starts and ends, and a vector of counts to fill with read data, 1 slot per %ID break
def prepare_matrices(contig_file_name, width, bin_height, id_lower):
	#Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
	id_breaks = []
	current_break = 100
	
	while current_break > id_lower:
		id_breaks.append(current_break)
		current_break -= bin_height

	id_breaks = id_breaks[::-1]
	
	
	
	#Each row will look like so.
	zeroes = []
	for i in id_breaks:
		zeroes.append(0)
	
	#current_contig = ""
	
	matrices = {}
	#begin reading contigs and determining their lengths.
	fh = open(contig_file_name)
	
	#the first line should always be a contig, and this grabs it. We only want the unique ID before any spaces, as this is what read aligners grab.
	current_contig = fh.readline()[1:].strip().split()[0]
	contig_length = 0
	
	for line in fh:
		if line[0] == ">":
			#create starts, ends, and counts vectors for the current contig
			starts = []
			ends = []
			pct_id_counts = []
			
			#get number of genome pos bins and their approx. width
			
			num_bins = int(contig_length / width)
			bin_width = (contig_length / num_bins)-1	
			
			cur_bin_start = 1
			#intentionally skip an iteration 
			for i in range(1, num_bins):
				starts.append(int(cur_bin_start))
				ends.append(int((cur_bin_start+bin_width)))
				pct_id_counts.append(zeroes[:])
				cur_bin_start+=bin_width+1
			
			#fill the last end with the actual contig end, in case of rounding errors.
			starts.append(int(cur_bin_start))
			ends.append(contig_length)
			pct_id_counts.append(zeroes[:])
			
			#last thing to do is to fill the matrix with the finished start-end-count set
			matrices[current_contig] = [starts, ends, pct_id_counts]
			
			#set to new contig. One final loop of starts ends counts is needed
			current_contig = line[1:].strip().split()[0]
			contig_length = 0
		else :
			contig_length += len(line.strip())
	
	fh.close()
	
	#need to fill the last contig slot; the loop terminates upon the last DNA line, and does not add this for the final contig otherwise
	starts = []
	ends = []
	pct_id_counts = []
		

	#print(contig_length, width, contig_length/width)
		
	num_bins = int(contig_length/width)
		
	#print(int(contig_length/width))
		
	bin_width = (contig_length / num_bins)-1	
			
	cur_bin_start = 1
	#intentionally skip an iteration 
	for i in range(1, num_bins):
		starts.append(int(cur_bin_start))
		ends.append(int((cur_bin_start+bin_width)))
		pct_id_counts.append(zeroes[:])
		cur_bin_start+=bin_width+1
			
	#fill the last end with the actual contig end, in case of rounding errors.
	starts.append(int(cur_bin_start))
	ends.append(contig_length)
	pct_id_counts.append(zeroes[:])
			
	#last thing to do is to fill the matrix with the finished start-end-count set
	matrices[current_contig] = [starts, ends, pct_id_counts]	

		
	return(matrices, id_breaks)

def numpy_matrices(contig_file_name, width, bin_height, id_lower):
	#Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
	id_breaks = []
	current_break = 100
	
	while current_break > id_lower:
		id_breaks.append(current_break)
		current_break -= bin_height

	id_breaks = id_breaks[::-1]

	
	#TODO
	#In the final version, this will be accessec through the database. I have this here only to develop code.
	matrices = {}
	#begin reading contigs and determining their lengths.
	fh = open(contig_file_name)
	
	#the first line should always be a contig, and this grabs it. We only want the unique ID before any spaces, as this is what read aligners grab.
	current_contig = fh.readline()[1:].strip().split()[0]
	contig_length = 0
	
	for line in fh:
		if line[0] == ">":
			#create starts, ends, and counts vectors for the current contig
			
			#get number of genome pos bins and their approx. width
			
			num_bins = int(contig_length / width)
						
			starts = np.linspace(1, contig_length, num = num_bins, dtype = np.uint64, endpoint = False)
			ends = np.append(starts[1:]-1, contig_length)
			
			numpy_arr = np.zeros((len(starts), len(id_breaks)), dtype = np.uint32)
			
			matrices[current_contig] = [starts, ends, numpy_arr]
			
			current_contig = line[1:].strip().split()[0]
			contig_length = 0
		else :
			contig_length += len(line.strip())
			
	fh.close()
	
	#Last iter. never happens in the loop - happens here, now
	num_bins = int(contig_length / width)
	
	starts = np.linspace(1, contig_length, num = num_bins, dtype = np.uint64, endpoint = False)
	ends = np.append(starts[1:]-1, contig_length)
			
	numpy_arr = np.zeros((len(starts), len(id_breaks)), dtype = np.uint32)
			
	matrices[current_contig] = [starts, ends, numpy_arr]
			
	return(matrices, id_breaks)			
	
def numpy_fill(fillable_matrix, breaks, rec_file):
	print("fill nummies")
	rec = open(rec_file)
	
	for line in rec:
		segment = line.split()
				
		start = int(segment[1])
		end = int(segment[2])
		pct_id = float(segment[3])
		ref = segment[4]
		total_count = end-start+1
		
		which_id = bisect.bisect_right(breaks, pct_id)-1
		start_idx = bisect.bisect_left(fillable_matrix[ref][1], start)
		end_idx = bisect.bisect_left(fillable_matrix[ref][1], end)
		
		#which_id = np.searchsorted(breaks, pct_id, side = 'right')-1
		#start_idx = np.searchsorted(fillable_matrix[ref][1], start)
		#end_idx = np.searchsorted(fillable_matrix[ref][1], end)
		
		if end_idx == len(fillable_matrix[ref][1]):
			continue
		
		#If the read just falls into 1 bin, throw it right on in
		if start_idx==end_idx:
			fillable_matrix[ref][2][start_idx][which_id] += total_count
		#if the read crosses bin boundaries, add the appropriate amount to each successive bin until the read is in the final bin to fill, then throw it in as above.
		else:
			for j in range(start_idx, end_idx+1):
				overflow = end - fillable_matrix[ref][1][j]
				if overflow > 0:
					fillable_matrix[ref][2][j][which_id] += (total_count-overflow)
					total_count = overflow
				else :
					fillable_matrix[ref][2][j][which_id] += total_count
	
	rec.close()
	
	print("nummies full")
	
	return(fillable_matrix)

#This function is designed to take sam-format lines piped from a samtools view command 
#and fill a recruitment matrix with them, terminating further processing there.		
def receive_sam(matrix, breaks, export, reads):
	
	if reads == "":

		for line in sys.stdin:
		#DON'T pop the line out to console again in this version
			if export:
				print(line, end='') 
			#We still want to pass the sam header to samtools, but don't want to work on it for L/R processing
			if "MD:Z:" not in line:
				continue
			else :
				segment = line.split()
				
				ref = segment[2]
				
				if ref not in matrix:
					continue
				
				#Often the MD:Z: field will be the last one in a magicblast output, but not always. Therefore, start from the end and work in.
				iter = len(segment)-1
				mdz_seg = segment[iter]
				
				#If it's not the correct field, proceed until it is.
				while not mdz_seg.startswith("MD:Z:"):
					iter -= 1
					mdz_seg = segment[iter]
				
				#Remove the MD:Z: flag from the start
				mdz_seg = mdz_seg[5:]
				
				match_count = re.findall('[0-9]+', mdz_seg)
				
				sum=0
				
				for num in match_count:
					sum+=int(num)
				
				total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
				
				pct_id = (sum/(total_count))*100
				which_id = bisect.bisect_right(breaks, pct_id)-1
				
				start = int(segment[3])
				end = start+total_count-1
				
				
				#Find the first bin end >= to the read's start and end points. This is the set of bins covered by each read
				start_idx = bisect.bisect_left(matrix[ref][1], start)
				end_idx = bisect.bisect_left(matrix[ref][1], end)
				
				#there are (a very small number of) reads which align past the end of the contig. This removes them.
				if end_idx == len(matrix[ref][1]):
					continue
				
				#If the read just falls into 1 bin, throw it right on in
				if start_idx==end_idx:
					matrix[ref][2][start_idx][which_id] += total_count
				#if the read crosses bin boundaries, add the appropriate amount to each successive bin until the read is in the final bin to fill, then throw it in as above.
				else:
					for j in range(start_idx, end_idx+1):
						overflow = end - matrix[ref][1][j]
						if overflow > 0:
							matrix[ref][2][j][which_id] += (total_count-overflow)
							total_count = overflow
						else :
							matrix[ref][2][j][which_id] += total_count
	else:
	
		rd = open(reads, "r")
	
		for line in rd:
		
			#We still want to pass the sam header to samtools, but don't want to work on it for L/R processing
			if "MD:Z:" not in line:
				continue
			else :
				segment = line.split()
				
				ref = segment[2]
				
				if ref not in matrix:
					continue
				
				#Often the MD:Z: field will be the last one in a magicblast output, but not always. Therefore, start from the end and work in.
				iter = len(segment)-1
				mdz_seg = segment[iter]
				
				#If it's not the correct field, proceed until it is.
				while not mdz_seg.startswith("MD:Z:"):
					iter -= 1
					mdz_seg = segment[iter]
				
				#Remove the MD:Z: flag from the start
				mdz_seg = mdz_seg[5:]
				
				match_count = re.findall('[0-9]+', mdz_seg)
				
				sum=0
				
				for num in match_count:
					sum+=int(num)
				
				total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
				
				pct_id = (sum/(total_count))*100
				which_id = bisect.bisect_right(breaks, pct_id)-1
				
				start = int(segment[3])
				end = start+total_count-1
				
				
				#Find the first bin end >= to the read's start and end points. This is the set of bins covered by each read
				start_idx = bisect.bisect_left(matrix[ref][1], start)
				end_idx = bisect.bisect_left(matrix[ref][1], end)
				
				#there are (a very small number of) reads which align past the end of the contig. This removes them.
				if end_idx == len(matrix[ref][1]):
					continue
				
				#If the read just falls into 1 bin, throw it right on in
				if start_idx==end_idx:
					matrix[ref][2][start_idx][which_id] += total_count
				#if the read crosses bin boundaries, add the appropriate amount to each successive bin until the read is in the final bin to fill, then throw it in as above.
				else:
					for j in range(start_idx, end_idx+1):
						overflow = end - matrix[ref][1][j]
						if overflow > 0:
							matrix[ref][2][j][which_id] += (total_count-overflow)
							total_count = overflow
						else :
							matrix[ref][2][j][which_id] += total_count
	
		rd.close()
		
	return(matrix)

#formats and prints a TSV with the correct output matrices.
def print_super_rec(matrix, breaks, step, prefix, mags):
	#Produces names for %ID bins. Upon readin by plotting code, this provides a consistent way to identify exactly the boundaries of the %ID bins
	names = []
	for i in range(0, len(breaks)):
		names.append(str(round(breaks[i]-step, 3))+"-"+str(round(breaks[i], 3)))
	
	rec = open(prefix+"_recruitment_matrices.tsv", "w")
	
	#Header line in output with colnames
	print("MAG_name", "Contig_name", "Start", "End", *names, sep = "\t", file=rec)
	
	#The data structure is { contig_name: { Pct_ID_bin: [ [starts], [stops], [counts] ] } }
	#the printing format is contig_name \t start_pos \t end_pos \t [tab sep. counts from low %ID to high %ID]
	#This requires wrangling.
	for key in matrix:
		#collect the start/stop bins we only need one copy of these- 
		#there's always a (100-height) to 100% ID bin no matter what; this value is not allowed to be absent
		
		which_mag = mags[key]
		
		for i in range(0, len(matrix[key][0])):
			#print(key, matrix[key][0][i], matrix[key][1][i], *matrix[key][2][i][::-1], sep = "\t", file=rec)
			print(which_mag, key, matrix[key][0][i], matrix[key][1][i], *matrix[key][2][i], sep = "\t", file=rec)
	
	rec.close()
	
#formats and prints a TSV with the correct output matrices.
def print_super_rec_anir(matrix, breaks, step, prefix, anir_mat, mags):
	#Produces names for %ID bins. Upon readin by plotting code, this provides a consistent way to identify exactly the boundaries of the %ID bins
	names = []
	for i in range(0, len(breaks)):
		names.append(str(breaks[i]-step)+"-"+str(breaks[i]))
	
	rec = open(prefix+"_recruitment_matrices.tsv", "w")
	
	#because we don't know what contigs go with which MAG at this point, outputting the contribution to ANIr at each bin value is the correct way to go here.
	#The contributions over all contigs belonging a MAG with pct_id >= some cutoff can be simply summed to calculate MAG ANIr.
	anir = open(prefix+"_ANIr_raw.tsv", "w")
	
	
	#Header line in output with colnames
	print("MAG_name", "Contig_name", "Start", "End", *names, sep = "\t", file=rec)
	print("MAG_name", "contig_name", "component", *names, file = anir)
	
	#Should this be two files, numerators and denominators, with resp. bin mins as cols and all values per contig as data? Probably
	for contig in anir_mat:
		which_mag = mags[contig]
		for i in [0, 1]:
			if i == 0:
				print(which_mag, contig, "numerator", *anir_mat[contig][0], file = anir)
			if i == 1:
				print(which_mag, contig, "denominator", *anir_mat[contig][1], file = anir)
			
	anir.close()
	
	#The data structure is { contig_name: { Pct_ID_bin: [ [starts], [stops], [counts] ] } }
	#the printing format is contig_name \t start_pos \t end_pos \t [tab sep. counts from low %ID to high %ID]
	#This requires wrangling.
	for key in matrix:
		#collect the start/stop bins we only need one copy of these- 
		#there's always a (100-height) to 100% ID bin no matter what; this value is not allowed to be absent

		which_mag = mags[key]
		
		for i in range(0, len(matrix[key][0])):
			print(which_mag, key, matrix[key][0][i], matrix[key][1][i], *matrix[key][2][i], sep = "\t", file=rec)
	
	rec.close()

def receive_blast_like(matrix, breaks, export, reads):
	
	if reads == "":
	
		for line in sys.stdin:
		
			if line.startswith("#"):
				pass
				
			if export:
				print(line, end='') 
		
			segment = line.split("\t")
			
			ref = segment[1]
			
			if ref not in matrix:
				continue

			pct_id = float(segment[2])
			
			pos1 = int(segment[8])
			pos2 = int(segment[9])
			
			start = min(pos1, pos2)
			end = start+(max(pos1, pos2)-min(pos1, pos2))
			
			total_count = end-start+1
			which_id = bisect.bisect_right(breaks, pct_id)-1
			
			#Find the first bin end >= to the read's start and end points. This is the set of bins covered by each read
			start_idx = bisect.bisect_left(matrix[ref][1], start)
			end_idx = bisect.bisect_left(matrix[ref][1], end)
			
			#there are (a very small number of) reads which align past the end of the contig. This removes them.
			if end_idx == len(matrix[ref][1]):
				continue
			
			#If the read just falls into 1 bin, throw it right on in
			if start_idx==end_idx:
				matrix[ref][2][start_idx][which_id] += total_count
			#if the read crosses bin boundaries, add the appropriate amount to each successive bin until the read is in the final bin to fill, then throw it in as above.
			else:
				for j in range(start_idx, end_idx+1):
					overflow = end - matrix[ref][1][j]
					if overflow > 0:
						matrix[ref][2][j][which_id] += (total_count-overflow)
						total_count = overflow
					else :
						matrix[ref][2][j][which_id] += total_count
						
	else:
	
		rd = open(reads, "r")
		
		for line in rd:
		
			if line.startswith("#"):
				pass
		
			segment = line.split("\t")
			
			ref = segment[1]
			
			if ref not in matrix:
				continue

			pct_id = float(segment[2])
			
			pos1 = int(segment[8])
			pos2 = int(segment[9])
			
			start = min(pos1, pos2)
			end = start+(max(pos1, pos2)-min(pos1, pos2))
			
			total_count = end-start+1
			which_id = bisect.bisect_right(breaks, pct_id)-1
			
			#Find the first bin end >= to the read's start and end points. This is the set of bins covered by each read
			start_idx = bisect.bisect_left(matrix[ref][1], start)
			end_idx = bisect.bisect_left(matrix[ref][1], end)
			
			#there are (a very small number of) reads which align past the end of the contig. This removes them.
			if end_idx == len(matrix[ref][1]):
				continue
			
			#If the read just falls into 1 bin, throw it right on in
			if start_idx==end_idx:
				matrix[ref][2][start_idx][which_id] += total_count
			#if the read crosses bin boundaries, add the appropriate amount to each successive bin until the read is in the final bin to fill, then throw it in as above.
			else:
				for j in range(start_idx, end_idx+1):
					overflow = end - matrix[ref][1][j]
					if overflow > 0:
						matrix[ref][2][j][which_id] += (total_count-overflow)
						total_count = overflow
					else :
						matrix[ref][2][j][which_id] += total_count
		rd.close()

	
	return(matrix)
		
#A function for reading args
def opts():
	parser = OptionParser()
	parser.add_option("-c", "--contigs", dest="contigs", help = "This should be a FASTA file containing all and only the contigs that you would like to be part of your recruitment plot.")
	parser.add_option("-m", "--mags", dest="mags", default="", help = "A tab separated file containing the names of MAGs in the first column and contigs in the second column. Every contig should have its parent MAG listed in this file.")	
	parser.add_option("-r", "--reads", dest="reads", default = "", help = "This should be a file with reads aligned to your contigs in any of the following formats: tabular BLAST(outfmt 6), SAM, or Magic-BLAST")
	parser.add_option("-f", "--format", dest="format", default="blast", help="The format of the reads file (write 'blast' or 'sam'). Defaults to tabular BLAST.")
	#parser.add_option("-g", "--genes", dest="genes", default = "", help = "Optional GFF3 file containing gene starts and stops to be use in the recruitment plot.")
	parser.add_option("-i", "--ID-step", dest="id_step", default = 0.5, help = "Percent identity bin width. Default 0.5.")
	parser.add_option("-w", "--bin-width", dest="bin_width", default = 1000, help = "Approximate genome bin width in bp. Default 1000.")
	parser.add_option("-o", "--output", dest="out_file_name", default = "recruitment_plot", help = "Prefix of results to be output. Default: 'recruitment_plot'")
	parser.add_option("-e", "--export", dest="output_line", action='store_true', help = "Output sam lines to stdout?")
	parser.add_option("-s", "--stats", dest="stats", action='store_true', help = "Write ANIr prep file?")
	parser.add_option("--interact", dest="lim_rec", action='store_true', help = "Create lim/rec files for each MAG; do NOT make recruitment matrices.")
	return(parser.parse_args())
	
#Main should be at bottom.
def main():
	options, args = opts()
	
	contigs = options.contigs
	reads = options.reads
	MAGs = options.mags
	format = options.format
	#genes = options.genes
	step = float(options.id_step)
	width = int(options.bin_width)
	prefix = options.out_file_name
	export_lines = options.output_line
	#do_stats = options.stats
	interact = options.lim_rec
	
	mags = {}
	
	#empty_mats, ids = numpy_matrices(contigs, width, step, 70)
	empty_mats, ids = prepare_matrices(contigs, width, step, 70)
	
	
	rec_f = "recplot_mass.rec"
	
	full_mats = numpy_fill(empty_mats, ids, rec_f)
	
	
	fh = open("numpy_experimental.tsv", "w")
	
	for key in full_mats:
		for i in range(0, len(full_mats[key][0])):
			print(key, full_mats[key][0][i], full_mats[key][1][i], *full_mats[key][2][i], file = fh)
	
	fh.close()
	
	
	
#Just runs main.
if __name__ == "__main__":main()

