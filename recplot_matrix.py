#!/usr/bin/env python3

import sys
import re
import bisect
import sqlite3
import argparse
from sys import argv
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np


def sqldb_creation(contigs, mags, sample_reads, format, database):
    """[summary]
    
    Arguments:
        contigs {[type]} -- [description]
        mags {[type]} -- [description]
        sample_reads {[type]} -- [description]
        format {[type]} -- [description]
        database {[type]} -- [description]
    """
    
    # ===== Database and table creation =====
    # Create or open database
    print("Creating databases...")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    # Create lookup table (always creates a new one)
    cursor.execute('DROP TABLE IF EXISTS lookup_table')
    cursor.execute('CREATE TABLE lookup_table \
        (mag_name TEXT, mag_id INTEGER, contig_name TEXT, contig_id INTEGER)')
    # Create sample_info, mag_info and gene_info tables
    cursor.execute('DROP TABLE IF EXISTS mag_info')
    cursor.execute('DROP TABLE IF EXISTS gene_info')
    cursor.execute('DROP TABLE IF EXISTS sample_info')
    cursor.execute('CREATE TABLE mag_info \
        (mag_id INTEGER, contig_id INTEGER, contig_len INTEGER)')
    cursor.execute('CREATE TABLE gene_info \
        (mag_id INTEGER, contig_id INTEGER, gene TEXT, gene_start INTEGER, gene_stop INTEGER)')
    cursor.execute('CREATE TABLE sample_info \
        (sample_name TEXT, sample_id TEXT, sample_number INTEGER)')
    # ========

    # === Extract sample information and save in into DB ===
    # Rename samples provided to avoid illegal names on files
    sampleid_to_sample = {}
    samples_to_db = []
    sample_number = 1
    for sample_name in sample_reads:
        sample_id = "sample_" + str(sample_number)
        samples_to_db.append((sample_name, sample_id, sample_number))
        sampleid_to_sample[sample_id] = sample_name
        sample_number += 1
    # Enter information into table
    cursor.execute("begin")
    cursor.executemany('INSERT INTO sample_info VALUES(?, ?, ?)', samples_to_db)
    cursor.execute('CREATE UNIQUE INDEX sample_index ON sample_info (sample_name)')
    cursor.execute("commit")
    # ========

    # === Extract contig information and MAG correspondence. Save into DB. ===
    # Get contig sizes
    contig_sizes = read_contigs(contigs)
    # Get contig - MAG information
    contig_mag_corresp = get_mags(mags)
    # Initialize variables
    contig_identifiers = []
    mag_ids = {}
    mag_id = 0
    contig_id = 1
    # The dictionary contig_information is important for speed when filling tables
    contig_information = {}
    # Iterate through contig - MAG pairs
    for contig_name, mag_name in contig_mag_corresp.items():
        # Store MAG (and contig) names and ids
        if mag_name in mag_ids:
            contig_identifiers.append((mag_name, mag_ids[mag_name], contig_name, contig_id))
            contig_information[contig_name] = [contig_id, mag_name, mag_ids[mag_name]]
            contig_id += 1
        else:
            mag_id += 1
            mag_ids[mag_name] = mag_id
            contig_identifiers.append((mag_name, mag_ids[mag_name], contig_name, contig_id))
            contig_information[contig_name] = [contig_id, mag_name, mag_ids[mag_name]]
            contig_id += 1
    cursor.executemany('INSERT INTO lookup_table VALUES(?, ?, ?, ?)', contig_identifiers)
    cursor.execute('CREATE INDEX mag_name_index ON lookup_table (mag_name)')
    conn.commit()
    # ========

    # === Fill contig length table ===
    contig_lengths = []
    for contig, contig_len in contig_sizes.items():
        # Get mag_id and contig_id
        sql_command = 'SELECT mag_id, contig_id from lookup_table WHERE contig_name = ?'
        cursor.execute(sql_command, (contig,))
        mag_contig_id = cursor.fetchone()
        contig_lengths.append((mag_contig_id[0], mag_contig_id[1], contig_len))
    cursor.executemany('INSERT INTO mag_info VALUES(?, ?, ?)', contig_lengths)
    cursor.execute('CREATE INDEX mag_id_index ON mag_info (mag_id)')
    conn.commit()
    # ========
    print("Create mag_info_table {}".format(datetime.datetime.now()))
    # === Create one table with information per sample ===
    for sample_name in sampleid_to_sample.keys():
        # Drop if they exist
        cursor.execute('DROP TABLE IF EXISTS ' + sample_name)
        # Create tables once again
        cursor.execute('CREATE TABLE ' + sample_name + \
            ' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
    # === Retrieve information from read mapping and store it ===
    # Read read mapping file for each sample and fill corresponding table
    for sample_name, mapping_file in sampleid_to_sample.items():
        print("Parsing {}... ".format(mapping_file))
        with open(mapping_file) as input_reads:
            record_counter = 0
            records = []
            if format == "blast":
                for line in input_reads:
                    # Commit changes after 500000 records
                    if record_counter == 500000:
                        cursor.execute("begin")
                        cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                        cursor.execute("commit")
                        record_counter = 0
                        records = []
                    if line.startswith("#"):
                        pass
                    else:
                        segment = line.split("\t")
                        try:
                            contig_ref = segment[1]
                        except:
                            print(segment)
                        if contig_ref not in contig_mag_corresp:
                            continue
                        else:
                            pct_id = float(segment[2])
                            pos1 = int(segment[8])
                            pos2 = int(segment[9])
                            start = min(pos1, pos2)
                            end = start+(max(pos1, pos2)-min(pos1, pos2))
                            mag_id = contig_information[contig_ref][2]
                            contig_id = contig_information[contig_ref][0]
                            records.append((mag_id, contig_id, pct_id, start, end))
                            record_counter += 1
                # Commit remaining records
                if record_counter > 0:
                    cursor.execute("begin")
                    cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                    cursor.execute("commit")
                # Create index for faster access
                cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
            if format == "sam":
                for line in input_reads:
                    if record_counter == 500000:
                        cursor.execute("begin")
                        cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                        cursor.execute("commit")
                        record_counter = 0
                        records = []
                    if "MD:Z:" not in line:
                        continue
                    else :
                        segment = line.split()
                        contig_ref = segment[2]
                        if contig_ref not in contig_mag_corresp:
                            continue
                        else:
                            # Often the MD:Z: field will be the last one in a magicblast output, but not always.
                            # Therefore, start from the end and work in.
                            iter = len(segment)-1
                            mdz_seg = segment[iter]
                            # If it's not the correct field, proceed until it is.
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
                            start = int(segment[3])
                            end = start+total_count-1
                            # Get mag_id and contig_id
                            mag_id = contig_information[contig_ref][2]
                            contig_id = contig_information[contig_ref][0]
                            records.append((mag_id, contig_id, pct_id, start, end))
                            record_counter += 1
                # Commit remaining records
                if record_counter > 0:
                    cursor.execute("begin")
                    cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                    cursor.execute("commit")
                # Create index for faster access
                cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
        print("Done")
    conn.commit()
    conn.close()
    # ========


def read_contigs(contig_file_name):
    """ Reads a FastA file and returns
        sequence ids and sizes
    
    Arguments:
        contig_file_name {[str]} -- FastA file location
    Returns:
        [dict] -- Dictionary with ids and sizes
    """
    print("Reading contigs... ", end="", flush=True)
    contig_sizes = {}
    with open(contig_file_name, 'r') as contigs:
        for identifier, sequence in SimpleFastaParser(contigs):
            contig_sizes[identifier] = len(sequence)
    print("done!")
    return contig_sizes


def get_mags(mag_file):
    """ Reads a file with columns:
        Contig_name MAG_name
        and returns the corresponding MAG per contig
    
    Arguments:
        mag_file {[str]} -- MAG correspondence file location
    Returns:
        [dict] -- Dictionary with contigs and corresponding MAG
    """
    mag_dict = {}
    with open(mag_file, 'r') as mags:
        for line in mags:
            mag_contig = line.split()
            mag_dict[mag_contig[0]] = mag_contig[1]
    return mag_dict


def numpy_matrices(database, mag_name, width, bin_height, id_lower):
    """ This function reads a contigs file a prepares the numpy
        matrices required to be filled.
    
    Arguments:
        contig_file_name {[type]} -- [description]
        width {[type]} -- [description]
        bin_height {[type]} -- [description]
        id_lower {[type]} -- [description]
    """
    print("Preparing recruitment matrices...", end="", flush=True)
    # Prep percent identity breaks - always starts at 100 and proceeds 
    # down by bin_height steps until it cannot do so again without passing id_lower
    id_breaks = []
    current_break = 100
    while current_break > id_lower:
        id_breaks.append(current_break)
        current_break -= bin_height
    id_breaks = np.array(id_breaks[::-1])
    

    # Retrieve mag_id from provided mag_name
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT mag_id from lookup_table WHERE mag_name = ?'
    cursor.execute(sql_command, (mag_name,))
    mag_id = cursor.fetchone()[0]
    # Retrieve all contigs from mag_name and their sizes
    sql_command = 'SELECT contig_id, contig_len from mag_info WHERE mag_id = ?'
    cursor.execute(sql_command, (mag_id,))
    contig_sizes = cursor.fetchall()
    # Create matrices for each contig in the mag_name provided
    matrix = {}
    for id_len in contig_sizes:
        contig_length = id_len[1]
        num_bins = int(contig_length / width)
        if num_bins < 1:
            num_bins = 1
        starts = np.linspace(1, contig_length, num = num_bins, dtype = np.uint64, endpoint = False)
        ends = np.append(starts[1:]-1, contig_length)
        numpy_arr = np.zeros((len(id_breaks), len(starts)), dtype = np.uint32)
        matrix[id_len[0]] = [starts, ends, numpy_arr]
    print("Done")
    return mag_id, matrix, id_breaks


def numpy_fill(database, mag_id, sample_name, matrices, id_breaks):
    print("Filling matrices...")
    # Retrieve sample_id from sample_name provided
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT sample_id from sample_info WHERE sample_name = ?'
    cursor.execute(sql_command, (sample_name,))
    sample_id = cursor.fetchone()[0]
    # Retrieve all read information from mag_name and sample_name provided
    sql_command = 'SELECT * from ' + sample_id + ' WHERE mag_id = ?'
    cursor.execute(sql_command, (mag_id,))
    read_information = cursor.fetchall()
    # read_information is (mag_id contig_id perc_id read_start read_stop)
    for read_mapped in read_information:
        if read_mapped[1] in matrices:
            contig_id = read_mapped[1]
            read_start = read_mapped[3]
            read_stop = read_mapped[4]
            read_len = read_stop - read_start + 1
            read_id_index = bisect.bisect_right(id_breaks, read_mapped[2]) - 1
            # print(read_mapped)
            # print(id_breaks)
            # print(read_id_index, id_breaks[read_id_index])
            read_start_loc = bisect.bisect_left(matrices[contig_id][1], read_start)
            read_stop_loc = bisect.bisect_left(matrices[contig_id][1], read_stop)
            # If the read falls entirely on a bin add all bases to the bin
            if read_start_loc == read_stop_loc:
                matrices[contig_id][2][read_id_index][read_start_loc] += read_len
            # On the contrary split bases between two or more bins
            else:
                for j in range(read_start_loc, read_stop_loc + 1):
                    overflow = read_stop - matrices[contig_id][1][j]
                    if overflow > 0:
                        matrices[contig_id][2][read_id_index][j] += (read_len - overflow)
                        read_len = overflow
                    else :
                        matrices[contig_id][2][read_id_index][j] += read_len
    print("Done")
    return matrices

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
        names.append(str(breaks[i]-step)+"-"+str(breaks[i]))
    
    
    
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

#This function is designed to take sam-format lines piped from a samtools view command 
#and fill a recruitment matrix with them, terminating further processing there.		
def receive_sam_anir(matrix, breaks, export):
    
    #contig:min_pct_id:contribution, divisor
    ANIr_collections = {}
    
    zeroes = []
    
    for i in breaks:
        zeroes.append(0)
    
    #we need the read
    for contig in matrix :
        ANIr_collections[contig] = [zeroes[:], zeroes[:]]
    
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
            
            ANIr_collections[ref][0][which_id] += (total_count*pct_id)/100
            ANIr_collections[ref][1][which_id] += total_count
            
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
                        
    return(matrix, ANIr_collections)

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
            
def receive_blast_like_anir(matrix, breaks):
    #contig:min_pct_id:contribution, divisor
    ANIr_collections = {}
    
    zeroes = []
    
    for i in breaks:
        zeroes.append(0)
    
    #we need the read
    for contig in matrix :
        ANIr_collections[contig] = [zeroes[:], zeroes[:]]

    for line in sys.stdin:
        segment = line.split("\t")
        
        if line.startswith("#"):
            pass
        
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
        
        ANIr_collections[ref][0][which_id] += (total_count*pct_id)/100
        ANIr_collections[ref][1][which_id] += total_count		

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
                        
                        
    return(matrix, ANIr_collections)		
    
def blast_rec_file(reads, MAGS, prefix):
    print("Reading BLAST reads... ", end="", flush=True)	
    
    if reads == "":
        for line in sys.stdin:
            
            if line.startswith("#"):
                pass
        
            segment = line.split("\t")
            
            ref = segment[1]
            
            if ref not in MAGS:
                continue

            pct_id = float(segment[2])
            
            pos1 = int(segment[8])
            pos2 = int(segment[9])
            
            start = min(pos1, pos2)
            end = start+(max(pos1, pos2)-min(pos1, pos2))
            
            print(MAGS[ref], str(start), str(end), pct_id, "", sep = "\t", file = rec)
            
    else:
        reads = open(reads, "r")
        rec = open(prefix+".rec", "w")
            
        for line in reads:
        
            if line.startswith("#"):
                pass
        
            segment = line.split("\t")
            
            ref = segment[1]
            
            if ref not in MAGS:
                continue

            pct_id = segment[2]
            
            pos1 = int(segment[8])
            pos2 = int(segment[9])
            
            start = min(pos1, pos2)
            end = start+(max(pos1, pos2)-min(pos1, pos2))
            
            print(MAGS[ref], str(start), str(end), pct_id, ref, sep = "\t", file = rec)
            
        reads.close()
        rec.close()
        
    print("done!")
        
def sam_rec_file(reads, MAGS, prefix):
    print("Reading SAM reads... ", end="", flush=True)	
    if reads == "":

        rec = open(prefix+".rec", "w")
        
        for line in sys.stdin:
        #DON'T pop the line out to console again in this version
            #print(line, end='') 
            #We still want to pass the sam header to samtools, but don't want to work on it for L/R processing
            if "MD:Z:" not in line:
                continue
            else :
            
                segment = line.split()
                
                ref = segment[2]
                
                if ref not in MAGS:
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
                
                start = int(segment[3])
                end = start+total_count-1
                
                print(MAGS[ref], str(start), str(end), str(pct_id), ref, sep = "\t", file = rec)	

        rec.close()
                
    else: 
    
        reads = open(reads, "r")
        rec = open(prefix+".rec", "w")

        for line in reads:
        #DON'T pop the line out to console again in this version
            #print(line, end='') 
            #We still want to pass the sam header to samtools, but don't want to work on it for L/R processing
            if "MD:Z:" not in line:
                continue
            else :
                segment = line.split()
                
                ref = segment[2]
                
                if ref not in MAGS:
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
                
                start = int(segment[3])
                end = start+total_count-1	

                print(MAGS[ref], str(start), str(end), str(pct_id), ref, sep = "\t", file = rec)
                
        reads.close()
        rec.close()
    
    print("done!")
    
def print_lim(contig_length, mag_dict, prefix):
    lim = open(prefix+".lim", "w")
    
    for key in contig_length:
        print(mag_dict[key], key, "1", contig_length[key], sep = "\t", file = lim)
    
    lim.close()
    
#A function for reading args
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds recruitment plots (COMPLETE DESCRIPTION HERE)\n'''
            '''Usage: ''' + argv[0] + ''' COMPLETE\n'''
            '''Global mandatory parameters: -g [Genome Files] OR -p [Protein Files] OR -s [SCG HMM Results] -o [AAI Table Output]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')

    parser.add_argument("-c", "--contigs", dest="contigs", 
    help = "This should be a FASTA file containing all and only the contigs that you would like to be part of your recruitment plot.")
    parser.add_argument("-m", "--mags", dest="mags", default="", help = "A tab separated file containing the names of MAGs in the first column and contigs in the second column. Every contig should have its parent MAG listed in this file.")	
    parser.add_argument("-r", "--reads", dest="reads", nargs='+', help = "This should be a file with reads aligned to your contigs in any of the following formats: tabular BLAST(outfmt 6), SAM, or Magic-BLAST")
    parser.add_argument("-f", "--format", dest="map_format", default="blast", help="The format of the reads file (write 'blast' or 'sam'). Defaults to tabular BLAST.")
    #parser.add_argument("-g", "--genes", dest="genes", default = "", help = "Optional GFF3 file containing gene starts and stops to be use in the recruitment plot.")
    parser.add_argument("-i", "--ID-step", dest="id_step", default = 0.5, help = "Percent identity bin width. Default 0.5.")
    parser.add_argument("-w", "--bin-width", dest="bin_width", default = 1000, help = "Approximate genome bin width in bp. Default 1000.")
    parser.add_argument("-o", "--output", dest="out_file_name", default = "recruitment_plot", help = "Prefix of results to be output. Default: 'recruitment_plot'")
    parser.add_argument("-e", "--export", dest="output_line", action='store_true', help = "Output sam lines to stdout?")
    parser.add_argument("-d", "--database", dest="sql_database", action='store', help = "SQLite database to create or update")
    parser.add_argument("-s", "--stats", dest="stats", action='store_true', help = "Write ANIr prep file?")
    parser.add_argument("--interact", dest="lim_rec", action='store_true', help = "Create lim/rec files for each MAG; do NOT make recruitment matrices.")
    
    args = parser.parse_args()
    
   
    contigs = args.contigs
    reads = args.reads
    mags = args.mags
    map_format = args.map_format
    # genes = args.genes
    step = float(args.id_step)
    width = int(args.bin_width)
    prefix = args.out_file_name
    export_lines = args.output_line
    do_stats = args.stats
    interact = args.lim_rec
    sql_database = args.sql_database
    

    sqldb_creation(contigs, mags, reads, map_format, sql_database)
    mag_id, matrix, id_breaks = numpy_matrices("TEST_DB", "IIa.A_ENTP2013_S02_SV82_300m_MAG_01", width, step, 70)
    matrix = numpy_fill("TEST_DB", mag_id, "03.All_SAR11--ETNP_2013_S02_SV89_300m.blast.bh", matrix, id_breaks)

    mags = {}
    
    # if interact :
    
    #     c_len = read_contigs(contigs)
    
    #     if MAGs == "":
    #         for key in c_len:
    #             mags[key] = key
        
    #         print_lim(c_len, mags, prefix)
            
    #         if format == "blast":
    #             blast_rec_file(reads, mags, prefix)
    #         else:
    #             sam_rec_file(reads, mags, prefix)
            
    #     else:
    #         mags = get_mags(MAGs)
            
    #         removed_contigs = []
            
    #         for key in c_len:
    #             if key not in mags:
    #                 removed_contigs.append(key)
                    
    #         for c in removed_contigs:
    #             del c_len[c]
                
    #         print_lim(c_len, mags, prefix)	
            
    #         if format == "blast":
    #             blast_rec_file(reads, mags, prefix)
    #         else:
    #             sam_rec_file(reads, mags, prefix)
            
    # else:
    
    #     mat, breaks = prepare_matrices(contigs, width, step, 70)
    
    #     if MAGs == "":
    #         #If MAGs aren't supplied, add a column specifying this
    #         for key in mat:
    #             mags[key] = key
            
    #         if format == "blast":
                
    #             mat = receive_blast_like(mat, breaks, export_lines, reads)
    #             print_super_rec(mat, breaks, step, prefix, mags)
                
    #         else:
            
    #             mat = receive_sam(mat, breaks, export_lines, reads)
    #             print_super_rec(mat, breaks, step, prefix, mags)
            
    #     else: 	
    #         mags = get_mags(MAGs)
            
    #         removed_contigs = []
            
    #         for key in mat:
    #             if key not in mags:
    #                 removed_contigs.append(key)
                    
    #         for c in removed_contigs:
    #             del mat[c]
                
    #         if format == "blast":
    #             mat = receive_blast_like(mat, breaks, export_lines, reads)
    #             print_super_rec(mat, breaks, step, prefix, mags)
                    
    #         else:
    #             mat = receive_sam(mat, breaks, export_lines, reads)
    #             print_super_rec(mat, breaks, step, prefix, mags)

#Just runs main.
if __name__ == "__main__":main()

