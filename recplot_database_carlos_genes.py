#!/usr/bin/env python3

import sys
import re
import bisect
import sqlite3
import argparse
from sys import argv
import platform

if platform.system() == "Windows":
    #print("Pysam cannot be loaded on windows. You will be unable to process BAM files directly.")
else:
    import pysam

def get_sys():
    return(platform.system())

def sqldb_creation(contigs, mags, sample_reads, map_format, database):
    """ Read information provided by user and creates SQLite3 database
    
    Arguments:
        contigs {str} -- Location of fasta file with contigs of interest.
        mags {str} -- Location of tab separated file with contigs and their corresponding mags.
        sample_reads {list} -- Location of one or more read mapping results file(s).
        format {str} -- Format of read mapping (blast or sam).
        database {str} -- Name (location) of database to create.
    """
    
    # ===== Database and table creation =====
    # Create or open database
    #print("Creating database...")
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    # Create lookup table (always creates a new one)
    cursor.execute('DROP TABLE IF EXISTS lookup_table')
    cursor.execute('CREATE TABLE lookup_table \
        (mag_name TEXT, mag_id INTEGER, contig_name TEXT, contig_id INTEGER)')
    # Create sample_info, mag_info, mags_per_sample, and gene_info tables
    cursor.execute('DROP TABLE IF EXISTS mag_info')
    cursor.execute('DROP TABLE IF EXISTS sample_info')
    cursor.execute('DROP TABLE IF EXISTS mags_per_sample')
    cursor.execute('CREATE TABLE mag_info \
        (mag_id INTEGER, contig_id INTEGER, contig_len INTEGER)')
    cursor.execute('CREATE TABLE sample_info \
        (sample_name TEXT, sample_id TEXT, sample_number INTEGER)')
    cursor.execute('CREATE TABLE mags_per_sample \
        (sample_name TEXT, mag_name TEXT)')
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
            contig_information[contig_name] = [mag_name, mag_ids[mag_name], contig_id]
            contig_id += 1
        else:
            mag_id += 1
            mag_ids[mag_name] = mag_id
            contig_identifiers.append((mag_name, mag_ids[mag_name], contig_name, contig_id))
            contig_information[contig_name] = [mag_name, mag_ids[mag_name], contig_id]
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
        mags_in_sample = []
        #print("Parsing {}... ".format(mapping_file), end = "",flush = True)
        contigs_in_sample = save_reads_mapped(mapping_file, sample_name, map_format, cursor, conn)
        cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
        all_contigs = cursor.fetchall()
        for element in all_contigs:
            if element[0] in contigs_in_sample:
                if element[1] not in mags_in_sample:
                    mags_in_sample.append(element[1])
                else:
                    continue
            else:
                continue
        mags_in_sample = [(mapping_file, x) for x in mags_in_sample]
        cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
        conn.commit()
        #print("Database creation finished!")
    cursor.close()
    conn.commit()
    conn.close()

def save_reads_mapped(mapping_file, sample_name, map_format, cursor, conn):
    """ This script reads a read mapping file, extracts the contig to which each read maps,
        the percent id, the start and stop, and stores it in a table per sample.
    
    Arguments:
        mapping_file {str} -- Location of read mapping file.
        sample_name {str} -- Name of database sample name (form sample_#)
        map_format {str} -- Format of read mapping results (blast or sam)
        cursor {obj} -- Cursor to execute db instructions
        conn {obj} -- Connection handler to db.
    """
    assert (map_format == "sam" or map_format == "bam" or map_format == "blast"), "Mapping format not recognized. Must be one of 'sam' 'bam' or 'blast'"
    contig_mag_corresp = {}
    contigs_in_sample = []
    # Retrieve mag information as mag_name, mag_id, contig_name, contig_id
    sql_command = 'SELECT * from lookup_table'
    cursor.execute(sql_command)
    contig_correspondence = cursor.fetchall()
    for contig_mag in contig_correspondence:
        contig_mag_corresp[contig_mag[2]] = [contig_mag[0], contig_mag[1], contig_mag[3]]
    
    # Read mapping files and fill sample tables
    if map_format == "blast":
        #print("Parsing tabular BLAST format reads... ", end = "", flush = True)
        with open(mapping_file) as input_reads:
            record_counter = 0
            records = []
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
                    contig_ref = segment[1]
                    # Exclude reads not associated with MAGs of interest
                    if contig_ref not in contig_mag_corresp:
                        continue
                    else:
                        if contig_ref not in contigs_in_sample:
                            contigs_in_sample.append(contig_ref)
                        pct_id = float(segment[2])
                        pos1 = int(segment[8])
                        pos2 = int(segment[9])
                        start = min(pos1, pos2)
                        end = start+(max(pos1, pos2)-min(pos1, pos2))
                        mag_id = contig_mag_corresp[contig_ref][1]
                        contig_id = contig_mag_corresp[contig_ref][2]
                        records.append((mag_id, contig_id, pct_id, start, end))
                        record_counter += 1
            # Commit remaining records
            if record_counter > 0:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
            # Create index for faster access
            cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
            
    if map_format == "sam":
        #print("Parsing SAM format reads... ", end = "", flush = True)
        record_counter = 0
        records = []
        with open(mapping_file) as input_reads:
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
                    # Exclude reads not associated with MAGs of interest
                    if contig_ref not in contig_mag_corresp:
                        continue
                    else:
                        if contig_ref not in contigs_in_sample:
                            contigs_in_sample.append(contig_ref)
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
                        mag_id = contig_mag_corresp[contig_ref][1]
                        contig_id = contig_mag_corresp[contig_ref][2]
                        records.append((mag_id, contig_id, pct_id, start, end))
                        record_counter += 1
            # Commit remaining records
            if record_counter > 0:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
            # Create index for faster access
            cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
    
    if map_format == "bam":
        #print("Parsing BAM format reads... ", end = "", flush = True)
        record_counter = 0
        records = []
        
        #This reader has some odd properties - most of it exists as a C interface
        #As a result, the entries are NOT the individual lines of the file and cannot be accessed as such
        #Instead, the iterator 'entry' returns a pointer to a location in memory based on the file
        #This iterator has a set of builtin functions that are called to access pos, ref name, MD:Z:
        input_reads = pysam.AlignmentFile(mapping_file, "rb")
        for entry in input_reads:
            #This line could allow processing to work like in SAM fmt, but is slower.
            #line = entry.to_string()
            if record_counter == 500000:
                cursor.execute("begin")
                cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
                cursor.execute("commit")
                record_counter = 0
                records = []
            #has_tag returns true if the entry has a %ID relevant field
            if not entry.has_tag("MD"):
                continue
            else :
                #No longer needed because of pysam accesses
                #segment = line.split()
                
                #The individual read has a reference ID, and the file has a list of names via IDs.
                #The entry.ref_ID gets the ID number, and the .get_ref returns the actual name from the number
                contig_ref = input_reads.get_reference_name(entry.reference_id)
                
                                
                # Exclude reads not associated with MAGs of interest
                if contig_ref not in contig_mag_corresp:
                    continue
                else:
                    if contig_ref not in contigs_in_sample:
                        contigs_in_sample.append(contig_ref)
                    
                    #Returns the MD:Z: segment
                    mdz_seg = entry.get_tag("MD")
                    match_count = re.findall('[0-9]+', mdz_seg)
                    sum=0
                    for num in match_count:
                        sum+=int(num)
                    total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
                    pct_id = (sum/(total_count))*100
                    
                    
                    #BAM files, unlike SAM files, are zero indexed. This +1 adjustment ensures SAM/BAM/R consistency
                    start = entry.reference_start+1
                    
                    
                    end = start+total_count-1
                    # Get mag_id and contig_id
                    mag_id = contig_mag_corresp[contig_ref][1]
                    contig_id = contig_mag_corresp[contig_ref][2]
                    
                    ##print(mag_id, contig_id)
                    
                    records.append((mag_id, contig_id, pct_id, start, end))
                    
                    ##print(*records)
                    
                    record_counter += 1
        # Commit remaining records
        if record_counter > 0:
            cursor.execute("begin")
            cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
            cursor.execute("commit")
        # Create index for faster access
        cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
    conn.commit()
    #print("done!")
    return contigs_in_sample


def add_sample(database, new_mapping_files, map_format):
    contig_mag_corresp = {}
    samples_dict = {}
    last_sample = 0
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    # Retrieve all sample information
    sql_command = 'SELECT * from sample_info'
    cursor.execute(sql_command)
    sample_information = cursor.fetchall()
    for sample in sample_information:
        samples_dict[sample[0]] = sample[1]
        if sample[2] > last_sample:
            last_sample = sample[2]
	
    # Retrieve contig - mag correspondence
    sql_command = 'SELECT * from lookup_table'
    cursor.execute(sql_command)
    contig_correspondence = cursor.fetchall()
    for contig_mag in contig_correspondence:
        contig_mag_corresp[contig_mag[2]] = [contig_mag[0], contig_mag[1], contig_mag[3]]
    for new_sample in new_mapping_files:
        # Check if new sample exists
        mags_in_sample = []
        if new_sample in samples_dict:
            #print("Dropping {} table and re-building it".format(new_sample))
            sample_name = samples_dict[new_sample]
            # If it does, drop the reads that table from that sample and re-build it
            cursor.execute('DROP TABLE IF EXISTS ' + sample_name)
            cursor.execute('CREATE TABLE ' + sample_name + \
                ' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
            cursor.execute('DELETE FROM mags_per_sample WHERE sample_name = ?', (new_sample,))
            conn.commit()
            #print("Adding {}... ".format(new_sample))
            contigs_in_sample = save_reads_mapped(new_sample, sample_name, map_format, cursor, conn)
            cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
            all_contigs = cursor.fetchall()
            for element in all_contigs:
                if element[0] in contigs_in_sample:
                    if element[1] not in mags_in_sample:
                        mags_in_sample.append(element[1])
                    else:
                        continue
                else:
                    continue
            mags_in_sample = [(new_sample, x) for x in mags_in_sample]
            cursor.execute("begin")
            cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
            cursor.execute("commit")

        else:
            # Otherwise create the new table and add the read information
            #print("Adding {} to existing database {}... ".format(new_sample, database))
            sample_name = "sample_" + str(last_sample + 1)
            last_sample += 1
            cursor.execute('CREATE TABLE ' + sample_name + \
                ' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
            contigs_in_sample = save_reads_mapped(new_sample, sample_name, map_format, cursor, conn)
            cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
            all_contigs = cursor.fetchall()
            for element in all_contigs:
                if element[0] in contigs_in_sample:
                    if element[1] not in mags_in_sample:
                        mags_in_sample.append(element[1])
                    else:
                        continue
                else:
                    continue
            mags_in_sample = [(new_sample, x) for x in mags_in_sample]
            cursor.execute("begin")
            cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
            cursor.execute("commit")
            # Add to sample_info table
            new_record = (new_sample, sample_name, last_sample)
            cursor.execute('INSERT INTO sample_info VALUES(?, ?, ?)', new_record)
            conn.commit()
        #print("Done")
        conn.commit()
    conn.close()


def parse_prodigal_genes(prodigal_gff):
    gene_info = []
    with open(prodigal_gff, 'r') as prodigal_genes:
        for line in prodigal_genes:
            if line.startswith("#"):
                continue
            else:
                line = line.strip().split()
                contig = line[0]
                start = min(int(line[3]), int(line[4]))
                end = max(int(line[3]), int(line[4]))
                strand = line[6]
                annotation = line[8]
                gene_id = annotation.split(";")[0].split("_")[1]
                gene_id = contig + "_" + gene_id
                gene_info.append((contig, gene_id, start, end, strand, annotation))
    return gene_info


def add_gene_information(database, gene_info):
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute('DROP TABLE IF EXISTS gene_info')
    cursor.execute('CREATE TABLE gene_info \
        (contig_name TEXT, gene_name TEXT, gene_start INTEGER, gene_stop INTEGER, strand TEXT, annotation TEXT)')
    cursor.execute("begin")
    cursor.executemany('INSERT INTO gene_info VALUES(?, ?, ?, ?, ?, ?)', gene_info)
    cursor.execute("commit")


def add_gene_annotation(database, annotation):
    annotations = []
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    cursor.execute('DROP TABLE IF EXISTS gene_annotation')
    cursor.execute('CREATE TABLE gene_annotation \
        (gene_name TEXT, annotation TEXT)')
    with open(annotation) as gene_annot:
        for line in gene_annot:
            line = line.strip().split("\t")
            annotations.append((line[0], line[1]))
    cursor.execute("begin")
    cursor.executemany('INSERT INTO gene_annotation VALUES(?, ?)', annotations)
    cursor.execute("commit")

def read_contigs(contig_file_name):
    """ Reads a FastA file and returns
        sequence ids and sizes
    
    Arguments:
        contig_file_name {[str]} -- FastA file location
    Returns:
        [dict] -- Dictionary with ids and sizes
    """
    #print("Reading contigs... ", end="", flush=True)
    
    contig_sizes = {}
    contig_length = 0
    contigs =  open(contig_file_name, 'r')
    
    #The ensuing loop commits a contig to the contig lengths dict every time a new contig is observed, i.e. whenever a current sequence has terminated.
    #This works both for single and splitline multi fastas.
    
    #The first line is manually read in so that its count can be gathered before it is committed - this basically skips the first iteration of the loop.
    current_contig = contigs.readline()[1:].strip().split()[0]
    
    for line in contigs:
        if line[0] == ">":
            #Add the contig that had 
            contig_sizes[current_contig] = contig_length
            
            #set to new contig. One final loop of starts ends counts is needed
            current_contig = line[1:].strip().split()[0]
            contig_length = 0
        else :
            contig_length += len(line.strip())
    
    contigs.close()
    
    #The loop never gets to commit on the final iteration, so this statement adds the last contig.
    contig_sizes[current_contig] = contig_length
    
    #print("done!")
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

#The purpose of this function is to take an empty recplot matrix object associated with one MAG, query the database for the sample and MAG in question,
#And fill the database with the returned information.
def fill_matrices(database, mag_id, sample_name, matrices, id_breaks):
    """[summary]
    
    Arguments:
        database {str} -- Name of database to use (location).
        mag_id {int} -- ID of mag of interest.
        sample_name {str} -- Name of sample to retrieve reads from.
        matrices {dict} -- Empty matrices to fill.
        id_breaks {list} -- List of identity percentages to include.
    
    Returns:
        matrix [dict] -- Dictionary with list of arrays of start and stop positions
                         and filled matrix to plot.
    """
    #print("Filling matrices...", end = "",flush = True)
    # Retrieve sample_id from sample_name provided
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT sample_id from sample_info WHERE sample_name = ?'
    cursor.execute(sql_command, (sample_name,))
    sample_id = cursor.fetchone()[0]
    # Retrieve all read information from mag_name and sample_name provided
    sql_command = 'SELECT * from ' + sample_id + ' WHERE mag_id = ?'
    cursor.execute(sql_command, (mag_id,))
	
	#TODO: Consider keeping/removing this in final...
    #read_counter = 0
    
    #TODO: We shouldn't need to fetch all reads. We can iterate on the cursor without doing this.
    #read_information = cursor.fetchall()
    #read_information is (mag_id contig_id perc_id read_start read_stop)
    #for read_mapped in read_information:
    for read_mapped in cursor:
        if read_mapped[1] in matrices:
            #read_counter += 1
            contig_id = read_mapped[1]
            read_start = read_mapped[3]
            read_stop = read_mapped[4]
            read_len = read_stop - read_start + 1
            read_id_index = bisect.bisect_right(id_breaks, read_mapped[2]) - 1
            # #print(read_mapped)
            # #print(id_breaks)
            # #print(read_id_index, id_breaks[read_id_index])
            read_start_loc = bisect.bisect_left(matrices[contig_id][1], read_start)
            read_stop_loc = bisect.bisect_left(matrices[contig_id][1], read_stop)
            # If the read falls entirely on a bin add all bases to the bin
			
            if read_start_loc == read_stop_loc:
                matrices[contig_id][2][read_start_loc][read_id_index] += read_len
            # On the contrary split bases between two or more bins
            else:
                for j in range(read_start_loc, read_stop_loc + 1):
                    overflow = read_stop - matrices[contig_id][1][j]
                    if overflow > 0:
                        matrices[contig_id][2][j][read_id_index] += (read_len - overflow)
                        read_len = overflow
                    else :
                        matrices[contig_id][2][j][read_id_index] += read_len
    #print("done!")
	##print("done!", read_counter, "reads collected!")
    return matrices

#The purpose of this function is to prepare an empty recplot matrix from a set of contig names and lengths associated with one MAG
def prepare_matrices(database, mag_name, width, bin_height, id_lower):
    #Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
    #print("Preparing recruitment matrices...", end="", flush=True)
    # Prep percent identity breaks - always starts at 100 and proceeds 
    # down by bin_height steps until it cannot do so again without passing id_lower
    id_breaks = []
    current_break = 100
    while current_break > id_lower:
        id_breaks.append(current_break)
        current_break -= bin_height
    id_breaks = id_breaks[::-1]
    
    zeroes = []
    for i in id_breaks:
        zeroes.append(0)
    
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
    #begin reading contigs and determining their lengths.
    
    #id_len is a list of contig name, contig_length
    for id_len in contig_sizes:
        
        starts = []
        ends = []
        pct_id_counts = []
        
        contig_length = id_len[1]
        
        num_bins = int(contig_length / width)
        if num_bins < 1:
            num_bins = 1
        
        bin_width = (contig_length / num_bins)-1

        cur_bin_start = 1
        
        for i in range(1, num_bins):
            starts.append(int(cur_bin_start))
            ends.append(int((cur_bin_start+bin_width)))
            pct_id_counts.append(zeroes[:])
            cur_bin_start+=bin_width+1
		
		#Append final bin; guarantees the final bin = contig length
        starts.append(int(cur_bin_start))
        ends.append(contig_length)
        pct_id_counts.append(zeroes[:])
        
        matrix[id_len[0]] = [starts, ends, pct_id_counts]
        
    #print("done!")

    return(mag_id, matrix, id_breaks)

def prepare_matrices_genes(database, mag_name, bin_height, id_lower):
    #Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
    #print("Preparing recruitment matrices...", end="", flush=True)
    # Prep percent identity breaks - always starts at 100 and proceeds 
    # down by bin_height steps until it cannot do so again without passing id_lower
    id_breaks = []
    current_break = 100
    while current_break > id_lower:
        id_breaks.append(current_break)
        current_break -= bin_height
    id_breaks = id_breaks[::-1]
    
    zeroes = []
    for i in id_breaks:
        zeroes.append(0)
    
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
	
	#Contig ID is used above, but contig names will be needed to translate from the genes table and the contig data
    contig_curs = get_contig_names(database, mag_name)
    contig_names = []
	
    for item in contig_curs:
        contig_names.append(item[0])
			
	#Retrieve all gene info - order listed
    sql_command = 'SELECT contig_name, gene_name, gene_start, gene_stop, strand, annotation from gene_info WHERE contig_name IN (SELECT contig_name from lookup_table WHERE mag_id = ?)'
    cursor.execute(sql_command, (mag_id,))
    gene_sizes = cursor.fetchall()
	
    gene_matrix = {}
	
	#Separates gene data into contigs to allow for access by contig during matrix creation.
    for item in gene_sizes:
        if item[0] not in gene_matrix:
            gene_matrix[item[0]] = [[item[1]], [item[2]], [item[3]], [item[4]], [item[5]]]
        else:
            gene_matrix[item[0]][0].append(item[1])
            gene_matrix[item[0]][1].append(item[2])
            gene_matrix[item[0]][2].append(item[3])
            gene_matrix[item[0]][3].append(item[4])
            gene_matrix[item[0]][4].append(item[5])
	
	#Final data structures initialization
    matrix = {}
    annotation_matrix = {}
	
    #Since genes are separated, id_len can be iterated through to give a length of contig for each and all the assoc. gene data can be accessed by name
    
    #id_len is a list of contig name, contig_length
    for id_len in contig_sizes:
        
		#These are the starts and ends of each gene on this contig
        starts = gene_matrix[str(contig_names[id_len[0]-1])][1][:]
        ends = gene_matrix[str(contig_names[id_len[0]-1])][2][:]
		
        gene_names = gene_matrix[str(contig_names[id_len[0]-1])][0][:]
        strand = gene_matrix[str(contig_names[id_len[0]-1])][3][:]
        annotation = gene_matrix[str(contig_names[id_len[0]-1])][4][:]
		
        final_gene_names = []
        final_gene_strands = []
        final_gene_annots = []
		#I shove around the starts/ends of overlapping genes to make sure that they no longer do, but that also means that 
		#the start/stop used for plotting is not guaranteed to be the real start/stop. I save the true ones here.
        final_gene_starts = []
        final_gene_ends = []

        pct_id_counts = []
        
        contig_length = id_len[1]
		
        final_starts = [1]
        final_ends = []
		
		#We want to see the genes in the context of the contig, so we have to add bins for all intergenic regions
		#This loop starts with 1 and creates a bin end based on the next gene start, then adds the gene's start and end, then adds a new start based on the gene's end
        for i in range(0, len(starts)):
			#We have to pad out the annotation information with blanks for the intergenic bins - 
			#determining their locations is easiest done alongside the initial creation of the bins rather than awkwardly afterwards
			
			#blanks
            final_gene_names.append("N/A")
            final_gene_strands.append("N/A")
            final_gene_annots.append("N/A")
            final_gene_starts.append("N/A")
            final_gene_ends.append("N/A")
			#actual gene info
            final_gene_names.append(gene_names[i])
            final_gene_strands.append(strand[i])
            final_gene_annots.append(annotation[i])
            final_gene_starts.append(str(starts[i]))
            final_gene_ends.append(str(ends[i]))
			#starts and stops added here
            final_ends.append(starts[i]-1)
            final_starts.append(starts[i])
            final_ends.append(ends[i])
            final_starts.append(ends[i]+1)
		
		#Caps the ends with the end of the contig
        final_ends.append(contig_length)
		
		#Finishes out the annotaion list - the last bin is always declared as non-genic, and it gets removed if that's untrue.
        final_gene_names.append("N/A")
        final_gene_strands.append("N/A")
        final_gene_annots.append("N/A")
        final_gene_starts.append("N/A")
        final_gene_ends.append("N/A")
		
		#If a gene starts at 1, then there should be a 1, 0 pair in the first position of the final starts/ends and the bin is unnecessary. This removes those
        if final_ends[0] < 1:
            del final_starts[0]
            del final_ends[0]
            del final_gene_annots[0]
            del final_gene_names[0]
            del final_gene_strands[0]
            del final_gene_starts[0]
            del final_gene_ends[0]		
			
		#If a gene terminates at the end of the contig, then the final start would be contig length + 1 and the bin is unnecessary. This removes those
        if final_starts[len(final_starts)-1] > final_ends[len(final_starts)-1]:
            del final_starts[len(final_starts)-1]
            del final_ends[len(final_starts)-1]
            del final_gene_annots[len(final_starts)-1]
            del final_gene_names[len(final_starts)-1]
            del final_gene_strands[len(final_starts)-1]
            del final_gene_starts[len(final_starts)-1]
            del final_gene_ends[len(final_starts)-1]

			
		
        starts = []
        ends = []
        gene_names = []
        strand = []
        annotation = []
        annot_start = []
        annot_end = []
		#If genes overlap, two problems are created:
		# (1) Intergenic region will have end >= start
		# (2) The end of the first overlapping gene will be <= the start of the second
		#This removes the intergenic regions and shoves the start/ends back roughly equal distances when they overlap - marginally inaccurate, but necessary to maintain histogram behavior
		#Has to be added to new lists because I cannot modify length of final starts/final ends inside of loop
        for i in range(0, len(final_starts)):
            if final_ends[i] >= final_starts[i]:
                starts.append(final_starts[i])
                ends.append(final_ends[i])
                gene_names.append(final_gene_names[i])
                strand.append(final_gene_strands[i])
                annotation.append(final_gene_annots[i])
                annot_start.append(final_gene_starts[i])
                annot_end.append(final_gene_ends[i])				
            else:
                midpt = int((final_ends[i-1] + final_starts[i+1])/2)
                #last element added was final_ends[i-1], and this is what needs updated
                ends[len(ends)-1] = midpt
				#This will naturally be added next loop, which should always happen because there will always be 1 more index to iterate over
                final_starts[i+1] = midpt + 1
				
        
				

		#Now add zeroes.
        pct_id_counts = []
        for index in starts:
            pct_id_counts.append(zeroes[:])
        
		#Add them to the data item
        matrix[id_len[0]] = [starts, ends, pct_id_counts]
        annotation_matrix[id_len[0]] = [gene_names, annot_start, annot_end, strand, annotation]

        
    #print("done!")
    return(mag_id, matrix, id_breaks, annotation_matrix)
	
#This function orchestrates calls to prepare_matrices and fill_matrices. 
#Translations between R and python through reticulate are inefficient and may cause errors with data typing.
#Constraining the transfers to arguments passed from R and a return passed from python alleviates these issues.
def extract_MAG_for_R(database, sample, mag_name, width, bin_height, id_lower):
    #print("Making recruitment matrix for:", mag_name, "in sample:", sample)
    mag_id, matrix, id_breaks = prepare_matrices(database, mag_name, width, bin_height, id_lower)
    
    matrix = fill_matrices(database, mag_id, sample, matrix, id_breaks)
    
    return(matrix, id_breaks)
    
def extract_genes_MAG_for_R(database, sample, mag_name, bin_height, id_lower):
    #print("Making recruitment matrix for:", mag_name, "in sample:", sample)
    mag_id, matrix, id_breaks, gene_table = prepare_matrices_genes(database, mag_name, bin_height, id_lower)
    
    matrix = fill_matrices(database, mag_id, sample, matrix, id_breaks)
    
    return(matrix, id_breaks, gene_table)
    
	
#This function queries the database and returns the names of all of the samples present within it
def assess_samples(database):
    #print("Acquiring samples in "+ database)
	# Retrieve sample_id from sample_name provided
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT sample_name from sample_info'
    cursor.execute(sql_command)
    
    samples = []
	
    for samp in cursor:
        samples.append(samp)
		
    cursor.close()
	
    return(samples)

#This function queries the database and returns the MAGs covered by a given sample.    
def assess_MAGs(database, sample):
    #print("Acquiring MAGs in "+sample)
	# Retrieve sample_id from sample_name provided
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT mag_name from mags_per_sample WHERE sample_name = ?'
    cursor.execute(sql_command, (sample,))
    
    mags = []
	
    for mag in cursor:
        mags.append(mag)
		
    cursor.close()
	
    return(mags)
	
def get_contig_names(database, mag_name):
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = 'SELECT contig_name from lookup_table WHERE mag_name = ?'
    cursor.execute(sql_command, (mag_name,))
    
    names = []
	
    for contig in cursor:
        names.append(contig)
		
    cursor.close()
	
    return(names)
	
def add_genes_to_db(database, genes_file, gene_format):
    if(gene_format == "prodigal"):
        gene_information = parse_prodigal_genes(genes_file)
    else:
        #print("I don't do that yet")
	
    add_gene_information(database, gene_information)
#add_gene_annotation(database, annotation)

def check_presence_of_genes(database):
    conn = sqlite3.connect(database)
    cursor = conn.cursor()
    sql_command = "SELECT name FROM sqlite_master WHERE type='table'"
    cursor.execute(sql_command,)
	
    checker = False
	
    for name in cursor:
        if str(name) == "('gene_info',)":
            checker = True
	
    cursor.close()
    return(checker)



	
#A function for reading args
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
            description='''This script builds recruitment plots (COMPLETE DESCRIPTION HERE)\n'''
            '''Usage: ''' + argv[0] + ''' COMPLETE\n'''
            '''Global mandatory parameters: -g [Genome Files] OR -p [Protein Files] OR -s [SCG HMM Results] -o [AAI Table Output]\n'''
            '''Optional Database Parameters: See ''' + argv[0] + ' -h')

    parser.add_argument("-c", "--contigs", dest="contigs", 
    help = "This should be a FASTA file containing all and only the contigs that you would like to be part of your recruitment plot.")
    parser.add_argument("-m", "--mags", dest="mags", default="", help = "A tab separated file containing the names of contigs in the first column and MAGs in the second column. Every contig should have its parent MAG listed in this file.")    
    parser.add_argument("-r", "--reads", dest="reads", nargs='+', help = "This should be a file with reads aligned to your contigs in any of the following formats: tabular BLAST(outfmt 6), SAM, or Magic-BLAST")
    parser.add_argument("-f", "--format", dest="map_format", default="blast", help="The format of the reads file (write 'blast' or 'sam'). Defaults to tabular BLAST.")
    parser.add_argument("-g", "--genes", dest="genes", default = "", help = "Optional GFF3 file containing gene starts and stops to be use in the recruitment plot.")
    parser.add_argument("-a", "--annot", dest="annotation", default = "", help = "Optional file with gene name in the first column and annotations in the second column.")
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
    genes = args.genes
    annotation = args.annotation
    step = float(args.id_step)
    width = int(args.bin_width)
    prefix = args.out_file_name
    export_lines = args.output_line
    do_stats = args.stats
    interact = args.lim_rec
    sql_database = args.sql_database
    
    # Create databases
    #sqldb_creation(contigs, mags, reads, map_format, sql_database)

    # Prepare user requested information
    # mag_id, matrix, id_breaks = prepare_matrices(sql_database, "IIa.A_ENTP2013_S02_SV82_300m_MAG_01", width, step, 70)
    # matrix = fill_matrices(sql_database, mag_id, "03.All_SAR11--ETNP_2013_S02_SV89_300m.blast.bh", matrix, id_breaks)

    # Add new sample to database
    # add_sample(sql_database, ["03.First_Mapping.blast.bh", "TEST"], map_format)

    # Add gene information
    #gene_information = parse_prodigal_genes(genes)
    #add_gene_information(sql_database, gene_information)
    #add_gene_annotation(sql_database, annotation)

    
   

#Just runs main.
#if __name__ == "__main__":main()

