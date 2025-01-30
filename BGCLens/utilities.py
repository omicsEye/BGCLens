# from __future__ import print_function
"""
BGCLens utilities
"""
import re
import sys
import os
import shutil
import csv
import pdb
import tarfile
from urllib.request import urlretrieve

import Bio
import gzip
import logging
import subprocess
import tempfile
import traceback
import bz2
from datetime import datetime
from matplotlib import pyplot as plt

from Bio import Seq
import pandas as pd
import config

import datetime
import time

# name global logging instance
logging.basicConfig()
logger = logging.getLogger(__name__)

# constants
# ---------------------------------------------------------------

c_strat_delim = "|"
c_taxon_delim = "."
c_name_delim = ": "
c_multiname_delim = ";"
c_str_unknown = "NO_NAME"
c_ungrouped = "UNGROUPED"
c_unmapped = "UNMAPPED"
c_unintegrated = "UNINTEGRATED"
c_many_bytes = 1e8
c_zip_multiplier = 10
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE = "#"

# the extension used for biom files
BIOM_FILE_EXTENSION = ".biom"

c_topsort = {
    c_unmapped: 0,
    c_ungrouped: 1,
    c_unintegrated: 2,
    "UniRef50_unknown": 3,
    "UniRef90_unknown": 4,
}


# from zopy.utils import iter_rows, tprint, warn

# ==== Helper function from zopy by Eric Franzosa ========
def warn(*args):
    script = "?"
    if sys.argv[0] != "":
        script = os.path.split(sys.argv[0])[1].upper()
    args = ["WARNING ({}):".format(script)] + list(args)
    print >> sys.stderr, " ".join(map(str, args))


def iter_rows(path):
    """ easy table loading """
    lens = []
    with try_open(path) as fh:
        for row in reader(fh):
            lens.append(len(row))
            yield row
    if len(set(lens)) != 1:
        warn("rows didn't all have equal lengths:", set(lens))


def tprint(*args, **kwargs):
    """ coerce list of items to strings then print with tabs between """
    print >> kwargs.get("file", sys.stdout), "\t".join(map(str, args))


def try_open(path, *args):
    """ open an uncompressed or gzipped file; fail gracefully """
    fh = None
    try:
        if re.search(r".gz$", path):
            print >> sys.stderr, "Treating", path, "as gzipped file"
            fh = gzip.GzipFile(path, *args)
        else:
            fh = open(path, *args)
    except:
        sys.exit("Problem opening:" + path)
    return fh


# ---------------------------------------------------------------
# text manipulation
# ---------------------------------------------------------------

def reader(file_handle):
    """ my favorite options for csv reader """
    for aItems in csv.reader(file_handle, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE):
        yield aItems


def make_directory(output_dir, force=False):
    if not os.path.isdir(output_dir):
        try:
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
    elif force:
        try:
            print("Removing the old output directory: " + output_dir)
            shutil.rmtree(output_dir)
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
    else:
        print(("Output directory %s exists and files will be overwrite...") % (output_dir))

    if not os.access(output_dir, os.W_OK):
        sys.exit("CRITICAL ERROR: The output directory is not " +
                 "writeable. This software needs to write files to this directory.\n" +
                 "Please select another directory.")


def read_ppanini_imp_genes_table_dead(filename):
    gene_table = pd.read_csv(filename, sep='\t', index_col=0)
    # print gene_table.columns.values
    # print gene_table['abundance']
    ppanini_table = {'genes': list(gene_table.index), 'abundance_rank': gene_table['abundance'],
                     'prevalence_rank': gene_table['prevalence'], 'ppanini_score': gene_table['ppanini_score']}
    return ppanini_table


def read_fasta(fasta_filename):
    '''Reads a fasta_file and returns a fasta dict
    Input: fasta_filename = path_to_fasta_file
    Output: fasta_seq = {sequence_header: sequence, ...}'''
    logger.debug('read_fasta ' + fasta_filename)

    fasta_file = open(fasta_filename)
    fasta_seq = {}
    name = ''

    for line in fasta_file:
        if not line.startswith('#'):
            if line.startswith('>'):
                name = re.sub('>', '', line.split(' ')[0].strip())
            else:
                if name not in fasta_seq:
                    fasta_seq[name] = re.sub('[\r\t\n]', '', line)
                else:
                    fasta_seq[name] += re.sub('[\r\t\n]', '', line)
    # print fasta_seq
    return fasta_seq


def parse_table(filename, fasta_filename):
    '''Parse the BLAST results to give gene hits to genomes
    Input:
    filename = filename of blast results
    fasta_filename = filename of corresponding fasta file

    Output:
    table = {gene: [List of genomes]}'''

    fasta_dict = read_fasta(fasta_filename)

    for seq in fasta_dict:
        fasta_dict[seq] = float(len(fasta_dict[seq]))
    table = {}
    foo = open(filename)
    for line in foo:
        split_i = line.split('\t')
        try:
            threshold = float(split_i[2]) * float(split_i[3]) / fasta_dict[split_i[0]]
        except:
            raise Exception('Gene ' + split_i[0] + ' not found in fasta: ' + fasta_filename)
        if threshold > 90.0:
            split_sp = split_i[1].split('|')
            sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
            if split_i[0] not in table:
                table[split_i[0]] = sp
            elif sp not in table[split_i[0]]:
                table[split_i[0]] += sp
    return table


def count_genomes(filename):
    table = pd.DataFrame.from_csv(filename, sep='\t', index_col=None, header=None)
    table.columns = ['gene', 'genome']
    gene_count_genome = table['gene'].value_counts()
    return gene_count_genome


def gene2genomes(filename):
    '''Read parsed table for {gene: genomes}
    Input:
    filename = filename of blast results

    Output:
    table = {gene: [List of genomes]}'''

    table = {}

    foo = open(filename)

    for line in foo:
        split_i = [i.strip() for i in line.split('\t')]
        try:
            table[split_i[0]] += [split_i[1]]
        except:
            table[split_i[0]] = [split_i[1]]
    return table


def number_of_unique_genomes(mg_file):
    metagenomic_table = pd.DataFrame.from_csv(mg_file, sep='\t', index_col=None, header=None)
    metagenomic_table.columns = ['gene', 'genome']
    uniq_genomes = []
    for gene in metagenomic_table:
        for genome in metagenomic_table[gene]:
            if genome not in uniq_genomes:
                uniq_genomes += [genome]
    no_uniq_genomes = len(uniq_genomes)
    # print 'No. of unique genomes: '+str(no_uniq_genomes)
    # ppanini_output = pd.DataFrame.from_csv(ppanini_output_file, sep='\t', index_col=0, header =0)#read_ppanini_imp_genes_table(ppanini_output_file)
    return no_uniq_genomes


def read_data(mg_file, ppanini_output_file):
    metagenomic_table = pd.DataFrame.from_csv(mg_file, sep='\t', index_col=None, header=None)
    metagenomic_table.columns = ['gene', 'genome']
    uniq_genomes = []
    for gene in metagenomic_table:
        for genome in metagenomic_table[gene]:
            if genome not in uniq_genomes:
                uniq_genomes += [genome]
    no_uniq_genomes = len(uniq_genomes)
    print('No. of unique genomes: ' + str(no_uniq_genomes))
    ppanini_output = pd.DataFrame.from_csv(ppanini_output_file, sep='\t', index_col=0,
                                           header=0)  # read_ppanini_imp_genes_table(ppanini_output_file)
    return metagenomic_table, ppanini_output, no_uniq_genomes


def read_abund_prev(filename):
    keys = {'abund': 0, 'alpha': 0, 'beta': 0}
    abund = []
    alpha = []
    genes = []
    with open(filename) as foo:
        for line in foo:
            split_line = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')]
            if line.startswith('Centroids'):
                for i, val in enumerate(split_line):
                    if 'abund' in val:
                        keys['abund'] = i
                    elif 'alpha' in val:
                        keys['alpha'] = i
                    elif 'beta' in val:
                        keys['beta'] = i
            else:
                # pdb.set_trace()
                genes += [split_line[0]]
                abund += [float(split_line[keys['abund']])]
                alpha += [float(split_line[keys['alpha']])]
    abund_prev = {'genes': genes, 'abundance': abund, 'prevalence': alpha}
    return abund_prev


def pullgenes_fromcontigs(contig_file, gff3_file, fna_file, faa_file):
    '''Pulls genes from contigs using the coordinates from GFF3 file provided

    Input: contig_file: FASTA file containing contigs in FNA format.
           gff3_file: File containing coordinates of genes in GFF3 format
           fna_file: filepath for genes written in nucleotides sequences
           faa_file: filepath for genes written in amino-acid sequences'''
    logger.debug('pullgenes_fromcontigs ' + contig_file + ' ' + gff3_file)

    gene_contig_mapper = {}
    contig_gene_mapper = {}
    gene_start_stop = {}

    [gene_contig_mapper, gene_start_stop, contig_gene_mapper] = read_gff3(gff3_file)

    genes_fasta = {}
    contigs_fasta_dict = read_fasta(contig_file)
    for contig in contigs_fasta_dict:
        if contig in contig_gene_mapper:
            for gene in contig_gene_mapper[contig]:
                [start_x, stop_x, strand] = gene_start_stop[gene]
                try:
                    if strand == '+':
                        genes_fasta[gene] = contigs_fasta_dict[contig][start_x - 1:stop_x + 1]
                    else:
                        contig_len = len(contigs_fasta_dict[contig])
                        start_minus = -1 * stop_x
                        stop_minus = -1 * (start_x - 1)
                        genes_fasta[gene] = Bio.Seq.reverse_complement(contigs_fasta_dict[contig])[
                                            start_minus:stop_minus]
                except:
                    raise Exception('Circular DNA Detected')
                    new_stop_x = stop_x - len(contigs_fasta_dict[contig]) + 1
                    genes_fasta[gene] = contigs_fasta_dict[contig] + contigs_fasta_dict[contig][:new_stop_x]

    write_fasta(genes_fasta, fna_file, False)  # FNA
    write_fasta(genes_fasta, faa_file, True)  # FAA


def read_gff3(filename):
    '''Reads GFF3 files and returns the relevant information
    Input: filename = path_to_gff3_file
    Output: [gene_contig_mapper, gene_start_stop, contig_gene_mapper]
    gene_contig_mapper = {gene: contig, ...}
    gene_start_stop = {gene: [start, stop, strand], ...}
    contig_gene_mapper = {contig: [List of genes], ...}'''

    logger.debug('read_gff3 ' + filename)

    gene_contig_mapper = {}
    contig_gene_mapper = {}
    gene_start_stop = {}

    with open(filename, 'r') as foo:
        for line in foo:
            if re.match('(\w+)\t(\w+)\tgene\t(\w+)', line):
                split_line = [re.sub('[\r\t\n]', '', i).strip() for i in line.split('\t')]
                gid = split_line[-1].split('=')[-1]
                gene_contig_mapper[gid] = split_line[0]
                gene_start_stop[gid] = [int(split_line[3]), int(split_line[4]), split_line[6]]
                if split_line[0] in contig_gene_mapper:
                    contig_gene_mapper[split_line[0]] += [gid]
                else:
                    contig_gene_mapper[split_line[0]] = [gid]
    return [gene_contig_mapper, gene_start_stop, contig_gene_mapper]


def write_fasta(seqs_dict, filename, to_translate):
    '''Writes dictionary of fasta sequences into text file
    Input: seqs_dict = {geneID: sequence}
           filename = path_to_output_genes_fastafile'''
    logger.debug('write_fasta to ' + filename)

    with open(filename, 'w') as foo:

        test = ''
        for i in seqs_dict:
            test = seqs_dict[i]
            break
        format = is_protein(test)
        to_translate = to_translate and not format

        if to_translate:  # if not FAA already and to be translated
            for seq in seqs_dict:
                t_seq = Bio.Seq.translate(seqs_dict[seq], to_stop=True)
                len_seq = str(len(t_seq) * 3)
                foo.writelines(['>' + seq + '|' + len_seq + '\n'])
                foo.writelines([t_seq + '\n'])
        else:
            ind = 1
            if format:
                ind = 3  # amino acids * 3 nucleotides
            for seq in seqs_dict:
                len_seq = str(len(seqs_dict[seq] * ind))
                foo.writelines(['>' + seq + '|' + len_seq + '\n'])
                foo.writelines([seqs_dict[seq] + '\n'])


def is_protein(sequence):
    '''Returns True if the sequence is protein.
    Input: (str) format sequence
    Output: (bool) True if amino acids; False if nucleotides sequence'''

    logger.debug('is_protein')

    format = False  # 'FNA'
    try:
        Bio.Seq.translate(sequence)
    except:
        format = True  # 'FAA'
    return format


def create_folders(list_folders):
    '''Creates the list of folders if they dont exist already
    Input: list_folders = [List of folders to be created]

    Output: Folders created'''

    logger.debug('create_folders ' + '\t'.join(list_folders))

    for fname in list_folders:
        try:
            os.stat(fname)
        except:
            os.mkdir(fname)


def read_dict(gene_annotations_file):
    '''Reads tabulated file into a dictionary

    Input: gene_annotations_file = path_to_output_gene_annotations_table

    Output: dictX = {geneID: annotation}'''

    logger.debug('read_dict ' + gene_annotations_file)

    dictX = {}

    with open(gene_annotations_file) as foo:
        for line in foo:
            if not line.startswith('#'):
                split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
                dictX[split_line[0]] = split_line[1]
    return dictX


def read_dict_num(gene_annotations_file):
    '''Reads tabulated file into a dictionary

    Input: gene_annotations_file = path_to_output_gene_annotations_table

    Output: dictX = {geneID: annotation}'''

    logger.debug('read_dict)num ' + gene_annotations_file)

    dictX = {}

    with open(gene_annotations_file) as foo:
        for line in foo:
            if not line.startswith('#'):
                split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split('\t')]
                if len(split_line) < 2:
                    split_line = [re.sub('[\t\r\n]', '', i).strip() for i in line.split(' ')]
                try:
                    dictX[split_line[0]] = float(split_line[1])
                except:
                    pdb.set_trace()
    return dictX


def write_dict(dictX, gene_annotations_file):
    '''Writes dictionary of genes and their annotations into text file
    Input: dictX = {geneID: annotation}
           gene_annotations_file = path_to_output_gene_annotations_table'''

    logger.debug('write_dict to ' + gene_annotations_file)

    with open(gene_annotations_file, 'w') as foo:
        foo.writelines('#GENEID\tANNOTATION\n')
        for i in dictX:
            foo.writelines(['\t'.join([i, dictX[i]]) + '\n'])


def is_present(metadata, meta_type):
    '''Returns True if meta_type is present in metadata extracted from mappert_file

    Input: metadata = [metadata strings]; Rows with # as first character in table
           meta_type = Type of metadata that you are querying e.g. FASTAS, NICHE etc.

    Output: [line, ind]
            line = The corresponding line from metadata, [] if not present
            ind = The index of the line in the metadata sequence, [] if not present'''

    logger.debug('is_present ' + meta_type)

    line = []
    ind = []
    for i, val in enumerate(metadata):
        if val.upper().startswith(meta_type):
            line = val
            ind = i
            break
    return [line, ind]


def generate_gene_table(abundance_dict, annotations_dict, niche_flag, mapper, output_table, samples):
    '''Input: abundance_dict: {sample:{gene:abundance,...},...}
              annotation_dict: {gene: annotation},
              niche_flag: True if NICHE exists
              mapper'''
    logger.debug('generate_gene_table ' + output_table)
    #	samples = abundance_dict.keys()
    # fasta_row = [mapper[i]['FAAS'] for i in samples]

    with open(output_table, 'w') as foo:
        if niche_flag:
            niche_row = [mapper[i]['NICHE'] for i in samples]
            foo.writelines([str.join('\t', ['#NICHE'] + niche_row) + '\n'])  # header
        # foo.writelines([str.join('\t', ['#FAAS']+fasta_row)+'\n']) #header
        foo.writelines([str.join('\t', ['#SAMPLES'] + samples) + '\n'])  # header

        # for i, sample in enumerate(samples):
        for gene in abundance_dict:
            # abund_x_i = abundance_dict[sample][gene]
            # data_row = numpy.zeros(len(samples))
            # data_row[i] = abund_x_i
            str_data_row = [str(ele) for ele in abundance_dict[gene]]
            if gene in annotations_dict:
                annot_x_i = annotations_dict[gene]
                if annot_x_i.startswith('UniRef90'):
                    umap_i = 'UniRef50_unknown'
                    annot_x = annot_x_i + '|' + umap_i
                else:
                    annot_x = 'UniRef90_unknown|' + annot_x_i
            else:
                annot_x = 'UniRef90_unknown|UniRef50_unknown'
            foo.writelines([str.join('\t', [gene + '|' + annot_x] + str_data_row) + '\n'])


def extract_mapping():
    '''Uses centroids from first file and
    produces mapped centroids to genes through mapper'''
    foo = open(sys.argv[1])
    centroids = [re.sub('[\t\r\n]', '', i).strip() for i in foo]

    mapper = open(sys.argv[2])
    mapping = {}
    for line in mapper:
        split_line = [re.sub('[\r\n\t]', '', i) for i in line.split('\t')]
        try:
            mapping[split_line[0]] += [split_line[1]]
        except:
            mapping[split_line[0]] = [split_line[1]]

    for i in centroids:
        if i in mapping:
            for j in mapping[i]:
                print('\t'.join([i, j]))
        else:
            print(i + '\tNA')


# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------
def read_map(map_obj):
    csv_map = csv.reader(open(map_obj), csv.excel_tab)
    uniref_go = {}
    for line in csv_map:
        for uid in line[1:]:
            uniref_go[uid.strip()] = line[0]
    return uniref_go


def size_warn(path):
    m = 1 if ".gz" not in path else c_zip_multiplier
    if m * os.path.getsize(path) > c_many_bytes:
        print("  This is a large file, one moment please...", file=sys.stderr)


# def try_zip_open(path, write=None):
#     """
#     open an uncompressed or gzipped file; fail gracefully
#     """
#     fh = None
#
#     # set the open mode
#     if write:
#         open_mode = "w"
#     elif path.endswith(".bz2"):
#         open_mode = "r"
#     else:
#         open_mode = "rt"
#
#     try:
#         if path.endswith(".gz"):
#             fh = gzip.open(path, open_mode)
#         elif path.endswith(".bz2"):
#             fh = bz2.BZ2File(path, open_mode)
#         else:
#             fh = open(path, open_mode)
#     except EnvironmentError:
#         sys.exit("Problem opening file: " + path)
#     return fh


def read_biom_table(path):
    """
    return the lines in the biom file
    """

    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software." +
                 " This software is required since the input file is a biom file.")

    try:
        tsv_table = biom.load_table(path).to_tsv().split("\n")
    except (EnvironmentError, TypeError):
        sys.exit("ERROR: Unable to read biom input file.")

    return tsv_table


def try_zip_open( path, write=None ):
    """
    open an uncompressed or gzipped file; fail gracefully
    """
    fh = None

    # set the open mode
    if write:
        open_mode = "w"
    elif path.endswith(".bz2"):
        open_mode = "r"
    else:
        open_mode = "rt"

    try:
        if path.endswith(".gz"):
            fh = gzip.open( path, open_mode )
        elif path.endswith(".bz2"):
            fh = bz2.BZ2File( path, open_mode )
        else:
            fh = open( path, open_mode )
    except EnvironmentError:
        sys.exit( "Problem opening file: " + path)
    return fh

def gzip_bzip2_biom_open_readlines(path):
    """
    return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
    """

    # if the file is biom, convert to text and return lines
    if path.endswith(BIOM_FILE_EXTENSION):
        for line in read_biom_table(path):
            yield line
    else:
        with try_zip_open(path) as file_handle:
            for line in file_handle:
                if path.endswith(".bz2"):
                    # convert the line to text from binary
                    yield line.decode('utf-8').rstrip()
                else:
                    yield line.rstrip()


def unnamed_temp_file(prefix=None):
    """
    Return the full path to an unnamed temp file
    stored in the unnamed temp folder
    """

    if not prefix:
        prefix = "tmp"

    try:
        file_out, new_file = tempfile.mkstemp(dir=config.unnamed_temp_dir, prefix=prefix)
        os.close(file_out)
    except EnvironmentError:
        sys.exit("ERROR: Unable to create temp file in directory: " + config.unnamed_temp_dir)

    return (new_file)


def append_filename2contignames(fna_file):
    """
    Bind the file name to the contig names
    """

    # create a named temp file
    new_file = name_temp_file('.fna')
    print("new file",new_file)
    exe = "ppanini_rename_contigs"
    args = ["-i", fna_file, "-o", new_file]

    message = "Running " + exe + " ........"
    logger.info(message)
    print("\n" + message + "\n")

    execute_command(exe, args, [fna_file], [new_file])
    return new_file


def fasta_or_fastq(file):
    """
    Check to see if a file is of fasta or fastq format

    Fastq format short example:
    @SEQ_ID
    GATCTGG
    +
    !''****

    Fasta format short example:
    >SEQ_INFO
    GATCTGG

    Returns error if not of fasta or fastq format
    """

    format = "error"

    # check file exists
    file_exists_readable(file)

    # read in first 2 lines of file to check format
    file_handle = open(file, "rt")

    first_line = file_handle.readline()
    second_line = file_handle.readline()

    # check that second line is only nucleotides or amino acids
    if re.search("^[A-Z|a-z]+$", second_line):
        # check first line to determine fasta or fastq format
        if re.search("^@", first_line):
            format = "fastq"
        if re.search("^>", first_line):
            format = "fasta"

    file_handle.close()

    return format


def abundance(genes_file, alignment_file):
    """
    gene abundance table with featureCounts
    """
    # featureCounts -T 8 -g ID -t CDS  -a ../../prodigal_output/hmp_sub_nares/renamed_contigs_SRS015051.gff -o counts.txt SRS015051.sam

    # name the abundance file
    abundance_file = name_temp_file('.txt')

    exe = "featureCounts"
    opts = config.featureCounts_opts

    args = ["-a", genes_file, "-o", abundance_file, "-T", config.threads, alignment_file]

    message = "Running " + exe + " ........"
    logger.info(message)
    print("\n" + message + "\n")

    args = opts + args

    # create temp file for stdout and stderr
    # tmpfile=unnamed_temp_file("featureCounts_stdout_")
    # tmpfile2=unnamed_temp_file("featureCounts_stderr_")
    print(exe, args, [genes_file, alignment_file], [abundance_file])
    execute_command(exe, args, [genes_file, alignment_file], [abundance_file])  # ,
    # stdout_file=tmpfile, stderr_file=tmpfile2)

    return abundance_file


def index(custom_database):
    """
    Index database and run alignment with bowtie2
    """

    '''$ bowtie2-build -f renamed_contigs_SRS015051.fna  renamed_contigs_SRS015051_bowtie2_index_db 
	* for all samples $ sh ~/ppanini_stuff/scripts/mkbowtie2_dbs.sh '''
    # name the index
    index_name = name_temp_file(config.bowtie2_index_name)

    exe = "bowtie2-build"
    opts = config.bowtie2_build_opts

    args = ["-f", custom_database, index_name]

    outfiles = [index_name + ext for ext in config.bowtie2_index_ext_list]

    # if custom_database is large (>4G) then use the --large-index flag
    if os.path.getsize(custom_database) > config.bowtie2_large_index_threshold:
        args += ["--large-index"]
        outfiles = [index_name + config.bowtie2_large_index_ext]

    # index the database
    message = "Running " + exe + " ........"
    logger.info(message)
    print("\n" + message + "\n")

    args += opts

    # create temp file for stdout and stderr
    # tmpfile=unnamed_temp_file("bowtie2_stdout_")
    # tmpfile2=unnamed_temp_file("bowtie2_stderr_")

    execute_command(exe, args, [custom_database], outfiles)
    # stdout_file=tmpfile, stderr_file=tmpfile2)

    return index_name


def alignment(user_fastq, index_name):
    """
    Run alignment with bowtie2
    """
    '''
	$ bowtie2 -q -p 8 -x SRS015051_bowtie2_index_db -U SRS015051.fastq.gz -S SRS015051.sam
	* for all samples sh ~/ppanini_stuff/scripts/align_bowtie2.sh '''
    # name the alignment file
    alignment_file = name_temp_file('.sam')

    # align user input to database
    exe = "bowtie2"
    opts = config.bowtie2_align_opts

    # determine input type as fastq or fasta
    input_type = fasta_or_fastq(user_fastq)

    logger.debug("Nucleotide input file is of type: %s", input_type)

    # determine input type flag
    # default flag to fastq
    input_type_flag = "-q"
    if input_type == "fasta":
        input_type_flag = "-f"

    args = [input_type_flag, "-x", index_name, "-U", user_fastq, "-S", alignment_file]

    # add threads
    if config.threads > 1:
        args += ["-p", config.threads]

    # run the bowtie2 alignment
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    args += opts

    execute_command(exe, args, [user_fastq], [alignment_file])

    return alignment_file


def genecall(contig_file):
    """
    Run gene call with Prodigal
    """

    # name the genes file
    genes_file_gff = name_temp_file('.gff')
    genes_file_fna = name_temp_file('.fna')
    genes_file_faa = name_temp_file('.faa')

    # align user input to database
    exe = "prodigal"
    opts = config.prodigal_opts

    args = ["-i", contig_file, "-o", genes_file_gff, "-f", "gff", '-d', genes_file_fna,
            '-a', genes_file_faa, '-p', 'meta']

    # run the prodigal gene caller
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    args += opts

    execute_command(exe, args, [contig_file], [genes_file_gff, genes_file_fna, genes_file_faa])

    return genes_file_gff, genes_file_fna, genes_file_faa


def diamond_alignment(genes_file, uniref_db):
    """
    Run diamond alignment on database formatted for diamond
    """

    exe = "diamond"
    # $ diamond blastp --quiet --query hmp_sub_nares.faa
    # --db /n/huttenhower_lab/data/humann2_databases/uniref_annotated/uniref90/v1.1_uniref90/uniref90_annotated.1.1.dmnd --threads 8
    # --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen --out hmp_sub_nares.uniref90hits

    # Select the command based on a protein or nucleotide database search
    args = []
    opts = []
    '''
    if config.pick_frames_toggle == "on":
        args=[config.diamond_cmmd_protein_search]
    else:
        args=[config.diamond_cmmd_nucleotide_search]

    opts=config.diamond_opts'''
    alignment_file = name_temp_file('genes.uniref90hits')
    bypass = check_outfiles([alignment_file])
    args += ["blastp", "--quiet", "--query", genes_file,  # "--evalue",config.evalue_threshold,
             "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
             "send", "evalue", "bitscore", "qlen", "slen",
             '--db', uniref_db,
             '--out', alignment_file, "--threads", config.threads]

    message = "Running " + exe + " ........"
    logger.info(message)
    print("\n" + message + "\n")

    if not bypass:
        args += opts
        execute_command(exe, args, [genes_file], [alignment_file])
    else:
        message = "Bypass"
        logger.info(message)
        print(message)
    return alignment_file


def Infer_aligmnets(alignment_file, output):
    """
    Run infer_abundance to get sufficnet maaped genes (hits) and insufficient genes (no_hits)
    """
    # name the hits and no hits file
    hits = output + '/hits.txt'  # name_temp_file('_hits.txt')
    no_hits = output + '/no_hits.txt'  # no_hits_map = name_temp_file('_no_hits.txt')
    uniref_gene_map = output + '/uniref_gene_map.txt'
    # print (hits_map)
    # align user input to database
    exe = "ppanini_infer_gene"
    opts = ''

    args = [alignment_file, "--min-percid", .9, '--min-qcover', .8, "--min-scover", .8, '--output', output]

    # run the prodigal gene caller
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    args += opts

    execute_command(exe, args, [alignment_file], [hits, no_hits])

    return hits, no_hits, uniref_gene_map


def select_sequnces(fasta_file, sequnces_names, output_name='selected_genes.faa'):
    """
    Gets a fasta file and a list of names and returen sequnces for the list
    """
    # name the  output
    genes_faa = name_temp_file(output_name)

    # align user input to database
    exe = "ppanini_fasta_select"

    args = ["-i", fasta_file, "-f", sequnces_names, "-o", genes_faa]

    # run the prodigal gene caller
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    execute_command(exe, args, [fasta_file, sequnces_names], [genes_faa])

    return genes_faa


def cluster_genes(genes_fasta_file):
    """
    Run CD-Hit to cluster genes
    """
    # name the genes cluster output
    cluster_alignments = name_temp_file('no_hits_reads.clust90')
    cluster_gene_file = name_temp_file('no_hits_reads.clust90.clstr')

    # align user input to database
    exe = "cd-hit"
    # cd-hit -d 0 -c .9 -aL .8 -G 0 -T 2 -i ${infer_output}/no_hits_reads.faa -o ${infer_output}/no_hits_reads.clust90
    args = ["-i", genes_fasta_file, "-o", cluster_alignments, "-T", config.threads, "-M", config.cd_hit_memory]
    args += config.cd_hit_opts
    # run the prodigal gene caller
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    execute_command(exe, args, [genes_fasta_file], [cluster_gene_file, cluster_alignments])

    return cluster_gene_file, cluster_alignments


def mapping_clusters_genes(cluster_gene_file):
    """
    generate a map file for clusters (generated by CD-Hit) to genes
    """
    # name the genes cluster output
    # mapper_cluster_dir = name_temp_file('')
    mapper_cluster_genes = name_temp_file('cd-hit_gene_map.txt')

    # align user input to database
    exe = "ppanini_cluster2genes"
    args = ["-i", cluster_gene_file, "-o", config.temp_dir]

    # run the ppanini_cluster2genes
    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    execute_command(exe, args, [cluster_gene_file], [mapper_cluster_genes])

    return mapper_cluster_genes


def gene2genefamilies(tables, mapping_cluster, mapping_uniref, scale):
    """
   #ppanini_join_tables -i tables/
    #-o hmp_sub_nares_genefamilies_abund.tsv
    #--mapping-cluster cd_hit_clust_temp/cluster_gene_map.txt
    #--mapping-uniref ./hmp_nares_map_uniref_gene.txt
    """
    # name the genes cluster output
    gene_families_table = name_temp_file('gene_families.txt')

    # align user input to database
    exe = "ppanini_join_tables"
    opts = ''
    args = ["-i", tables, "-o", gene_families_table, "--mapping-cluster", mapping_cluster, "--mapping-uniref",
            mapping_uniref, "--scale", scale]

    message = "Running " + exe + " ........"
    print("\n" + message + "\n")

    args += opts

    execute_command(exe, args, [mapping_cluster, mapping_uniref], [gene_families_table])

    return gene_families_table


def name_temp_file(file_name):
    """
    Return the full path to a new temp file
    using the sample name and temp dir location
    """
    if config.file_basename != '':
        print("name_temp_file结果1",os.path.join(config.temp_dir, config.file_basename + file_name))
        return os.path.join(config.temp_dir,
                            config.file_basename + file_name)
    else:
        print("name_temp_file结果2", os.path.join(config.temp_dir, file_name))
        return os.path.join(config.temp_dir, file_name)


def return_exe_path(exe):
    """
    Return the location of the exe in $PATH
    """
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path, exe)
        if os.path.exists(fullexe) and os.path.isfile(fullexe):
            if os.access(fullexe, os.X_OK):
                return path
    return ""


def remove_file(file):
    """
    If file exists, then remove
    """

    try:
        if os.path.isfile(file):
            os.unlink(file)
            logger.debug("Remove file: %s", file)
    except OSError:
        message = "Unable to remove file"
        logger.error(message)


def remove_directory(dir):
    """
    Remove directory if exists
    """
    if os.path.isdir(dir):
        try:
            shutil.rmtree(dir)
            logger.debug("Remove directory: " + dir)
        except EnvironmentError:
            logger.error("Unable to remove directory: " + dir)
    else:
        logger.debug("Request to remove directory that does not exist: " + dir)


def file_exists_readable(file, raise_IOError=None):
    """
    Exit with error if file does not exist or is not readable
    Or raise an IOerror if selected
    """

    if not os.path.isfile(file):
        message = "Can not find file " + file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)
            raise IOError
        else:
            sys.exit("CRITICAL ERROR: " + message)

    if not os.access(file, os.R_OK):
        message = "Not able to read file " + file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)
            raise IOError
        else:
            sys.exit("CRITICAL ERROR: " + message)


def check_outfiles(outfiles):
    """
    If outfiles already_exist, then remove or bypass
    """
    bypass = []
    for file in outfiles:
        if os.path.isfile(file):
            if config.resume and os.path.getsize(file) > 0:
                bypass.append(True)
            else:
                bypass.append(False)
        else:
            bypass.append(False)

    if False in bypass or not bypass:
        # remove any existing files
        for file in outfiles:
            remove_file(file)
        return False
    else:
        return True

def byte_to_gigabyte(byte):
	"""
	Convert byte value to gigabyte
	"""

	return byte / (1024.0 ** 3)

def log_system_status():
    """
    Print the status of the system
    """

    module_available = True
    try:
        import psutil
    except ImportError:
        module_available = False

    if module_available:
        try:
            # record the memory used
            memory = psutil.virtual_memory()
            logger.info("Total memory = " + str(byte_to_gigabyte(memory.total)) + " GB")
            logger.info("Available memory = " + str(byte_to_gigabyte(memory.available)) + " GB")
            logger.info("Free memory = " + str(byte_to_gigabyte(memory.free)) + " GB")
            logger.info("Percent memory used = " + str(memory.percent) + " %")

            # record the cpu info
            logger.info("CPU percent = " + str(psutil.cpu_percent()) + " %")
            logger.info("Total cores count = " + str(psutil.cpu_count()))

            # record the disk usage
            disk = psutil.disk_usage('/')
            logger.info("Total disk = " + str(byte_to_gigabyte(disk.total)) + " GB")
            logger.info("Used disk = " + str(byte_to_gigabyte(disk.used)) + " GB")
            logger.info("Percent disk used = " + str(disk.percent) + " %")

            # record information about this current process
            process = psutil.Process()
            process_memory = process.memory_info()
            process_create_time = datetime.datetime.fromtimestamp(
                process.create_time()).strftime("%Y-%m-%d %H:%M:%S")
            process_cpu_times = process.cpu_times()
            # two calls required to cpu percent for non-blocking as per documentation
            process_cpu_percent = process.cpu_percent()
            process_cpu_percent = process.cpu_percent()

            logger.info("Process create time = " + process_create_time)
            logger.info("Process user time = " + str(process_cpu_times.user) + " seconds")
            logger.info("Process system time = " + str(process_cpu_times.system) + " seconds")
            logger.info("Process CPU percent = " + str(process_cpu_percent) + " %")
            logger.info("Process memory RSS = " + str(byte_to_gigabyte(process_memory.rss)) + " GB")
            logger.info("Process memory VMS = " + str(byte_to_gigabyte(process_memory.vms)) + " GB")
            logger.info("Process memory percent = " + str(process.memory_percent()) + " %")

        except (AttributeError, OSError, TypeError, psutil.Error):
            pass


def return_module_path(module):
    """
    Return the full path to the python module
    """

    # if this is the full path to the module then return the parent directory
    if os.path.isabs(module):
        if os.path.exists(module):
            return os.path.abspath(os.path.join(module, os.pardir))
    else:
        # search for the path to the module
        for path in sys.path:
            full_module = os.path.join(path, module)
            if os.path.exists(full_module):
                return path
    return ""

def execute_command(exe, args, infiles, outfiles, stdout_file=None,
                    stdin_file=None, raise_error=None, stderr_file=None):
    """
    Execute third party software or shell command with files
    """

    if exe == sys.executable:
        # check that the python module can be found
        module_path = return_module_path(args[0])
        if not module_path:
            message = "Can not find python module " + args[0]
            logger.critical(message)
            if raise_error:
                raise EnvironmentError
            else:
                sys.exit("CRITICAL ERROR: " + message)
        # update the module to the full path if not already the full path
        elif not os.path.isabs(args[0]):
            args[0] = os.path.join(module_path, args[0])

        logger.debug("Using python module : " + args[0])
    else:
        # check that the executable can be found
        exe_path = return_exe_path(exe)
        if not exe_path:
            message = "Can not find executable " + exe
            logger.critical(message)
            if raise_error:
                raise EnvironmentError
            else:
                sys.exit("CRITICAL ERROR: " + message)
        # update the executable to the full path
        else:
            exe = os.path.join(exe_path, exe)

        logger.debug("Using software: " + exe)

    # check that the input files exist and are readable
    for file in infiles:
        file_exists_readable(file, raise_IOError=raise_error)

    # check if outfiles already exist
    bypass = check_outfiles(outfiles)

    # convert numbers to strings
    args = [str(i) for i in args]

    if not bypass:

        cmd = [exe] + args

        message = " ".join(cmd)
        logger.info("Execute command: " + message)
        if config.verbose:
            print("\n" + message + "\n")

        # Open the input and output files (stdin, stdout, stderr)
        stdin = None
        stdout = None
        stderr = None

        if stdin_file:
            try:
                stdin = open(stdin_file, "rt")
            except EnvironmentError:
                message = "Unable to open file: " + stdin_file
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)

        if stdout_file:
            # check for file open mode
            try:
                stdout_file_name, mode = stdout_file
            except ValueError:
                stdout_file_name = stdout_file
                mode = "w"

            try:
                stdout = open(stdout_file_name, mode)
            except EnvironmentError:
                message = "Unable to open file: " + stdout_file_name
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)

        if stderr_file:
            try:
                stderr = open(stderr_file, "w")
            except EnvironmentError:
                message = "Unable to open file: " + stderr_file
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)

        try:
            if stdin_file or stdout_file or stderr_file:
                # run command, raise CalledProcessError if return code is non-zero
                p = subprocess.check_call(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
            else:
                p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                logger.debug(p_out)
        except (EnvironmentError, subprocess.CalledProcessError) as e:
            message = "Error executing: " + " ".join(cmd) + "\n"
            if hasattr(e, 'output') and e.output:
                message += "\nError message returned from " + os.path.basename(exe) + " :\n" + e.output.decode("utf-8")
            logger.critical(message)
            logger.critical("TRACEBACK: \n" + traceback.format_exc())
            log_system_status()
            if raise_error:
                raise
            else:
                sys.exit("CRITICAL ERROR: " + message)

        # check that the output files exist and are readable
        for file in outfiles:
            file_exists_readable(file, raise_IOError=raise_error)

    else:
        if config.verbose:
            print("Bypass: \n" + exe + " " + " ".join(args) + "\n")
        else:
            print("Bypass\n")


def process_gene_table_with_header(gene_table, allow_for_missing_header=None):
    """
    Process through the header portion of the gene table file
    """

    # try to open the file
    try:
        lines = gzip_bzip2_biom_open_readlines(gene_table)
    except EnvironmentError:
        sys.exit("Unable to read file: " + gene_table)

    # find the headers
    header = ""
    first_data_line = ""
    for line in lines:
        if line[0] == GENE_TABLE_COMMENT_LINE:
            header = line
        else:
            first_data_line = line
            break

    if not header and not allow_for_missing_header:
        sys.exit("File does not have a required header: " + gene_table +
                 " . Please add a header which includes the indicator: " +
                 GENE_TABLE_COMMENT_LINE)

    # provide the header, if one was found
    if header:
        yield header

    # provide the first data line
    yield first_data_line

    # now provide the remaining lines
    for line in lines:
        yield line


def rev_load_polymap(path_in='', path_out='', start=0, skip=None, allowed_keys=None, allowed_values=None,
                     write_output=True, sep='\t'):
    """
    Load a file like:
    A 1 2
    B 3
    C 4 5
    To a nested dict structure:
    {1:{A}, 2:{A}, 3:{B:1}, 4:{C:1}, 5:{C:1}}
    Inner values are not important (set to 1)
    """
    # if it used from command line
    if path_in == '':
        if len(sys.argv) < 3:
            sys.exit('Please provide an input and output path')
        path_in = sys.argv[1]
        path_out = sys.argv[2]
    polymap_all = {}
    print("--- Loading mapping file from:", path_in, file=sys.stderr)
    size_warn(path_in)
    for line in gzip_bzip2_biom_open_readlines(path_in):
        row = line.split("\t")
        key = row[start]
        # if the row input format is like: A\t1;2
        if sep == ';':
            row = row[1].split(";")
        if allowed_keys is None or key in allowed_keys:
            for i, value in enumerate(row):
                # value = value.replace(".","_")
                if i != start and (skip is None or i not in skip):
                    if allowed_values is None or value in allowed_values:
                        polymap_all.setdefault(value, {})[key] = 1  # polymap.setdefault( key, {} )[value] = 1

    if write_output:
        with gzip.open(path_out + '_dict.txt.gz', 'wt') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            for key, values in polymap_all.items():
                writer.writerow([key, ";".join(values)])
        print("--- Mappping is done")
    else:
        # print (polymap_all)
        return polymap_all


def load_polymap_dic(path, key_ind=0, value_ind=1, skip=None, allowed_keys=None, allowed_values=None):
    """
    Loads a file to a dictionry
    INPUT:
    for each row in file format is:
    UniRefID    GOID1;GOID2
    OUTPUT:
    Dictionary : UniRefID is key and value is GOID1;GOID2
    """
    polymap_all = {}
    print("--- Loading mapping file from:", path, file=sys.stderr)
    size_warn(path)
    for line in gzip_bzip2_biom_open_readlines(path):
        row = line.split("\t")
        key = row[key_ind]
        polymap_all[key] = row[value_ind]
    '''with open(path) as csv_file:
        reader = csv.reader(csv_file)
        polymap_all = dict(reader)'''
    return polymap_all


def load_polymap(path, start=0, skip=None, allowed_keys=None, allowed_values=None):
    """
    Load a file like:
    A 1 2
    B 1
    B 3
    C 1 2 4
    To a nested dict structure:
    {A:{1:1, 2:1}, B:{1:1, 3:1}, C:{1:1, 2:2, 4:1}
    Inner values are not important (set to 1)
    """
    polymap = {}
    print("--- Loading mapping file from:", path, file=sys.stderr)
    size_warn(path)
    for line in gzip_bzip2_biom_open_readlines(path):
        row = line.split("\t")
        key = row[start]
        if allowed_keys is None or key in allowed_keys:
            for i, value in enumerate(row):
                if i != start and (skip is None or i not in skip):
                    if allowed_values is None or value in allowed_values:
                        polymap.setdefault(key, {})[value] = 1
    return polymap


def uniref2go(ppanini_table, uniref_go_path):
    # Load mapping UniRef90--Go term dictionary
    go1000_uniref90_dic = load_polymap_dic(uniref_go_path)

    # add the GO terms to the end of each line (last column) of the ppanini_table
    ppanini_table['GO'] = ppanini_table.index.to_series().apply(go1000_uniref90_dic.get)


def check_cmd(cmd):
    # check if cmd is installed
    cmd_wh = "where" if sys.platform.system() == "Windows" else "which"
    try:
        subprocess.call([cmd_wh, cmd])
        return True
    except:
        print("No executable %s! Please install it! ", cmd)


def ncolors(n, colormap="jet"):
    """utility for defining N evenly spaced colors across a color map"""
    cmap = plt.get_cmap(colormap)
    cmap_max = cmap.N
    return [cmap(int(k * cmap_max / (n - 1))) for k in range(n)]


def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(set(a) & set(b)))


class ReportHook():
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(byte_to_gigabyte(total_size)) + " GB\n")
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} GB ".format(byte_to_gigabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += "{:3.2f}".format(percent_downloaded) + " %  " + \
                          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                          "{:2.0f}".format(estimated_minutes) + " min " + \
                          "{:2.0f}".format(estimated_seconds) + " sec "
            status += "        \r"
            sys.stdout.write(status)


def download_tar_and_extract_with_progress_messages(url, filename, folder):
    """
    Download the file at the url
    """

    # check for local file
    local_file = False
    if os.path.isfile(url):
        local_file = True

    if not local_file:
        print("Download URL: " + url)

    try:
        if not local_file:
            url_handle = urlretrieve(url, filename, reporthook=ReportHook().report)
        else:
            filename = url

        print("\nExtracting: " + filename)
        # comment before
        tarfile_handle = tarfile.open(filename)
        tarfile_handle.extractall(path=folder)

    except (EnvironmentError, tarfile.ReadError):
        if local_file:
            sys.exit("CRITICAL ERROR: Unable to extract from local file: " + url)
        else:
            sys.exit("CRITICAL ERROR: Unable to download and extract from URL: " + url)


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / (1024.0 ** 2)

# runs MAP and ID steps
def map_id(seq_file1, seq_file2, ref_database, out_dir, out_align, threads, index_dir, tag):
    '''
    make use of the function!!!
    '''
    if seq_file2 is None:
        cmd = f"BGCLens MAP -U {seq_file1} -numThreads {threads} -indexDir {index_dir}" \
              f" -targetIndexPrefixes {ref_database} -outDir {out_dir} -outAlign {out_align}.sam"
    else:
        cmd = f"BGCLens MAP -1 {seq_file1} -2 {seq_file2} -numThreads {threads}" \
              f" -indexDir {index_dir} -targetIndexPrefixes {ref_database} -outDir {out_dir} -outAlign {out_align}.sam"

    os.system(cmd)

    cmd = f"BGCLens ID -alignFile {os.path.join(out_dir, out_align + '.sam')}" \
          f" -fileType sam -outDir results_ID -expTag {tag}"
    os.system(cmd)
    return 'DONE!!'