# default values for user options
import configparser
import os
import sys

version = '1.0.0'
input_table = ''
output_folder = '~/Documents/GitHub/BGCLens/BGCLens/Test'
basename = 'BGCLens'
log_level = 'DEBUG'
verbose = 'DEBUG'
nprocesses = 1

# prioritizing thresholds
tshld_abund = 75
tshld_prev = 75
abundance_detection_level = 0.0001
beta = 0.5
tshld = None
summary_table = None
niches = []
niche_flag = None
ppanini_niche_score_labels = []

# temp_folder = ''
centroids_list = ''
essantial_genes_uniref90_id = None
genomic_score = False
uniref2go = ''

# uniref formatting
uniref_delimiter = "|"
uniref_gene_index = -2
uniref_length_index = -1

# bowtie2 options and threshold
bowtie2_large_index_threshold = 4000000000
bowtie2_index_ext_list = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                          ".rev.1.bt2", ".rev.2.bt2"]
bowtie2_large_index_ext = ".1.bt2l"
bowtie2_version = {
    "flag": "--version",
    "major": 2,
    "minor": 2,
    "line": 0,
    "column": 2}

# prodigal options
prodigal_opts = ["-q"]

one_contig = False

bowtie2_build_opts = ["-q"]  # "--threads "+str(threads)
bowtie2_align_opts = ["--sensitive"]  # "--threads "+str(threads)
bowtie2_index_name = "_bowtie2_index"

threads = 1
cd_hit_memory = 200000

cd_hit_opts = ["-d", 0, "-c", .9, "-aL", .8, "-G", 0]
featureCounts_opts = ["-g", "ID", "-t", "CDS"]

# translated alignment options
translated_alignment_choices = ["usearch", "rapsearch", "diamond", "vsearch"]
translated_alignment_selected = translated_alignment_choices[2]

# file naming
temp_dir = ''
unnamed_temp_dir = temp_dir
file_basename = ''

# diamond options
diamond_database_extension = ".dmnd"
diamond_opts_uniref50 = ["--max-target-seqs", 20, "--sensitive", "--outfmt", 6]
diamond_opts_uniref90 = ["--max-target-seqs", 20, "--outfmt", 6]
diamond_cmmd_protein_search = "blastp"
diamond_cmmd_nucleotide_search = "blastx"
diamond_version = {
    "flag": "--version",
    "major": 0,
    "minor": 8,
    "second minor": 22,
    "line": 0,
    "column": 2}
pick_frames_toggle = 'on'

resume = False

# Modify on 07-25-2022 - Xinyang

user_edit_config_file = "BGCLens.cfg"

# print(os.path.abspath(__file__))
#
# full_path_user_edit_config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                                                user_edit_config_file)
# print(full_path_user_edit_config_file)

full_path_user_edit_config_file = os.path.join(os.path.dirname(os.path.abspath(".")),"BGCLens/BGCLens",
                                               user_edit_config_file)
print("Full path",full_path_user_edit_config_file)

def update_user_edit_config_file_single_item(section, name, value):
    """
    Update the settings to the user editable config file for one item
    """
    print("section",section)
    print("name",name)
    print("value",value)

    new_config_items = {section: {name: value}}
    print("new_config_items",new_config_items)

    update_user_edit_config_file(new_config_items)

    print("BGCLens configuration file updated: " + section + " : " + name + " = " + str(value))


def update_user_edit_config_file(new_config_items):
    """
    Update the settings to the user editable config file
    """

    config = configparser.RawConfigParser()

    # start with the current config settings
    config_items = read_user_edit_config_file()
    print("new_config_items1111111",new_config_items)

    # update with the new config items
    for section in new_config_items:
        for name, value in new_config_items[section].items():
            if section in config_items:
                if name in config_items[section]:
                    config_items[section][name] = value
                else:
                    sys.exit("ERROR: Unable to add new name ( " + name +
                             " ) to existing section ( " + section + " ) to " +
                             " config file: " + full_path_user_edit_config_file)
            else:
                sys.exit("ERROR: Unable to add new section ( " + section +
                         " ) to config file: " + full_path_user_edit_config_file)

    for section in config_items:
        config.add_section(section)
        for name, value in config_items[section].items():
            value = str(value)
            if "file" in section or "folder" in section:
                # convert to absolute path if needed
                if not os.path.isabs(value):
                    value = os.path.abspath(value)
            config.set(section, name, value)

    try:
        print("start reading")
        file_handle = open(full_path_user_edit_config_file, "wt")
        print(file_handle)
        config.write(file_handle)
        file_handle.close()
    except EnvironmentError:
        sys.exit("Unable to write to the BGCLens config file.")


def read_user_edit_config_file():
    """
    Read the settings from the config file
    """

    config = configparser.ConfigParser()
    # print("config",config)
    # print("full_path_user_edit_config_file",full_path_user_edit_config_file)

    try:
        # print("full_path_user_edit_config_file",full_path_user_edit_config_file)
        config.read(full_path_user_edit_config_file)
        # config.read("/Users/xinyangzhang/Documents/GitHub/BGCLens/BGCLens/BGCLens.cfg")

    except EnvironmentError:
        sys.exit("Unable to read from the config file: " + full_path_user_edit_config_file)

    # print("config.sections()",config.sections())
    # read through all of the sections
    config_items = {}
    for section in config.sections():
        # print("section", section)
        config_list = config.items(section)
        # print("config_list", config_list)
        config_items[section] = {}
        for name, value in config_list:
            if "file" in section or "folder" in section:
                # if not absolute path, then return absolute path relative to this folder
                if not os.path.isabs(value):
                    value = os.path.abspath(os.path.join(os.path.dirname(full_path_user_edit_config_file), value))
            config_items[section][name] = value

    # print("config_items",config_items)
    return config_items

# if __name__=='__main__':
# 	pass
