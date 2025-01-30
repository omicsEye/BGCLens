import sys

# # Try to load one of the BGCLens src modules to check the installation
# try:
#     from .. import check
# except ImportError:
#     sys.exit("CRITICAL ERROR: Unable to find the BGCLens python package." +
#              " Please check your install.")
#
# # Check the python version
# check.python_version()

import argparse

from .. import config


def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="BGCLens Configuration\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--print",
        dest="print_config",
        action="store_true",
        help="print the configuration\n")
    parser.add_argument(
        "--update",
        nargs=3,
        metavar=("<section>", "<name>", "<value>"),
        help="update the section : name to the value provided\n")

    return parser.parse_args()


def main():
    # Parse arguments from the command line
    args = parse_arguments(sys.argv)

    if args.update:
        # update the config file
        config.update_user_edit_config_file_single_item(args.update[0], args.update[1], args.update[2])

    if args.print_config or not args.update:
        # print the current configuration
        current_config_items = config.read_user_edit_config_file()
        print("BGCLens Configuration ( Section : Name = Value )")
        for section in current_config_items:
            for name, value in current_config_items[section].items():
                print(section + " : " + name + " = " + str(value))