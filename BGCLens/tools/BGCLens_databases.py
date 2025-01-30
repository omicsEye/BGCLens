#!/usr/bin/env python

"""
BGCLens: BGCLens_databases module
Download databases an update config settings
Dependencies: None
To Run: BGCLens_databases --download <database> <build> <install_location>
Copyright (c) 2014 Harvard School of Public Health
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import argparse
import os
import sys

from .. import config
from .. import utilities

# the locations of the current databases to download
current_downloads = {
    "bgc":
        {
            "MIBiG": "https://gwu.box.com/shared/static/s9g2v012kyy97p50juy4boybeojrfhda.gz"
        },
    "fungi":
        {
            "fungiBowtie2Index": "https://gwu.box.com/shared/static/mm6oxn4ovfcc5rhaly27lxh40c3filro.gz"
        },
    "human":
        {
            "GRCh38Bt2": "https://gwu.box.com/shared/static/b0uu1sqyk1y6rb4wzllssxv1vox84205.gz"
        },
    "viral":
        {
            "viralBowtie2Index": "https://gwu.box.com/shared/static/pho5jc7nwbr00hxqboxew7zrug4p6bvf.gz"
        },
    "bacterial":
        {
            "bacterialBowtie2Index": "https://gwu.box.com/shared/static/ibzpdz9fac2twxssquunp99rcld9aaf5.gz"
        },
    "archaeal":
        {
            "archaealBowtie2Index": "https://gwu.box.com/shared/static/afj4s72zi0hbjgbt1ga6j35icf93nypp.gz"
        },
}

database_type = {
    "bgc": "bgc",
    "fungi": "fungi",
    "human": "human",
    "viral": "viral",
    "bacterial": "bacterial",
    "archaeal": "archaeal"
}


def check_user_database(original, user):
    """ Check that the user database is of the expected version """

    if original.split("/")[-1] == user.split("/")[-1]:
        return True
    else:
        sys.exit("The user database selected does not match that expected: " + original.split("/")[-1])


def download_database(database, build, location, database_location):
    """
    Download and decompress the selected database
    """

    install_location = ""
    if database in current_downloads:
        if build in current_downloads[database]:
            # create a subfolder to hold the contents of the database
            install_location = os.path.join(location, database)
            if not os.path.isdir(install_location):
                try:
                    print("Creating subdirectory to install database: " + install_location)
                    os.mkdir(install_location)
                except EnvironmentError:
                    sys.exit("CRITICAL ERROR: Unable to create directory: " + install_location)

            # download the database
            downloaded_file = os.path.join(location, current_downloads[database][build].split('/')[-1])

            if database_location:
                check_user_database(current_downloads[database][build], database_location)
                print("download_tar_and_extract_with_progress_messages11", database_location)
                print("downloaded_file", downloaded_file)
                print("install_location", install_location)
                utilities.download_tar_and_extract_with_progress_messages(database_location,
                                                                          downloaded_file, install_location)
            else:
                print("download_tar_and_extract_with_progress_messages22", database_location)
                print("current_downloads22", current_downloads[database][build])
                print("downloaded_file22", downloaded_file)
                print("install_location22", install_location)
                utilities.download_tar_and_extract_with_progress_messages(current_downloads[database][build],
                                                                          downloaded_file, install_location)

            # remove the download (if not a local file provided by the user)
            if not database_location or (database_location and not os.path.isfile(database_location)):
                try:
                    os.unlink(downloaded_file)
                except EnvironmentError:
                    print("Unable to remove file: " + downloaded_file)

            print("\nDatabase installed: " + install_location + "\n")
        else:
            sys.exit("ERROR: Please select an available build.")
    else:
        sys.exit("ERROR: Please select an available database.")

    return install_location


def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="BGCLens Databases\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--available",
        action="store_true",
        help="Print the available databases\n")
    parser.add_argument(
        "--download",
        nargs=3,
        metavar=("<database>", "<build>", "<install_location>"),
        help="Download the selected database to the install location\n")
    parser.add_argument(
        "--update-config",
        default="yes",
        choices=["yes", "no"],
        help="Update the config file to set the new database as the default [DEFAULT: yes]\n")
    parser.add_argument(
        "--database-location",
        help="Location (local or remote) to pull the database")

    return parser.parse_args()




def main():
    # Parse arguments from the command line
    args = parse_arguments(sys.argv)

    if args.download:
        # download the database
        database = args.download[0]
        build = args.download[1]
        location = os.path.abspath(args.download[2])

        print("database", database)
        print("build", build)
        print("location", location)

        # create the installation location if it does not already exist
        print("os.path.isdir(location)", os.path.isdir(location))
        if not os.path.isdir(location):
            try:
                print("Creating directory to install database: " + location)
                os.mkdir(location)
            except EnvironmentError:
                sys.exit("CRITICAL ERROR: Unable to create directory: " + location)

        install_location = download_database(database, build, location, args.database_location)
        print("install_location", install_location)
        if args.update_config == "yes":
            # update the config file with the installed location
            print("database_type[database]", database_type[database])
            config.update_user_edit_config_file_single_item("database_folders",
                                                            database_type[database], install_location)

    if args.available or not args.download:
        # print the available databases
        current_config_items = config.read_user_edit_config_file()
        print("BGCLens Databases ( database : build = location )")
        for database in current_downloads:
            for build, location in current_downloads[database].items():
                print(database + " : " + build + " = " + location)



if __name__ == "__main__":
    main()