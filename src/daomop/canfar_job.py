"""
Build an input file to be run via calling condor on canfar.
"""
import argparse
import sys
import os
import stat
from vos import Client
from . import storage
from .params import qrunid_end_date, qrunid_start_date


def bad_exposure_list():
    """
    Retreieve the list of bad exposrures from S. Gwyn's VOSpace listing.
    
    :return: list of bad exposures
     :rtype list
    """

    Client().copy(u'vos:sgwyn/tkBAD', u'tkBAD')

    bad_expnum_list = []
    for lines in open('tkBAD').readlines():
        # noinspection PyBroadException
        try:
            bad_expnum_list.append(lines.split()[0])
        except:
            sys.stderr.write(lines)
    return bad_expnum_list


def create_build_cat_command(command):

    with open(command, 'w') as fout:
        fout.write("#!/bin/bash -i\n")
        fout.write("export HOME=`cd ~ ; pwd`\n")
        fout.write("source ${HOME}/.bashrc\n")
        fout.write("source activate ossos\n")
        fout.write("cp cadcproxy.pem ${HOME}/.ssl/\n")
        fout.write("""echo "Processing ",$1\n""")
        fout.write("""daomop_populate $1 --verbose\n""")
        fout.write("""daomop_cat $1 --verbose\n""")

    os.chmod(command, stat.S_IXGRP | stat.S_IRWXU | stat.S_IXOTH | stat.S_IROTH | stat.S_IRGRP)
    return command


def create_build_cat(qrunid, force=False, job_filename="build_cat_job.in", command_filename="build_cat.sh"):
    """
    Create a CANFAR job file to run the build_cat script on all available CFIS 'r' band data.
    
    :param command_filename: name of the command the will run build_cat on th VM in CANFAR.
     :param qrunid:  qrun are we building catalog for.
     :param force: compute catalog even if previously successful.
     :param job_filename:
    :return: name of the job submission file and the build_cat command file.
    """
    command_filename = create_build_cat_command(command_filename)
    bad_expnum_list = bad_exposure_list()
    force = force and " --force " or ""

    with open(job_filename, 'w') as fout:
        fout.write("""Universe   = vanilla\n""")
        fout.write("""should_transfer_files = YES\n""")
        fout.write("""when_to_transfer_output = ON_EXIT_OR_EVICT\n""")
        fout.write("""RunAsOwner = True\n""")
        fout.write("""transfer_output_files = /dev/null\n\n""")
        fout.write("Executable = {}\n".format(command_filename))

        for row in storage.get_cfis_exposure_table(start_date=qrunid_start_date(qrunid),
                                                   end_date=qrunid_end_date(qrunid)):
            if row['observationID'] in bad_expnum_list:
                continue
            fout.write("Arguments = {} {}\n".format(row['observationID'], force))
            fout.write("Output = {}.out\n".format(row['observationID']))
            fout.write("Log = {}.log\n".format(row['observationID']))
            fout.write("Error = {}.err\n".format(row['observationID']))
            fout.write("Queue\n\n")

    return job_filename, command_filename

      
def create_stationary_command(command_filename="stationary.sh"):
    """
    Create the shell script that is the processing step.

    :param command_filename:  target_name of file to contain the processing shell script.
    :return:
    """
    with open(command_filename, 'w') as fout:
        fout.write("""#!/bin/bash -i\n""")
        fout.write("""export HOME=`cd ~ ; pwd`\n""")
        fout.write("""source ${HOME}/.bashrc\n""")
        fout.write("""source activate ossos\n\n""")
        fout.write("""cp cadcproxy.pem ${HOME}/.ssl/\n""")
        fout.write("""echo "Building Catalog ",$1 $2\n""")
        fout.write("""daomop_stationary $1 $2 --verbose """)

    os.chmod(command_filename, stat.S_IXGRP | stat.S_IRWXU | stat.S_IXOTH | stat.S_IROTH | stat.S_IRGRP)
    return command_filename


def create_stationary_job(qrunid, force=False, job_filename="stationary_job.in", command_filename="stationary.sh", runids=['17BC99']):
    """Build the stationary catalog builder job submission script.

    :param qrunid: Which qrunid should this stationary catalog be for
    :param force: build stationary catalogs even for images that are marked as completed.
    :param job_filename: Name of the executable script that will be run.
    :param command_filename: name of file that will hold the executable command.
    :return:
    """
    command_filename = create_stationary_command(command_filename)

    with open(job_filename, 'w') as fout:
        fout.write("""Universe   = vanilla\n""")
        fout.write("""should_transfer_files = YES\n""")
        fout.write("""when_to_transfer_output = ON_EXIT_OR_EVICT\n""")
        fout.write("""RunAsOwner = True\n""")
        fout.write("""transfer_output_files = /dev/null\n""")
        fout.write("""\nExecutable = {}\n""".format(command_filename))

        force = force and "--force" or ""
        for healpix in storage.list_healpix(start_date=qrunid_start_date(qrunid),
                                            end_date=qrunid_end_date(qrunid), runids=runids):
            params = {"Arguments": "{} {} {}".format(healpix, qrunid, force),
                      "Log": "{}.log".format(healpix),
                      "Output": "{}.out".format(healpix),
                      "Error": "{}.err".format(healpix)
                      }

            for param in params:
                fout.write("{} = {}\n".format(param, params[param]))
            fout.write("Queue\n\n")

    return job_filename, command_filename


def main():

    parser = argparse.ArgumentParser(description="Build the job and input files needed to launch a canfar job.")
    parser.add_argument("command", choices=['stationary', 'build_cat'], help="Which command to build job for")
    parser.add_argument("qrunid", help="which CFHT Queue run is this job for.", action="store", default="%")
    parser.add_argument("runids", nargs='*', help="which CFHT Queue run is this job for.", action="store", default="%")
    parser.add_argument("--force", help="Force job to run on previously processed images.", action='store_true',
                        default=False)

    args = parser.parse_args()
    if args.command == "stationary":
        return create_stationary_job(args.qrunid, force=args.force, runids=args.runids)
    elif args.command == "build_cat":
        return create_build_cat(args.qrunid, force=args.force)


if __name__ == "__main__":
    sys.exit(main())
