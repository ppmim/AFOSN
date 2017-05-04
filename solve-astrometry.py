#! /usr/bin/env python
# encoding:UTF-8

# Copyright (c) 2013 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of OSN CCD astrometry procedure
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



from __future__ import division

import sys
import os
import os.path
import shutil
import fileinput
import logging
import subprocess
import multiprocessing
import glob
from optparse import OptionParser

import pyfits as fits
import logging as log
import time
import datetime

# TODO
#
# - Tiempo límite de resolución de un fichero --cpulimit (default to 300s)
# - distorsion promedio (encontre un mail de Dustin donde hablaba de eso)


# Project modules
import clfits


def readHeader(filename, check_filename=False):
    """
    Read from the FITS header values required for astrometric calibration
    """

    try:
        fits = clfits.ClFits(filename)
    except Exception as e:
        msg = 'Error reading FITS file: ' + filename
        log.error(msg)
        log.error(str(e))
        raise e
    else:
        # Return values
        scale = fits.pix_scale
        ra = fits.ra
        dec = fits.dec
        instrument = fits.getInstrument()
        is_science = fits.isScience()

        # Sometimes, mainly for flat-fields files, the IMAGETYP has not the right value, and it has the value 'LIGHT'
        # instead of 'FLAT' o 'SKY_FLAT'. So, we optionally check the filename for not science frames to detect
        # sky flats.
        if check_filename and ('flat' or 'bias' or 'dark') in filename.lower() and is_science:
            log.info("Discarding file by filename %s:" % filename)
            is_science = False
            if 'flat' in filename.lower():
                with fits.open(filename, mode='update') as ff:
                    ff[0].header['IMAGETYP'] = 'FLAT'

        return scale, ra, dec, instrument, is_science


def solveField(filename, tmp_dir, pix_scale=None, blind=False, patch_calibs=True):
    """
    Do astrometric calibration to the given filename using Astrometry.net 
    function 'solve-field'
    
    Parameters
    ----------
    filename : str
        File name of the FITS file to resolve.

    tmp_dir: str
        Directory where output (temporal and results) files are saved
    
    pix_scale: float
        Default pixel scale to use in case it cannot be find out from header
    
    blind: bool
        When True, a blind calibration with no reference coordinates is done.

    patch_calibs: bool
        When True, the calibration files (bias, dark, flat) are also patched,
        both header and filename.      
        
    Returns
    -------
    New filename of just solved file (or Exception if not solved).
    
    """

    #
    # Create temporal directory
    #    
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    #
    # Read header parameters
    #    
    (scale, ra, dec, instrument, is_science) = readHeader(filename, check_filename=True)

    # Whether no scale was found out and some was given as default, we use it
    if scale == -1 and pix_scale != None:
        scale = pix_scale

    if not is_science:
        log.info("Frame %s is not a science frame" % filename)
        if patch_calibs:
            # patch the header
            with fits.open(filename, mode='update') as new_fits:
                logging.info("Patching the header")
                patch_header(new_fits[0].header)
            # patch the filename
            new_fn = patch_filename(filename, rename=False, ext='fits')

            # Copy to output dir (as the science frames solved)
            try:
                shutil.copyfile(filename, tmp_dir + '/' + os.path.basename(new_fn))
            except OSError as ex:
                logging.error("Cannot rename file %s" % new_fn)
                return new_fn + '.not_science'
            else:
                return tmp_dir + '/' + os.path.basename(new_fn) + '.not_science'
        else:
            return "{0}.not_science".format(filename)

    logging.debug("Starting to solve-field for: %s  Scale=%s  RA= %s Dec= %s \
    INSTRUMENT= %s" % (filename, scale, ra, dec, instrument))

    # Hardcoded the path !
    path_astrometry = "/usr/local/astrometry/bin"

    #
    # We must distinguish different cases
    #

    # 0) blind calibration
    if blind:
        log.debug("Blind calibration, nothing is supposed")
        str_cmd = "%s/solve-field -O -p -D %s -m %s %s\
        " % (path_astrometry, tmp_dir, tmp_dir, filename)

    # 1) RA, Dec and Scale are known
    elif ra != -1 and dec != -1 and scale != -1:
        log.debug("RA, Dec and Scale are known")
        # To avoid problems with wrong RA,Dec coordinates guessed, a wide 
        # radius is used (0.5 degrees)
        # Although --downsample is used, scale does not need to be modified
        str_cmd = "%s/solve-field --cpulimit 60 -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s --ra %s --dec %s --radius 0.5 -D %s -m %s %s --downsample 2\
        " % (path_astrometry, scale - 0.05, scale + 0.05, ra, dec, tmp_dir, tmp_dir, filename)
    # 2) RA, Dec are unknown but scale is
    elif (ra == -1 or dec == -1) and scale != -1:
        log.debug("RA, Dec are unknown but scale is")
        str_cmd = "%s/solve-field -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s -D %s -m %s %s\
        " % (path_astrometry, scale - 0.1, scale + 0.1, tmp_dir, tmp_dir, filename)
    elif (ra != -1 and dec != -1):
        str_cmd = "%s/solve-field --cpulimit 60 -O -p  \
        --ra %s --dec %s --radius 0.5 -D %s -m %s %s --downsample 2\
        " % (path_astrometry, ra, dec, tmp_dir, tmp_dir, filename)
    # 4) None is known -- blind calibration
    else:
        log.debug("Nothing is known")
        str_cmd = "%s/solve-field -O -p -D %s -m %s %s\
        " % (path_astrometry, tmp_dir, tmp_dir, filename)

    clean_output = True
    if clean_output:
        str_cmd += " --corr none --rdls none --match none --wcs none"

    log.debug("CMD=" + str_cmd)

    try:
        p = subprocess.Popen(str_cmd, bufsize=0, shell=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True)
    except Exception as ex:
        log.error("Some error while running subprocess: " + str_cmd)
        log.error(str(ex))
        raise ex

    # Warning:
    # We use communicate() rather than .stdin.write, .stdout.read or .stderr.read 
    # to avoid deadlocks due to any of the other OS pipe buffers filling up and 
    # blocking the child process.(Python Ref.doc)

    (stdoutdata, stderrdata) = p.communicate()
    solve_out = stdoutdata + "\n" + stderrdata

    # print "STDOUTDATA=",stdoutdata
    # print "STDERRDATA=",stderrdata

    if len(solve_out) > 1:
        logging.info("Solve-field output:")
        print solve_out
    #
    # Look for filename.solved to know if field was solved
    #
    solved_file = tmp_dir + "/" + os.path.splitext(os.path.basename(filename))[0] + ".solved"

    # either succeeded or failed, clean up
    try:
        os.remove(solved_file.split('.solved')[0] + '.axy')
    except OSError:
        pass
    ## If solved
    if os.path.exists(solved_file):
        logging.info("Field solved !")
        # clean up        
        try:
            os.remove(solved_file.split('.solved')[0] + '-indx.xyls')
            os.remove(solved_file.split('.solved')[0] + '.solved')
        except OSError:
            pass
        # patch the header
        with fits.open(solved_file.split('.solved')[0] + '.new', mode='update') as new_fits:
            logging.info("Patching the header")
            patch_header(new_fits[0].header)

        # patch the filename        
        new_fn = patch_filename(solved_file.split('.solved')[0] + '.new', ext='fits')

        return new_fn
    ## If not solved, then second try blind solve (without coordinates)
    else:
        if not blind:
            log.error("First try to solve failed. Lets try again...")
            try:
                return solveField(filename, tmp_dir,
                                  pix_scale=None, blind=True)
            except Exception as ex:
                raise ex
        else:
            raise Exception("Field was not solved")


def calc(args):
    """
        Method used only to use with Pool.map_asycn() function
        """
    return solveField(*args)


def runMultiSolver(files, tmp_dir, pix_scale=None):
    """
    Run a parallel processing to solve astrometry for the input files taking
    advantage of multi-core CPUs.
    """

    # use all CPUs available in the computer
    n_cpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=n_cpus)

    log.debug("N_CPUS :" + str(n_cpus))
    log.debug("FILES_TO_SOLVE :" + str(files))
    log.debug("TMP_DIR :" + str(tmp_dir))

    results = []
    solved = []
    not_solved = []
    for file in files:
        red_parameters = (file, tmp_dir, pix_scale)
        try:
            # Instead of pool.map() that blocks until
            # the result is ready, we use pool.map_async()
            results += [pool.map_async(calc, [red_parameters])]
        except Exception as ex:
            not_solved.append(file)
            log.error("Error processing/solving file: %s\n" % file)
            log.error(str(ex))

    for result in results:
        try:
            result.wait()
            # the 0 index is *ONLY* required if map_async is used !!!
            solved.append(result.get()[0])
        except Exception as e:
            log.error("Cannot get result: %s\n" %(str(e)))
            

    # Prevents any more tasks from being submitted to the pool. 
    # Once all the tasks have been completed the worker 
    # processes will exit.
    pool.close()

    # Wait for the worker processes to exit. One must call 
    # close() or terminate() before using join().
    pool.join()

    log.info("Finished parallel calibration")

    return solved, not_solved


def patch_header(header):
    """
    Patch the FITS header.

    Parameters
    ----------
    header : object (FITS header class. This class exposes both a dict-like 
                    interface and a list-like interface to FITS headers.)
        Header to be patched. 

    """

    # remove useless keyword
    useless_keywords = ['WCSDIM', 'LTM1_1', 'LTM2_2', 'WAT0_001', 'WAT1_001',
                        'WAT2_001', 'RADECSYS', 'SWCREATE', 'SWOWNER']

    for key in useless_keywords:
        # If key is not in header, return None.
        header.pop(key, None)

    # Now, remove all keywords starting with '_'
    for key in header:
        if key.startswith('_'):
            # If key is not in header, return None.
            header.pop(key, None)

    # Remove all comments starting with 'Original'
    new_cards = [card for card in header.cards
                 if not (card[0] == 'COMMENT' and card[1].startswith('Original'))]
    header.cards._header._cards = new_cards

    # patch OBSERVER keyword (remove @iaa.es)
    try:
        if header['OBSERVER']:
            header['OBSERVER'] = header['OBSERVER'].split('@')[0]
        else:
            header['OBSERVER'] = 'osn_astronomer'
    except Exception:
        pass


def patch_filename(filename, rename=True, ext=None, check_data_obs=True):
    """
    Patch the FITS filename and optionally rename the file on disk.

    Parameters
    ----------
    filename : str 
        Full file path of the file to be checked and patched if needed.
    
    rename: bool
        If  True, proceed to rename the original filename on disk
    
    ext: str
       Filename extension to be set replacing the original filename extension.

    check_data_obs: bool
        True if required a prefix on filename based on DATE-OBS keyword; it 
        compare to DATE-OBS, and in case of mismatch, it is added as prefix.

    Returns
    -------
    New patched filename.   

    """
    original = filename
        
    # Replace non-standar or characters known to the GAVO staff to be 
    # hazardous in URLs.   
    import re
    for c in re.findall(r'[^A-Za-z0-9_\-\\\/.]', filename):
        filename = filename.replace(c, "_")
    
    # Check the filename prefix is YYYYMMDD.hhmmss obteined from DATE-OBS 
    # keyword (YYYY-MM-DDThh:mm:ss).
    # Due to Astrometry.net output solved files are saved on the same directory,
    # it needed that each file has an unique filename to avoid overwriting.
    base_filename = os.path.basename(filename) 
    if check_data_obs:
        if len(base_filename) > 15:
            prefix = base_filename[:15]
        else:
            prefix = ""
        # read the DATE-OBS keyword        
        try:
            logging.info("Original filename: %s" % original)
            with fits.open(original) as f:
                date_obs = f[0].header['DATE-OBS'][:19]
       
            date_obs = date_obs.replace('-','')
            date_obs = date_obs.replace('T','.')
            date_obs = date_obs.replace(':','')
            if date_obs != prefix:
                print "DATE_OBS=",date_obs
                print "prefix=",prefix
                # prefix does not match, then **add** new prefix based on date_obs
                base_filename = date_obs + "." + base_filename
                filename = os.path.dirname(filename) + "/" + base_filename
                logging.info("New filename -> %s" %filename)                                   
 
        except Exception, e:
            logging.error("Cannot read DATA-OBS from %s" % original)
            raise e
                                                
    # if required, rename extension
    if ext:
        filename = filename.replace(filename.split(".")[-1], ext)

    if rename:
        try:
            os.rename(original, filename)
        except Exception as e:
            logging.error("Cannot rename file %s" % original)
            return original

    return filename


###############################################################################
# main
###############################################################################
#
# Procedure:
#  1. Check type of file (SCIENCE or CALIB)
#  2. IF CALIB:
#     2.1 Patch header and filename of original file
#     2.2 Copy to dst directory
#  3. IF SCIENCE:
#     3.1 Solve astrometrically the SCIENCE files
#     3.2 Patch header and filename or solved file 
#     3.3 Move solved and patched file to dst directory
#  4. Generate the report (txt files)
#  
# !! WARNING !! : the basename of source files must be unique !! it is because
# the (solved) output files are all copied to the **same** output directory.
# In principle, input files have the format YYYYMMDD.hhmmss.name_id.fit[s]  
# (Telescope id should be added !)

if __name__ == "__main__":

    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2 ..."
    desc = """Performs the astrometric calibration of a set of images,
in principle previously reduced, but not mandatory.

"""
    parser = OptionParser(usage, description=desc)

    parser.add_option("-s", "--source",
                      action="store", dest="source_file",
                      help="Source file list of data frames. "
                           "It can be either a file or a directory name.")

    parser.add_option("-o", "--output_dir",
                      action="store", dest="output_dir", default="/tmp",
                      help="Place all output files in the specified directory [default=%default]")

    parser.add_option("-p", "--pixel_scale",
                      action="store", dest="pixel_scale", type=float,
                      help="Pixel scale of the images")

    parser.add_option("-r", "--recursive",
                      action="store_true", dest="recursive", default=False,
                      help="Recursive subdirectories (only first level)")

    (options, args) = parser.parse_args()

    # Logging setup
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='/tmp/field-solver.log',
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # logging = logging.getLogger('field-solver')
    # logging.setLevel(logging.DEBUG)
    # # create file handler which logs even debug messages
    # fh = logging.FileHandler('field-solver.log')
    # fh.setLevel(logging.DEBUG)
    # # create console handler with a higher log level
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.DEBUG)
    # # create formatter and add it to the handlers
    # FORMAT =  "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    # formatter = logging.Formatter(FORMAT)
    # ch.setFormatter(formatter)
    # fh.setFormatter(formatter)
    # # add the handlers to logger
    # logging.addHandler(ch)
    # logging.addHandler(fh)

    logging.debug("Logging setup done !")

    files_solved = []
    files_not_solved = []

    # args is the leftover positional arguments after all options have been processed
    if not options.source_file or len(args) != 0:
        parser.print_help()
        parser.error("incorrect number of arguments ")

    # Check if source_file is a FITS file or a text file listing a set of files
    tic = time.time()
    if os.path.exists(options.source_file):
        if os.path.isfile(options.source_file):
            try:
                hdulist = fits.open(options.source_file)
                filelist = [options.source_file]
            except:
                filelist = [line.replace("\n", "")
                            for line in fileinput.input(options.source_file)]
        elif os.path.isdir(options.source_file):
            filelist = glob.glob(options.source_file + "/*.fit")
            filelist += glob.glob(options.source_file + "/*.fits")
            # Look for subdirectories
            if options.recursive:
                subdirectories = [name for name in os.listdir(options.source_file) if
                                  os.path.isdir(os.path.join(options.source_file, name))]
                for subdir in subdirectories:
                    filelist += glob.glob(os.path.join(options.source_file, subdir) + "/*.fit")
                    filelist += glob.glob(os.path.join(options.source_file, subdir) + "/*.fits")

        # Parallel approach        
        files_solved, files_not_solved = runMultiSolver(filelist, options.output_dir,
                                      options.pixel_scale)
        """for f in filelist:
            new_f = patch_filename(f, rename=False, ext='fits')
            new_f = options.output_dir + "/" + os.path.basename(new_f)
            if new_f not in files_solved and new_f + ".not_science" not in files_solved:
                files_not_solved.append(new_f)
        """
        
    else:
        logging.error("Source file %s does not exists", options.source_file)

    # Clean up the tmp files
    # It seems that astrometry.net (some of its components) does not delete all the temporal
    # files created (ie., tmp.[ppm|wcs|rdls|uncompressed].XXXXXX) created in the tmp directory (-m).
    # So, here we clean up all these files.
    # Note: we are using for the same directory for ouput and tmp files.
    for tmp_file in glob.glob(options.output_dir + "/tmp.*.*"):
        os.remove(tmp_file)

    # Print reports
    toc = time.time()
    log.info("No. files = %s" % len(filelist))
    # calibration files (bias, dark, flats, ...) are considered as solved files
    log.info("No. files solved = %s" % (len(filelist) - len(files_not_solved)))
    log.info("----------------")
    log.info(files_solved)
    log.info("No. files NOT solved = %s", len(files_not_solved))
    log.info("--------------------")
    log.info(files_not_solved)
    log.info("Time : %s" % (toc - tic))

    # Write report files (solved/not_solved) to txt file
    f_solved_txt = "solved-{:%Y-%m-%dT%H:%M:%S}.txt".format(datetime.datetime.now())
    f_not_solved_txt = "not_solved-{:%Y-%m-%dT%H:%M:%S}.txt".format(datetime.datetime.now())
    with open(f_solved_txt, "w") as solved_txt:
        for item in files_solved:
            solved_txt.write("%s\n" % item)
    with open(f_not_solved_txt, "w") as not_solved_txt:
        for item in files_not_solved:
            not_solved_txt.write("%s\n" % item)

    # Finally, exit
    sys.exit()
