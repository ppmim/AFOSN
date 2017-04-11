#! /usr/bin/env python
#encoding:UTF-8

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
import fileinput
import logging
import subprocess
import multiprocessing
import glob
from optparse import OptionParser

import pyfits
import logging as log
import time
import datetime


# TODO
#
# - Comprobacion tipo de imagen no es bias, dark, flat, test  ---DONE
# - Contabilidad de ficheros resueltos y no  ---DONE
# - Tiempo límite de resolución de un fichero --cpulimit (default to 300s)
# - Opcion de añadir header wcs a la cabecera (image.new) 
# - Multiprocessing     --- DONE
# - Estadísicas de errores de calibracion
# - distorsion promedio (encontre un mail de Dustin donde hablaba de eso)
# - limpiar de ficheros temporales/salida creados excepto el .wcs
# - calibracion "fuerza bruta" 
# - log file


# Project modules
import clfits

def readHeader(filename):
    """
    Read from the FITS header values required for astrometric calibration
    """
    
    try:
        fits = clfits.ClFits(filename)
    except Exception, e:
        msg = "Error reading FITS file: " + filename
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
        
        return (scale, ra, dec, instrument, is_science)
        
    
def solveField(filename, tmp_dir, pix_scale=None, blind=False):
    """
    Do astrometric calibration to the given filename using Astrometry.net 
    function 'solve-field'
    
    Parameters
    ----------
    
    pix_scale: float
        Default pixel scale to use in case it cannot be find out from header
    
    tmp_dir: str
        Directory where output files are saved
        
    Returns
    -------
    Filename of solved file (filename.solved) 
    
    """
    
    #
    #Create temporal directory
    #    
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    
    #
    # Read header parameters
    #    
    (scale, ra, dec, instrument, is_science) = readHeader(filename)

    # Whether no scale was found out and some was given as default, we use it
    if scale == -1 and pix_scale != None:
        scale = pix_scale
        
    if not is_science:
        log.info("Frame %s is not a science frame"%filename)
        return filename + ".not_science"
        
    logging.debug("Starting to solve-field for: %s  Scale=%s  RA= %s Dec= %s \
    INSTRUMENT= %s"%(filename, scale, ra , dec, instrument))
    
    # Hardcoded the path !
    path_astrometry = "/usr/local/astrometry/bin"  

    #
    # We must distinguish different cases
    #

    # 0) blind calibration
    if blind:
        log.debug("Blind calibration, nothing is supposed")
        str_cmd = "%s/solve-field -O -p -D %s -m %s %s\
        "% (path_astrometry, tmp_dir, tmp_dir, filename)
        
    # 1) RA, Dec and Scale are known
    elif ra != -1 and dec != -1 and scale != -1:
        log.debug("RA, Dec and Scale are known")
        # To avoid problems with wrong RA,Dec coordinates guessed, a wide 
        # radius is used (0.5 degrees)
        # Although --downsample is used, scale does not need to be modified
        str_cmd = "%s/solve-field --cpulimit 60 -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s --ra %s --dec %s --radius 0.5 -D %s -m %s %s --downsample 2\
        "%(path_astrometry, scale-0.05, scale+0.05, ra, dec, tmp_dir, tmp_dir, filename)
    # 2) RA, Dec are unknown but scale is
    elif (ra == -1 or dec == -1) and scale!=-1:
        log.debug("RA, Dec are unknown but scale is")
        str_cmd = "%s/solve-field -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s -D %s -m %s %s\
        "% (path_astrometry, scale-0.1, scale+0.1, tmp_dir, tmp_dir, filename)
    elif (ra != -1 and dec != -1):
        str_cmd = "%s/solve-field --cpulimit 60 -O -p  \
        --ra %s --dec %s --radius 0.5 -D %s -m %s %s --downsample 2\
        "%(path_astrometry, ra, dec, tmp_dir, tmp_dir, filename)
    # 4) None is known -- blind calibration
    else:
        log.debug("Nothing is known")
        str_cmd = "%s/solve-field -O -p -D %s -m %s %s\
        "% (path_astrometry, tmp_dir, tmp_dir, filename)
    
    log.debug("CMD=" + str_cmd)
    
    clean_output = True
    if clean_output:
        str_cmd += " --corr none --rdls none --match none --wcs none"

    try:
        p = subprocess.Popen(str_cmd, bufsize=0, shell=True, 
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True)
    except Exception, ex:
        log.error("Some error while running subprocess: " + str_cmd)
        log.error(str(ex))
        raise ex

    # Warning:
    # We use communicate() rather than .stdin.write, .stdout.read or .stderr.read 
    # to avoid deadlocks due to any of the other OS pipe buffers filling up and 
    # blocking the child process.(Python Ref.doc)

    (stdoutdata, stderrdata) = p.communicate()
    solve_out =  stdoutdata + "\n" + stderrdata

    #print "STDOUTDATA=",stdoutdata
    #print "STDERRDATA=",stderrdata
    
    if len(solve_out) > 1:
        logging.info("Solve-field output:")
        print solve_out
    #
    # Look for filename.solved to know if field was solved
    #
    solved_file = tmp_dir + "/" + os.path.splitext(os.path.basename(filename))[0] + ".solved"

    # whether we succeeded or failed, clean up
    try:
        os.remove(solved_file.split('.solved')[0] + '.axy')
    except OSError:
        pass    
    # If solved      
    if os.path.exists(solved_file):
        logging.info("Field solved !")
        # clean up        
        try:
            os.remove(solved_file.split('.solved')[0] + '-indx.xyls')
            os.remove(solved_file.split('.solved')[0] + '.solved')
        except OSError:
            pass
        return filename
    # If not solved, then second try blind search (without coordinates)
    else:
        if not blind:
            log.error("First try to solve failed. Lets try again...")
            try:
                return solveField(filename, tmp_dir, 
                                  pix_scale=None, blind=True)
            except Exception, ex:
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
    Run a parallel proceesing to solve astrometry for the input files taking
    advantege of multi-core CPUs
    """

    # use all CPUs available in the computer
    n_cpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=n_cpus)
    
    log.debug("N_CPUS :" + str(n_cpus))
    log.debug("FILES_TO_SOLVE :" + str(files))
    log.debug("TMP_DIR :" + str(tmp_dir))
    
    results = []
    solved = []
    for file in files:
        red_parameters = (file, tmp_dir, pix_scale)
        try:
            # Instead of pool.map() that blocks until
            # the result is ready, we use pool.map_async()
            results += [pool.map_async(calc, [red_parameters])]
        except Exception, ex:
            log.error("Error processing file: " + file)
            log.error(str(ex))
            
    for result in results:
        try:
            result.wait()
            # the 0 index is *ONLY* required if map_async is used !!!
            solved.append(result.get()[0])
        except Exception,e:
            log.error("Cannot process file \n" + str(e))
            
    
    # Here we could try again to solve fields that were not solved but
    # using other parameters or index files


    # Prevents any more tasks from being submitted to the pool. 
    # Once all the tasks have been completed the worker 
    # processes will exit.
    pool.close()

    # Wait for the worker processes to exit. One must call 
    #close() or terminate() before using join().
    pool.join()
    
    log.info("Finished parallel calibration")
    
    return solved

def patch_header(header):
    """
    Add minimal information to the FITS headers.

    Parameters
    ----------
    dir : str, optional
        Directory containing the files to be patched. Default is the current
        directory, ``.``

    new_file_ext : str, optional
        Name added to the FITS files with updated header information. It is
        added to the base name of the input file, between the old file name
        and the `.fit` or `.fits` extension. Default is 'new'.

    save_location : str, optional
        Directory to which the patched files should be written, if not `dir`.

    overwrite : bool, optional
        Set to `True` to replace the original files.

    purge_bad : bool, optional
        Remove "bad" keywords form header before any other processing. See
        :func:`purge_bad_keywords` for details.

    add_time : bool, optional
        If ``True``, add time information (e.g. JD, LST); see
        :func:`add_time_info` for details.

    add_apparent_pos : bool, optional
        If ``True``, add apparent position (e.g. alt/az) to headers. See
        :func:`add_object_pos_airmass` for details.

    add_overscan : bool, optional
        If ``True``, add overscan keywords to the headers. See
        :func:`add_overscan_header` for details.

    fix_imagetype : bool, optional
        If ``True``, change image types to IRAF-style. See
        :func:`change_imagetype_to_IRAF` for details.

    add_unit : bool, optional
        If ``True``, add image unit to FITS header.
    """

    # remove useless keyword
    useless_keyowrds = ['WCSDIM', 'LTM1_1','LTM2_2', 'WAT0_001', 'WAT1_001', 
                       'WAT2_001', 'RADECSYS','SWCREATE','SWOWNER']

    for key in useless_keywords:
        header.pop(key, None)

    # Now, remove all keword starting with '_'
    for key in header:
         if key.startswith('_'):
            header.pop(key, None)

    # Remove all comments starting with 'Original'
    for card in header.cards:
        if card[0]=='COMMENT' and card[1].startswith('Original'): 
            del card
    
###############################################################################
# main
###############################################################################
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
                  "It can be a file or directory name.")
    
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
    if not options.source_file  or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    # Check if source_file is a FITS file or a text file listing a set of files
    tic = time.time()
    if os.path.exists(options.source_file):
        if os.path.isfile(options.source_file):
            try:
                hdulist = pyfits.open(options.source_file)
                filelist = [options.source_file]
            except:
                filelist = [line.replace( "\n", "") 
                            for line in fileinput.input(options.source_file)]
        elif os.path.isdir(options.source_file):
            filelist = glob.glob(options.source_file + "/*.fit")
            filelist += glob.glob(options.source_file + "/*.fits")
            # Look for subdirectories
            if options.recursive:
                subdirectories = [ name for name in os.listdir(options.source_file) if os.path.isdir(os.path.join(options.source_file, name)) ]
                for subdir in subdirectories:
                    filelist += glob.glob(os.path.join(options.source_file, subdir) + "/*.fit")
                    filelist += glob.glob(os.path.join(options.source_file, subdir) + "/*.fits")
                    
                
        # Parallel approach        
        files_solved = runMultiSolver(filelist, options.output_dir, 
                                      options.pixel_scale)
        for file in filelist:
            if file not in files_solved and file + ".not_science" not in files_solved:
                files_not_solved.append(file)
        
        # Serial approach
                    
        #for file in filelist:
        #    try:
        #        solveField(file, options.output_dir, options.pixel_scale)
        #        files_solved.append(file)                
        #    except Exception,e:
        #        files_not_solved.append(file)
        #        logging.error("Error solving file %s  [%s] "%(file,str(e)))
                    
    else:
        logging.error("Source file %s does not exists",options.source_file)

    toc = time.time()

    # print "\n"
    # print "No. files solved = ", len(files_solved)
    # print "------------------------"    
    # print files_solved
    # print "\n"
    # print "No. files not solved = ", len(files_not_solved)
    # print "------------------------"    
    # print files_not_solved
    
    log.info("No. files = %s"%len(filelist))
    # calibracion files (bias, dark, flats, ...) are considered as solved files
    log.info("No. files solved = %s"%(len(filelist)-len(files_not_solved)))
    log.info("----------------")
    log.info(files_solved)
    log.info("No. files NOT solved = %s", len(files_not_solved))
    log.info("--------------------")
    log.info(files_not_solved)
    log.info("Time : %s"%(toc-tic))
    
    # Write files (solved/not_solved) to txt file
    f_solved_txt = 'solved-{:%Y-%m-%dT%H:%M:%S}.txt'.format(datetime.datetime.now())
    f_not_solved_txt = 'not_solved-{:%Y-%m-%dT%H:%M:%S}.txt'.format(datetime.datetime.now())
    with open(f_solved_txt, "w") as solved_txt:
        for item in files_solved:
            solved_txt.write("%s\n" % item)
    with open(f_not_solved_txt, "w") as not_solved_txt:
        for item in files_not_solved:
            not_solved_txt.write("%s\n" % item)
    
    sys.exit()
        
    
