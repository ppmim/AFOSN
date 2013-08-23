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
import shutil
import tempfile
from optparse import OptionParser
import fileinput
import logging
import subprocess
import glob

import numpy
import pywcs
import pyfits
import sys
import logging as log


# TODO
#
# - Comprobacion tipo de imagen no es bias, dark, flat, test  ---DONE
# - Contabilidad de ficheros resueltos y no  ---DONE
# - Tiempo límite de resolución de un fichero --cpulimit (default to 300s)
# - Opcion de añadir header wcs a la cabecera (image.new) 
# - Multiprocessing
# - Estadísicas de errores de calibracion
# - distorsion promedio (encontre un mail de Dustin donde hablaba de eso)
# -


# Project modules
import clfits

def readHeader(filename):
    """
    Read from the FITS header values required for astrometric calibration
    """
    
    try:
        fits = clfits.ClFits(filename)
    except Exception,e:
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
        


def initWCS( input_image, pixel_scale):
    """
    Call this routine to write rough WCS into FITS header and update RA,DEC
    coordinates to J2000.0 equinox; and this way allow SCAMP/Astrometry.net 
    make astrometry with external Catalogs (in J2000??).
    
    Warning: The original file header (input_image) will be modified
    """
    
    try:
        f = datahandler.ClFits ( input_image )
    except Exception,e:
        raise e
    
    fits_file = pyfits.open(input_image, 'update', ignore_missing_end=True)

    if f.isMEF(): # is a MEF
        raise Exception("Sorry, currently this function only works with simple "
                        "FITS files with no extensions")
    else:  # is a simple FITS
        header = fits_file[0].header
        try:
            checkWCS(header)
            log.debug("FITS looks having a right WCS header")
        except Exception,e:
            log.debug("No WCS compliant header, trying to create one ...")
            try:
                # Read some basic values
                naxis1 = f.getNaxis1()
                naxis2 = f.getNaxis2()
                ra = f.ra
                dec = f.dec
                equinox0 = f.getEquinox()
                # 
                # Transform RA,Dec to J2000 -->fk5prec(epoch0, 2000.0, &ra, &dec);
                # EQUINOX precessing is DONE by SCAMP !!!
                WCS_J2000 = 1  #J2000(FK5) right ascension and declination
                WCS_B1950 = 2  #B1950(FK4) right ascension and declination
                #[new_ra, new_dec]=wcscon.wcscon(WCS_J2000, WCS_J2000, equinox0, 2000.0, ra, dec, 0)
                # Find out PIXSCALE
                if "PIXSCALE" in header:
                    scale = header['PIXSCALE']
                    degscale = scale/3600.0
                else:
                    scale = pixel_scale
                    degscale = scale/3600.0
                    log.warning("Cannot find out the PIXSCALE for the image.")
                    log.warning("Using Pixel scale = %s"%pixel_scale)
                    #fits_file.close()
                    #raise Exception("Cannot find out PIXSCALE for the image")

                create_wcs = True                
                if create_wcs:
                    #Create initial WCS
                    #
                    header.update("CRPIX1", naxis1/2.0, "Ref. pixel in <axis direction>")
                    header.update("CRPIX2", naxis2/2.0, "Ref. pixel in <axis direction>")
                    header.update("CRVAL1", ra, "Coordinate value of ref. pixel")
                    header.update("CRVAL2", dec, "Coordinate value of ref. pixel")
                    #header.update("RA", new_ra, "Coordinate value of ref. pixel")
                    #header.update("DEC", new_dec, "Coordinate value of ref. pixel")
                    header.update("CTYPE1", "RA---TAN", "Pixel coordinate system")
                    header.update("CTYPE2", "DEC--TAN", "Pixel coordinate system")
                    #header.update("RADECSYS","FK5","Coordinate reference frame")
                    # CD matrix (the CDi_j elements) encode the sky position angle,
                    # the pixel scale, and a possible flipping.
                    # CD1_1 is <0 because East is supposed at Left = flipX
                    # CD2_2 is >0 because North is supposed at Up
                    # In addition, it must be noted that:
                    # CD1_1 = cos(r), CD1_2 = sin(r), CD2_1 = -sin(r), CD2_2 = cos(r)
                    # r = clockwise rotation_angle  
                    header.update("CD1_1", -degscale, "Translation matrix element")
                    header.update("CD1_2", 0.0, "Translation matrix element")
                    header.update("CD2_1", 0.0, "Translation matrix element")
                    header.update("CD2_2", degscale, "Translation matrix element")
                    header.update("SCALE", scale, "Image scale")
                    #header.update("EQUINOX", 2000.0, "Standard FK5(years)")
                else:
                    header.update("RA", ra, "Right Ascension (degree)")
                    header.update("DEC", dec, "Declination (degree)")
                    
                # clean incompatible CDi_j and CDELT matrices
                if "CDELT1" in header:
                    del header["CDELT1"]
                if "CDELT2" in header:
                    del header["CDELT2"]
                
                log.debug("Successful WCS header created !")
                
            except Exception,e:
                log.error("Some error while creating initial WCS header: %s"%str(e))
                fits_file.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                raise e
        
        fits_file.close(output_verify='ignore')
        log.debug("Right WCS info")
            
def checkWCS( header ):
    """
    Checks for a variety of WCS keywords and raise an Exception if the header 
    lacks a proper combination of them.  This is needed because wcstools will 
    not raise any sort of error if a WCS isn't present or is malformed, and 
    SCAMP (E.Bertin) need a initial WCS information. This is probably 90% 
    complete in terms of its checking for the types of FITS files that we are 
    likely to be using.
    
    If you find any WCS keywords that cause wcstools to behave in an erratic
    manner without signaling errors, add them to this method.  Experience has
    shown that the astrophysical community has an uncanny ability to produce
    data sets that cause FITS readers and WCS projections to break.  It is
    important that we check for irregular cases and flag them before the code
    runs and produces confusing results.  Our only defense against this is
    experience with unusual data sets, so the more checks here the better.
    
    TODO(): Implement stricter checking under CDELT case, in
    particular for full PC matrices (both kinds), as well as LATPOLE and
    LONPOLE.  It would also be good to check for illegal values, but that's
    a lot of work.
    """
    
    keywords_to_check=['NAXIS1','NAXIS2','CTYPE1','CTYPE2','CRVAL1','CRVAL2',
                       'CRPIX1','CRPIX2']


    raise Exception("Forzamos la creacion de un nuevo header")

    # Every header must have these keywords.
    for kw in keywords_to_check:
        if kw not in header:
            log.debug("Keyword %s not found",kw)
            raise Exception("Keyword %s not found"%kw)
    
    # Check for the equinox, which can be specified in more than 1 way.
    if 'EPOCH' not in header and 'EQUINOX' not in header:
        log.debug("Missing keyword EPOCH or EQUINOX")
        raise Exception("Missing keyword EPOCH or EQUINOX")
        
    # Check some values
    if header['CTYPE1']=='PIXEL' or header['CTYPE2']=='PIXEL':
        log.debug("Wrong CTYPE value (PIXEL-Cartesian Coordinates) for WCS header")
        raise Exception ("Wrong CTYPE value (PIXEL-Cartesian Coordinates) for WCS header")
        
    # Check for CDi_j or CDELT matrix
    # CDELT matrix : Here we should probably be more rigorous and check
    # for a full PC matrix or CROTA value, but for now
    # this is pretty good.
    if 'CD1_1' not in header or 'CD1_2' not in header \
        or 'CD2_1' not in header or 'CD2_2' not in header:
            if 'CDELT1' not in header or 'CDELT2' not in header:
                log.debug("Couldn't find a complete set of CDi_j matrix or CDELT")
                raise Exception("Couldn't find a complete set of CDi_j matrix or CDELT")
                

    
def solveField(filename, tmp_dir, pix_scale=None):
    """
    Do astrometric calibration to the given filename using Astrometry.net 
    function 'solve-field'
    
    Parameters
    ----------
    
    pix_scale: float
        Default pixel scale to use in case it cannot be find out from header
    
    tmp_dir: str
        Directory where output files are saved
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
    if scale==-1 and pix_scale!=None:
        scale = pix_scale
        
    if not is_science:
        log.info("Frame %s is not a science frame"%filename)
        return 1
        
    logging.debug("Starting to solve-field for: %s  Scale=%s  RA= %s Dec= %s \
    INSTRUMENT= %s"%(filename, scale, ra , dec, instrument))
    
    path_astrometry = "/usr/local/astrometry/bin"  

    #
    # We must distinguish different cases
    #

    # 1) RA, Dec and Scale are known
    if ra!=-1 and dec!=-1 and scale!=-1:
        log.debug("RA, Dec and Scale are known")
        # To avoid problems with wrong RA,Dec coordinates guessed, a wide 
        # radius is used (0.5 degrees)
        str_cmd = "%s/solve-field -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s --ra %s --dec %s --radius 0.5 -D %s %s\
        "%(path_astrometry, scale-0.05, scale+0.05, ra, dec, tmp_dir, filename)
    # 2) RA, Dec are unknown but scale is
    elif ra==-1 or dec==-1:
        log.debug("RA, Dec are unknown but scale is")
        str_cmd = "%s/solve-field -O -p --scale-units arcsecperpix --scale-low %s \
        --scale-high %s -D %s %s\
        "%(path_astrometry, scale-0.1, scale+0.1, tmp_dir, filename)
    # 3) None is known -- blind calibration
    if (ra==-1 or dec==-1) and scale==-1:
        log.debug("Nothing is known")
        str_cmd = "%s/solve-field -O -p -D %s %s\
        "%(path_astrometry, tmp_dir, filename)
    
    log.debug("CMD="+str_cmd)
    
    try:
        p = subprocess.Popen(str_cmd, bufsize=0, shell=True, 
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True)
    except Exception, e:
        log.error("Some error while running subprocess: " + str_cmd)
        log.error(str(e))
        raise e

    # Warning:
    # We use communicate() rather than .stdin.write, .stdout.read or .stderr.read 
    # to avoid deadlocks due to any of the other OS pipe buffers filling up and 
    # blocking the child process.(Python Ref.doc)

    (stdoutdata, stderrdata) = p.communicate()
    solve_out =  stdoutdata + "\n" + stderrdata

    #print "STDOUTDATA=",stdoutdata
    #print "STDERRDATA=",stderrdata
    
    if len(solve_out)>1:
        logging.info("Solve-field output:")
        print solve_out
    #
    # Look for filename.solved to know if field was solved
    #
    solved_file = tmp_dir + "/" + os.path.splitext(os.path.basename(filename))[0] + ".solved"
    #print "FILE=",solved_file
    if os.path.exists(solved_file):
        logging.info("Field solved !")
        return 1
    else:
        log.error("Field was not solved.")
        raise Exception("Field was not solved")
            

def runMultiSolver(files, tmp_dir, pix_scale=None):
    """
    Run a parallel proceesing to solve astrometry for the input files taking
    advantege of multi-core CPUs
    """
    
    
                  
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
    
                                
    (options, args) = parser.parse_args()
    
    
    ## Logging setup
    FORMAT =  "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging_level = logging.DEBUG
    logging.basicConfig(format = FORMAT, level = logging_level)
    logging.debug("Logging setup done !")
    
    files_solved = []
    files_not_solved = []
    
    # args is the leftover positional arguments after all options have been processed
    if not options.source_file  or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    # Check if source_file is a FITS file or a text file listing a set of files
    if os.path.exists(options.source_file):
        if os.path.isfile(options.source_file):
            try:
                hdulist = pyfits.open(options.source_file)
                filelist = [options.source_file]
            except:
                filelist = [line.replace( "\n", "") 
                            for line in fileinput.input(options.source_file)]
        elif os.path.isdir(options.source_file):
            filelist = glob.glob(options.source_file+"/*.fit*")
                        
        for file in filelist:
            try:
                solveField(file, options.output_dir, options.pixel_scale)
                files_solved.append(file)                
            except Exception,e:
                files_not_solved.append(file)
                logging.error("Error solving file %s  [%s] "%(file,str(e)))
                    
    else:
        logging.error("Source file %s does not exists",options.source_file)
        
    print "\n"
    print "No. files solved = ", len(files_solved)
    print "------------------------"    
    print files_solved
    print "\n"
    print "No. files not solved = ", len(files_not_solved)
    print "------------------------"    
    print files_not_solved
    
    
    sys.exit()
        
    
