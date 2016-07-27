import datetime as dt
from davitpy.pydarn.sdio.fetchUtils import fetch_local_files
from davitpy.pydarn.sdio import radDataOpen, radDataReadRec
from davitpy.pydarn.sdio import radDataPtr
import davitpy
import logging
import os
import string
import matplotlib.pyplot as plt

# input parameters
ctr_time = dt.datetime(2010,1,15)
rad = "bks"
channel = None
bmnum = 7
ftype = "fitacf"
filtered = True
scr = "local"
localdirfmt = "/sd-data/{year}/{ftype}/{radar}/"
localdict = {"ftype" : "fitacf", "radar" : "bks", "channel" : None}
tmpdir = "/tmp/sd/"
fnamefmt = ['{date}.{hour}......{radar}.{channel}.{ftype}', '{date}.{hour}......{radar}.{ftype}']
#davitpy.rcParams['verbosity'] = "debug"


def fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt):
    
    # expend the time to three days
    stime = ctr_time - dt.timedelta(days=1)
    etime = ctr_time + dt.timedelta(days=1)
    radcode = localdict["radar"]
    ftype = localdict["ftype"]
    channel = localdict["channel"]

    # fetch the data for three days
    file_list = fetch_local_files(stime, etime, localdirfmt, localdict, tmpdir, fnamefmt)

    # check if we have found files
    if len(file_list) != 0:
        # concatenate the files into a single file
        logging.info('Concatenating all the files in to one')
        # choose a temp file name with time span info for cacheing
        if (channel is None):
            tmp_name = '%s%s.%s.%s.%s.%s.%s' % \
                      (tmpdir, stime.strftime("%Y%m%d"),
                       stime.strftime("%H%M%S"),
                       etime.strftime("%Y%m%d"),
                       etime.strftime("%H%M%S"), radcode, ftype)
        else:
            tmp_name = '%s%s.%s.%s.%s.%s.%s.%s' % \
                      (tmpdir, stime.strftime("%Y%m%d"),
                       stime.strftime("%H%M%S"),
                       etime.strftime("%Y%m%d"),
                       etime.strftime("%H%M%S"),
                       radcode, channel, ftype)
        logging.debug('cat ' + string.join(file_list) + ' > ' + tmp_name)
        os.system('cat ' + string.join(file_list) + ' > ' + tmp_name)

        # remove the unneeded files from the tmpdir
        for file_name in file_list:
            logging.debug('rm ' + file_name)
            os.system('rm ' + file_name)
            os.system('rm ' + file_name+".bz2")
            logging.info("removed unneeded files")
    else:
        tmp_name = None
        
    fname = tmp_name
    return fname

def boxcar_filter(fname):

    # do the boxcar filter
    if fname is not None:
        ftype = fname.split(".")[-1]
        if not ftype+'f' in fname:
            try:
                ffname = fname + 'f'
                command = '/davit/lib/rst/rst/bin/fitexfilter ' + fname + ' > ' + ffname
                logging.debug("performing: {:s}".format(command))
                os.system(command)
                logging.info("done filtering")
            except Exception, e:
                estr = 'problem filtering file, using unfiltered'
                logging.warning(estr)
                #ffname = fname
        else:
            ffname = fname

        os.system("mv /tmp/sd/" + ffname.split("/")[-1] + " ./data/")
        ffname = "./data/" + ffname.split("/")[-1]


    return ffname

def cluster_data(myPtr, bmnum, params=["velocity"], tbands=None):
    """Reads data from the file pointed to by myPtr

    Parameter
    ---------
    myPtr :
        a davitpy file pointer object
    myBeam : 
        a davitpy beam object
    bmnum : int
        beam number of data to read in
    params : list
        a list of the parameters to read
    tbands : list
        a list of the frequency bands to separate data into

    Returns
    -------
    A dictionary of the data. Data is stored in lists and separated in
    to tbands.

    Example
    -------
        from davitpy import pydarn
        from datetime import datetime
        myPtr = pydarn.sdio.radDataOpen(datetime(2012,11,24),'sas')
        myBeam = myPtr.readRec()
        data_dict = read_data(myPtr, myBeam, 7, ['velocity'], [8000,20000])

    Written by Muhammad 20160722

    """

    if tbands is None:
        tbands = [8000, 20000]

    # Initialize some things.
    data = dict()
    data_keys = ['vel', 'pow', 'wid', 'elev', 'phi0', 'times', 'freq', 'cpid',
                 'nave', 'nsky', 'nsch', 'slist', 'mode', 'rsep', 'nrang',
                 'frang', 'gsflg', 'velocity_error', "clust_num"]
    for d in data_keys:
        data[d] = []

    # Read the parameters of interest.
    myPtr.rewind()
    myBeam = myPtr.readRec()
    while(myBeam is not None):
        if(myBeam.time > myPtr.eTime): break
        if(myBeam.bmnum == bmnum and (myPtr.sTime <= myBeam.time)):
            if (myBeam.prm.tfreq >= tbands[0] and
                myBeam.prm.tfreq <= tbands[1]):
                data['times'].append(myBeam.time)
                data['cpid'].append(myBeam.cp)
                data['nave'].append(myBeam.prm.nave)
                data['nsky'].append(myBeam.prm.noisesky)
                data['rsep'].append(myBeam.prm.rsep)
                data['nrang'].append(myBeam.prm.nrang)
                data['frang'].append(myBeam.prm.frang)
                data['nsch'].append(myBeam.prm.noisesearch)
                data['freq'].append(myBeam.prm.tfreq / 1e3)
                data['mode'].append(myBeam.prm.ifmode)
                data['gsflg'].append(myBeam.fit.gflg)
                data['slist'].append(myBeam.fit.slist)
                # To save time and RAM, only keep the data specified
                # in params.
                if('velocity' in params):
                    data['vel'].append(myBeam.fit.v)
                if('power' in params):
                    data['pow'].append(myBeam.fit.p_l)
                if('width' in params):
                    data['wid'].append(myBeam.fit.w_l)
                if('elevation' in params):
                    data['elev'].append(myBeam.fit.elv)
                if('phi0' in params):
                    data['phi0'].append(myBeam.fit.phi0)
                if('velocity_error' in params):
                    data['velocity_error'].append(myBeam.fit.v_e)

        myBeam = myPtr.readRec()

    return data

def dopsearch(ctr_time, bmnum, localdirfmt, localdict, tmpdir, fnamefmt):

    # fetch and concatenate the three consecutive days of data centered on the target date 
    #concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)

    # box car fiter the data
    #ffname = boxcar_filter(concated_file)

    ffname = "./data/20100114.000000.20100116.000000.bks.fitacff"

    import sys
    sys.path.append("/home/muhammad/softwares/davitpy_MuhammadVT/davitpy/pydarn/plotting/")
    from rti import plot_rti

    stm = ctr_time - dt.timedelta(days=1)
    etm = ctr_time + dt.timedelta(days=1)

    #plot_rti(stm, "bks", eTime=etm, bmnum=7)
    #plot_rti(stm, "bks", eTime=etm, bmnum=7, fileName=ffname)
    #plot_rti(dt.datetime(2013,3,16), 'bks', eTime=dt.datetime(2013,3,16,14,30), bmnum=12)

    ftype = "fitacf"
    scales = [[-150, 150]]
    #ftype = "fitex"
    myPtr = radDataOpen(stm, "bks", eTime=etm, bmnum=bmnum, fileName=ffname, fileType=ftype)
    #myPtr = radDataOpen(stm, "bks", eTime=etm, bmnum=7, fileType="fitacf")
    #plot_rti(stm, 'bks', eTime=etm, bmnum=7, fileType=ftype, myFile=myPtr,
    #         fileName=ffname, params=["velocity"], coords='rng', gsct=True, scales=scales)

    #plot_rti(stm, 'bks', eTime=etm, bmnum=7, fileType=ftype, params=["velocity"],
    #          coords='rng', gsct=True, scales=scales, filtered=True)

    list_of_clusters = cluster_data(myPtr, bmnum, params=["velocity"], tbands=None)
    N = 1
    return list_of_clusters


data = dopsearch(ctr_time, bmnum, localdirfmt, localdict, tmpdir, fnamefmt)

#    #myPtr = radDataOpen(stime, rad, etime, channel=channel, bmnum=bmnum, fileType=ftype,
#    #                    filtered=filtered, local_dirfmt=localdirfmt)
#
#    myPtr = radDataPtr(stime, rad, etime, channel=channel, bmnum=bmnum, fileType=ftype,
#                        filtered=filtered, local_dirfmt=localdirfmt)
#

#concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)
#ffname = boxcar_filter(concated_file)

#from davitpy.pydarn.plotting import plot_rti



