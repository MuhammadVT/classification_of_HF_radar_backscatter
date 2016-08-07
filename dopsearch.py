import datetime as dt
from davitpy.pydarn.sdio.fetchUtils import fetch_local_files
from davitpy.pydarn.sdio import radDataOpen, radDataReadRec
from davitpy.pydarn.sdio import radDataPtr
import davitpy
import logging
import os
import string
import matplotlib.pyplot as plt
from glob import glob
import pdb
import numpy as np

def fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt):

    """ fetches files for three days centered at ctr_time.day, then unzips and concatenates
    them into a single file """
    
    # expend the time to three days
    stime = ctr_time - dt.timedelta(days=1)
    etime = ctr_time + dt.timedelta(days=2)
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
            #os.system('rm ' + file_name+".gz")
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
            print "file " + fname + " exists"
            ffname = fname

    return ffname

def prepare_file(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt):
    """ A wrapper for file fetching and boxcar filtering"""

    # fetch and concatenate the three consecutive days of data centered on the target date 
    concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)

    # box car fiter the data
    ffname = boxcar_filter(concated_file)

    return ffname


def read_data_for_rtiplot(myPtr, bmnum, params=["velocity"], tbands=None, coords="geo"):

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
    coords : string 
        converts the range-time cell position (clat, clon) into the value 
        given by coords. Has to be one of ["mag", "geo", "mlt"]
        (Note: only works for "geo" so far due to speed issue)

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
        data_dict = read_data(myPtr, bmnum, params=['velocity'], [8000,20000])

    Written by Muhammad 20160722

    """

    from davitpy import pydarn

    if tbands is None:
        tbands = [8000, 20000]

    # Initialize some things.
    data = dict()

    # use the following data_keys if you want to plot the data using rtiplot function
    data_keys = ['vel', 'pow', 'wid', 'elev', 'phi0', 'times', 'freq', 'cpid',
                 'nave', 'nsky', 'nsch', 'slist', 'mode', 'rsep', 'nrang',
                 'frang', 'gsflg', 'velocity_error', 'bmazm']
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
                data['bmazm'].append(myBeam.prm.bmazm)
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

#                # save the center lat, lon position of each scatter. This takes too much time
#                site = pydarn.radar.network().getRadarById(myBeam.stid) \
#                    .getSiteByDate(myBeam.time)
#                myFov = pydarn.radar.radFov.fov(site=site, ngates=myBeam.prm.nrang,
#                                                nbeams=site.maxbeam,
#                                                rsep=myBeam.prm.rsep,
#                                                coords=coords,
#                                                date_time=myBeam.time)
#                data['clats'].append(myFov.latCenter[bmnum][myBeam.fit.slist])
#                data['clons'].append(myFov.lonCenter[bmnum][myBeam.fit.slist])
                
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

    return {bmnum:data}



def read_data(myPtr, params=["velocity"], tbands=None, coords="geo"):

    """Reads data from the file pointed to by myPtr

    Parameter
    ---------
    myPtr :
        a davitpy file pointer object
    myBeam : 
        a davitpy beam object
    params : list
        a list of the parameters to read
    tbands : list
        a list of the frequency bands to separate data into
    coords : string 
        converts the range-time cell position (clat, clon) into the value 
        given by coords. Has to be one of ["mag", "geo", "mlt"]
        (Note: only works for "geo" so far due to speed issue)

    Returns
    -------
    A list of dicts of the data. Each dict element stores data for a certain beam.
    Data in each dict is stored in lists and separated in to tbands.

    Example
    -------
        from davitpy import pydarn
        from datetime import datetime
        myPtr = pydarn.sdio.radDataOpen(datetime(2012,11,24),'sas')
        myBeam = myPtr.readRec()
        data_dict = read_data(myPtr, params=['velocity'], [8000,20000])

    Written by Muhammad 20160722

    """

    from davitpy import pydarn
    import copy

    if tbands is None:
        tbands = [8000, 20000]

    # Initialize some things.
    data = dict()

    plotrti=False

    if not plotrti:
        # use the following data_keys if you do not want to plot the data using rtiplot function
        data_keys = ['vel', 'times', 'slist', 'rsep', 'nrang', 'frang', 'gsflg', 'bmazm']

    else:

        # use the following data_keys if you want to plot the data using rtiplot function
        data_keys = ['vel', 'pow', 'wid', 'elev', 'phi0', 'times', 'freq', 'cpid',
                     'nave', 'nsky', 'nsch', 'slist', 'mode', 'rsep', 'nrang',
                     'frang', 'gsflg', 'velocity_error', 'bmazm', 'clats', 'clons']

    for d in data_keys:
        data[d] = []

    # list of dicts. each dict stores data for a certain bmnum
    max_nbms = 30
    all_beams = dict()
    all_beams = {bm:copy.deepcopy(data) for bm in xrange(max_nbms)}      

    # Read the parameters of interest.
    myPtr.rewind()
    myBeam = myPtr.readRec()
    while(myBeam is not None):
        if(myBeam.time > myPtr.eTime): break
        if(myPtr.sTime <= myBeam.time):
            if (myBeam.prm.tfreq >= tbands[0] and
                myBeam.prm.tfreq <= tbands[1]):
                bmnum = myBeam.bmnum
                all_beams[bmnum]['times'].append(myBeam.time)
                all_beams[bmnum]['bmazm'].append(round(myBeam.prm.bmazm,2))
                all_beams[bmnum]['rsep'].append(myBeam.prm.rsep)
                all_beams[bmnum]['nrang'].append(myBeam.prm.nrang)
                all_beams[bmnum]['frang'].append(myBeam.prm.frang)
                all_beams[bmnum]['gsflg'].append(myBeam.fit.gflg)
                all_beams[bmnum]['slist'].append(myBeam.fit.slist)

####################################################################

                if plotrti:
                    all_beams[bmnum]['cpid'].append(myBeam.cp)
                    all_beams[bmnum]['nave'].append(myBeam.prm.nave)
                    all_beams[bmnum]['nsky'].append(myBeam.prm.noisesky)
                    all_beams[bmnum]['nsch'].append(myBeam.prm.noisesearch)
                    all_beams[bmnum]['freq'].append(myBeam.prm.tfreq / 1e3)
                    all_beams[bmnum]['mode'].append(myBeam.prm.ifmode)

####################################################################


                # To save time and RAM, only keep the data specified
                # in params.
                if('velocity' in params):
                    all_beams[bmnum]['vel'].append(myBeam.fit.v)

        myBeam = myPtr.readRec()

    return all_beams

def read_file(ffname, rad, stm, etm, params, ftype="fitacf"):

    """ A wrapper for reading a file. It reads all beams at once 
    """

    myPtr = radDataOpen(stm, rad, eTime=etm, fileName=ffname, fileType=ftype)
    data_dict = read_data(myPtr, params=params, tbands=None)

    return data_dict 


def read_file_for_rtiplot(ffname, rad, stm, etm, bmnum, params, ftype="fitacf"):

    """ A wrapper for reading a file. this will read data for a certain beam 
        so that rtiplot function can plot the data.
    """

    myPtr = radDataOpen(stm, rad, eTime=etm, fileName=ffname, fileType=ftype)
    data_dict = read_data_for_rtiplot(myPtr, bmnum, params=params, tbands=None)

    return data_dict 

def create_nodes(data_dict):
    """ Create nodes using time indices and gate numbers from the data_dict. 
    Nondes are list of lists. Each list element is a collection of nodes for a given time_index
    A node is a range-gate cell in the rti plot.
    """

    nodes = [[(x,y) for y in data_dict['slist'][x]] for x in xrange(len(data_dict['times']))]

    return nodes

def find_start_node(nodes, visited_nodes=None):

    if visited_nodes is None:
        visited_nodes = set()
    start_node = None
    for sublist in nodes:
        for itm in sublist:
            if itm[1] >=7 and (itm not in visited_nodes):
                start_node = itm
                break
        if start_node is not None:
            break
    return start_node


def create_graph(vertex, data):
    """ all nodes should be fed to data argument """

    tm_indx = vertex[0]
    if tm_indx == 0:
        xx = [0, 1]   
    elif tm_indx == len(data)-1:
        xx = [tm_indx, tm_indx-1]  # time indices centered at the time index of the vertex 
    else:   
        xx = [vertex[0] + yi for yi in [-1, 0, 1]]  # time indices centered at the time index of the vertex 
    yy = [vertex[1] + yi for yi in [-1, 0, 1]]  #range-gate indices centered at the range-gate index of the vertex 

    # create a tuple that inclues actual neighbering vortice
    xy = [(i,j) for i in xx for j in set([x[1] for x in data[i]]).intersection(set(yy))] 
    xy.remove(vertex)    # remove the vertex from xy 
    adjacent_nodes = set()
    for tpl in xy:
        adjacent_nodes.add(tpl)
    G = {vertex:adjacent_nodes}
    return G


def search_tree(data, start, data_dict):
    """ all nodes should be fed to data argument """

    # create a grash as a start
    G = create_graph(start, data)
    # do the breath_firt_search
    visited = set()
    queue = [start]
    while queue:
        vertex = queue.pop(0)

        if (vertex[1] >= 7) and (vertex not in visited):
            visited.add(vertex)
            queue.extend(G[vertex] - visited)
        try:
            next_node = queue[0]
            G = create_graph(next_node, data)
        except:

            end_indx = max([x[0] for x in visited])
            gates_tmp = [x[1] for x in visited if x[0]==end_indx]
            possible_gates = set([x for y in gates_tmp \
                                for x in [y-1, y, y+1]])

            k = 2
            end_tm = data_dict['times'][end_indx]
            del_time = 0

            while del_time<6.0:
            #while k<=6:
                # use try-except to avoid hitting the end of time indices 
                try:
                    next_gates_tmp = [x[1] for x in data[end_indx+k]]
                    actual_gates = possible_gates.intersection(set(next_gates_tmp))
                    if (len(actual_gates)>0) and (max(actual_gates)>=7):
                        queue = [x for x in data[end_indx+k] if x[1] in actual_gates]
                        next_node = queue[0]
                        G = create_graph(next_node, data)
                        break

                    tm_next = data_dict['times'][end_indx + k]
                    del_time = round((tm_next - end_tm).total_seconds() / 60.)

                    k += 1
                except:
                    break
                    
    return visited

def push_stm_etm(cluster, data_dict, vel_threshold=15.):
    import datetime as dt
    # write cluster as list of lists. Each list element stors the data for a given time
    tm_indices = sorted(list(set([x[0] for x in cluster])))

    # change cluster to a list of lists
    cluster_lol = [[x for x in cluster if y==x[0]] for y in tm_indices]

    # initialize the time indices
    stm_indx = tm_indices[0]
    etm_indx = tm_indices[-1]
    stm = data_dict['times'][stm_indx]
    etm = data_dict['times'][etm_indx]

    # check the time duration
    tm_del = etm - stm
    if tm_del <= dt.timedelta(hours=1):
        return cluster

    gates_width = 4
    update_stm = True
    update_etm = True

    # initialize cluster_lol indices
    sindx = 0
    eindx = len(cluster_lol)-1 
    for ii in xrange(len(tm_indices)):
        # determine the starting time of the cluster 
        if update_stm:
            cluster_left = cluster_lol[0:ii+gates_width+1] if ii < gates_width \
                    else cluster_lol[ii-gates_width:ii+gates_width+1]
            # flatten cluster_left
            cluster_left = [x for y in cluster_left for x in y]
            
            # get the velocities
            vels_left = [data_dict['vel'][item[0]][(data_dict['slist'][item[0]]).index(item[1])] \
                    for item in cluster_left]

            high_vels_num_left = len([x for x in vels_left if abs(x) > vel_threshold])
            low_vels_num_left = len(vels_left) - high_vels_num_left

            # exclude the case where low_vels_num is 0
            try:
                high_to_low_ratio_left = (high_vels_num_left *1.0) / low_vels_num_left
            except:
                high_to_low_ratio_left = 10

            # update the indices
            if high_to_low_ratio_left <= 0.75:
                #sleft_indx = tm_indices[0] if (ii<gates_width) else tm_indices[ii-gates_width]
                sindx = ii+1
                stm_indx = tm_indices[sindx] 
                #sright_indx = tm_indices[ii+1+gates_width]
            else:
                update_stm = False 

        # determine the ending time of the cluster 
        if update_etm:
            cluster_right = cluster_lol[-ii-gates_width-1:] if ii < gates_width+1 \
                    else cluster_lol[-ii-gates_width-1:-ii+gates_width]

            # flatten cluster_right
            cluster_right = [x for y in cluster_right for x in y]

            vels_right = [data_dict['vel'][item[0]][(data_dict['slist'][item[0]]).index(item[1])] \
                    for item in cluster_right]

            high_vels_num_right = len([x for x in vels_right if abs(x) > vel_threshold])
            low_vels_num_right = len(vels_right) - high_vels_num_right

            # exclude the case where low_vels_num is 0
            try:
                high_to_low_ratio_right = (high_vels_num_right *1.0) / low_vels_num_right
            except:
                high_to_low_ratio_right = 10

            # update the indices
            if high_to_low_ratio_right <= 0.75:
                #eright_indx = tm_indices[-1] if (ii<gates_width+1) else tm_indices[-ii+gates_width]
                eindx = (len(cluster_lol)-1) - ii - 1
                etm_indx = tm_indices[eindx - ii - 1] 
                #eleft_indx = tm_indices[-ii-gates_width-1]
            else:
                update_etm = False 


        # check the time duration of the cluster
        stm = data_dict['times'][stm_indx]
        etm = data_dict['times'][etm_indx]
        tm_del = etm - stm
        if tm_del <= dt.timedelta(hours=1):
            break

    # update cluster_lol
    cluster_lol = cluster_lol[sindx:eindx+1]
    # flatten cluster_lol and convert to a set
    cluster = set([x for y in cluster_lol for x in y])

    return cluster

def isevent(cluster, data_dict, vel_threshold=15.):
    import datetime as dt
    # find time indices
    tm_indices = sorted(list(set([x[0] for x in cluster])))

    stm_indx = tm_indices[0]
    etm_indx = tm_indices[-1]
    stm = data_dict['times'][stm_indx]
    etm = data_dict['times'][etm_indx]
    tm_del = etm - stm
    result = False
    if tm_del <= dt.timedelta(hours=1):
        pass

    elif tm_del >= dt.timedelta(hours=14):
        pass
    else:
        cluster_vels = [data_dict['vel'][item[0]][(data_dict['slist'][item[0]]).index(item[1])] \
                for item in cluster]

        high_vels_num = len([x for x in cluster_vels if abs(x) > vel_threshold])
        low_vels_num = len(cluster_vels) - high_vels_num

        # exclude the case where low_vels_num is 0
        try:
            high_to_low_ratio = (high_vels_num *1.0) / low_vels_num
        except:
            high_to_low_ratio = 10

        if tm_del <= dt.timedelta(hours=2):
            if high_to_low_ratio > 0.475:
                result = True
        elif tm_del <= dt.timedelta(hours=3):
            if high_to_low_ratio > 0.33:
                result = True
        #elif tm_del < dt.timedelta(hours=14):
        else :
            if high_to_low_ratio > 0.2:
                result = True

    return result
         

def is_low_vel_event(cluster, data_dict, vel_lim=120):
    import datetime as dt

    result = False

    # low_vel_event has to be an iscat event
    if isevent(cluster, data_dict, vel_threshold=15.):

        cluster_vels = sorted([data_dict['vel'][item[0]][(data_dict['slist'][item[0]]).\
                index(item[1])] for item in cluster])

        # find indices that corresponds to 3 and 97 percentile 
        third_percentile = np.percentile(cluster_vels, 3)
        nightyseventh_percentile = np.percentile(cluster_vels, 97)

        if (third_percentile>-vel_lim and nightyseventh_percentile<vel_lim):
            result = True

    return result
       

def change_gsflg(cluster, data_dict, gscat_value=0):
    for tpl in cluster:
        x1, x2 = tpl 
        indx = data_dict['slist'][x1].index(x2)
        data_dict['gsflg'][x1][indx] = gscat_value 


def remove_gscat(all_iscat, data_dict):


    """ removes the gscat and leave only iscat

    all_iscat : set
        a set of tuples. each tuple is in the form of (tm_indx, gate_num)
    data_dict : dict
        holds parameters from a certain beam. 

    Returns
    iscat_dict : dict
        a dict of dicts similiar to data_dict but only holdes the iscat data

    
    """
    t1 = dt.datetime.now()

    # write all_iscat set as list of lists. Each list element stors iscat cells for a given time
    tm_indices = sorted(list(set([x[0] for x in all_iscat])))
    all_iscat_lol = [[x for x in all_iscat if y==x[0]] for y in tm_indices]

    kys_tmp = data_dict.keys()
    iscat_dict = dict()
    # initialize iscat_dict
    for d in kys_tmp:
        iscat_dict[d] = []

    kys_a = []    # stores parameters like "slist", "vel" 
    kys_b = []    # stores parameters like "bmazm", "times"

    for ky in kys_tmp:
        if data_dict[ky] == []:
            iscat_dict.pop(ky)
            continue
        if isinstance(data_dict[ky][tm_indices[0]], list):
            kys_a.append(ky)
        else:
            kys_b.append(ky)

    for i, tm in enumerate(tm_indices):

        for ky in kys_tmp:
            if ky in kys_a:
                iscat_dict[ky].append([])
                indx_tmp = [data_dict['slist'][tpl[0]].index(tpl[1]) for tpl in all_iscat_lol[i]]
                iscat_dict[ky][i] = [data_dict[ky][tm][ix] for ix in indx_tmp]
            if ky in kys_b:
                iscat_dict[ky].append(data_dict[ky][tm])

    t2 = dt.datetime.now()
    print ("remove_gscat takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    return iscat_dict

def select_target_interval(data_dict, stm, etm):
    
    """ selects data for a time interval between stm and etm """
    stm_indx = np.argmin(np.array([abs((x-stm).total_seconds()) for x in data_dict['times']]))
    etm_indx = np.argmin(np.array([abs((x-etm).total_seconds()) for x in data_dict['times']]))

    kys_tmp = data_dict.keys()
    for ky in kys_tmp:
        data_dict[ky] = data_dict[ky][stm_indx:etm_indx+1]

    return data_dict

def search_iscat_event(data_dict, ctr_time, bmnum, params, 
                       low_vel_iscat_event_only=True, no_gscat=False):
    """ do the classification for 3 days data one beam at a time.
    
    low_vel_iscat_event_only : bool
        It set to True, returns low velocity inospheric scatter event only
        
    Returns : dict
        A dict of dicts in the form of {bmnum:dict}.
        returns data of a beam with all its points' gsflg marked as 1 (gscat) except for iscat
    """

    # create nodes, whic is a list of lists, from data_dict.
    # Each node is represented by (time_index, gate_number)

    data_dict = data_dict[bmnum]
    nodes = create_nodes(data_dict)

    # cluster the data using depth_first_search algorithm

    # cluster the data
    clusters = []
    visited_nodes_all = set() 
    start_node = find_start_node(nodes, visited_nodes=None)
    while start_node:
        visited_nodes = search_tree(nodes, start_node, data_dict)    # returns a set
        clusters.append(visited_nodes)
        visited_nodes_all.update(visited_nodes) 
        #visited_nodes_all = set([x for y in clusters for x in y])
        start_node = find_start_node(nodes, visited_nodes=visited_nodes_all)


    # pul all the clusters classified as events 
    all_iscat = set([])
    #clusters = clusters[:5]
    for cluster in clusters:
        
        # find the starting and ending times of a cluster
        cluster = push_stm_etm(cluster, data_dict, vel_threshold=15.)

        # classify the cluster
        if low_vel_iscat_event_only:
            event_logic = is_low_vel_event(cluster, data_dict)
        else:
            event_logic = isevent(cluster, data_dict)
        if event_logic:
            # change the gsflg values to 0(isact)
            change_gsflg(cluster, data_dict, gscat_value=0)
            all_iscat.update(cluster)    

    nodes_flat = set([x for y in nodes for x in y])
    if no_gscat:
        # remove gscat points
        data_dict = remove_gscat(all_iscat, data_dict)
    else:
        # change the gsflg values of non-events to 1(gsact)
        all_gscat = set(nodes_flat) - all_iscat
        change_gsflg(all_gscat, data_dict, gscat_value=1)

    # limit the data to to the day of ctr_time(center date)
    stm_target = ctr_time
    etm_target = ctr_time + dt.timedelta(days=1)
    data_dict = select_target_interval(data_dict, stm_target, etm_target)

    return {bmnum:data_dict}

#def remove_gscat(data_dict):
#    """ removes the gscat and leave only iscat
#
#    data_dict : dict
#        holds parameters from a certain beam. 
#    
#    """
#    kys_tmp = data_dict.keys()
#    kys_a = []    # stores parameters like "slist", "vel" 
#    kys_b = []    # stores parameters like "bmazm", "times"
#    for ky in kys_tmp:
#        if data_dict[ky] == []:
#            data_dict.pop(ky)
#            continue
#        if isinstance(data_dict[ky][0], list):
#            kys_a.append(ky)
#        else:
#            kys_b.append(ky)
#
#    kys_a.remove('gsflg')
#    for ky in kys_a:
#        data_tmp = [[x for j, x in enumerate(l) if data_dict['gsflg'][i][j]==0]\
#                for i, l in enumerate(data_dict[ky]) if l]
#        data_dict[ky] = data_tmp
#    for ky in kys_b:
#        data_tmp = [x for i, x in enumerate(data_dict[ky])\
#                if (not np.all(data_dict['gsflg'][i])) and (data_dict['gsflg'][i] is not None)]
#        data_dict[ky] = data_tmp
#
#    data_tmp = [[x for j, x in enumerate(l) if x == 0]\
#            for i, l in enumerate(data_dict['gsflg']) if l]
#    data_dict['gsflg'] = data_tmp
#
#
#    return

def iscat_event_searcher(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt,
                   params=["velocity"], low_vel_iscat_event_only=False,
                   search_allbeams=True, bmnum=7, no_gscat=False, ffname=None):
    """ A wrapper that does all of file prepareting, file reading, and 
        searching for iscat events.
        
    params : list
        works for params=["velocity"] only
    search_allbeams : bool
        if set to true, iscat event searching will be performed on all the 
        beams, and ignores the bmnum argument. 
    bmnum : int
        bmnum argument only works in search_allbeams is set to False
    no_gscat : removes all the gscat
    ffname : string
        if ffname is not set to None, ffname will be be read

    Returns : dict
        A dict of dicts in the form of {bmnum:dict}.
        if no_gscat==False, returns data all its points'
        gsflg marked as 1 (gscat) except for iscat.
        if no_gscat==True, returns only the iscat (gsflg=0)
        
    """
    stm = ctr_time - dt.timedelta(days=1)
    etm = ctr_time + dt.timedelta(days=2)
    rad = localdict["radar"]
    ftype = localdict["ftype"]

    # prepare the data
    t1 = dt.datetime.now()
    if not ffname:
        ffname = prepare_file(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt) 
    t2 = dt.datetime.now()
    print ("prepare_file takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    t1 = dt.datetime.now()
    # read the file. Returns a dict of dicts with bmnums as key words.
    all_beams = read_file(ffname, rad, stm, etm, params, ftype=ftype)
    t2 = dt.datetime.now()
    print ("read_file takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    # search for iscat events
    t1 = dt.datetime.now()
    real_bmnums = [int(x) for x in all_beams.keys() if len(all_beams[x]['times'])>0]
    events = {}
    if search_allbeams:
        for b in real_bmnums:
            # search for event. Returns a dict of dicts with a single bmnum as keyword.
            # this is because search_iscat_event workds on a single beam at a time
            events.update(search_iscat_event(all_beams, ctr_time, b, params,
                low_vel_iscat_event_only=low_vel_iscat_event_only, no_gscat=no_gscat))
    else:
        events.update(search_iscat_event(all_beams, ctr_time, bmnum, params,
            low_vel_iscat_event_only=low_vel_iscat_event_only, no_gscat=no_gscat))
    t2 = dt.datetime.now()
    print ("iscat event searching process takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    return events


def rtiplot(rad, stm, etm, bmnum, params, data_dict=None, fileType="fitacf"):

    from myrti import plot_rti

    data_dict = data_dict[bmnum]
    scales = [[-120, 120]]
    yrng = [0, 70]
    filtered=False
    #fig = plot_rti(stm, "bks", eTime=etm, bmnum=7, gsct=True,
    #        params=["velocity"], scales=[[-120, 120]], colors="aj", yrng=[0, 70])
    fig = plot_rti(stm, rad, eTime=etm, bmnum=bmnum, data_dict=data_dict, gsct=True,
            params=params, scales=scales, colors="aj", yrng=yrng, fileType=fileType,
            filtered=filtered)


    plt.show()
    return fig 

# run the code
def test_code(plotRti=False):

    # input parameters
    ctr_time = dt.datetime(2010,1,15)
    #ctr_time = dt.datetime(2008,9,17)
    rad = "bks"
    channel = None
    bmnum = 7
    params=['velocity']
    ftype = "fitacf"
    #ftype = "fitex"
    filtered = True
    scr = "local"
    localdirfmt = "/sd-data/{year}/{ftype}/{radar}/"
    localdict = {"ftype" : ftype, "radar" : rad, "channel" : channel}
    #tmpdir = "/tmp/sd/"
    tmpdir = "../data/"
    fnamefmt = ['{date}.{hour}......{radar}.{channel}.{ftype}', '{date}.{hour}......{radar}.{ftype}']
    #davitpy.rcParams['verbosity'] = "debug"

    stm = ctr_time - dt.timedelta(days=0)
    etm = ctr_time + dt.timedelta(days=1)
    #etm = ctr_time + dt.timedelta(hours=12)


    # prepare the data
    #ffname = prepare_file(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)
    ffname = tmpdir + "20100114.000000.20100117.000000.bks.fitacff"
    #ffname = tmpdir + "20080916.000000.20080918.000000.bks.fitacff"
    #ffname = tmpdir + "20080916.000000.20080918.000000.bks.fitexf"

#    # make an rti plot
#    if plotRti:
#        # read the file
#        t1 = dt.datetime.now()
#        data_dict = read_file_for_rtiplot(ffname, rad, stm, etm, bmnum, params, ftype=ftype)
#        t2 = dt.datetime.now()
#        print ("read_file takes " + str((t2-t1).total_seconds() / 60.)) + " mins"
#    else:
#        # read the file
#        t1 = dt.datetime.now()
#        data_dict = read_file(ffname, rad, stm, etm, params, ftype=ftype)
#        t2 = dt.datetime.now()
#        print ("read_file takes " + str((t2-t1).total_seconds() / 60.)) + " mins"
#
#    t1 = dt.datetime.now()
#    #data_dict, clusters = search_iscat_event(data_dict, ctr_time, bmnum, params)
#    data_dict = search_iscat_event(data_dict, ctr_time, bmnum, params,
#                                   low_vel_iscat_event_only=False)
#    t2 = dt.datetime.now()
#    print ("search_iscat_event takes " + str((t2-t1).total_seconds() / 60.)) + " mins"

    if plotRti:

        events = iscat_event_searcher(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt,
                       params=params, low_vel_iscat_event_only=False,
                       search_allbeams=False, bmnum=bmnum, no_gscat=False, ffname=ffname)
        data_dict = events

        fig = rtiplot(rad, stm, etm, bmnum, params, data_dict=data_dict, fileType=ftype)

    return data_dict

if __name__ == "__main__":
    #data_dict = test_code(plotRti=False)
    data_dict = test_code(plotRti=True)

