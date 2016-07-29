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
params=['velocity']
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

def read_data(myPtr, bmnum, params=["velocity"], tbands=None):
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
                 'frang', 'gsflg', 'velocity_error']
    for d in data_keys:
        data[d] = []


    
    # Read the parameters of interest.
    myPtr.rewind()
    myBeam = myPtr.readRec()


#    xx = 0
#    stack = []
#    visited = set()
#    if (myBeam is not None):
#        stack.append((xx, myBeam.fit.slist[0]))

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

def create_graph(vertex, data):

    if vertex[0] == 0:
        xx = [0, 1]  # time indices centered at the time index of the vertex 
    else:
        xx = [vertex[0] + yi for yi in [-1, 0, 1]]  # time indices centered at the time index of the vertex 

    tm_indx = vertex[0]
    if tm_indx == 0:
        xx = [0, 1]  # time indices centered at the time index of the vertex 
    elif tm_indx == len(data)-1:
        xx = [tm_indx, tm_indx-1]  # time indices centered at the time index of the vertex 
    else:   
        xx = [vertex[0] + yi for yi in [-1, 0, 1]]  # time indices centered at the time index of the vertex 
    yy = [vertex[1] + yi for yi in [-1, 0, 1]]  #range-gate indices centered at the range-gate index of the vertex 

    # create a tuple that inclues actual neighbering vortice
    xy = [(i,j) for i in xx for j in set([x[1] for x in data[i]]).intersection(set(yy))] 
    xy.remove(vertex)    # remove the vertex from xy 
    nodes = set()
    for tpl in xy:
        nodes.add(tpl)
    G = {vertex:nodes}
    return G


def search_tree(data, start, visited = None):

    # create a grash as a start
    G = create_graph(start, data)
    # do the breath_firt_search
    if visited is None:
        visited = set()
    queue = [start]
    while queue:
        vertex = queue.pop(0)
        if vertex[1] < 7: 
            next_node = queue[0]
            G = create_graph(next_node, data)
            continue
        if vertex not in visited:
            visited.add(vertex)
            #try:
            #    pdb.set_trace()
            queue.extend(G[vertex] - visited)
            #except:
                #return visited
        try:
            next_node = queue[0]
            G = create_graph(next_node, data)
        except:
            pass
    return visited

def connect_clusters(clusters):
    pass

def change_gsflg(cluster, data_dict, gscat_value=0):
    for tpl in cluster:
        x1, x2 = tpl 
        indx = data_dict['slist'][x1].index(x2)
        data_dict['gsflg'][x1][indx] = gscat_value 

def isevent(cluster, data_dict, vel_threshold=15.):
    import datetime as dt
    # convert cluster to list of lists
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

#        cluster_vels = [] 
#        for item in cluster:
#            try:
#                vel_tmp = data_dict['vel'][item[0]][(data_dict['vel'][item[0]]).index(item[1])]
#                cluster_vels.append(vel_tmp)
#            except:
#                print item
        
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

def classify_event(cluster, data_dict):

    # change the gsflg values
    event_logic = isevent(cluster, data_dict)
    if event_logic:
        change_gsflg(cluster, data_dict, gscat_value=0)
    else:
        change_gsflg(cluster, data_dict, gscat_value=1)
         
    
    # find the starting and ending time of a cluster
    #while True:

    #cluster = [[x for x in cluster if x[0]==y] for y in tm_indices]
#    return cluster
    

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


def dopsearch(ctr_time, bmnum, params, localdirfmt, localdict, tmpdir, fnamefmt):

    # fetch and concatenate the three consecutive days of data centered on the target date 
    #concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)

    # box car fiter the data
    #ffname = boxcar_filter(concated_file)

    ffname = "./data/20100114.000000.20100116.000000.bks.fitacff"

    #import sys
    #sys.path.append("/home/muhammad/softwares/davitpy_MuhammadVT/davitpy/pydarn/plotting/")
    #from rti import plot_rti
    from myrti import plot_rti

    #stm = ctr_time - dt.timedelta(days=1)
    stm = ctr_time
    etm = ctr_time + dt.timedelta(days=1)

    #plot_rti(stm, "bks", eTime=etm, bmnum=7)
    #plot_rti(stm, "bks", eTime=etm, bmnum=7, fileName=ffname)
    #plot_rti(dt.datetime(2013,3,16), 'bks', eTime=dt.datetime(2013,3,16,14,30), bmnum=12)

    ftype = "fitacf"
    scales = [[-120, 120]]
    #ftype = "fitex"
    myPtr = radDataOpen(stm, "bks", eTime=etm, bmnum=bmnum, fileName=ffname, fileType=ftype)
    #myPtr = radDataOpen(stm, "bks", eTime=etm, bmnum=7, fileType="fitacf")
    #plot_rti(stm, 'bks', eTime=etm, bmnum=7, fileType=ftype, myFile=myPtr,
    #         fileName=ffname, params=["velocity"], coords='rng', gsct=True, scales=scales)

    #plot_rti(stm, 'bks', eTime=etm, bmnum=7, fileType=ftype, params=["velocity"],
    #          coords='rng', gsct=True, scales=scales, filtered=True)

    data_dict = read_data(myPtr, bmnum, params=params, tbands=None)

    # create nodes from data_dict. Each node is represented by (time_index, gate_number)
    nodes = create_nodes(data_dict)

    # cluster the data using depth_first_search algorithm

    # cluster the data using breath_first_search algorithm




    cluster_list = []
    start_node = find_start_node(nodes, visited_nodes=None)
    while start_node:
        visited_nodes = search_tree(nodes, start_node)
        cluster_list.append(visited_nodes)
        start_node = find_start_node(nodes,
                                     visited_nodes=set([x for y in cluster_list for x in y]))

    #visited_nodes=set([x for y in cluster_list for x in y])


    for cluster in cluster_list:
        classify_event(cluster, data_dict)

    #start_node = nodes[0][0]
    #start_node = (1349, 7)
    #visited_nodes = search_tree(nodes, start_node)

#    nodes_flat = list(set([x for y in nodes for x in y]) - visited_nodes)
#    tm_indices = list(set([x[0] for x in nodes_flat]))
#    nodes_tmp = [[x for x in nodes_flat if y==x[0]] for y in tm_indices]


#    for sublist in nodes:
#        for itm in sublist:
#            if itm[1] >=7 and (itm not in visited_nodes):
#                start_node = itm
#            break
#    print start_node in visited_nodes
#    visited_nodes_tmp = search_tree(nodes, start_node, visited = visited_nodes)
#
#    #visited_nodes = visited_nodes.union(search_tree(nodes, nodes[0][0]))
#    #print visited_nodes - visited_nodes_tmp
#    visited_nodes = visited_nodes.union(visited_nodes_tmp)

#    # change the gsflg 
#    for tpl in visited_nodes:
#        x1, x2 = tpl 
#        indx = data_dict['slist'][x1].index(x2)
#        data_dict['gsflg'][x1][indx] = 0 

#    for x1 in range(len(data_dict["times"])):
#        try:
#            for x2 in range(len(data_dict["gsflg"][x1])):
#                data_dict['gsflg'][x1][x2] = 0 
#        except:
#            continue

    stm = ctr_time
    etm = ctr_time + dt.timedelta(days=1)
    fig = plot_rti(stm, "bks", eTime=etm, bmnum=7, data_dict=data_dict, gsct=True,
            params=["velocity"], scales=[[-120, 120]], colors="aj")


    plt.show()
    return data_dict, visited_nodes, nodes 



import pdb
data_dict, visited_nodes, nodes = dopsearch(ctr_time, bmnum, params, localdirfmt, localdict, tmpdir, fnamefmt)


#    #myPtr = radDataOpen(stime, rad, etime, channel=channel, bmnum=bmnum, fileType=ftype,
#    #                    filtered=filtered, local_dirfmt=localdirfmt)
#
#    myPtr = radDataPtr(stime, rad, etime, channel=channel, bmnum=bmnum, fileType=ftype,
#                        filtered=filtered, local_dirfmt=localdirfmt)
#

#concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)
#ffname = boxcar_filter(concated_file)

#from davitpy.pydarn.plotting import plot_rti



