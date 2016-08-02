import datetime as dt
from davitpy.pydarn.sdio.fetchUtils import fetch_local_files
from davitpy.pydarn.sdio import radDataOpen, radDataReadRec
from davitpy.pydarn.sdio import radDataPtr
import davitpy
import logging
import os
import string
import matplotlib.pyplot as plt

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

def prepare_file(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt):

    # fetch and concatenate the three consecutive days of data centered on the target date 
    concated_file = fetch_concat(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)

    # box car fiter the data
    ffname = boxcar_filter(concated_file)

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

def read_file(ffname, rad, ctr_time, bmnum, params, ftype="fitacf"):

    #stm = ctr_time - dt.timedelta(days=1)
    stm = ctr_time
    etm = ctr_time + dt.timedelta(days=1)

    #myPtr = radDataOpen(stm, "bks", eTime=etm, bmnum=bmnum, cp=153,
    myPtr = radDataOpen(stm, rad, eTime=etm, bmnum=bmnum,
                        fileName=ffname, fileType=ftype)
    data_dict = read_data(myPtr, bmnum, params=params, tbands=None)

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
                    del_time = round((tm_next - end_tm).seconds / 60.)

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
         
    
    # find the starting and ending time of a cluster
    #while True:

    #cluster = [[x for x in cluster if x[0]==y] for y in tm_indices]
#    return cluster
    

def change_gsflg(cluster, data_dict, gscat_value=0):
    for tpl in cluster:
        x1, x2 = tpl 
        indx = data_dict['slist'][x1].index(x2)
        data_dict['gsflg'][x1][indx] = gscat_value 



def dopsearch(data_dict, ctr_time, bmnum, params):


    # create nodes, whic is a list of lists, from data_dict.
    # Each node is represented by (time_index, gate_number)
    nodes = create_nodes(data_dict)

    # cluster the data using depth_first_search algorithm

    # cluster the data
    clusters = []
    visited_nodes_all = set() 
    start_node = find_start_node(nodes, visited_nodes=None)
    while start_node:
        #print start_node
        #if start_node == (830, 37):
        #    pdb.set_trace()
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
        event_logic = isevent(cluster, data_dict)
        if event_logic:
            # change the gsflg values to 0(isact)
            change_gsflg(cluster, data_dict, gscat_value=0)
            all_iscat.update(cluster)

    nodes_flat = set([x for y in nodes for x in y])
    all_gscat = set(nodes_flat) - all_iscat
    # change the gsflg values to 1(gsact)
    change_gsflg(all_gscat, data_dict, gscat_value=1)

    #return data_dict, clusters 
    #print "cluster_num = ", len(clusters)
    return data_dict 


def rtiplot(rad, stm, etm, bmnum, params, data_dict=None, fileType="fitacf"):

    from myrti import plot_rti

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
import pdb

# input parameters
#ctr_time = dt.datetime(2010,1,15)
ctr_time = dt.datetime(2008,9,17)
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
tmpdir = "/tmp/sd/"
fnamefmt = ['{date}.{hour}......{radar}.{channel}.{ftype}', '{date}.{hour}......{radar}.{ftype}']
#davitpy.rcParams['verbosity'] = "debug"

# prepare the data
#ffname = prepare_file(ctr_time, localdirfmt, localdict, tmpdir, fnamefmt)
#ffname = "./data/20100114.000000.20100116.000000.bks.fitacff"
ffname = "./data/20080916.000000.20080918.000000.bks.fitacff"
#ffname = "./data/20080916.000000.20080918.000000.bks.fitexf"

# read the file
data_dict = read_file(ffname, rad, ctr_time, bmnum, params, ftype=ftype)

#data_dict, clusters = dopsearch(data_dict, ctr_time, bmnum, params)
data_dict = dopsearch(data_dict, ctr_time, bmnum, params)

# make an rti plot
#stm = dt.datetime(2010,1,15, 12)
#etm = dt.datetime(2010,1,15, 14)
stm = ctr_time
#etm = ctr_time + dt.timedelta(days=12)
etm = ctr_time + dt.timedelta(hours=12)
fig = rtiplot(rad, stm, etm, bmnum, params, data_dict=data_dict, fileType=ftype)


