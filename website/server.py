#!/usr/bin/env python

import os, ujson, pickle, urllib, re
from collections import Counter
import networkx as nx
from networkx.readwrite import json_graph
import tornado.httpserver
import tornado.ioloop
import tornado.options
import tornado.web
import tornado.websocket
from tornado.escape import json_encode, native_str, xhtml_escape

interaction_network = None
proteins_details = None

pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/interaction_network.pickle'
with open(pickle_file) as h:
    interaction_network = pickle.load(h)
pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/protein_details.pickle'
with open(pickle_file) as h:
    proteins_details = pickle.load(h)

def extract_ontology_network(network_name, go_ids):
    """
    Extract a subnetwork from a list of go_ids and for a given network name (biological_process, molecular_function or cellular_component)
    """
    pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/%s_network.pickle'%network_name
    with open(pickle_file) as h:
        network = pickle.load(h)
    nodes = [] #the nodes making all the shortest paths
    for go_id in go_ids:
        nodes += nx.shortest_path(network, go_id).keys() #will search for the shortest path between the go_id and all other nodes
    #we count how many times we found each go id
    c = Counter(nodes)
    nodes = list(set(nodes)) #remove redundancy
    sub_network = network.subgraph(nodes) #extract the sub_network made with the nodes of the shortest paths
    for k in c.keys():
        sub_network.node[k]['counts'] = c[k] #add the counts attribute to each node (i.e. how many times we found each go id)
    return sub_network

def extract_interaction_network(proteins_selected, all_proteins, keep_full_paths = False):
    """
    Extract a subnetwork from interaction network and a list of protein ids
    """
    all_paths = []
    all_ids = []
    #first we extract all the shortest path between each pairs of protein. We want to find a shortest path (meaning a link of interaction) for each pair of proteins
    for i in range(0, len(proteins_selected)-1):
        protein_without_path = True
        for j in range(i+1, len(proteins_selected)):
            if proteins_selected[i] == proteins_selected[j]: #it is not excluded that a protein is recorded several times in the CSV file
                continue
            try:
                path = nx.shortest_path(interaction_network, proteins_selected[i], proteins_selected[j])
                all_paths.append(path)
                all_ids += path
                protein_without_path = False
            except nx.exception.NetworkXNoPath as e:
                pass
            except AttributeError as e:
                pass
        if protein_without_path and not proteins_selected[i] in all_ids: #the second condition is due to the fact that the current protein could have been in a path with a protein previously tested
            all_paths.append([proteins_selected[i]])
    #the last one (if any) is in its own shortest path (if it was not already in a path, that's why we have the second condition)
    if i < len(proteins_selected)-1 and not proteins_selected[i] in all_ids:
        all_paths.append([proteins_selected[i+1]])
    core = []
    counter = Counter(all_ids)
    max_counter = 0
    if len(counter):
        max_counter = max(counter.values())
    #now we need to remove the redundancy...
    merged_paths = ["not empty"]
    uniq_paths = []
    if not keep_full_paths:
        while len(merged_paths) != 0: #we stop the process when we have no more new merged path
            merged_paths = []
            #for each path linking two proteins, we search for another path with the maximum in common
            for i in range(0, len(all_paths)-1):
                merged_path = None
                for j in range(i+1, len(all_paths)):
                    p = set(all_paths[i]) & set(all_paths[j])
                    if merged_path == None and len(p):
                        merged_path = p
                    elif merged_path != None and len(p) > len(merged_path):
                        merged_path = p
                #if we found the best path with the maximum in common, we extract their similarities into a new path
                if merged_path != None: #a common path with at least 1 element has been found
                    merged_paths.append(merged_path)
                    #print 'found best with length %i'%max_common_length
                #otherwise, for the current path tested, it is now uniq
                else:
                    #print 'no best found'
                    if not all_paths[i] in uniq_paths:
                        uniq_paths.append(all_paths[i])
            #the last one (if any)
            if i < len(all_paths)-1 and not all_paths[i+1] in uniq_paths:
                uniq_paths.append(all_paths[i+1])
            #we relaunch a step of merging with the new merged path
            all_paths = merged_paths
    else:
        uniq_paths = all_paths
    #print uniq_paths
    #we gather all the nodes of the core...
    for p in uniq_paths:
        core += p
    core = set(core)
    sub_network = interaction_network.subgraph(core)
    for node_id in sub_network.nodes():
        sub_network.node[node_id]['core'] = True
        sub_network.node[node_id]['selection'] = False
        sub_network.node[node_id]['selected'] = False
        sub_network.node[node_id]['sample'] = False
        sub_network.node[node_id]['connexity'] = counter[node_id]
        if max_counter:
            sub_network.node[node_id]['normalized connexity'] = float(counter[node_id])/float(max_counter)
        else:
            sub_network.node[node_id]['normalized connexity'] = 0.0
    for protein_selected in proteins_selected:
        try:
            sub_network.node[protein_selected]['selected'] = True
        except KeyError as e:
            pass
    for protein in all_proteins:
        try:
            sub_network.node[protein]['sample'] = True
        except KeyError as e:
            pass
    return sub_network

def extend_interaction_network(network, protein_id, protein_ids_selected, all_proteins):
    """
    Using a protein id, extend a subnetwork extracted from the interaction network
    """
    #we get its nodes
    core = network.nodes()
    core_ids = []
    for node_id in core:
        core_ids.append(node_id)
    #we get the nodes for all the shortest paths between the protein_id selected by the user and the nodes of the current network displayed
    all_paths = []
    for c in core:
        try:
            path = nx.shortest_path(interaction_network, c, protein_id)
            all_paths.append(path)
        except nx.exception.NetworkXNoPath as e:
            #print e
            all_paths.append([protein_id])
        except AttributeError as e:
            #print e
            pass
    for p in all_paths:
        core += p
    core = set(core)
    sub_network = interaction_network.subgraph(core)
    for core_id in sub_network.nodes():
        sub_network.node[core_id]['core'] = False
        sub_network.node[core_id]['selection'] = False
        sub_network.node[core_id]['selected'] = False
    for core_id in core_ids:
        sub_network.node[core_id]['core'] = True
    for protein_id_selected in protein_ids_selected:
        try:
            sub_network.node[protein_id_selected]['selected'] = True
        except KeyError as e:
            pass
    for protein in all_proteins:
        try:
            sub_network.node[protein]['sample'] = True
        except KeyError as e:
            pass
    try:
        sub_network.node[protein_id]['selection'] = True
    except:
        pass
    return sub_network

app = None
static_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')
pages_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates')

class Index(tornado.web.RequestHandler):
    def get(self):
        self.render('index.html')

class WebSocket(tornado.websocket.WebSocketHandler):

    def open(self, *args):
        print "New client connected"

    def on_message(self, message):
        message = ujson.loads(message)
        if message['header'] == 'construct networks':
            proteins_selected = message['protein_selected']
            all_proteins = message['all_proteins']
            table_columns = message['table_columns'] #the data to fill the columns besite 'ID' and 'name'
            keep_full_paths = message['keep_full_paths']
            answer = {}
            proteins = {} # a dict of protein that will be used to construct the DataTable in the index.html file
            network = None
            go_ids = []
            sub_interaction_network = extract_interaction_network(proteins_selected, all_proteins, keep_full_paths = keep_full_paths)
            for protein_selected in proteins_selected:
                #first we construct an array with the details of the proteins selected. These data will be listed in the table
                protein_details = proteins_details[protein_selected]
                proteins[protein_selected] = []
                for column in table_columns:
                    if column == 'connexity':
                        #we extract the connexity from the sub_interaction_network just computed (if the subnetwork has still the protein as a node)
                        if protein_selected in sub_interaction_network.nodes():
                            proteins[protein_selected].append(sub_interaction_network.node[protein_selected]['connexity'])
                        else:
                            proteins[protein_selected].append(0)
                    else:
                        proteins[protein_selected].append(protein_details[column])

                #gather go_ids
                #if message['network type'] == 'cellular_component':
                #    go_ids += protein_details['cellular_component'].keys()
                #elif message['network type'] == 'molecular_function':
                #    go_ids += protein_details['molecular_function'].keys()
                #elif message['network type'] == 'biological_process':
                #    go_ids += protein_details['biological_process'].keys()
            answer = {
                'header': 'networks constructed',
                'interactions-network': json_graph.node_link_data(sub_interaction_network),
                #'ontology-network': json_graph.node_link_data(ontology_network),
                'proteins': proteins.values()
            }
            self.write_message(answer, binary = False)
        elif message['header'] == 'extend network':
            answer = {
                'header': 'network extended',
                'network': json_graph.node_link_data(extend_interaction_network(json_graph.node_link_graph(message['network']), message['protein_id'], message['protein_ids_selected'], message['all_proteins']))
            }
            self.write_message(answer, binary = False)

    def on_close(self):
        print "Client disconnected"

class Application(tornado.web.Application):
    def __init__(self):

        handlers = [
            (r'/', Index),
            (r'/websocket', WebSocket)
        ]

        settings = {
            'template_path': pages_dir,
            'static_path': static_dir,
            'debug': True
        }

        tornado.web.Application.__init__(self, handlers, **settings)

if __name__ == '__main__':
    webserver_port = 8080
    app = Application()
    server = tornado.httpserver.HTTPServer(app)
    server.listen(webserver_port)
    main_loop = tornado.ioloop.IOLoop.instance()
    main_loop.start()
