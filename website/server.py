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

def extract_network(network_name, go_ids):
    """
    Extract the subgraph from a list of go_ids and for a given network name (biological_process, molecular_function or cellular_component)
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
    return json_graph.node_link_data(sub_network)

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
            protein_ids = message['protein_ids']
            answer = {}
            proteins = {} # a dict of protein that will be used to construct a javascript table
            network = None
            #first we extract the go terms for the protein ids
            # biological_process_go_terms = []
            # molecular_function_go_terms = []
            # cellular_component_go_terms = []
            # pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/%s_terms.pickle'%message['network type']
            # terms = None
            # with open(pickle_file) as h:
            #     terms = pickle.load(h)
            for protein_id in protein_ids:
                 details = proteins_details[protein_id]
            #     go_terms = []
            #     #we recover the name for each go term
            #     for go_term_id in protein_details[message['network type']]:
            #         go_terms.append(terms[go_term_id])
            #     #proteins[protein_id] = [protein_id, protein_details['name'] go_terms]
                 proteins[protein_id] = [protein_id, details['name']]
            #     #print go_terms
            #     biological_process_go_terms += protein_details['biological_process']
            #     molecular_function_go_terms += protein_details['molecular_function']
            #     cellular_component_go_terms += protein_details['cellular_component']
            #if message['network type'] == 'cellular_component':
            #    network = extract_network('cellular_component', cellular_component_go_terms)
            #elif message['network type'] == 'molecular_function':
            #    network = extract_network('molecular_function', molecular_function_go_terms)
            #elif message['network type'] == 'biological_process':
            #    network = extract_network('biological_process', biological_process_go_terms)
            all_paths = []
            proteins_without_paths = []
            #first we extract all the shortest path between each pairs of protein. We want to find a shortest path (meaning a link of interaction) for each pair of proteins
            for i in range(0, len(protein_ids)-1):
                protein_without_path = True
                for j in range(i+1, len(protein_ids)):
                    try:
                        path = nx.shortest_path(interaction_network, protein_ids[i], protein_ids[j])
                        all_paths.append(path)
                        protein_without_path = False
                    except nx.exception.NetworkXNoPath as e:
                        #print e
                        pass
                    except AttributeError as e:
                        #print e
                        pass
                if protein_without_path:
                    proteins_without_paths.append(protein_ids[i])
            core = []
            #now we need to remove the redundancy...
            merged_path = ["not empty"]
            uniq_paths = []
            while len(merged_path) != 0: #we stop the process when we have no more new merged path
                merged_path = []
                #for each path linking two proteins, we search for another path with the maximum in common
                for i in range(0, len(all_paths)-1):
                    max_common_length = 0
                    best_common_path = None
                    for j in range(i+1, len(all_paths)):
                        l = len(set(all_paths[i]) & set(all_paths[j]))
                        if l > max_common_length:
                            best_common_path = all_paths[j]
                            max_common_length = l
                    #if we found the best path with the maximum in common, we extract their similarities into a new path
                    if best_common_path != None: #a common path with at least 1 element has been found
                        merged_path.append(set(all_paths[i]) & set(best_common_path))
                        #print 'found best with length %i'%max_common_length
                    #otherwise, for the current path tested, it is now uniq
                    else:
                        #print 'no best found'
                        if not all_paths[i] in uniq_paths:
                            uniq_paths.append(all_paths[i])
                #we relaunch a step of merging with the new merged path
                all_paths = merged_path
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
            for protein_id_selected in protein_ids:
                try:
                    sub_network.node[protein_id_selected]['selected'] = True
                except KeyError as e:
                    pass
            answer = {
                'header': 'networks constructed',
                'interactions-network': json_graph.node_link_data(sub_network),
                'proteins': proteins.values()
            }
            self.write_message(answer, binary = False)
        elif message['header'] == 'extend network':
            #we recover the current interactions network displayed
            subnetwork = json_graph.node_link_graph(message['network'])
            #we keep the ids of the nodes that made the current network displayed (the core)
            #we get its nodes
            core = subnetwork.nodes()
            core_ids = []
            for node_id in core:
                core_ids.append(node_id)
            #we get the nodes for all the shortest paths between the protein_id selected by the user and the nodes of the current network displayed
            all_paths = []
            for c in core:
                try:
                    path = nx.shortest_path(interaction_network, c, message['protein_id'])
                    all_paths.append(path)
                except nx.exception.NetworkXNoPath as e:
                    #print e
                    all_paths.append([message['protein_id']])
                except AttributeError as e:
                    #print e
                    pass
            for p in all_paths:
                core += p
            core = set(core)
            sub_network = interaction_network.subgraph(core)
            for core_id in sub_network.nodes():
                sub_network.node[core_id]['core'] = True
                sub_network.node[core_id]['selection'] = False
                sub_network.node[core_id]['selected'] = False
            for protein_id_selected in message['protein_ids_selected']:
                try:
                    sub_network.node[protein_id_selected]['selected'] = True
                except KeyError as e:
                    pass
            try:
                sub_network.node[message['protein_id']]['selection'] = True
            except:
                pass
            answer = {
                'header': 'network extended',
                'network': json_graph.node_link_data(sub_network)
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
