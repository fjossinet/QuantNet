#!/usr/bin/env python

"""
This script allows to:
- reconstruct the network of ontology terms using the OBO file go-basic.obo
"""

import urllib2, os, pickle, os.path, json, sys, re
from collections import Counter
import networkx as nx
from networkx.readwrite import json_graph
import pandas as pd
from pandas import read_csv

def find_protein_partners(network, protein_id):
    """
    Search for protein partners in a Uniprot entry. Store this partner and the interaction in a network.
    """
    if not network.has_node(protein_id):
        network.add_node(protein_id, id = protein_id)
    print "Fetch http://www.uniprot.org/uniprot/%s.txt"%protein_id
    try:
        response = urllib2.urlopen("http://www.uniprot.org/uniprot/%s.txt"%protein_id)
        content = str(response.read())
        load_interactions = False
        partner_ids = []
        for l in content.split('\n'):
            if l.startswith('CC   -!- INTERACTION:'):
                load_interactions = True
                continue
            elif l.startswith('CC   -!-'):
                load_interactions = False
            if load_interactions and l.startswith('CC       '):
                line = re.split('^CC', l.split(';')[0])[-1].strip()
                if line == 'Self':
                    continue
                partner_id = line.split(':')[0]
                if not network.has_node(partner_id):
                    network.add_node(partner_id, id = protein_id)
                    partner_ids.append(partner_id)
                if not network.has_edge(protein_id, partner_id) and not network.has_edge(partner_id, protein_id): # the graph in directed, so we need to check both options
                    network.add_edge(protein_id, partner_id)
        print "new protein to process: %i"%len(partner_ids)
        for partner_id in partner_ids:
            find_protein_partners(network, partner_id)
    except:
        print "Log: Pb with address http://www.uniprot.org/uniprot/%s.txt"%protein_id

def get_details(protein_id):
    """
    Extract informations from a Uniprot entry
    """
    response = urllib2.urlopen("http://www.uniprot.org/uniprot/%s.txt"%protein_id)
    content = str(response.read())
    biological_process_go_terms = []
    molecular_function_go_terms = []
    cellular_component_go_terms = []
    name = None
    for l in content.split('\n'):
        if l.startswith('DR'):
            l = l.split('DR')[-1].strip()
            if l.startswith('GO;'):
                a = map(str.strip, l.split(';'))
                if a[2].startswith('C:'):
                    cellular_component_go_terms.append(a[1])
                elif a[2].startswith('F:'):
                    molecular_function_go_terms.append(a[1])
                elif a[2].startswith('P:'):
                    biological_process_go_terms.append(a[1])
        elif l.startswith('DE'):
            name = l.split('=')[-1].split('{')[0].replace(";", "")
    return {
        'id': protein_id,
        'name': name,
        'biological_process': biological_process_go_terms,
        'cellular_component': cellular_component_go_terms,
        'molecular_function': molecular_function_go_terms
    }

def get_protein_details(csv_file):
    protein_ids = []
    protein_details = {}
    df = read_csv(csv_file)
    for index, row in df.iterrows():
        protein_ids.append(row['ID'])
    for protein_id in protein_ids:
         protein_details[protein_id] = get_details(protein_id)
    pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/protein_details.pickle'
    h = open(pickle_file, 'w')
    pickle.dump(protein_details, h)
    h.close()

def reconstruct_interaction_network(csv_file):
    """
    Using a list of uniprot ids, this function reconstructs their interaction network based on their Uniprot entry.
    """
    protein_ids = []
    df = read_csv(csv_file)
    for index, row in df.iterrows():
        protein_ids.append(row['ID'])
    network = nx.DiGraph() #a DiGraph will compute the simple paths (networkx.all_simple_paths()) between two proteins faster than with a undirected Graph
    for protein_id in protein_ids:
        find_protein_partners(network, protein_id)
    print "Log: %i proteins"%len(network.nodes())
    print "Log: %i interactions"%len(network.edges())
    pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/interaction_network.pickle'
    h = open(pickle_file, 'w')
    pickle.dump(network, h)
    h.close()

def reconstruct_ontologies():
    """
    Parse the OBO file provided by the Gene Ontology to produce the ontologies (as NetworkX graphs) for 3 domains: biological_process, molecular_function or cellular_component.
    This ontologies will be dumped as pickle files in website/static/data.
    """
    obo_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/go-basic.obo'
    reconstruct_ontology(obo_file, 'biological_process')
    reconstruct_ontology(obo_file, 'molecular_function')
    reconstruct_ontology(obo_file, 'cellular_component')

def reconstruct_ontology(obo_file, domain = 'biological_process'):
    """
    Reconstruct and dump an ontology (as a NetworkX graph) for a given domain: biological_process, molecular_function or cellular_component
    """
    print "Ontology reconstruction from OBO file for %s"%domain
    network = nx.DiGraph()
    ids_2_names = {}
    with open(obo_file) as h:
        current_go_id = None
        name = None
        for line in h.readlines():
            if line.startswith('id: GO:'):
                current_go_id = line.split('id:')[-1].strip()
            if line.startswith('name:'):
                name = line.split('name:')[-1].strip()
                ids_2_names[current_go_id] = name
            if 'namespace: %s'%domain in line:
                network.add_node(current_go_id, name = name)
    with open(obo_file) as h:
        current_go_id = None
        new_node = False
        for line in h.readlines():
            if line.startswith('id: GO:'):
                current_go_id = line.split('id:')[-1].strip()
            if line.startswith('namespace:'):
                if domain in line:
                    new_node = True
                else:
                    new_node = False
            if line.startswith('is_a:') and new_node:
                is_a_go_id = line.split('is_a:')[-1].split('!')[0].strip()
                network.add_edge(current_go_id, is_a_go_id, type='is_a')
            if line.startswith('relationship: part_of') and new_node:
                part_of_go_id = line.split('relationship: part_of')[-1].split('!')[0].strip()
                network.add_edge(current_go_id, is_a_go_id, type='part_of')
    pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/%s_network.pickle'%domain
    h = open(pickle_file, 'w')
    pickle.dump(network, h)
    h.close()

    pickle_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/%s_terms.pickle'%domain
    h = open(pickle_file, 'w')
    pickle.dump(ids_2_names, h)
    h.close()

if __name__ == '__main__':
    csv_file = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/website/static/data/sample.csv'
    get_protein_details(csv_file)
    reconstruct_ontologies()
    reconstruct_interaction_network(csv_file)
