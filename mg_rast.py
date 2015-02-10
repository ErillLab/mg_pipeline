# -*- coding: utf-8 -*-
"""
This script handles API calls and download queries to MG-RAST.

See http://api.metagenomics.anl.gov/api.html for info.

Created on Wed Oct 29 14:31:16 2014

@author: Talmo
"""

import requests
import json
import urllib2
import os
import json

""" Configuration """
# Base API URL
API_URL = 'http://api.metagenomics.anl.gov'


def validate_mg_id(mg_id):
    """
    This function ensures that the metagenome ID follows a standard MG-RAST
    format, including the mgm prefix.
    
    Args:
        mg_id: MG-RAST metagenome ID to be tested
        
    Returns:
        mg_id: Validated ID of the format 'mgm#######.#'
    """
    if type(mg_id) != str:
        mg_id = str(mg_id)
    if mg_id[0:3] != 'mgm':
        mg_id = 'mgm' + mg_id
    return mg_id


def get_download_links(mg_id, query=True):
    """
    Returns the links to download files associated with a metagenome.
    
    Args:
        mg_id: The MG-RAST ID of the metagenome
        query: If True, queries the API for the links, otherwise returns the
            links without querying. If False, the returned links may not
            actually exist (default = True).
            
    Returns:
        links: dictionary of 'file_id': 'url' pairs.
        file_names: dictionary of 'file_id': 'file_name' pairs.
            
    Example:
        >>> links, file_names = get_download_links('mgm4508947.3')
        >>> links
        {'050.1': 'http://api.metagenomics.anl.gov/download/mgm4508947.3?file=050.1',
         ...
         '700.1': 'http://api.metagenomics.anl.gov/download/mgm4508947.3?file=700.1'}
        >>> file_names
        {'050.1': 'mgm4508947.3.050.upload.fna',
        ...
        '700.1': 'mgm4508947.3.700.annotation.sims.filter.seq'}


    """
    
    # Make sure we have the correct format of the ID ('mgm...')
    mg_id = validate_mg_id(mg_id)
    
    if query:
        # Query the API for the links
        query_url = API_URL + '/download/' + mg_id;
        r = requests.get(query_url)
        data = json.loads(r.content)['data']
        links = {str(link['file_id']): str(link['url']) for link in data}
        file_names = {str(link['file_id']): str(link['file_name']) for link in data}

    else:
        # Assume the links follow a format
        file_ids = ['050.1', '100.1', '100.2', '150.1', '150.2', '299.1', '350.1', '425.1', '440.1', '440.2', '450.1', '550.1', '550.2', '650.1', '700.1']
        file_names = file_ids
        links = {file_id: API_URL + '/download/' + mg_id + '?file=' + file_id for file_id in file_ids}
        
    return links, file_names


def download_url(url, path):
    req = urllib2.Request(url)
    with open(path, "w") as local_file:
        local_file.write(urllib2.urlopen(req).read())

def download(mg_id, stage, folder_path='.'):
    """ Downloads the data file associated with the given MG-RAST stage. """
    # Make sure we have the correct format of the ID ('mgm...')
    mg_id = validate_mg_id(mg_id)
    
    # Append MG ID to the path
    folder_path = folder_path + '/' + mg_id
    
    # Check if folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    # Get links for metagenome
    links, file_names = get_download_links(mg_id, query=True)
    url = links[stage]
    #file_name = file_names[stage].replace(mg_id + '.', '')
    file_name = file_names[stage]
    
    # Download
    download_url(url, folder_path + '/' + file_name)

def get_metadata(mg_id, parse=True):
    """ Gets all metadata associated with the specified metagenome. """
    
    # Make sure we have the correct format of the ID ('mgm...')
    mg_id = validate_mg_id(mg_id)
    
    # Query API
    r = requests.get(API_URL + "/metagenome/" + mg_id + "?verbosity=full")
    metadata = r.content
    
    # Parse JSON response
    if parse:
        metadata = json.loads(metadata)
    
    return metadata
