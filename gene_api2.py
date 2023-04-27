#!/usr/bin/env python3

#Anna Victoria Lavelle, Kamilla Madera, Dieu-Quyen Nguyen
#April 26, 2023

from flask import Flask, request, send_file
from collections import Counter, OrderedDict
import requests
import redis
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import jobs

app = Flask(__name__)

def get_redis0():
    redis_ip = os.environ.get('REDIS_IP')
    if not redis_ip:
        raise Exception()
    red = redis.Redis(host=redis_ip, port=6379, db=0, decode_responses=True)
    return red

def get_redis1():
    redis_ip1 = os.environ.get('REDIS_IP')
    if not redis_ip1:
        raise Exception()
    red1 = redis.Redis(host=redis_ip1, port=6379, db=1)
    return red1

def get_redis2():
    redis_ip2 = os.environ.get('REDIS_IP')
    if not redis_ip2:
        raise Exception()
    red2 = redis.Redis(host=redis_ip2, port=6379, db=2, decode_responses=True)
    return red2

rd = get_redis0()
rd1 = get_redis1()
rd2 = get_redis2()

def get_data():
    response = requests.get(url = 'https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/json/hgnc_complete_set.json')
    response = response.json()['response']['docs']
    return response


@app.route('/data', methods = ['POST', 'GET', 'DELETE'])
def handle_data():
    """
    POST: Posts the data to the database and returns confirmation of this to the 
        user.
    GET: Retrieves the data from the database and returns it to the user.
    DELETE: Deletes the data from the database and returns confirmaiton of this to
        the user.
    Args:
        POST: None.
        GET: None.
        DELETE: None.
    Returns:
        POST (str): Returns "Data loaded" message.
        GET (output_list): Returns the data from the database.
        DELETE (str): Returns "Data deleted, there are 0 keys in the db" message.
    """

    if request.method == 'GET':
        output_list = []
        if len(rd.keys()) < 1:
            return ("No data in the database to retrieve. Please use a POST route first.\n")
        for item in rd.keys():
            output_list.append(json.loads(rd.get(item)))
        return output_list
    elif request.method == 'POST':
        data = get_data()
        for item in data:
            key = f'{item["hgnc_id"]}'
            rd.set(key, json.dumps(item))
        return f'Data loaded.\n'
    elif request.method == 'DELETE':
        if len(rd.keys()) < 1:
            return ("No data in the database to delete.\n")
        rd.flushdb()
        return f'Data deleted, there are {len(rd.keys())} keys in the db.\n'
    else:
        return 'The method you tried does not work.\n'

@app.route('/genes', methods = ['GET'])
def get_genes() -> list:
    """
    A route that returns a json-formatted list of all hgnc_ids.
    Args:
        None.
    Returns:
        output (list): List of all hgnc_ids.
    """

    output = []
    if len(rd.keys()) < 1:
        return ("No data available in the database. Please use a POST route first.\n")
    return rd.keys()

@app.route('/locusdata', methods = ['GET'])
def get_locusdata() -> dict:
    """
    A route that returns the amount of each locus group in the data.
    Args:
        None.
    Returns:
        locus (dict): Dictionary with how many of each locus group there are.    
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    counts = []
    groups = []
    for item in rd.keys():
        gene = json.loads(rd.get(item))
        group = gene['locus_group']
        groups.append(group)
    groupd = dict(Counter(groups))
    groupd = dict(sorted(groupd.items()))
    title = {"Locus Group": "Number of Entries"}
    title.update(groupd)
    return title

@app.route('/genes/<string:hgnc_id>', methods = ['GET'])
def get_gene(hgnc_id: str) -> dict:
    """
    A route that returns all data associated with a specified hgnc_id.
    Args:
        hgnc_id (str): The specified hgnc_id.
    Returns:
        items (dict): Dictionary with all data for the hgnc_id.
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    items = json.loads(rd.get(hgnc_id))
    return items

@app.route('/locus/<string:hgnc_id>', methods = ['GET'])
def get_locus(hgnc_id: str) -> dict:
    """
    A route that returns the locus_group of a specified gene.
    Args:
        hgnc_id (str): The specified hgnc_id.
    Returns:
        group (dict): Dictionary with locus group of the hgnc_id.
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    items = json.loads(rd.get(hgnc_id))

    group = {}
    for item in items:
        if item == "locus_group":
            group["locus group"] = items[item]

    return group

@app.route('/when/<string:hgnc_id>', methods = ['GET'])
def get_date(hgnc_id: str) -> dict:
    """
    A route that returns when a specified HGNC ID was first approved, last
    modified, or had their gene symbol or name changed. Not all of these attributes
    are available for each gene, but every attribute available is returned.
    Args:
        hgnc_id (str): The specified hgnc_id.
    Returns:
        dates (dict): Dictionary with all of the relevant dates.
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    items = json.loads(rd.get(hgnc_id))

    dates = {}
    for item in items:
        if item == "date_approved_reserved":
            dates["date first approved"] = items[item]
        elif item == "date_modified":
            dates["date last modified"] = items[item]
        elif item ==  "date_symbol_changed":
            dates["date symbol last changed"] = items[item]
        elif item == "date_name_changed":
            dates["date name last changed"] = items[item]

    return dates

@app.route('/image', methods = ['POST','GET', 'DELETE'])
def get_image():
    """
    POST: Reads data from the database, creates an image of a plot, and posts it 
        into another database. Returns confirmation of this to the user.
    GET: Retrieves the image from the database and returns it to the user.
    DELETE: Deletes the image from the database and returns confirmation of this 
        to the user.
    Args:
        POST: start (int): Starting year for the plot.
              end (int): Ending year for the plot.
        GET: None.
        DELETE: None.
    Returns:
        POST (str): Returns "Image created" message.
        GET (file): Returns the image to the user which becomes accessible using 
            scp.
        DELETE (str): Returns "Image deleted, there are 0 images in the db" message.
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    counts = []
    years = []
    start = request.args.get('start', 1986)
    end = request.args.get('end', 2023)

    if start:
        try:
            start = int(start)
        except ValueError:
            return ("Enter a positive integer for start year.", 400)

    if end:
        try:
            end = int(end)
        except ValueError:
            return ("Enter a positive integer for end year.", 400)

    if start and end:
        if start >= end:
            return ("Your start year must be less than your end year.")

    if request.method == 'POST':
        for item in rd.keys():
            gene = json.loads(rd.get(item))
            year = gene['date_approved_reserved'][0:4]
            years.append(year)
        yeard = dict(Counter(years))
        yeard = dict(sorted(yeard.items()))
        y = []
        c = []
        for item in yeard:
            y.append(item)
            c.append(yeard[item])
        sind = y.index(str(start))
        eind = y.index(str(end))+1
        y = y[sind:eind]
        c = c[sind:eind]
        plt.figure(figsize=(28,6))
        plt.bar(y, c, width = 0.35)
        plt.xlabel("Years")
        plt.ylabel("Number of Entries Approved")
        plt.title("Genes Approved Each Year")
        plt.savefig('approvalyears.png')
        file_bytes = open('./approvalyears.png', 'rb').read()
        dset = {}
        for x in range(len(y)):
            dset.update({str(y[x]):x[x]})
        rd1.set('genes_approved', file_bytes)
        #rd1.set('image_data', json.dumps(yeard))
        rd1.set('image_data', json.dumps(dset))
        return ("Image created\n")
    elif request.method == 'GET':
        path = './myapprovalyears.png'
        with open(path, 'wb') as f:
            try:
                f.write(rd1.get('genes_approved'))
            except TypeError:
                return ("No image in the database. Please use a POST route first.\n")
            f.write(rd1.get('genes_approved'))
        return send_file(path, mimetype='image/png', as_attachment=True)
    elif request.method == 'DELETE':
        if len(rd1.keys()) < 1:
            return ("No image in the database to delete. Please use a POST route first.\n")
        rd1.flushdb()
        return f'Image deleted, there are {len(rd1.keys())} images in the db.\n'
    else:
        return 'The method you tried does not work.\n'

@app.route('/imagedata', methods = ['GET'])
def get_imagedata() -> dict:
    """
    A route that returns the data used to create the plot in the /image POST route.
    Args:
        None.
    Returns:
        approved (dict): Dictionary with how many genes were approved each year.    
    """

    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    if len(rd1.keys()) < 1:
        return("No image data populated yet. Please use a POST route first.\n")

    approved = json.loads(rd1.get("image_data"))
    title = {"Years": "Number of Entries Approved"}
    title.update(approved)
    return title
              
@app.route('/help', methods = ['GET'])
def get_help() -> str:
    """
    A route that provides help test fro the user that describes each route.
    Args:
        None
    Returns:
        help (str): Help text for the user.
    """

    intro = "\nThese are the routes for the gene_api.py.\n"
    head1 = "\nretrieve elements from the data in the database\n"
    head2 = "\npost data to the database\n"
    head3 = "\ndelete data from the database\n"
    head4 = "\nget help\n"

    one ="   /data (GET)                                    Return all the data in the database\n"
    two ="   /data (POST)                                   Post the data to the database\n"
    thr ="   /data (DELETE)                                 Delete the data from the database\n"
    fou ="   /genes (GET)                                   Return a list of all HGNC IDs\n"
    fiv ="   /genes/<hgnc_id> (GET)                         Returns all of the information for a specified HGNC ID\n"
    six ="   /help (GET)                                    Return help text for the user\n"
    sev ="   /image?start=year&end=year (POST)              Generate a plot of the specified years and post it to the database\n"
    eig ="   /image (DELETE)                                Delete image from the database\n"
    nin ="   /image (GET)                                   Return image to the user\n"
    ten ="   /when/<hgnc_id> (GET)                          Return dates of approval or modification for a specified HGNC ID\n"
    ele ="   /imagedata (GET)                               Return the data used to generate the image from the /image route\n"
    twe ="   /locusdata (GET)                               Return the number of entries in each locus group\n"
    thi ="   /locus/<hgnc_id> (GET)                         Return the locus group of a specified HGNC ID\n"
    fot ="   /jobs (POST)                                   Create a new job to do some analysis of the data\n"
    return intro + head2 + two + sev + head1 + one + fou + fiv + nin + ele+ ten +twe + thi+ head3 + thr + eig + head4 + six

@app.route('/jobs', methods=['POST'])
def jobs_api():
    """
    API route for creating a new job to do some analysis. This route accepts a JSON payload
    describing the job to be created.
    """
    try:
        job = request.get_json(force=True)
    except Exception as e:
        return True, json.dumps({'status': "Error", 'message': 'Invalid JSON: {}.'.format(e)})
    return json.dumps(jobs.add_job(job['start'], job['end']))

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
