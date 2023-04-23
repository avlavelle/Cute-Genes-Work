#!/usr/bin/env python3

#Anna Victoria Lavelle
#April 20, 2023

from flask import Flask, request, send_file
from collections import Counter, OrderedDict
import requests
import redis
import json
import os
import matplotlib.pyplot as plt
import numpy as np

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
    A route that posts the data, deletes the data, or outputs the data depending on 
    the user's query.

    Args:
        None.

    Returns:
        output_list (list): Returns the redis data.
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

@app.route('/genes/<string:hgnc_id>', methods = ['GET'])
def get_gene(hgnc_id: str) -> dict:
    """
    A route that returns all data associated with a hgnc_id.

    Args:
        hgnc_id (str): The specified hgnc_id.

    Returns:
        items (dict): Dictionary with all data for the hgnc_id.
    """
    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    items = json.loads(rd.get(hgnc_id))
    return items

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

@app.route('/image', methods = ['POST','GET', 'DELETE'])
def get_image():
    """
    A route that reads data from the database, creates an image of a plot, posts it in    to the database, gets the image from the database, and deletes the image from the
    database.

    Args:
        None.

    Returns:
        (file_bytes): Returns the file_bytes of the image in the POST route.
        send_file(): Returns the image to the user in the GET route.
    """
    if len(rd.keys()) < 1:
        return ("No data in the database. Please use a POST route first.\n")
    years = []
    counts = []
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
        plt.figure(figsize=(28,6))
        plt.bar(y, c, width = 0.35)
        plt.xlabel("Years")
        plt.ylabel("Number of Entries Approved")
        plt.title("Genes Approved Each Year")
        plt.savefig('approvalyears.png')
        file_bytes = open('./approvalyears.png', 'rb').read()
        rd1.set('genes_approved', file_bytes)
        return ("Image created.\n")
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
        return f'Image deleted, there are {len(rd1.keys())} keys in the db.\n'
    else:
        return 'The method you tried does not work.\n'
    
    






if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
