#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wikidataintegrator import wdi_core
import urllib.request
from bs4 import BeautifulSoup
import xlsxwriter

__author__ = "Riemer van der Vliet"
__copyright__ = "Copyright 2020, Laboratory of Systems and Synthetic Biology"
__credits__ = ["Riemer van der Vliet", "Jasper Koehorst"]
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Riemer van der Vliet"
__email__ = "riemer.vandervliet@wur.nl"
__status__ = "Development"

"""
Example script using WDI triple query. 
Parses data and retrieves wikibase IDs and values
Output excel file.
"""
# Wikibase SPARQL endpoint
endpoint_url = "https://bioparts.wiki.opencura.com/query/sparql"

# Example data
data = [[{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.18.3', 'KEGG': 'R09518'}},
         {'reaction': {'EC': '5.1.2.1', 'KEGG': 'R01450'}}],
        [{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.18.3', 'KEGG': 'R09518'}},
         {'reaction': {'EC': '1.1.1.27', 'KEGG': 'R00703'}}],
        [{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.18.3', 'KEGG': 'R09518'}},
         {'reaction': {'EC': '1.2.1.22', 'KEGG': 'R01446'}}],
        [{'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01142'}},
         {'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.1.1.27', 'KEGG': 'R00703'}}],
        [{'reaction': None},
         {'reaction': {'EC': '1.11.1.7', 'KEGG': 'R00602'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01142'}},
         {'reaction': None},
         {'reaction': {'EC': '5.1.2.1', 'KEGG': 'R01450'}}],
        [{'reaction': None},
         {'reaction': {'EC': '1.11.1.7', 'KEGG': 'R00602'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01142'}},
         {'reaction': {'EC': '1.1.1.27', 'KEGG': 'R00703'}},
         {'reaction': None}],
        [{'reaction': None},
         {'reaction': {'EC': '1.11.1.7', 'KEGG': 'R00602'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01142'}},
         {'reaction': {'EC': '1.2.1.22', 'KEGG': 'R01446'}},
         {'reaction': None}],
        [{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01143'}},
         {'reaction': {'EC': '1.1.1.27', 'KEGG': 'R00703'}}],
        [{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01143'}},
         {'reaction': {'EC': '5.1.2.1', 'KEGG': 'R01450'}}],
        [{'reaction': {'EC': '1.1.1.244', 'KEGG': 'R00605'}},
         {'reaction': {'EC': '1.14.13.25', 'KEGG': 'R01143'}},
         {'reaction': {'EC': '1.2.1.22', 'KEGG': 'R01446'}}]]


def make_query(ID: str) -> str:
    """Creates triple query using correct prefixes

    :param ID: Enzyme Commision number

    :return: query with EC number
    """

    query = """
PREFIX wd: <http://bioparts.wiki.opencura.com/entity/>
PREFIX wdt: <http://bioparts.wiki.opencura.com/prop/direct/>
PREFIX wikibase: <http://wikiba.se/ontology#>
PREFIX p: <http://bioparts.wiki.opencura.com/prop/>
PREFIX ps: <http://bioparts.wiki.opencura.com/prop/statement/>
PREFIX pq: <http://bioparts.wiki.opencura.com/prop/qualifier/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX bd: <http://www.bigdata.com/rdf#>

SELECT ?ID ?item ?type ?Seq 
Where {
    ?item p:P47 ?alignment.
    ?alignment pq:P10 \"""" + ID + """\".
    ?item wdt:P23 ?Seq.
    ?item wdt:P27 ?x.
    ?x rdfs:label ?type.
    ?item wdt:P38 ?ID.
} LIMIT 100
"""
    return query


def get_results(endpoint_url: str, query: str) -> dict:
    """Gets results

    :param endpoint_url: Wikibase endpoint
    :param query: query string with EC number

    :return results: output of sparql query
    """
    results = wdi_core.WDItemEngine.execute_sparql_query(query=query, endpoint=endpoint_url)
    return results


def parse_results(result: dict, vari: list) -> dict:
    """Parses the results and creates dictionary

    :param result: output of sparql query
    :param vari: list of variables used in query.

    :return: BB_dict with key variables and value corresponding result
    """

    # placeholder dictionary
    BB_dict = {}

    # iterates variables in list
    for item in vari:

        # retrieves sequence
        if item == "Seq":
            BB_dict[item] = get_sequence(result[item]['value'])

        # fills dictionary
        else:
            BB_dict[item] = result[item]['value']
    return BB_dict


def get_sequence(URL: str) -> str:
    """Retrieves nucleotide sequence given URL

    :param URL: URL to iGEM biobrick registry biobrick page

    :return seq: nucleotide sequence
    """

    # opens html
    with urllib.request.urlopen(URL) as response:
        html = response.read()

    # retrieves sequence from HTML body
    seq = str(BeautifulSoup(html, 'html.parser'))

    # strips line endings
    seq.rstrip()
    seq = seq[1:]

    return seq


def queries(EC: str) -> list:
    """
    gets results and splits in an output and variable list

    :param EC: string EC number

    :return: list of output dictionary or None
    """
    # placeholder list
    output = []

    # retrieves results
    print("retrieving:", EC)
    results_tot = get_results(endpoint_url, make_query(ID=EC))

    # sets list of variables
    vari = results_tot['head']['vars']

    # iterates results_tot and fills list output
    for result in results_tot['results']['bindings']:
        output.append(parse_results(result, vari))

    if len(output) <= 1:
        print("No hit for:", EC)
        output = None

    return output


def from_data(data) ->list:
    """Retrieves and parses data and checks for duplicate EC numbers in input

    :param data: Formatted data from call_full_analysis

    :return EC_list: List of individual EC numbers
    """

    # placeholder for list
    EC_list = []

    # iterates data and parses
    for line in data:
        for reaction in line:
            try:
                EC = reaction['reaction']['EC']
            except TypeError:
                pass

            # Checks for duplicates
            if EC in EC_list:
                pass
            else:
                # creates list of EC numbers
                EC_list.append(EC)

    return EC_list


def query_bioparts(data):
    """Queries the EC number and adds to workbook

    :param data: Formatted data from call_full_analysis
    :param EC: Enzyme Commission number
    """

    retrieved = {}
    EC_list = from_data(data)

    for EC in EC_list:
        retrieved[EC] = queries(EC)

    return retrieved

if __name__ == "__main__":
    retrieved = query_bioparts(data)

    for key, value in retrieved.items():
        print(key, value)

