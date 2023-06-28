#!/usr/bin/env python3

from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

# Entrez email here
Entrez.email = "Your Email address here"

# retrieve records as list
handle = Entrez.esearch(db="assembly", retmax=30000, term='"Escherichia coli"[Organism] AND latest_refseq[filter]', idtype="acc")

record = Entrez.read(handle, validate=False)

ids = record['IdList']

handle.close()

# define function to get biosample attributes from Entrez

def get_biosample_attributes(identifier):

    # assign None to variables to catch missing data
    isolation_source = None
    host = None
    sample_type = None
    location = None
    date = None

    try:
        
        # get biosample id 
        elink = Entrez.read(Entrez.elink(dbfrom="assembly", db='biosample', id=identifier), validate=False)  
        biosample_id = elink[0]['LinkSetDb'][0]['Link'][0]['Id']

        # get refseq accession numbers 
        summary = Entrez.read(Entrez.esummary(db="biosample", id=biosample_id))
        get_accession = Entrez.read(Entrez.esummary(db="assembly", id=identifier), validate=False)
        refseq_accession = get_accession['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']

        # get named attributes
        attributes = summary['DocumentSummarySet']['DocumentSummary'][0]['SampleData']
        tree = ET.fromstring(attributes)
        for item in tree.iter():
            if item.tag == 'Attribute':
                if item.attrib['attribute_name'] == "isolation_source":
                    isolation_source = (item.text)
                if item.attrib['attribute_name'] == "host":
                    host = (item.text)
                if item.attrib['attribute_name'] == "sample_type":
                    sample_type = (item.text)
                if item.attrib['attribute_name'] == "geo_loc_name":
                    location = (item.text)
                if item.attrib['attribute_name'] == "collection_date":
                    date = (item.text)
        
        print(biosample_id, refseq_accession, isolation_source, host, sample_type, location, date)

        return biosample_id, refseq_accession, isolation_source, host, sample_type, location, date
    
    except:
        biosample_id = None
        return biosample_id, refseq_accession, isolation_source, host, sample_type, location, date
    
# create empty lists to append sample data to
UIDs = []
refseq_accessions = []
biosample_ids = []
isolation_sources = []
hosts = []
sample_types = []
locations = []
dates = []

# loop counter
n = 0

for identifier in ids:
    # loop counter
    n = n + 1
    print(n)
    # get biosample data using function
    biosample_id, refseq_accession, isolation_source, host, sample_type, location, date = get_biosample_attributes(identifier)
    # append results to list
    UIDs.append(identifier)
    refseq_accessions.append(refseq_accession)
    biosample_ids.append(biosample_id)
    isolation_sources.append(isolation_source)
    hosts.append(host)
    sample_types.append(sample_type)
    locations.append(location)
    dates.append(date)

# zip lists into dataframe
biosample_df = pd.DataFrame(list(zip(UIDs, refseq_accessions, biosample_ids, 
                                     isolation_sources, hosts, sample_types,
                                     locations, dates)), 
            columns =['UID', 'refseq_id', 'biosample_id',
                      'isolation_source', 'host', 'sample_type',
                      'location', 'date'])

final_biosample_df = biosample_df.fillna("missing") 

final_biosample_df.to_csv('../../results/biosample_attributes_EC.csv')
