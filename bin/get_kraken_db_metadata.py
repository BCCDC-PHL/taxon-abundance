#!/usr/bin/env python3

import argparse
import json
import os

import yaml

def main(args):
    db_metadata_path = os.path.join(args.db, 'metadata.json')
    db_metadata = {}
    if os.path.exists(db_metadata_path):
        with open(db_metadata_path, 'r') as f:
            db_metadata = json.load(f)

    if not db_metadata:
        exit()

    provenance = {}
    if args.provenance:
        with open(args.provenance, 'r') as f:
            provenance = yaml.safe_load(f)

    database_provenance = {}
    database_name = db_metadata.get('dbname', None)
    if database_name:
        database_provenance['database_name'] = database_name
    database_version = db_metadata.get('version', None)
    if database_version:
        database_provenance['database_version'] = database_version

    for provenance_record in provenance:
        if provenance_record.get('process_name') == 'kraken2':
            provenance_record['databases'] = []
            provenance_record['databases'].append(database_provenance)

    with open(args.provenance, 'w') as f:
        yaml.dump(provenance, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get metadata from Kraken database')
    parser.add_argument('--provenance', help='Path to provenance file')
    parser.add_argument('--db', help='Path to Kraken database')
    args = parser.parse_args()
    main(args)
