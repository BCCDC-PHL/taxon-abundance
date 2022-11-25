#!/usr/bin/env python

import argparse
import csv
import sys
import re
import json


def parse_kraken_report(kraken_report_path):
    kraken_report_lines = []
    with open(kraken_report_path, 'r') as f:
        for line in f:
            line_split = line.strip().split(maxsplit=5)
            k = {}
            k['fraction_total_reads'] = float(line_split[0])
            k['kraken_assigned_reads'] = int(line_split[1])
            # k['num_reads_at_this_level'] = int(line_split[2])
            k['taxonomy_lvl'] = line_split[3]
            k['taxonomy_id'] = line_split[4]
            k['name'] = line_split[5]

            kraken_report_lines.append(k)

    return kraken_report_lines
        

def main(args):
    kraken_report = parse_kraken_report(args.kraken_report)
    
    kraken_report_unclassified = list(filter(lambda x: x['name'] == 'unclassified', kraken_report))[0]
    kraken_report_classified = list(filter(lambda x: x['name'] != 'unclassified', kraken_report))
    kraken_report_at_specified_level = list(filter(lambda x: x['taxonomy_lvl'] == args.taxonomy_lvl, kraken_report_classified))

    output_fields = [
        'sample_id',
        'name',
        'taxonomy_id',
        'taxonomy_lvl',
        'kraken_assigned_reads',
        'fraction_total_reads',
    ]

    output_lines = []
    kraken_report_unclassified['sample_id'] = args.sample_id
    output_lines.append(kraken_report_unclassified)
    for kraken_record in kraken_report_at_specified_level:
        kraken_record['sample_id'] = args.sample_id
        output_lines.append(kraken_record)

    output_lines_sorted = sorted(output_lines, key=lambda x: x['kraken_assigned_reads'], reverse=True)
        

    csv.register_dialect('unix-csv-quote-minimal', delimiter=',', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix-csv-quote-minimal')
    writer.writeheader()

    for line in output_lines_sorted:
        writer.writerow(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('kraken_report')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-l', '--taxonomy-lvl', default='S')
    args = parser.parse_args()
    main(args)
