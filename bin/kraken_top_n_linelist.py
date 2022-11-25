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
            k['percent_reads_below_this_level'] = float(line_split[0])
            k['num_reads_below_this_level'] = int(line_split[1])
            k['num_reads_at_this_level'] = int(line_split[2])
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
    kraken_report_at_specified_level_sorted = sorted(kraken_report_at_specified_level, key=lambda x: x['percent_reads_below_this_level'], reverse=True)

    output_fields = ['sample_id', 'taxonomy_level']
    output_line = {
        'sample_id': args.sample_id,
        'taxonomy_level': args.taxonomy_lvl,
    }
    
    for n in range(args.top_n):
        num = str(n + 1)

        name_field = 'abundance_' + num + '_name'
        try:
            output_line[name_field] = kraken_report_at_specified_level_sorted[n]['name']
        except IndexError as e:
            output_line[name_field] = "None"
        output_fields.append(name_field)

        taxonomy_id_field = 'abundance_' + num + '_ncbi_taxonomy_id'
        try:
            output_line[taxonomy_id_field] = kraken_report_at_specified_level_sorted[n]['taxonomy_id']
        except IndexError as e:
            output_line[taxonomy_id_field] = "None"
        output_fields.append(taxonomy_id_field)

        num_assigned_reads_field = 'abundance_' + num + '_num_assigned_reads'
        try:
            output_line[num_assigned_reads_field] = kraken_report_at_specified_level_sorted[n]['num_reads_at_this_level']
        except IndexError as e:
            output_line[num_assigned_reads_field] = 0
        output_fields.append(num_assigned_reads_field)

        fraction_total_reads_field = 'abundance_' + num + '_fraction_total_reads'
        try:
            output_line[fraction_total_reads_field] = "{:.3f}".format(kraken_report_at_specified_level_sorted[n]['percent_reads_below_this_level'] / 100)
        except IndexError as e:
            output_line[fraction_total_reads_field] = 0.0
        output_fields.append(fraction_total_reads_field)

    unclassified_name_field = 'unclassified_name'
    output_line[unclassified_name_field] = 'unclassified'
    output_fields.append(unclassified_name_field)

    unclassified_taxonomy_id_field = 'unclassified_ncbi_taxonomy_id'
    output_line[unclassified_taxonomy_id_field] = '0'
    output_fields.append(unclassified_taxonomy_id_field)

    num_unclassified_reads_field = 'unclassified_num_assigned_reads'
    try:
        output_line[num_unclassified_reads_field] = kraken_report_unclassified['num_reads_at_this_level']
    except IndexError as e:
        output_line[num_unclassified_reads_field] = 0
    output_fields.append(num_unclassified_reads_field)

    fraction_unclassified_reads_field = 'unclassified_fraction_total_reads'
    try:
        output_line[fraction_unclassified_reads_field] = kraken_report_unclassified['percent_reads_below_this_level']
    except IndexError as e:
        output_line[fraction_unclassified_reads_field] = 0.0
    output_fields.append(fraction_unclassified_reads_field)
        

    csv.register_dialect('unix-csv-quote-minimal', delimiter=',', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix-csv-quote-minimal')
    writer.writeheader()
    writer.writerow(output_line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('kraken_report')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-l', '--taxonomy-lvl', default='S')
    parser.add_argument('-n', '--top-n', type=int, default=5)
    args = parser.parse_args()
    main(args)
