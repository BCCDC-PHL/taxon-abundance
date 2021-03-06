#!/usr/bin/env python

import argparse
import csv
import sys
import re
import json


def parse_bracken_report(bracken_report_path):
    bracken_report_lines = []
    with open(bracken_report_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            b = {}
            b['sample_id'] = row['sample_id']
            b['name'] = row['name']
            b['taxonomy_id'] = row['taxonomy_id']
            b['taxonomy_lvl'] = row['taxonomy_lvl']
            b['kraken_assigned_reads'] = int(row['kraken_assigned_reads'])
            b['added_reads'] = int(row['added_reads'])
            b['new_est_reads'] = int(row['new_est_reads'])
            b['fraction_total_reads'] = float(row['fraction_total_reads'])
            bracken_report_lines.append(b)

    return bracken_report_lines
        

def main(args):
    bracken_report = parse_bracken_report(args.bracken_report)

    bracken_report_unclassified = list(filter(lambda x: x['name'] == 'unclassified', bracken_report))[0]
    bracken_report_non_unclassified = list(filter(lambda x: x['name'] != 'unclassified', bracken_report))
    
    bracken_report_sorted = sorted(bracken_report_non_unclassified, key=lambda k: k['fraction_total_reads'], reverse=True)

    output_fields = ['sample_id', 'taxonomy_level']
    output_line = {
        'sample_id': args.sample_id,
        'taxonomy_level': bracken_report_sorted[0]['taxonomy_lvl']
    }
    
    for n in range(args.top_n):
        num = str(n + 1)

        name_field = 'abundance_' + num + '_name'
        try:
            output_line[name_field] = bracken_report_sorted[n]['name']
        except IndexError as e:
            output_line[name_field] = "None"
        output_fields.append(name_field)

        taxonomy_id_field = 'abundance_' + num + '_ncbi_taxonomy_id'
        try:
            output_line[taxonomy_id_field] = bracken_report_sorted[n]['taxonomy_id']
        except IndexError as e:
            output_line[taxonomy_id_field] = "None"
        output_fields.append(taxonomy_id_field)

        num_assigned_reads_field = 'abundance_' + num + '_num_assigned_reads'
        try:
            output_line[num_assigned_reads_field] = bracken_report_sorted[n]['new_est_reads']
        except IndexError as e:
            output_line[num_assigned_reads_field] = 0
        output_fields.append(num_assigned_reads_field)

        fraction_total_reads_field = 'abundance_' + num + '_fraction_total_reads'
        try:
            output_line[fraction_total_reads_field] = bracken_report_sorted[n]['fraction_total_reads']
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
        output_line[num_unclassified_reads_field] = bracken_report_unclassified['new_est_reads']
    except IndexError as e:
        output_line[num_unclassified_reads_field] = 0
    output_fields.append(num_unclassified_reads_field)

    fraction_unclassified_reads_field = 'unclassified_fraction_total_reads'
    try:
        output_line[fraction_unclassified_reads_field] = bracken_report_unclassified['fraction_total_reads']
    except IndexError as e:
        output_line[fraction_unclassified_reads_field] = 0.0
    output_fields.append(fraction_unclassified_reads_field)
        

    csv.register_dialect('unix-csv-quote-minimal', delimiter=',', doublequote=False, lineterminator='\n', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fields, dialect='unix-csv-quote-minimal')
    writer.writeheader()
    writer.writerow(output_line)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bracken_report')
    parser.add_argument('-s', '--sample-id')
    parser.add_argument('-n', '--top-n', type=int)
    args = parser.parse_args()
    main(args)
