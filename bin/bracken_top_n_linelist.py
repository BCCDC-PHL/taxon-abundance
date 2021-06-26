#!/usr/bin/env python

import argparse
import csv
import sys
import re
import json


def parse_bracken_report(bracken_report_path):
    bracken_report_lines = []
    with open(bracken_report_path, 'r') as f:
        reader = csv.DictReader(f, dialect='unix')
        for row in reader:
            bracken_report_lines.append(row)

    return bracken_report_lines
        

def main(args):
    bracken_report = parse_bracken_report(args.bracken_report)

    bracken_report_sorted = sorted(bracken_report, key=lambda k: k['fraction_total_reads'], reverse=True) 
    
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
        fraction_total_reads_field = 'abundance_' + num + '_fraction_total_reads'
        try:
            output_line[fraction_total_reads_field] = bracken_report_sorted[n]['fraction_total_reads']
        except IndexError as e:
            output_line[fraction_total_reads_field] = 0.0
        output_fields.append(fraction_total_reads_field)
        

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
