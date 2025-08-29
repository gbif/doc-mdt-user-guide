#!/usr/bin/env python3

import argparse
import csv
import gzip
import os

def open_maybe_gz(filename):
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

def read_sample_ids(sample_file):
    sample_ids = set()
    with open_maybe_gz(sample_file) as f:
        reader = csv.reader(f, delimiter='\t' if sample_file.endswith('.tsv') else ',')
        next(reader, None)  # Skip header
        for row in reader:
            if row:
                sample_ids.add(row[0])
    return sample_ids

def read_taxonomy_ids(tax_file):
    otu_ids = set()
    with open_maybe_gz(tax_file) as f:
        reader = csv.reader(f, delimiter='\t' if tax_file.endswith('.tsv') else ',')
        next(reader, None)  # Skip header
        for row in reader:
            if row:
                otu_ids.add(row[0])
    return otu_ids

def read_otu_table_ids(otu_file, transposed=False):
    with open_maybe_gz(otu_file) as f:
        reader = csv.reader(f, delimiter='\t' if otu_file.endswith('.tsv') else ',')
        header = next(reader)
        if transposed:
            otu_ids = set(header[1:])  # skip first column
            sample_ids = set()
            for row in reader:
                if row:
                    sample_ids.add(row[0])
        else:
            sample_ids = set(header[1:])  # skip first column
            otu_ids = set()
            for row in reader:
                if row:
                    otu_ids.add(row[0])
    return sample_ids, otu_ids

def write_id_list(ids, path):
    with open(path, 'w') as f:
        for _id in sorted(ids):
            f.write(f"{_id}\n")

def main(otu_file, sample_file, tax_file, transposed, outdir):
    print("Reading sample table...")
    sample_ids = read_sample_ids(sample_file)
    print("Reading taxonomy table...")
    tax_ids = read_taxonomy_ids(tax_file)
    print("Reading OTU table...")
    otu_sample_ids, otu_ids = read_otu_table_ids(otu_file, transposed=transposed)

    # Sample comparison
    only_in_sample = sample_ids - otu_sample_ids
    only_in_otu_samples = otu_sample_ids - sample_ids
    common_samples = sample_ids & otu_sample_ids

    # OTU comparison
    only_in_tax = tax_ids - otu_ids
    only_in_otu_otus = otu_ids - tax_ids
    common_otus = tax_ids & otu_ids

    # Report
    print(f"\nSample table contains {len(sample_ids)} sample IDs")
    print(f"OTU table contains {len(otu_sample_ids)} sample IDs")
    print(f"{len(common_samples)} sample IDs are common")
    print(f"{len(only_in_sample)} sample IDs only in sample table")
    print(f"{len(only_in_otu_samples)} sample IDs only in OTU table")

    print(f"\nTaxonomy table contains {len(tax_ids)} OTU IDs")
    print(f"OTU table contains {len(otu_ids)} OTU IDs")
    print(f"{len(common_otus)} OTU IDs are common")
    print(f"{len(only_in_tax)} OTU IDs only in taxonomy table")
    print(f"{len(only_in_otu_otus)} OTU IDs only in OTU table")

    # Optional output
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        if only_in_sample:
            write_id_list(only_in_sample, os.path.join(outdir, "only_in_sample_table.txt"))
        if only_in_otu_samples:
            write_id_list(only_in_otu_samples, os.path.join(outdir, "only_in_otu_table_samples.txt"))
        if only_in_tax:
            write_id_list(only_in_tax, os.path.join(outdir, "only_in_taxonomy_table.txt"))
        if only_in_otu_otus:
            write_id_list(only_in_otu_otus, os.path.join(outdir, "only_in_otu_table_otus.txt"))
        print(f"\nWrote non-matching ID lists to directory: {outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare OTU, sample, and taxonomy tables")
    parser.add_argument('--otu', required=True, help="OTU table file (TSV or CSV, optionally .gz)")
    parser.add_argument('--samples', required=True, help="Sample metadata file (TSV or CSV, optionally .gz)")
    parser.add_argument('--taxonomy', required=True, help="Taxonomy metadata file (TSV or CSV, optionally .gz)")
    parser.add_argument('--otu-transposed', action='store_true',
                        help="Set this if OTU table has OTU IDs as column headers and sample IDs as row names")
    parser.add_argument('--outdir', default=None, help="Output directory for non-matching ID lists (optional)")
    args = parser.parse_args()
    main(args.otu, args.samples, args.taxonomy, args.otu_transposed, args.outdir)
