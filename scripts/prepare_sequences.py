#!/usr/bin/env python2
import argparse
import gzip
import os
import sys
from collections import namedtuple

import pysam

Variant = namedtuple('Variant', ['chrom', 'pos', 'ref', 'alt'], verbose=False)
Region = namedtuple('Region', ['variant_idx', 'sample_idx', 'is_hom', 'window'], verbose=False)
variants = []

def parse_args():
    parser = argparse.ArgumentParser(description='')
    subparsers = parser.add_subparsers(help='', dest='command')

    samples_parser = subparsers.add_parser('samples', help='Extracts all unique sample IDs from the HET and HOM INFO fields.')
    samples_parser.add_argument('-i', '--in', metavar='file', dest='input_file', required=True, help='Input compressed VCF.')

    bam_file_parser = subparsers.add_parser('cram', help='Generates CRAM file with sequences from heterozygous/homozygous samples.')
    bam_file_parser.add_argument('-i', '--in', metavar='file', dest='input_file', required=True, help='Input compressed VCF with HET and HOM INFO fields. Multi-allelic variants must be split into bi-allelic entries.')
    bam_file_parser.add_argument('-c', '--crams', metavar='/file/path/sample2 /file/path/sample2', dest='input_crams', required=True, nargs='*', help='Input paths to sample files seperated by space.')
    bam_file_parser.add_argument('-w', '--window', metavar='base-pair', dest='window', type=int, required=False, default=100, help='Window size around each variant in base-pairs.')
    bam_file_parser.add_argument('-o', '--out', metavar='file', dest='output_file', required=True, help='Unsorted output CRAM file.')
    return parser.parse_args()

def process_sample(cram_path, regions, ocram):
    with pysam.AlignmentFile(cram_path, 'rc') as icram:
        for region in regions:
            process_region(icram, region, ocram)

def process_region(icram, region, ocram):
    global variants
    variant = variants[region.variant_idx]
    qnames = dict()
    for read in icram.fetch(variant.chrom, variant.pos - region.window if variant.pos > region.window else 0, variant.pos + region.window):
        qname = qnames.get(read.query_name, None)
        if qname is None:
            qname = str(len(qnames) + 1)
            qnames[read.query_name] = qname
        read.query_name = '{}:{}:{}:{}{}:{}'.format(variant.pos, variant.ref, variant.alt, '0' if region.is_hom else '', region.sample_idx, qname)
        for tag, _ in read.get_tags():
            read.set_tag(tag, None)
        ocram.write(read)

def main():
    args = parse_args()

    if args.command == 'samples':
        unique_samples = set()
        with gzip.GzipFile(args.input_file, 'r') as ifile:
            for line in ifile:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                info = dict(map(lambda x: (x[0], x[1]) if len(x) > 1 else (x[0], None), (x.split('=', 1) for x in fields[7].split(';'))))
                hets = info.get('HET', None)
                homs = info.get('HOM', None)
                if hets is not None:
                    for sample in hets.strip().split(','):
                        unique_samples.add(sample.strip())
                if homs is not None:
                    for sample in homs.strip().split(','):
                        unique_samples.add(sample.strip())
        for sample in unique_samples:
            sys.stdout.write('{}\n'.format(sample))
    elif args.command == 'cram':
        header = {'HD': {'SO': 'coordinate', 'VN': '1.3'}, 'SQ': [], 'RG': []}

        crams = dict((os.path.basename(n), n) for n in args.input_crams if os.path.isfile(n))

        if len(crams) == 0:
            sys.exit(1)

        samples = dict()

        with gzip.GzipFile(args.input_file, 'r') as ifile:
            for line in ifile:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = dict(map(lambda x: (x[0], x[1]) if len(x) > 1 else (x[0], None), (x.split('=', 1) for x in fields[7].split(';'))))
                hets = info.get('HET', None)
                homs = info.get('HOM', None)
                variant_idx = len(variants)
                variants.append(Variant(chrom, pos, ref, alt))
                sample_idx = 0
                if hets is not None:
                    for sample in (x.strip() for x in hets.strip().split(',')):
                        if sample in crams:
                            if sample not in samples:
                                samples[sample] = []
                            sample_idx += 1
                            samples[sample].append(Region(variant_idx, sample_idx, False, args.window))
                if homs is not None:
                    for sample in (x.strip() for x in homs.strip().split(',')):
                        if sample in crams:
                            if sample not in samples:
                                samples[sample] = []
                            sample_idx += 1
                            samples[sample].append(Region(variant_idx, sample_idx, True, args.window))

        with pysam.AlignmentFile(crams[next(iter(crams))], 'rc') as icram:
            for sq_line in icram.header['SQ']:
                header['SQ'].append(sq_line)

        i = 0
        with pysam.AlignmentFile(args.output_file, 'wc', header=header) as ocram:
            for sample, sample_variants in samples.iteritems():
                process_sample(crams[sample], sample_variants, ocram)
                i += 1
                if i % 100 == 0:
                    sys.stdout.write('Processed {}/{} sample(s).\n'.format(i, len(samples)))
        sys.stdout.write('Done ({}/{}).\n'.format(i, len(samples)))

if __name__ == '__main__':
    main()
