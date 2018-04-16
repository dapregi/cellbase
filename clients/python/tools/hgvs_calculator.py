#!/usr/bin/env python

import sys
import os
import argparse
import logging
from itertools import islice
from pycellbase.cbclient import ConfigClient
from pycellbase.cbclient import CellBaseClient

_DEFAULT_HOST = 'bioinfo.hpc.cam.ac.uk:80/cellbase'
_DEFAULT_API_VERSION = 'v4'
_DEFAULT_SPECIES = 'hsapiens'
_DEFAULT_ASSEMBLY = 'GRCh38'

# Reference sequence notation
# http://www.hgvs.org/mutnomen/standards.html
_HGVS_REF_SEQ_LETTER = {
    'genomic': 'g.', 'cdna': 'c.', 'rna': 'r.', 'protein': 'p.',
    'mitochondrial': 'm.', 'noncoding': 'n.'
}


def _parse_arguments():
    """Parse arguments"""

    desc = 'This tool returns HGVS notation for variants'
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '-v', '--verbose', dest='verbosity', action='store_true',
        help='increase output verbosity'
    )
    parser.add_argument(
        'input',
        help='input file or comma-separated string\n(e.g. "19:45411941:T:C")'
    )
    parser.add_argument(
        '-o', '--output', dest='output_fpath', default=sys.stdout,
        help='output file path (default: stdout)'
    )
    parser.add_argument(
        '-r', '--ref_seq_type', dest='ref_seq_type',
        choices=_HGVS_REF_SEQ_LETTER.keys(),
        help='reference sequence type'
    )
    parser.add_argument(
        '--assembly', dest='assembly', choices=['grch37', 'grch38'],
        help='reference genome assembly (default: ' + _DEFAULT_ASSEMBLY + ')'
    )
    parser.add_argument(
        '--species', dest='species',
        help=('species name (default:' + _DEFAULT_SPECIES +
              ').\nOverrides configuration provided with "--config" parameter')
    )
    parser.add_argument(
        '--host', dest='host',
        help=('web services host (default: ' + _DEFAULT_HOST +
              ').\nOverrides configuration provided with "--config" parameter')
    )
    parser.add_argument(
        '--version', dest='api_version',
        help=('api version (default: ' + _DEFAULT_API_VERSION +
              ').\nOverrides configuration provided with "--config" parameter')
    )
    parser.add_argument(
        '--config', dest='config',
        help='CellBase configuration. Overrides default values.'
    )

    args = parser.parse_args()
    return args


def _set_logger(verbosity=False):
    """Set logging system"""

    if verbosity:
        logging.basicConfig(format="[%(levelname)s] %(message)s",
                            level=logging.DEBUG)
    else:
        logging.basicConfig(format="[%(levelname)s] %(message)s",
                            level=logging.WARNING)


def _read_file_in_chunks(fhand, number_of_lines=100, remove_empty=True):
    """Read file and retrieve a specified number of lines"""

    while True:
        n_lines = map(lambda x: x.rstrip(),
                      list(islice(fhand, number_of_lines)))

        # If end of file
        if not n_lines:
            break

        # Remove empty lines
        if remove_empty:
            n_lines = list(filter(None, n_lines))
        if not n_lines:
            continue

        yield n_lines


def _get_hgvs(query, cellbase_client, ref_seq_type, assembly):

    # Setting up CellBase Variant client
    vc = cellbase_client.get_variant_client()

    response = vc.get_annotation(query, include='hgvs', assembly=assembly)

    for i, query_response in enumerate(response):
        hgvs_list = []

        # Getting list of all hgvs notations
        for hgvs in query_response['result'][0]['hgvs']:
            hgvs_list.append(hgvs)

        # Setting up list to filter hgvs notations
        filtering = [1]*len(hgvs_list)

        # Checking reference sequence type
        if ref_seq_type is not None:
            for j, hgvs in enumerate(hgvs_list):
                notation = hgvs.split(':')[1]
                if _HGVS_REF_SEQ_LETTER[ref_seq_type] not in notation:
                    filtering[j] = 0

        yield (query.split(',')[i],
               [hgvs for i, hgvs in enumerate(hgvs_list) if filtering[i]])


def calculate_hgvs(input_data, output_fpath, cbc, ref_seq_type, assembly):

    # Checking output
    if output_fpath is sys.stdout:
        output_fhand = output_fpath
    else:
        output_fhand = open(output_fpath, 'w')

    # Creating input generator
    number_of_items = 100
    input_fhand = None
    if os.path.isfile(input_data):
        input_fhand = open(input_data, 'r')
        input_gen = _read_file_in_chunks(input_fhand,
                                         number_of_lines=number_of_items)
    else:
        input_data = input_data.split(',')
        input_gen = [input_data[i:i+number_of_items]
                     for i in range(0, len(input_data), number_of_items)]

    for query in input_gen:
        query = ','.join(query)
        for variant, hgvs in _get_hgvs(query, cbc, ref_seq_type, assembly):
            if not hgvs:
                hgvs = ['.']
            output_fhand.write('\t'.join([variant, ';'.join(hgvs)]) + '\n')

    output_fhand.close()
    if input_fhand is not None:
        input_fhand.close()


def main():

    # Getting args
    args = _parse_arguments()

    # Setting up logging system
    _set_logger(args.verbosity)

    # Setting up PyCellBase clients
    cc = ConfigClient(
        {"species": _DEFAULT_SPECIES, "version": _DEFAULT_API_VERSION,
         "rest": {"hosts": [_DEFAULT_HOST]}}
    )

    if args.config is not None:
        cc = ConfigClient(args.config)
    if args.species is not None:
        cc.species = args.species
    if args.api_version is not None:
        cc.version = args.api_version
    if args.host is not None:
        cc.host = args.host
    if args.assembly is not None:
        assembly = args.assembly
    else:
        assembly = _DEFAULT_ASSEMBLY
    cbc = CellBaseClient(cc)

    calculate_hgvs(args.input, args.output_fpath, cbc, args.ref_seq_type,
                   assembly)


if __name__ == '__main__':
    sys.exit(main())
