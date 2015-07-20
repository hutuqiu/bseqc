import sys

import logging

# ------------------------------------
#logging object
# ------------------------------------


logging.basicConfig(level=20,
    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr,
    filemode="w"
)
info=logging.info
error=logging.error
warning=logging.warning


def get_ref(ref_file):
    '''
    get the reference sequences, and save it as a dict
    '''
    info("Get reference sequences...")
    ref, cr, seq = {}, '', ''
    for line in open(ref_file):
    #read reference fasta
        if line.startswith('>'):
            if cr:
                ref[cr] = seq.upper()
            cr, seq = line[1:-1].split()[0], ''
        else:
            seq += line.rstrip()
    ref[cr] = seq.upper()
    info("The reference has been got!")
    return ref