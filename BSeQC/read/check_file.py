#check the input mapping files
import os
import sys
import logging

logging.basicConfig(level=20,
    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    stream=sys.stderr,
    filemode="w"
)
info=logging.info
error=logging.error
warning=logging.warning

def check_mapping_file_header(mapping_file, s_path):
    if mapping_file[-4:].upper() == '.SAM':
        sam_format = 1
        read_inf = open(mapping_file)
    elif mapping_file[-4:].upper() == '.BAM':
        sam_format, read_inf = 1, os.popen('%ssamtools view -h %s' % (s_path, mapping_file))
    else:
        error("The input mapping file is not SAM format or BAM format")
    return sam_format, read_inf

def check_mapping_file(mapping_file, s_path):
    if mapping_file[-4:].upper() == '.SAM':
        sam_format,read_inf =1, os.popen("%ssamtools view -S %s" %(s_path,mapping_file))
    elif mapping_file[-4:].upper() == '.BAM':
        sam_format, read_inf = 1, os.popen('%ssamtools view %s' % (s_path, mapping_file))
    else:
        error("The input mapping file is not SAM format or BAM format")
    return sam_format, read_inf

def check_mapping_file_flag(mapping_file, s_path):
    if mapping_file[-4:].upper() == '.BAM':
        sam_format, read_inf =1, os.popen("%ssamtools view -X %s" %(s_path,mapping_file))
    elif mapping_file[-4:].upper() == '.SAM':
        sam_format, read_inf =1, os.popen("%ssamtools view -XS %s" %(s_path,mapping_file))
    else:
        error("The input mapping file is not SAM format or BAM format")
    return sam_format, read_inf



