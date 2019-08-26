
import re
import gzip

# open gzipped or regular text vcf files
def open_file_for_reading(filename):
    if re.search("\.gz$", filename):
        return gzip.open(filename, 'rt') # t needed for python3 to look like regular text
    else:
        return open(filename, 'r') # regular text file
