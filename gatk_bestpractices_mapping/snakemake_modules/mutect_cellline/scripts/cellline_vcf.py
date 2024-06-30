import argparse
import re
import os

parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", help = "Input VCF")
parser.add_argument("-c", "--cmdline", help = "commandline used to generate the vcf we're reheadering")
#parser.add_argument("-o", "--output", required=False)

args = parser.parse_args()

## We're editing the header
pattern = r'##'

# Open the text file for reading
## This is a hack but it's few enough files to be ok
lines = []
with open(args.vcf, 'r') as file:
    for line in file:
        # Strip newline characters from the line
        line = line.strip()
        # Check if the line matches the regular expression pattern
        if re.match(pattern, line):
            lines.append(line)  # You can process or store the matched line here
        else:
            break  # Stop reading further lines if the pattern doesn't match

## Now first handle the new filters
filt_ind = max([index for index, value in enumerate(lines) if re.match("##FILTER", value)])+1
lines.insert(filt_ind, r'##FILTER=<ID=shared_external,Description="variant found in unrelated daughters">')
lines.insert(filt_ind, r'##FILTER=<ID=shared_parental,Description="variant found in >50% of study parents">')
lines.insert(filt_ind, r'##FILTER=<ID=depth,Description="variant has <15 read depth in parent">')


format_ind = max([index for index, value in enumerate(lines) if re.match("##FORMAT", value)]) + 1
lines.insert(format_ind, r'##FORMAT=<ID=LAB,Number=1,Type=String,Description="Label describing sharedness of mutation among cell line parents and other daughters">')

## Add the command line to this
with open(args.cmdline, "r") as file:
    cmdline = file.readlines()
mtime = os.path.getmtime(args.cmdline)

if len(cmdline) != 1:
    raise ValueError("cmdline file does not have exactly one line!")

cmdline = cmdline[0].strip()

info_ind = min([index for index, value in enumerate(lines) if re.match("##INFO", value)]) + 1

cmdline_entry = f"FilteringCommand=<ID=cellline_filtering,CommandLine={cmdline},Version=\"1.0\",Date=\"{mtime}\""
cmdline = r'##' + cmdline_entry
lines.insert(info_ind, cmdline)

[print(line) for line in lines]
