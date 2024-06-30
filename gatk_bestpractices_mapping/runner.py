import argparse
import snakemake


parser = argparse.ArgumentParser(
    prog="runner.py",
    description="Executes pipelines for WGS analysis"
)
group = parser.add_mutually_exclusive_group(required=True)

group.add_argument('--id', help = "sample IDs", nargs = "+", action = "extend")
group.add_argument('--idfile', help = "File containing a list of IDs")
parser.add_argument("--endpoint", help = "Name of endpoint requested")

script_args, unknown = parser.parse_known_args()

#ID parsing
if script_args.id is not None:
    id = script_args.id
elif script_args.idfile is not None:
    with open(script_args.idfile, "r") as f:
        id = f.readlines()
    id = [line.strip() for line in id]

snakemake.main(argv = unknown)