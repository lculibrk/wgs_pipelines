import argparse
import genomicsapi
from genomicsapi import connection
import genomicsapi.select
import genomicsapi.inserts
import genomicsapi.translate
import genomicsapi.update
import os

parser = argparse.ArgumentParser(
    prog="mark_complete.py",
    description="Marks an analysis as complete, updates the output files entry and updates the timestamp in the genomicsdb"
)

parser.add_argument("-i", "--id", help = "analysis ID")
parser.add_argument("-d", "--db", help = "database name")
parser.add_argument("outputs", nargs = "+", help = "List of output files")

args = parser.parse_args()

analysis_id = genomicsapi.translate.stringtoid(args.id)

from datetime import datetime
now = datetime.now()
timestamp = now.strftime('%Y-%m-%d %H:%M:%S')
genomicsapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_time", timestamp)
genomicsapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_complete", "True")
genomicsapi.update.update(args.db, "analyses", {"id":analysis_id}, "analysis_dir", os.path.dirname(args.outputs[0]))