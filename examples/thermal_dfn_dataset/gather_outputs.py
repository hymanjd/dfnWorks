"""
Gather finished runs' output data into results/<jobname>/ and delete their dfnWorks
job folders (what driver_v2.py now does at the end of every run):

python gather_outputs.py lhs_hf_*                 # gather and delete
python gather_outputs.py lhs_hf_* --keep          # gather only
python gather_outputs.py lhs_hf_* --results other_results

Anything that is not a folder (e.g. lhs_hf_1.log, matched by lhs_hf_*) is skipped, and
a folder is deleted only if it is a dfnWorks job with a thermal.npz that was gathered.
"""
import os
import sys
import argparse

from helper_functions import gather_and_clean

parser = argparse.ArgumentParser(description="Gather run outputs into results/ and delete the job folders")
parser.add_argument("jobdirs", nargs="+", help="run folders, e.g. lhs_hf_*")
parser.add_argument("--results", default="results", help="folder to gather into (default results)")
parser.add_argument("--keep", action="store_true", help="gather only; keep the job folders")
args = parser.parse_args()

results = os.path.abspath(args.results)
n = 0
for d in args.jobdirs:
    # skip files (e.g. lhs_hf_1.log) and the results folder itself
    if not os.path.isdir(d) or os.path.abspath(d) == results:
        continue
    if not os.path.isfile(os.path.join(d, "thermal.npz")):
        print(f"--> Skipping {d}: no thermal.npz (run unfinished or failed)")
        continue
    # copy the outputs to results/<jobname>/, then delete the job folder unless --keep
    gather_and_clean(d, results, cleanup=not args.keep)
    n += 1
print(f"--> {n} runs gathered into {results}")
