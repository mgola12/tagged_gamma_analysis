import uproot
import awkward as ak
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--run_file", type=str, required=True)
parser.add_argument("--n_events", type=int, default=1000)
args = parser.parse_args()

with uproot.open(args.run_file) as f:
    tree = f["WCTEReadoutWindows"]
    print(f"Total entries in tree: {tree.num_entries}")
    print(f"Reading first {args.n_events} events\n")

    branches = ["beamline_id", "beamline_id_name"]
    chunk = tree.arrays(branches, entry_start=0, entry_stop=args.n_events, library="ak")

    seen = {}
    for event in chunk:
        ids   = ak.to_list(event["beamline_id"])
        names = ak.to_list(event["beamline_id_name"])
        for bid, bname in zip(ids, names):
            if bid not in seen:
                seen[bid] = bname

    print(f"{'beamline_id':>15}  beamline_id_name")
    print("-" * 50)
    for bid in sorted(seen.keys()):
        print(f"{bid:>15}  {seen[bid]}")
