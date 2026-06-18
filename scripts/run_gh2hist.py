import argparse
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Axi_Schwarzschild'))
from data_prep.hist_losvd import convert_gh_to_hist


def main():
    parser = argparse.ArgumentParser(
        description="Convert Gauss-Hermite kinematics to binned LOSVD histograms")
    parser.add_argument("gh_file",
                        help="Input ECSV file with GH kinematics (e.g. kinematics/gauss_hermite_kins_o.ecsv)")
    parser.add_argument("--nvbins", type=int, default=21,
                        help="Number of velocity bins in histogram (default: 21)")
    parser.add_argument("--vmax", type=float, default=1000,
                        help="Maximum velocity for histogram range [-vmax, vmax] (default: 1000)")
    parser.add_argument("--nmc", type=int, default=2000,
                        help="Number of Monte Carlo realizations for errors (default: 2000)")
    parser.add_argument("--min-err", type=float, default=1e-3,
                        help="Minimum LOSVD error to avoid zeros (default: 1e-3)")
    parser.add_argument("--output", default=None,
                        help="Output ECSV file (default: kinematics/hist_kinematics.ecsv)")
    args = parser.parse_args()

    output = args.output or "kinematics/hist_kinematics.ecsv"

    print(f"=== GH-to-LOSVD histogram conversion ===")
    print(f"  Input:  {args.gh_file}")
    print(f"  Output: {output}")
    print(f"  nvbins: {args.nvbins}, vmax: {args.vmax}, nmc: {args.nmc}")

    hist_table = convert_gh_to_hist(
        args.gh_file, output, args.nvbins, args.vmax, args.nmc, args.min_err)

    print(f"\n=== Done! {len(hist_table)} bins written to {output} ===")


if __name__ == "__main__":
    main()
