import sys
import argparse
import logging
import warnings
import pandas as pd
import cellphoneDB_utils as cpu
warnings.filterwarnings("ignore")

logger = logging.getLogger("build_contact_map")
logger.setLevel(logging.INFO)
handler=logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s:%(step)s - %(levelname)s - %(message)s"))
logger.addHandler(handler)

def build_contact_map(sc_coordinate,out_path, divergence_cutoff, band_width=20, n_bootstrap=20,n_threads=30):

    cell_type_dummy_df = cpu.process_metadata(sc_coordinate)
    coord = sc_coordinate[['x_noise','y_noise']]
    dis_boot_array, dis_cons, mst_cons = cpu.KL_JS_boot_mst(dummy_df=cell_type_dummy_df,
                                                            coord_df=coord, h=band_width,
                                                            boot_n=n_bootstrap,
                                                            n_threads=n_threads)
    cpu.network_microenv(mst_cons, out_path=out_path, cutoff=divergence_cutoff)
    cpu.divergence_clustermap(dis_cons, out_path=out_path)
    cpu.network_draw(mst_cons, out_path=out_path)

def main(argsv):
    parser = argparse.ArgumentParser(description="Build cell-cell contact map and detect microenvironment.")
    parser.add_argument("--sc_coordinate", type=str, help="Path to sc_coordinate.csv")
    parser.add_argument('--divergence_cutoff', type=float, default=0.5, help='Filter out divergence lower than cutoff percentile, ranging from 0 to 1. [default: 0.5]')
    parser.add_argument("--band_width", type=int, default=20, help='Bandwidths for x and y directions, more details in KernelDensity function in sklearn.neighbors package. [default: 20]')
    parser.add_argument("--n_bootstrap", type=int, default=20, help='Number of bootstrapping iterations. [default: 20]')
    parser.add_argument("--n_threads", type=int, default=30, help='Number of threads available to bootstrapping [default: 30]')
    parser.add_argument("--outDir", type=str, help="Output directory")

    args = parser.parse_args()

    #### Read data
    logger.info("Reading data...", extra={'step': 'Main'})
    sc_coordinate = pd.read_csv(args.sc_coordinate,header=0,index_col=False)

    logger.info("Building contact map...", extra={'step': 'Main'})
    build_contact_map(sc_coordinate=sc_coordinate,
                      out_path=args.outDir,
                      divergence_cutoff=args.divergence_cutoff,
                      band_width=args.band_width,
                      n_bootstrap=args.n_bootstrap,
                      n_threads=args.n_threads)

    logger.info("All Finished!", extra={'step': 'Main'})

if __name__ == '__main__':
    import sys

    main(sys.argv)