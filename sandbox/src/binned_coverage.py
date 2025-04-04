import multiprocessing as mp
import pysam
import tqdm
import pandas as pd
import argparse
import numpy as np

def count_region(regions_df, bam_fn, return_dict, process_id):
    """Process a chunk of cytoband_df and update coverage per cell"""
    cb_coverage = {}
    with pysam.AlignmentFile(bam_fn, "rb") as bam:
        for i, r in regions_df.iterrows():
            for read in bam.fetch(r['chr'], r['start'], r['end']):
                CB = read.get_tag('CB')
                try:
                    cb_coverage[(i, CB)] += 1
                except KeyError:
                    cb_coverage[(i, CB)] = 1

    return_dict[process_id] = cb_coverage  # Store results in shared dict

def parallel_count_region(regions_df, bam_fn, num_workers=24):

    manager = mp.Manager()
    return_dict = manager.dict()  # Shared dictionary for results
    processes = []

    num_chunks = int(np.ceil(regions_df.shape[0] / num_workers))
        
    print(f'Num chunks={num_chunks}')
    
    groups = regions_df.groupby(np.arange(len(regions_df.index))//num_chunks)
    print(f'Num groups={groups}')

    for i, group in groups:
        p = mp.Process(target=count_region, args=(group, bam_fn, return_dict, i))
        processes.append(p)
        print(f'Launching {i}')
        p.start()
    
    print('waiting')
    for p in processes:
        p.join()  # Wait for all processes to complete

    # Merge results from all processes
    final_cb_coverage = {k: v for d in return_dict.values() for k, v in d.items()}

    return final_cb_coverage
    
def main(
    bam_fn, output_fn, bin_size=1_000_000, num_workers=4
):
    with pysam.AlignmentFile(bam_fn, "rb") as bam:
        bam_header = bam.header.to_dict()
        
        bed_dfs = []
        for contig in bam_header['SQ']:
            contig_name = contig['SN']
            contig_len = contig['LN']

            bins = np.arange(0, contig_len + bin_size, bin_size)
            bed_df = pd.DataFrame({'chr': [contig_name] * (len(bins) - 1), 'start': bins[:-1], 'end': bins[1:]})
            bed_dfs.append(bed_df)
        all_bed_dfs = pd.concat(bed_dfs).reset_index(drop=True)
        
    all_beds_fn = f'{output_fn.rsplit(".", 1)[0]}.all_beds.tsv'
    print(f"Num regions: {all_bed_dfs.shape}. Saved to {all_beds_fn}")
    
    all_bed_dfs.to_csv(all_beds_fn, sep='\t')
    
    cb_coverage = parallel_count_region(all_bed_dfs.head(20), bam_fn, num_workers=num_workers)

    cb_coverage_s = pd.Series(cb_coverage)
    cb_coverage_s.to_csv(output_fn, sep='\t')
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog='cytoband coverage',
        description='get coverage over cytobands per cell barcode',
    )
    
    parser.add_argument('--bam', help='bam path')
    parser.add_argument('--bin_size', type=int, help='bin_size to extract from the bam header')
    parser.add_argument('--num_workers', type=int, help='number of workers')
    parser.add_argument('--output_fn', help='path to save counts')
    
    args = parser.parse_args()
    main(bam_fn=args.bam, bin_size=args.bin_size, num_workers=args.num_workers, output_fn=args.output_fn)
    

