import multiprocessing as mp
import pysam
import tqdm
import pandas as pd
import argparse

def process_cytoband_chunk(chunk, bam_fn, return_dict, process_id):
    """Process a chunk of cytoband_df and update coverage per cell"""
    cb_coverage = {}
    with pysam.AlignmentFile(bam_fn, "rb") as bam:
        for i, r in chunk.iterrows():
            for read in bam.fetch(r['chr'], r['pos1'], r['pos2']):
                CB = read.get_tag('CB')
                try:
                    cb_coverage[(i, CB)] += 1
                except KeyError:
                    cb_coverage[(i, CB)] = 1

    return_dict[process_id] = cb_coverage  # Store results in shared dict

def parallel_cytoband_processing(cytobands_df, bam_fn, num_workers=24):
    """Split cytobands_df into chunks and process in parallel"""
    manager = mp.Manager()
    return_dict = manager.dict()  # Shared dictionary for results
    processes = []

    chunk_size = np.ceil(cytobands_df['length'].sum() / num_workers)
    # chunk by total genomic coverage
    chunks = []
    for i in range(num_workers):
        chunk_df = cytobands_df[(cytobands_df['cumulative_length'] >= chunk_size * i) & (cytobands_df['cumulative_length'] < chunk_size * (i + 1))]
        chunks.append(chunk_df)
        
    print(f'Num chunks={len(chunks)}')
        
    assert np.sum([df.shape[0] for df in chunks]) == cytobands_df.shape[0]

    for i, chunk in enumerate(chunks):
        p = mp.Process(target=process_cytoband_chunk, args=(chunk, bam_fn, return_dict, i))
        processes.append(p)
        p.start()
    
    print('waiting')
    for p in processes:
        p.join()  # Wait for all processes to complete

    # Merge results from all processes
    final_cb_coverage = {k: v for d in return_dict.values() for k, v in d.items()}

    return final_cb_coverage
    
def main(
    bam_fn, cytoband_fn, num_workers, output_fn
):
    cytobands_df = pd.read_csv(cytobands_fn, sep='\t', compression='gzip', header=None)
    cytobands_df.columns = "chr pos1 pos2 cytoband name".split()
    cytobands_df['length'] = cytobands_df['pos2'] - cytobands_df['pos1']
    cytobands_df['cumulative_length'] = np.cumsum(cytobands_df['length'])
    
    cb_coverage = parallel_cytoband_processing(cytobands_df, bam_fn, num_workers=num_workers)

    cb_coverage_s = pd.Series(cb_coverage)
    cb_coverage_s.to_csv(output_fn, sep='\t')
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        prog='cytoband coverage',
        description='get coverage over cytobands per cell barcode',
    )
    
    parser.add_argument('--bam', help='bam path')
    parser.add_argument('--cytoband', help='cytoband path')
    parser.add_argument('--num_workers', help='number of workers')
    parser.add_argument('--output_fn', help='path to save counts')
    
    args = parser.parse_args()
    main(bam_fn=args.bam, cytoband_fn=args.cytoband, num_workers=args.num_workers, output_fn=args.output_fn)
    

