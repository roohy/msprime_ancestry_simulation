import tskit
import argparse,gzip

def main():
    parser=argparse.ArgumentParser()
    
    parser.add_argument('--tree-seq','-t',help='address to the tskit tree sequence file',type=str,required=True,dest='tree_seq')
    parser.add_argument('--out','-o',help='Address of the output match file',type=str,required=True)
    parser.add_argument('--recomb-rate','r',help='Recombination rate (the script currently only supports one chromosome sims)',type=float,default=1.1485597641285933e-08,dest='recomb_rate')
    parser.add_argument('--min-length','-m',help='Minimum Length of IBD segments to be considered',type=float,default=3.0,dest='min_length')
    parser.add_argument('--min-agg-length','-l',help='Minimum aggregate IBD shared between pairs to be included in the final network',type=float,default=6.0,dest='min_agg_length')
    # parser.add_argument('')
    # parser.add_argument('--')


    args = parser.parse_args()

    ts = tskit.load(args.tree_seq)
    print(f"The file 'data/basics.trees' exists, and contains {ts.num_trees} trees")
    
    rate = args.recomb_rate

    minimum_bp_length = args.min_length/(100*rate)

    segments = ts.ibd_segments(store_segments=True,min_span=minimum_bp_length,max_time=200)
    with gzip.open(args.out,'wt') as output_file:
        for pair,segment_list in segments.items():
            total_length = 0
            for segment in segment_list:
                segment_length = (segment.right-segment.left)*rate*100
                total_length += segment_length
            if total_length > args.min_agg_length:
                print(f'{pair[0]} {pair[1]} {total_length}',file=output_file)
        # print((segment.right-segment.left)*rate*100)
    

if __name__ == '__main__':
    main()