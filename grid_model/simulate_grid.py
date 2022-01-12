import msprime,argparse
import numpy as np
from functools import partial



def check_set_symmetric_migration_change(d1, d2, migration_rate, demography, height, width):
    if d2[0] >= 0 and d2[0] < height:
        if d2[1] >= 0 and d2[1] < width:
            #print(f'd{d1[0]+1}_{d1[1]+1} - {d2[0]+1}_{d2[1]+1}')
            demography.set_symmetric_migration_rate([name_generator(d1[0],d1[1]),name_generator(d2[0],d2[1])],migration_rate)
            return True
    return False

def name_generator(row,col):
    return f'd{row+1}_{col+1}'
def check_multi_value(flag,value,row,col):
    if flag:
        return value[row,col]
    else:
        return value


def main():
    parser=argparse.ArgumentParser()
    
    parser.add_argument('--sample_size','-s',help='Either a single sample size for all of the demes or a list of sample sizes separated by spaces.',type=int,required=True,nargs='*')
    parser.add_argument('--output_dir','-o',dest='outdir',help='Output file prefix',type=str,required=True)
    parser.add_argument('--chr_length','-c',help='Either a single length for one chromosome simulation or a list of chromosome lengths separated by space.',
        type=int,default=[1e7],nargs='*')
    parser.add_argument("--ne",'-n',help='Effective population size. Either a single number for all demes or a list of numbers separated by spaces.',
        type=int,nargs='*',default=[1e4])
    parser.add_argument('--rho','-r',help='Recombination rate for each chromosome. Either a list or a single value.',type=float,default=[1e-08],nargs='*')
    parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08).",type=float,default=1e-08)
    parser.add_argument('--migration_rate','-m',help='Migration rate among the demes',type=float,default=0.05)
    parser.add_argument('--migration_dir',help='direction of migration, it can be all possible paths (all) or downward paths only (down).',
        type=str,default='all',choices=['all','down'])
    parser.add_argument('--deme_rows','-x',help='How many rows of demes to be simulated',type=int,default=3)
    parser.add_argument('--deme_columns','-y',help='How many columns of demes to be simulated',type=int,default=3)
    parser.add_argument('--dtwf_duration','-d',help='Number of generations simulated using DTWF model.',type=int,default=50)
    parser.add_argument('--time_to_merge','-t',help='Time (in generations) to panmixia. nonpositive numbers will be treated as inifinity',type=int,default=100)
    parser.add_argument('--ancestral_size','-a',help='Effective size of the ancestral population. Defauls is set to the effective population size of the first deme.',type=int,default=-1)
    parser.add_argument('--random_seed',help='Random seed for randomized parts of the algorithm (MSPRIME)',type=int,default=1234)
    args=parser.parse_args()

    print(args)
    demography = msprime.Demography()
    width = args.deme_columns
    height = args.deme_rows
    deme_count = width * height
    sample_size = args.sample_size
    effective_size = args.ne

    ancestral_size = args.ancestral_size if args.ancestral_size > 0 else sample_size[0]

    migration_dir = args.migration_dir
    if len(effective_size) > 1 and ancestral_size < 1:
        raise ValueError("Ambiguity in ancestral population size since multiple population sizes are available.")
    if len(sample_size) > 1 and len(sample_size) != deme_count:
        raise ValueError("Mismatch in number of demes and available initial samples sizes.")
    if len(sample_size) != 1 and len(effective_size) != 1 and len(sample_size) != len(effective_size):
        raise ValueError('Discrepancy between the number of initial sample sizes and effect population sizes passed.')
    migration_rate = args.migration_rate
    multiple_ss = len(sample_size) > 1
    multiple_ne = len(effective_size) > 1
    if multiple_ss:
        sample_size = np.array(sample_size).reshape((height, width))
    else:
        sample_size = sample_size[0]
    if multiple_ne:
        effective_size = np.array(effective_size).reshape((height,width))
    else:
        effective_size = effective_size[0]
    
    migration_partial = partial(check_set_symmetric_migration_change, migration_rate=migration_rate, demography=demography, height=height, width=width)
    samples = {}
    for row in range(height):
        for col in range(width):
            name = name_generator(row,col)
            temp_ss = check_multi_value(multiple_ss, sample_size, row, col)
            temp_ne = check_multi_value(multiple_ne, effective_size, row, col)
            demography.add_population(name=name,initial_size=temp_ne)
            samples[name] = temp_ss
    for row in range(height):
        for col in range(width):
            if migration_dir == 'all':
                coordinates = zip([row, row, row-1, row+1],[col-1, col+1, col, col])
            elif migration_dir == 'down':
                coordinates = [(row+1,col)]
            d1 = (row, col)
            print(list(map(lambda d2: migration_partial(d1,d2),coordinates)))
    demography.add_population(name='pan',initial_size=ancestral_size,description='Ancestral population!')
    demography.add_population_split(time=args.time_to_merge,ancestral='pan',derived=np.arange(height*width))
    model = [msprime.DiscreteTimeWrightFisher(duration=args.dtwf_duration),msprime.StandardCoalescent()]

    chr_lengths = args.chr_length
    chr_count = len(chr_lengths)
    recombination_rates = args.rho
    if chr_count != len(recombination_rates) and len(recombination_rates) > 1:
        raise ValueError('Mismatch in values provided for recombination rate on each chromosome and length of chromosomes')
    
    
    border_list = []
    head = 0
    for index,chr_len in enumerate(chr_lengths):
        border_list.append(head)
        border_list.append(head+chr_len)
        head += chr_len + 1 
    recomb_list = None
    if len(recombination_rates) == 1 :
        recomb_list = (chr_count*[recombination_rates[0],0.5])[:-1]
    else:
        recomb_list = np.zeros((2*chr_count-1))
        recomb_list[1::2] = .5
        recomb_list[::2] = recombination_rates
    
    rate_map = msprime.RateMap(position=border_list,rate=recomb_list)
    

    mu = args.mu
    ts = msprime.sim_ancestry(samples=samples,demography=demography,model= model,random_seed=args.random_seed,recombination_rate=rate_map)
    mts = msprime.sim_mutations(ts,rate=mu,random_seed=args.random_seed)
    print(border_list)
    for chr_num in range(chr_count):
        start,end = border_list[chr_num*2:chr_num*2+2]
        chrom_ts = mts.keep_intervals([[start, end]], simplify=False).trim()
        chrom_ts.dump(f'{args.outdir}_chr{chr_num+1}.ts')
        

    
if __name__ == '__main__':
    main()