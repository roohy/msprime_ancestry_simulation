from secrets import choice
import msprime,argparse,pickle,tskit,os
import numpy as np
from functools import partial
from config import genome_data
import make_bed
from recomb_map import RecombinationMap



class CellDemography:
    
    def __init__(self,height,width,migration_rate) -> None:
        self.pop = msprime.Demography()
        self.height = height
        self.width = width
        self.migration_rate = migration_rate
        self.deme_count = height*width
        self.ne = None
        self.ss = None
        self.samples = {}
    @staticmethod
    def name_generator(row, col):
        return f'd{row+1}_{col+1}'
    def migration_options(self,row,col):
        return zip([row, row, row-1, row+1],[col-1, col+1, col, col])
    def check_migration_path(self,d1,d2):
        if d2[0] >= 0 and d2[0] < self.height:
            if d2[1] >= 0 and d2[1] < self.width:
                return True
                #print(f'd{d1[0]+1}_{d1[1]+1} - {d2[0]+1}_{d2[1]+1}')
                mig_function([name_generator(d1[0],d1[1]),name_generator(d2[0],d2[1])],migration_rate)
        return False
    def add_migration_path(self,d1,d2):
        if self.check_migration_path(d1,d2):
            self.pop.set_symmetric_migration_rate([self.name_generator(d1[0],d1[1]),self.name_generator(d2[0],d2[1])],self.migration_rate)
            return True
        return False
    def setup_migration(self):
        for i in range(self.height):
            for j in range(self.width):
                
                #These two lines either also check for the options to see if it is a legit path or/and generate a name which is not ideal  checking should be done before name generation or independent of it (durrently the second option is used)
                #d1 = CellDemography.name_generator(i,j)
                # migration_options = [CellDemography.name_generator(*item) for item in self.migration_options(i,j) if self.check_migration_path(d1,item)]
                # migration_options = list(map(lambda x: CellDemography.name_generator(*x), self.migration_options(i,j))) 
                d1 = (i,j)
                res = list(map(lambda d2: self.add_migration_path(d1,d2),self.migration_options(i,j) ))
                print(res)
    def set_size(self,val,key):
        if len(val) > 1:
            if len(val) != self.deme_count:
                raise ValueError(f'Mismatch between the number of sizes passed and number of available demes for {key}!')
            self.__dict__[key] = np.array(val).reshape((self.height,self.width))
            self.__dict__['multi_'+key] = True
        elif len(val) == 1 :
            self.__dict__[key] = val[0]
            self.__dict__['multi_'+key] = False
    def set_effective_size(self,ne):
        self.set_size(ne,'ne')
      
    def set_sample_size(self,ss):
        self.set_size(ss,'ss')
    def get_size(self,row,col,key):
        if self.__dict__['multi_'+key]:
            return self.__dict__[key][row,col]
        return self.__dict__[key]
    def get_sample_size(self,row,col):
        return self.get_size(row,col,'ss')
    def get_effective_size(self,row,col):
        return self.get_size(row,col,'ne')

    def setup_population(self,ne=None,ss=None):
        if self.ne is None:
            if ne is None:
                raise Exception('Effective size is not setup')
            else:
                self.set_effective_size(ne)
        if self.ss is None:
            if ss is None:
                raise Exception('Sample size is not setup')
            else:
                self.set_sample_size(ss)

        for row in range(self.height):
            for col in range(self.width):
                name = CellDemography.name_generator(row,col)
                self.pop.add_population(name=name,initial_size=self.get_effective_size(row,col))
                self.samples[name] = self.get_sample_size(row,col)
        
    def add_ancestral_pop(self,time_to_merge=200,ancestral_size=None):
        if ancestral_size == None or ancestral_size < 1:
            if self.multi_ne:
                raise ValueError('Ambiguity in the population effective size for ancestral population!')
            ancestral_size = self.ne
        self.pop.add_population(name='pan',initial_size=ancestral_size,
                                description='Ancestral population!')
        self.pop.add_population_split(time=time_to_merge,ancestral='pan',
                                      derived=np.arange(self.height*self.width))
    def add_bottleneck(self,row,col,bn_size,start_time,bn_duration):
        deme_name = CellDemography.name_generator(row-1,col-1)
        original_size = self.get_effective_size(row-1,col-1)
        print(f'Adding bottleneck for population {deme_name}. This population will'+
              f'experience a bottleneck at generation {start_time} that lasts for {bn_duration} generations.'+
              f'During this time, the effective population size for this population will be {bn_size},'+
              f'before going back to the original size of {original_size}.')
        self.pop.add_population_parameters_change(start_time,initial_size=bn_size,population=deme_name)
        self.pop.add_population_parameters_change(start_time+bn_duration,initial_size=original_size,population=deme_name)
    
    @staticmethod
    def bottleneck_parser(self,input_addr):
        bn_list = []
        with open(input_addr,'r') as bn_file:
            for line in bn_file:
                data = line.strip().split()
                bn_list.append([int(item) for item in data])
        return bn_list
    def add_bottlenecks_from_file(self,input_addr):
        bn_list = CellDemography.bottleneck_parser()
        for bn_item in bn_list:
            self.add_bottleneck(*bn_item)
        self.pop.sort_events()
                


class DownwardDemography(CellDemography):
    def __init__(self, height, width, migration_rate) -> None:
        super().__init__(height, width, migration_rate)
    def migration_options(self, row, col):
        return ([row+1,col])
    def add_migration_path(self, d1, d2):
        if self.check_migration_path(d1,d2):
            self.pop.set_migration_rate(source=d2,destination=d1,rate=self.migration_rate)
            return True
        return False
    
        
    
class GridSimulation():
    def __init__(self) -> None:
        pass
    def setup_demography(self,height,width,migration_rate,migration_dir,sample_size,effective_size,ancestral_size,time_to_merge):
        if migration_dir == 'all':
            self.demo = CellDemography(height,width,migration_rate)
        elif migration_dir == 'down':
            self.demo = DownwardDemography(height,width,migration_rate)
        else:
            raise ValueError(f'Migration direction "{migration_dir}" is not supported!')
        self.demo.set_effective_size(effective_size)
        self.demo.set_sample_size(sample_size)
        self.demo.setup_population()
        self.demo.setup_migration()
        self.demo.add_ancestral_pop(time_to_merge,ancestral_size)
    def setup_recombination(self,chr_lengths,recomb_rates):
        self.recomb = RecombinationMap(chr_lengths,recomb_rates)
    def setup_model(self,dtwf_duration):
        self.model = [msprime.DiscreteTimeWrightFisher(duration=dtwf_duration),msprime.StandardCoalescent()]
    def setup_bottlenecks(self,bn_file):
        self.demo.add_bottlenecks_from_file(bn_file)
    def simulate(self,mu,random_seed=1234):
        self.random_seed = random_seed
        self.ts = msprime.sim_ancestry(samples=self.demo.samples,demography=self.demo.pop,model= self.model,random_seed=random_seed,recombination_rate=self.recomb.rate_map)
        self.mts = msprime.sim_mutations(self.ts,rate=mu,random_seed=random_seed)
        self.recomb.chr_divider(self.mts)

    def write_to_file(self,output_prefix,mode='single'): 
        if mode == 'single':
            self.recomb.write_to_file(output_prefix,self.mts)
        elif mode == 'multiple':
            self.recomb.write_to_file(output_prefix)
        elif mode == 'none':
            print('Tree Sequence data discarded')
    def write_vcf(self,output_prefix):
        self.recomb.write_vcf(output_prefix)
    def write_bed(self,output_prefix,maf=0):
        self.recomb.write_bed(output_prefix,maf)
    



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
    parser.add_argument('--time_to_merge','-t',help='Time (in generations) to panmixia. nonpositive numbers will be treated as inifinity',type=int,default=150)
    parser.add_argument('--ancestral_size','-a',help='Effective size of the ancestral population. Defauls is set to the effective population size of the first deme.',type=int,default=-1)
    parser.add_argument('--random_seed',help='Random seed for randomized parts of the algorithm (MSPRIME)',type=int,default=1234)
    # parser.add_argument('--no_tskit',help='also saves the tskit tree sequence file',dest='no_tskit',action='store_false',default=False)
    parser.add_argument('--tskit_save',help='How do you want the tskit data to be saved',dest='tskit_mode',type=str,default='single',choices=['single','multiple','none'])
    parser.add_argument('--make_vcf',help='save the vcf files',dest='make_vcf',action='store_true',default=False)
    parser.add_argument('--make_bed',help='save a single bed file',dest='make_bed',action='store_true',default=False)
    parser.add_argument('--bed_maf',help='minor allele frequency filtering for the bed file',default=0,type=float)
    parser.add_argument('--bottlenecks','-b', dest='bnfile',help='File with bottleneck list', type=str,default='')
    args=parser.parse_args()

    print(args)
    simulator = GridSimulation()
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
    simulator.setup_demography(height,width,migration_rate,migration_dir,sample_size,effective_size,ancestral_size,args.time_to_merge)
    if args.bnfile:
        simulator.setup_bottlenecks(args.bnfile)
    simulator.setup_model(args.dtwf_duration)
    if args.chr_length[0] == -1:
        lengths = [genome_data[key]['length'] for key in range(1,23) ]
        rates = [genome_data[key]['rate'] for key in range(1,23) ]
        simulator.setup_recombination(lengths,rates)
    elif args.rho[0] == -1:
        rates = [genome_data[key]['rate'] for key in range(1,23) ]
        simulator.setup_recombination(args.chr_length,rates[:len(args.chr_length)])
    else:
        simulator.setup_recombination(args.chr_length,args.rho)
    simulator.simulate(args.mu,args.random_seed)
    dirname = os.path.dirname(args.outdir)
    print('Simulation finished.')
    print(f'Checking the output directory at {dirname}...')
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    print('Generating the TS file...')
    simulator.write_to_file(args.outdir,args.tskit_mode)
    print('Sequence file finished.')
    print('Generating the genotype files...')
    if args.make_vcf:
        simulator.write_vcf(args.outdir)
    if args.make_bed:
        print()
        simulator.write_bed(args.outdir,args.bed_maf)
        print('Bed file finished')
    
if __name__ == '__main__':
    main()