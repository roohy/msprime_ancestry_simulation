import msprime,argparse
import numpy as np
from functools import partial
from config import genome_data


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
    
        
class RecombinationMap():
    def __init__(self,chr_lengths,recomb_rates) -> None:
        self.chr_lengths = chr_lengths
        self.chr_count = len(chr_lengths)
        self.recomb_rates  = recomb_rates
        if self.chr_count != len(self.recomb_rates) and len(self.recomb_rates) > 1:
            raise ValueError('Mismatch in values provided for recombination rate on each chromosome and length of chromosomes')
        self.border_list = []
        head = 0
        for index,chr_len in enumerate(self.chr_lengths):
            self.border_list.append(head)
            self.border_list.append(head+chr_len)
            head += chr_len + 1 
        self.recomb_list = None
        if len(self.recomb_rates) == 1 :
            recomb_list = (self.chr_count*[self.recomb_rates[0],0.5])[:-1]
        else:
            recomb_list = np.zeros((2*self.chr_count-1))
            recomb_list[1::2] = .5
            recomb_list[::2] = self.recomb_rates
    
        self.rate_map = msprime.RateMap(position=self.border_list,rate=recomb_list)
    def chr_divider(self, ts):
        self.chrom_ts_list = []
        for chr_num in range(self.chr_count):
            start,end = self.border_list[chr_num*2:chr_num*2+2]
            chrom_ts = ts.keep_intervals([[start, end]], simplify=False).trim()
            self.chrom_ts_list.append(chrom_ts)
    def write_to_file(self,output_prefix):
        for chr_num in range(self.chr_count):
            self.chrom_ts_list[chr_num].dump(f'{output_prefix}_chr{chr_num+1}.ts')
    def write_vcf(self,output_prefix):
        n_dip_indv = int(self.chrom_ts_list[0].num_samples / 2)
        indv_names = [f"id_{str(i)}" for i in range(1,n_dip_indv+1)]
        for chr_num in range(self.chr_count):
            with open(f'{output_prefix}_chr{chr_num+1}.vcf', "w") as vcf_file:
                self.chrom_ts_list[chr_num].write_vcf(vcf_file, individual_names=indv_names,contig_id=chr_num+1)
            
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
    def simulate(self,mu,random_seed=1234):
        self.random_seed = random_seed
        self.ts = msprime.sim_ancestry(samples=self.demo.samples,demography=self.demo.pop,model= self.model,random_seed=random_seed,recombination_rate=self.recomb.rate_map)
        self.mts = msprime.sim_mutations(self.ts,rate=mu,random_seed=random_seed)
        self.recomb.chr_divider(self.mts)

    def write_to_file(self,output_prefix): 
        self.recomb.write_to_file(output_prefix)
    def write_vcf(self,output_prefix):
        self.recomb.write_vcf(output_prefix)




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
    parser.add_argument('--no_tskit',help='also saves the tskit tree sequence file',dest='no_tskit',action='store_false',default=False)
    parser.add_argument('--no_vcf',help='save the vcf file',dest='no_vcf',action='store_false',default=False)
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
    simulator.setup_model(args.dtwf_duration)
    if args.chr_length[0] == -1:
        lengths = [genome_data[key]['length'] for key in genome_data ]
        rates = [genome_data[key]['rate'] for key in genome_data ]
        simulator.setup_recombination(lengths,rates)
    elif args.rho[0] == -1:
        rates = [genome_data[key]['rate'] for key in genome_data ]
        simulator.setup_recombination(args.chr_length,rates[:len(args.chr_length)])
    else:
        simulator.setup_recombination(args.chr_length,args.rho)
    simulator.simulate(args.mu,args.random_seed)
    if not args.no_tskit:
        simulator.write_to_file(args.outdir)
    if not args.no_vcf:
        simulator.write_vcf(args.outdir)
    
if __name__ == '__main__':
    main()