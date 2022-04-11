from random import sample
import numpy as np 
import msprime,math
import argparse
from recomb_map import RecombinationMap
from config import genome_data


class OutofAfricaDemography:
    generation_time = 25
    mutation_rate = 2.36e-8
    N0 = 7310  # initial population size
    Thum = 5920  # time (gens) of advent of modern humans
    Naf = 14474  # size of african population
    Tooa = 2040  # number of generations back to Out of Africa
    Nb = 1861  # size of out of Africa population
    mafb = 1.5e-4  # migration rate Africa and Out-of-Africa
    Teu = 920  # number generations back to Asia-Europe split
    Neu = 1032  # bottleneck population sizes
    Nas = 554
    mafeu = 2.5e-5  # mig. rates
    mafas = 7.8e-6
    meuas = 3.11e-5
    reu = 0.0038  # growth rate per generation in Europe
    ras = 0.0048  # growth rate per generation in Asia
    Tadmix = 12  # time of admixture
    Nadmix = 30000  # initial size of admixed population
    radmix = 0.05
    pop_labels = ['AFR','EUR','ASI','ADM']
    def __init__(self) -> None:
        self.demo = msprime.Demography()
        self.samples = {}
        self.demo.add_population(name =self.pop_labels[0],initial_size=self.Naf,growth_rate=0.0)
        # samples['AFR'] = 20
        self.demo.add_population(name=self.pop_labels[1],initial_size=self.Neu*math.exp(self.reu*self.Teu),growth_rate=self.reu)
        # samples['EUR'] = 20
        self.demo.add_population(name=self.pop_labels[2],initial_size=self.Nas*math.exp(self.ras*self.Teu),growth_rate=self.ras)
        # samples['ASI'] = 20
        self.demo.add_population(name=self.pop_labels[3],initial_size=self.Nadmix*math.exp(self.radmix*self.Tadmix),growth_rate=self.radmix)
        # samples['ADM'] = 20
        self.migration_matrix = [
                [0, self.mafeu, self.mafas, 0],
                [self.mafeu, 0, self.meuas, 0],
                [self.mafas, self.meuas, 0, 0],
                [0, 0, 0, 0],
            ]
        self.demo.add_mass_migration(time=self.Tadmix,source=self.pop_labels[3],dest=self.pop_labels[0],proportion=1.0/6.0)
        self.demo.add_mass_migration(time=self.Tadmix,source=self.pop_labels[3],dest=self.pop_labels[1],proportion=2.0/5.0)
        self.demo.add_mass_migration(time=self.Tadmix,source=self.pop_labels[3],dest=self.pop_labels[2],proportion=1.0)
        self.demo.add_migration_rate_change(time=self.Teu,rate=0.0)
        self.demo.add_mass_migration(time=self.Teu, source=self.pop_labels[2],dest=self.pop_labels[1],proportion=1.0)
        self.demo.add_population_parameters_change(time=self.Teu,initial_size=self.Nb,growth_rate=0.0,population=self.pop_labels[1])
        self.demo.add_symmetric_migration_rate_change(time=self.Teu,rate=self.mafb,populations=[self.pop_labels[0],self.pop_labels[1]])
        self.demo.add_migration_rate_change(time=self.Tooa,rate=0.0)
        self.demo.add_mass_migration(time=self.Tooa,source=self.pop_labels[1],dest=self.pop_labels[0],proportion=1.0)
        self.demo.add_population_parameters_change(time=self.Thum,initial_size=self.N0,population=self.pop_labels[0])
    def setup_demography(self,sample_size):
        if len(sample_size) == 1:
            sample_size = sample_size*4
        for index,pop in enumerate(self.pop_labels):
            self.samples[pop] = sample_size[index]
    def setup_recombination(self,chr_lengths,recomb_rates):
        self.recomb = RecombinationMap(chr_lengths,recomb_rates)
    def setup_model(self,dtwf_duration):
        self.model = [msprime.DiscreteTimeWrightFisher(duration=dtwf_duration),msprime.StandardCoalescent()]
    def simulate(self,mu=-1,random_seed=1234):
        if mu == -1:
            mu = self.mutation_rate
        self.random_seed = random_seed
        self.ts = msprime.sim_ancestry(samples=self.samples,demography=self.demo,model= self.model,random_seed=random_seed,recombination_rate=self.recomb.rate_map)
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
    
    parser.add_argument('--sample_size','-s',help='Either a single sample size for all of the populations'+
        ' or a list of sample sizes separated by spaces ordered as AFR-EUR-ASIA-ADMIX.',type=int,required=True,nargs='*')
    parser.add_argument('--output_dir','-o',dest='outdir',help='Output file prefix',type=str,required=True)
    parser.add_argument('--chr_length','-c',help='Either a single length for one chromosome simulation or a list of chromosome lengths separated by space.',
        type=int,default=[1e7],nargs='*')
    # parser.add_argument("--ne",'-n',help='Effective population size. Either a single number for all demes or a list of numbers separated by spaces.',
    #     type=int,nargs='*',default=[1e4])
    parser.add_argument('--rho','-r',help='Recombination rate for each chromosome. Either a list or a single value.',type=float,default=[1e-08],nargs='*')
    # parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08).",type=float,default=1e-08)
    
    
    parser.add_argument('--dtwf_duration','-d',help='Number of generations simulated using DTWF model.',type=int,default=50)
    
    
    parser.add_argument('--random_seed',help='Random seed for randomized parts of the algorithm (MSPRIME)',type=int,default=1234)
    # parser.add_argument('--no_tskit',help='also saves the tskit tree sequence file',dest='no_tskit',action='store_false',default=False)
    parser.add_argument('--tskit_save',help='How do you want the tskit data to be saved',dest='tskit_mode',type=str,default='single',choices=['single','multiple','none'])
    parser.add_argument('--make_vcf',help='save the vcf files',dest='make_vcf',action='store_true',default=False)
    parser.add_argument('--make_bed',help='save a single bed file',dest='make_bed',action='store_true',default=False)
    parser.add_argument('--bed_maf',help='minor allele frequency filtering for the bed file',default=0.0,type=float)
    
    args=parser.parse_args()

    print(args)
    simulator = OutofAfricaDemography()
    sample_size = args.sample_size
    if len(sample_size) > 1 and len(sample_size) != 4:
        raise ValueError('Mismatch in number of sample sizes and number of populations.')
    simulator.setup_demography(sample_size)
    
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
    simulator.simulate(random_seed=args.random_seed)
    simulator.write_to_file(args.outdir,args.tskit_mode)
    
    if args.make_vcf:
        simulator.write_vcf(args.outdir)
    if args.make_bed:
        simulator.write_bed(args.outdir,args.bed_maf)
    
if __name__ == '__main__':
    main()