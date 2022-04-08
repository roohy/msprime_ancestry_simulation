import msprime,tskit,pickle
import make_bed
import numpy as np
class TreeIO():
    def __init__(self,ts=None,border_list=None,chr_count=0) -> None:
        self.ts = ts
        self.border_list = border_list
        self.chr_count = chr_count
        self.chrom_ts_list = None
    def read_from_file(self,ts_addr,border_list_addr):
        self._read_ts_from_file(ts_addr)
        self._read_border_list_file(border_list_addr)
        
    def _read_ts_from_file(self,ts_addr):
        self.ts = tskit.load(ts_addr)
    def _read_border_list_file(self,border_list_addr):
        self.border_list = pickle.load(open(border_list_addr,'rb'))
        self.chr_count = int(len(self.border_list)//2)
    def chr_divider(self,ts=None,border_list=None):
        if type(ts) != type(None):
            self.ts = ts
        if border_list != None:
            self.border_list = border_list
        assert type(self.ts) != type(None) and self.border_list != None
        self.chrom_ts_list = []
        for chr_num in range(int(len(self.border_list)//2)):
            start,end = self.border_list[chr_num*2:chr_num*2+2]
            chrom_ts = self.ts.keep_intervals([[start, end]], simplify=False).trim()
            self.chrom_ts_list.append(chrom_ts)
    def _write_single_file(self,output_addr):
        pickle.dump(self.border_list,open(output_addr+'.bls.pkl','wb'))
        self.ts.dump(output_addr+'.ts')
    def _load_single_file(self, addr): #TODO: this function should not be used. It is only here for the sake of compatibility 
        ts = tskit.load(addr+'.ts')
        self.border_list = pickle.load(open(addr+'.bls.pkl','rb'))
        self.temp_ts = ts
        return ts
    def write_vcf(self,output_prefix):
        if self.chrom_ts_list is None:
            print('No chromosome based Tree Sequence available! So we are making them')
            self.chr_divider()
        n_dip_indv = int(self.chrom_ts_list[0].num_samples / 2)
        indv_names = [f"id_{str(i)}" for i in range(1,n_dip_indv+1)]
        for chr_num in range(self.chr_count):
            with open(f'{output_prefix}_chr{chr_num+1}.vcf', "w") as vcf_file:
                self.chrom_ts_list[chr_num].write_vcf(vcf_file, individual_names=indv_names,contig_id=chr_num+1)
    def write_bed(self,output_prefix,maf=0):
        if self.chrom_ts_list is None:
            print('No chromosome based Tree Sequence available! So we are making them')
            self.chr_divider()
        n_dip_indv = int(self.chrom_ts_list[0].num_samples / 2)
        indv_names = [f"id_{str(i)}" for i in range(1,n_dip_indv+1)]
        bed_writer = make_bed.BedWriter(self.chrom_ts_list[0],individual_names=indv_names,contig_id=1)
        with open(output_prefix+'.fam','w') as fam_output:
            with open(output_prefix+'.bim','w') as bim_output:
                with open(output_prefix+'.bed','wb') as bed_output:
                    for chr_num in range(self.chr_count):
                        bed_writer.contig_id = chr_num+1
                        bed_writer.tree_sequence = self.chrom_ts_list[chr_num]
                        bed_writer.write(bed_output,bim_output,fam_output,chr_num != 0,maf)

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
        self.tio = None
    
    def tio_init(self,ts):
        self.tio = TreeIO(ts=ts,border_list=self.border_list,chr_count=self.chr_count)
    def chr_divider(self, ts):
        if self.tio is None:
            self.tio_init(ts)
        self.tio.chr_divider()
        # self.chrom_ts_list = [] TODO:remove these lines after testing
        # for chr_num in range(self.chr_count):
        #     start,end = self.border_list[chr_num*2:chr_num*2+2]
        #     chrom_ts = ts.keep_intervals([[start, end]], simplify=False).trim()
        #     self.chrom_ts_list.append(chrom_ts)
    def _write_single_file(self,output_addr,ts=None):
        if self.tio is None:
            assert ts is not None
            self.tio_init(ts)
        self.tio._write_single_file(output_addr)
        
    
    def write_to_file(self,output_prefix,single_file=None):
        if single_file is not None:
            self._write_single_file(ts=single_file,output_addr=output_prefix)
        else:
            for chr_num in range(self.chr_count):
                self.chrom_ts_list[chr_num].dump(f'{output_prefix}_chr{chr_num+1}.ts')
    def write_vcf(self,output_prefix):
        assert self.tio is not None
        self.tio.write_vcf(output_prefix)
        
    def write_bed(self,output_prefix,maf=0):
        assert self.tio is not None
        self.tio.write_bed(output_prefix,maf)

