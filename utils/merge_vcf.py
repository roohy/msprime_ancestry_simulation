from email.mime import base
import sys,os
import numpy as np

class BedFileHandler:
    def __init__(self,bed_addr,bim_addr) -> None:
        self.bed_addr = bed_addr
        self.bim_addr = bim_addr
    def __enter__(self):
        self.bed_file = open(self.bed_addr, "wb")
        self.bim_file = open(self.bim_addr,'w')
        
        self.bed_file.write(b'\x6c\x1b\x01')

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.bed_file.close()
        self.bim_file.close()
        
class VCFConverter:
    def __init__(self,base_name,output_addr) -> None:
        self.base = base_name
        
        self.bed_addr = output_addr+'.bed'
        self.bim_addr = output_addr+'.bim'
        self.fam_addr = output_addr+'.fam'
    def load_fam(self):
        ids = None
        with open(f'{self.base}1.vcf','r') as chr1vcf:
            for line in chr1vcf:
                if line.startswith('#CHROM'): 
                    data = line.strip().split()
                    ids = data[9:]
                    break
        with open(self.fam_addr,'w') as fam_file:
            for item in ids: 
                fam_file.write(f'{item} {item} 0 0 0 -9\n')
        self.sample_count = len(ids)
        
        
    def merge(self):
        geno_array = np.zeros((2*self.sample_count),dtype=np.uint8)
        with BedFileHandler(self.bed_addr,self.bim_addr) as bed_handler:
            for i in range(1,23):
                vcf_addr = f'{self.base}{i}.vcf'
                VCFConverter.read_vcf(vcf_addr,bed_handler,geno_array)
    
    @staticmethod
    def read_vcf(vcf_addr,bed_handler,geno_array):
        with open(vcf_addr) as vcf_file:
            for line in vcf_file:
                if line.startswith('#CHROM'): 
                    break
            for line in vcf_file:
                data = line.strip().split()
                if ',' in data[3]: 
                    data[3] = data[3].replace(',','')
                if ',' in data[4]:
                    data[4] = data[4].replace(',','')
                bed_handler.bim_file.write(f'{data[0]} {data[2]} 0 {data[1]} {data[3]} {data[4]}\n')
                data = data[9:]
                index = 0 
                for item in data:
                    if item[0] == '1':
                        geno_array[2*index+1] = 1
                        if item[2] == '1':
                            geno_array[2*index] = 1 
                    elif item[2] == '1':
                        geno_array[2*index+1] = 1 
                    index += 1 
                packed_bits = np.packbits(geno_array,bitorder='little')
                bed_handler.bed_file.write(packed_bits)
                del packed_bits
                geno_array[:] = 0 



            
# def list_files(var_dir):         
#     var_files = [(os.path.join(var_dir,f),int(f[:-3])) for f in os.listdir(var_dir) if f.startswith('base') and f.endswith('.vcf')]
#     var_files.sort(key=lambda x:x[1])
#     return var_files


def vcf_to_bed(vcf_prefix, output_addr):
    converter = VCFConverter(vcf_prefix,output_addr)
    converter.load_fam()
    converter.merge()

def main():
    vcf_pre = sys.argv[1]
    output_addr = sys.argv[2]
    
    vcf_to_bed(vcf_pre,output_addr)

if __name__ == '__main__':
    main()