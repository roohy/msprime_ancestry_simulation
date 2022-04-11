import numpy as np
import sys,gzip
from config import genome_data
header_line = 'chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n'
def main():
    
    rates = [genome_data[key]['rate'] for key in range(1,23) ]
    combined_rates = [item*1_000_000*100 for item in rates]
    input_addr = sys.argv[1]
    output_addr = sys.argv[2]
    with open(input_addr,'r') as input_bim:
        with gzip.open(output_addr,'wb') as output_map:
            output_map.write(header_line.encode())
            for line in input_bim:
                data = line.strip().split()
                chrom = int(data[0])
                position = int(data[3])
                output_map.write(f'{chrom} {position} {combined_rates[chrom-1]} {position*rates[chrom-1]*100}\n'.encode())


if __name__ == '__main__':
    main()