import gzip
import argparse
import pandas as pd

def parse_vcf(vcf_file):
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'rt') as f:
        samples = []
        genotype_counts = {}

        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                samples = header[9:]
                for sample in samples:
                    genotype_counts[sample] = {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0}
            else:
                columns = line.strip().split('\t')
                genotypes = columns[9:]

                for i, genotype_info in enumerate(genotypes):
                    genotype = genotype_info.split(':')[0]
                    if genotype in genotype_counts[samples[i]]:
                        genotype_counts[samples[i]][genotype] += 1
                    else:
                        genotype_counts[samples[i]]['./.'] += 1

        return genotype_counts

def calculate_proportions(genotype_counts):
    proportions = {}
    for sample, counts in genotype_counts.items():
        total = sum(counts.values())
        proportions[sample] = {}
        for genotype, count in counts.items():
            proportions[sample][genotype] = count / total if total > 0 else 0
    return proportions

def main(input_vcf, output_file):
    genotype_counts = parse_vcf(input_vcf)
    proportions = calculate_proportions(genotype_counts)

    with open(output_file, 'w') as f:
        f.write("Sample\t0/0\t0/1\t1/1\t./.\t0/0(%)\t0/1(%)\t1/1(%)\t./.(%)\n")
        for sample, counts in genotype_counts.items():
            f.write(f"{sample}\t"
                    f"{counts['0/0']}\t{counts['0/1']}\t{counts['1/1']}\t{counts['./.']}\t"
                    f"{proportions[sample]['0/0']*100:.2f}\t"
                    f"{proportions[sample]['0/1']*100:.2f}\t"
                    f"{proportions[sample]['1/1']*100:.2f}\t"
                    f"{proportions[sample]['./.']*100:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count and calculate proportions of genotypes in a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file (can be gzipped).")
    parser.add_argument("output_file", help="Output file to save the results.")
    
    args = parser.parse_args()
    main(args.input_vcf, args.output_file)

