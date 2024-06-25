import gzip
import argparse
import pandas as pd

def parse_vcf(vcf_file, sample_list):
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'rt') as f:
        samples = []
        sample_indices = []
        consistent_sites = 0
        inconsistent_sites = 0
        total_sites = 0
        genotype_counts = {sample: {'0/0': 0, '0/1': 0, '1/1': 0, './.': 0} for sample in sample_list}

        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                header = line.strip().split('\t')
                samples = header[9:]
                sample_indices = [samples.index(sample) for sample in sample_list if sample in samples]
            else:
                columns = line.strip().split('\t')
                genotypes = [columns[9 + i].split(':')[0] for i in sample_indices]
                
                if len(set(genotypes)) == 1:
                    consistent_sites += 1
                else:
                    inconsistent_sites += 1
                    for i, genotype in enumerate(genotypes):
                        if genotype in genotype_counts[sample_list[i]]:
                            genotype_counts[sample_list[i]][genotype] += 1
                        else:
                            genotype_counts[sample_list[i]]['./.'] += 1

                total_sites += 1

        return consistent_sites, inconsistent_sites, total_sites, genotype_counts

def calculate_proportions(genotype_counts, inconsistent_sites):
    proportions = {}
    for sample, counts in genotype_counts.items():
        proportions[sample] = {}
        for genotype, count in counts.items():
            proportions[sample][genotype] = count / inconsistent_sites if inconsistent_sites > 0 else 0
    return proportions

def main(input_vcf, output_file, sample_list_file):
    with open(sample_list_file, 'r') as f:
        sample_list = [line.strip() for line in f.readlines()]
    
    consistent_sites, inconsistent_sites, total_sites, genotype_counts = parse_vcf(input_vcf, sample_list)
    proportions = calculate_proportions(genotype_counts, inconsistent_sites)

    with open(output_file, 'w') as f:
        f.write("Sample List: {}\n".format(', '.join(sample_list)))
        f.write("Consistent Sites: {}\n".format(consistent_sites))
        f.write("Inconsistent Sites: {}\n".format(inconsistent_sites))
        f.write("Total Sites: {}\n".format(total_sites))
        f.write("Proportion of Consistent Sites: {:.2%}\n".format(consistent_sites / total_sites if total_sites > 0 else 0))
        f.write("Proportion of Inconsistent Sites: {:.2%}\n".format(inconsistent_sites / total_sites if total_sites > 0 else 0))
        f.write("\nGenotype counts and proportions for inconsistent sites:\n")
        f.write("Sample\t0/0\t0/1\t1/1\t./.\t0/0(%)\t0/1(%)\t1/1(%)\t./.(%)\n")
        for sample, counts in genotype_counts.items():
            f.write(f"{sample}\t"
                    f"{counts['0/0']}\t{counts['0/1']}\t{counts['1/1']}\t{counts['./.']}\t"
                    f"{proportions[sample]['0/0']*100:.2f}\t"
                    f"{proportions[sample]['0/1']*100:.2f}\t"
                    f"{proportions[sample]['1/1']*100:.2f}\t"
                    f"{proportions[sample]['./.']*100:.2f}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare genotypes of specified samples at each site in a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file (can be gzipped).")
    parser.add_argument("output_file", help="Output file to save the results.")
    parser.add_argument("sample_list_file", help="File containing list of samples to compare (one per line).")
    
    args = parser.parse_args()
    main(args.input_vcf, args.output_file, args.sample_list_file)

