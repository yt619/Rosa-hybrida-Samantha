import pandas as pd
import argparse

def read_bed_file(bed_file):
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    return bed_df

def read_gff_file(gff_file):
    gff_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gff_df = pd.read_csv(gff_file, sep='\t', comment='#', header=None, names=gff_columns)
    return gff_df

def extract_genes(bed_df, gff_df):
    result = []
    for _, bed_row in bed_df.iterrows():
        chrom = bed_row['chrom']
        start = bed_row['start']
        end = bed_row['end']
        
        genes_in_region = gff_df[(gff_df['seqname'] == chrom) & 
                                 (gff_df['start'] >= start) & 
                                 (gff_df['end'] <= end) & 
                                 (gff_df['feature'] == 'gene')]
        
        if not genes_in_region.empty:
            result.append(genes_in_region)
    
    if result:
        return pd.concat(result)
    else:
        return pd.DataFrame()

def main(bed_file, gff_file, output_file):
    bed_df = read_bed_file(bed_file)
    gff_df = read_gff_file(gff_file)
    
    genes_in_bed_regions = extract_genes(bed_df, gff_df)
    
    if not genes_in_bed_regions.empty:
        genes_in_bed_regions.to_csv(output_file, sep='\t', index=False)
        print(f"Extracted genes saved to {output_file}")
    else:
        print("No genes found in the specified regions.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract genes from GFF file based on BED file regions.')
    parser.add_argument('bed_file', type=str, help='Input BED file')
    parser.add_argument('gff_file', type=str, help='Input GFF file')
    parser.add_argument('output_file', type=str, help='Output file to save extracted genes')
    
    args = parser.parse_args()
    
    main(args.bed_file, args.gff_file, args.output_file)

