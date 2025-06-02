import vcf
import argparse

def filter_vcf(input_vcf, output_vcf):
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        if record.CHROM not in ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']:
            continue
        if not record.FILTER or 'PASS' in record.FILTER:
            if (record.samples[0]['GT'][0]!=record.samples[0]['GT'][-1]) \
                    or (record.samples[0]['GT']!='0/0'):
                        vcf_writer.write_record(record)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter VCF file.')
    parser.add_argument('--input_vcf', required=True, help='Input VCF file.')
    parser.add_argument('--output_vcf', required=True, help='Output VCF file.')
    args = parser.parse_args()

    filter_vcf(args.input_vcf, args.output_vcf)