import pysam
import os

def normalize_chrom_name(chrom):
    """Convert chromosome name between UCSC and Ensembl formats."""
    if chrom.startswith('chr'):
        # Convert from UCSC to Ensembl style
        return chrom[3:] if chrom != 'chrM' else 'MT'
    else:
        # Convert from Ensembl to UCSC style
        return f"chr{chrom}" if chrom != 'MT' else 'chrM'

def annotate_vcf_with_tabix(input_vcf, dbsnp_file, output_vcf):
    """
    Annotate VCF file with dbSNP information using tabix indices for efficient access.
    Handles chromosome name differences between files.
    """
    print("Starting VCF annotation using tabix indices...")
    
    try:
        vcf_in = pysam.VariantFile(input_vcf, 'r')
        dbsnp = pysam.VariantFile(dbsnp_file, 'r')
        vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
        
        # Add dbSNP info field to header
        vcf_in.header.info.add('dbSNP_RS', '1', 'String', 'dbSNP RS ID')
        
        # Check chromosome naming in both files
        vcf_has_chr = any(contig.startswith('chr') for contig in vcf_in.header.contigs)
        dbsnp_has_chr = any(contig.startswith('chr') for contig in dbsnp.header.contigs)
        
        print(f"VCF uses 'chr' prefix: {vcf_has_chr}")
        print(f"dbSNP uses 'chr' prefix: {dbsnp_has_chr}")
        
        # Process variants
        for variant in vcf_in:
            # Convert chromosome name if needed
            if vcf_has_chr != dbsnp_has_chr:
                query_chrom = normalize_chrom_name(variant.chrom)
            else:
                query_chrom = variant.chrom
                
            # Fetch overlapping variants from dbSNP
            try:
                for dbsnp_record in dbsnp.fetch(query_chrom, variant.pos-1, variant.pos):
                    # Check if variants match
                    if (dbsnp_record.pos == variant.pos and 
                        dbsnp_record.ref == variant.ref and 
                        any(a1 == a2 for a1 in variant.alts for a2 in dbsnp_record.alts)):
                        
                        # Add dbSNP ID to variant
                        variant.id = dbsnp_record.id
                        break
            
            except ValueError as e:
                print(f"Warning: Could not fetch region {query_chrom}:{variant.pos}")
                print(f"Available contigs in dbSNP: {list(dbsnp.header.contigs)}")
                # Continue with next variant
                continue
                
            # Write annotated variant
            vcf_out.write(variant)
            
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure both .tbi index files exist and you have pysam installed")
        return False
    
    finally:
        # Close files
        vcf_in.close()
        dbsnp.close()
        vcf_out.close()
    
    print(f"Annotation complete. Output written to {output_vcf}")
    return True

def main():
    # Check for required files
    input_vcf = "data/filtered.snp.vcf.gz"
    dbsnp_file = "data/00-All.vcf.gz"
    output_vcf = "data/annotated.snp.vcf.gz"
    
    # Check for index files
    if not os.path.exists(input_vcf + '.tbi'):
        print(f"Error: Index file {input_vcf}.tbi not found")
        return
    if not os.path.exists(dbsnp_file + '.tbi'):
        print(f"Error: Index file {dbsnp_file}.tbi not found")
        return
    
    # Run annotation
    annotate_vcf_with_tabix(input_vcf, dbsnp_file, output_vcf)

if __name__ == "__main__":
    main()