import argparse
import sqlite3
import pyfaidx
import vcf
import logging
import gffutils
import os
from Bio.Seq import Seq
import csv
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="To input files for iteration")
parser.add_argument('--vcf', required=True, help='VCF file')
parser.add_argument('--gff', required=True, help='GFF file')
parser.add_argument('--fasta', required=True, help='Fasta file')
parser.add_argument('--quality', default=20, help='Quality filter')
parser.add_argument('--output', required=True, help='output filenames')
args = parser.parse_args()

# logger print out to log file
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('{}.log'.format(args.output))
fh.setLevel(logging.INFO)
# formatter for clarity and to separate multiple runs
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)

# all arguments used as input in command line
logger.info('The files used in command line:\nVCF: {}\tGFF: {}\tfasta: {} \nQuality used: {}'.format(args.vcf, args.gff, args.fasta, args.quality))

# make database prefix same as gff file
db = args.gff.replace('.gff', '.db')

# opening the gff file to create the database
if not os.path.isfile(db):
    try:
        db = gffutils.create_db(args.gff, dbfn=db, keep_order=True)
    except sqlite3.OperationalError:
        logger.error('Unable to create database {}\n'.format(db))
        raise SystemExit(1)
    except ValueError:
        logger.error('File {} does not exist\n'.format(db))
        raise SystemExit(1)

# if database already exists then just access
else:
    try:
        db = gffutils.FeatureDB(db, keep_order=True)
    except ValueError:
        logger.error('Cannot access database {}\n'.format(db))
        raise SystemExit(1)
    except sqlite3.DatabaseError:
        logger.error('File {} is not a database\n'.format(db))
        raise SystemExit(1)


# open vcf file with try and except block
try:
    vcffile = vcf.Reader(filename=args.vcf)
except FileNotFoundError:
    logger.error('File {} does not exist\n'.format(args.vcf))
    raise SystemExit(1)


# empty list for storing of all filtered records
allrecords = []
# counts for the different types of variants
syn = 0
nonsyn = 0
noncode = 0
failfilter = 0
for entry in vcffile:
    try:
        # filter for quality based on input
        if entry.QUAL > float(args.quality):
            # get just the chromosome number for tsv file
            chrom = str(entry.CHROM)
            chrom = chrom.split('_')
            chrom = chrom[1]

            # checking if entry has CDS features, if not then it goes into non-coding type
            if len(list(db.region(seqid=entry.CHROM, start=entry.POS, end=entry.POS, featuretype='CDS'))) > 0:
                # selecting features that contain our SNPs
                for snps in db.region(seqid=entry.CHROM, start=entry.POS, end=entry.POS, featuretype='CDS'):
                    # going up to transcript level and selecting for mRNA, thus ignoring pseudogenes
                    for parent in db.parents(snps.id, featuretype='mRNA'):
                        snpseq = ''
                        # selecting all the CDS of transcript to recreate spliced mRNA sequence
                        for child in db.children(parent.id, featuretype='CDS'):
                            try:
                                # joining all CDS of transcript to create untranslated mRNA sequence
                                snpseq = snpseq + child.sequence(args.fasta, use_strand=True)
                            except pyfaidx.FastaNotFoundError:
                                logger.error('Fasta not found in {}\n'.format(args.fasta))
                                raise SystemExit(1)

                        # calculating position of SNP in untranslated mRNA sequence
                        mutpos = 0
                        for mutation in db.children(parent.id, featuretype='CDS'):
                            # if SNP in feature then calculate position of base then break the loop
                            if mutation in db.region(seqid=entry.CHROM, start=entry.POS, end=entry.POS, featuretype='CDS'):
                                if parent.strand == '+':
                                    mutpos = mutpos + (entry.POS - mutation.start + 1)
                                    break
                                elif parent.strand == '-':
                                    mutpos = mutpos + (mutation.end - entry.POS + 1)
                                    break
                            # if SNP not in feature then calculate length of feature and add it to mutpos
                            else:
                                mutpos = mutpos + (mutation.end - mutation.start + 1)
                                continue

                        # amino acid position has to be calculated to check if change has occurred
                        aapos = 0
                        # if SNP in last position of codon, then amino acid position is mutpos/3
                        if mutpos % 3 == 0:
                            aapos = int(mutpos / 3)
                        # if SNP not in last position, then mutpos/3 will be 1 less than the codon the SNP is in
                        else:
                            aapos = int((mutpos / 3) + 1)

                        # since entry.ALT is a list
                        for i in entry.ALT:
                            # turn the untranslated sequence into a list for easy mutation
                            mutlist = list(snpseq)
                            mutlist[mutpos - 1] = str(i)
                            # once mutation is done, turn the list back into a string
                            mutatedseq = ''.join(mutlist)
                            # convert original and mutated sequences in Seq objects
                            snpseq = Seq(snpseq)
                            mutatedseq = Seq(mutatedseq)
                            # once they're Seq objects, translate can be used
                            seqtrans = snpseq.translate()
                            muttrans = mutatedseq.translate()
                            # select the codon which contains the SNP in the reference protein
                            refaa = seqtrans[aapos - 1]

                            # list to store all data for each valid entry, it will later
                            # be appended to bigger list created at start
                            record = []
                            # check if codon in reference protein matches the one in mutated protein
                            if seqtrans[aapos-1] == muttrans[aapos-1]:
                                # codons are the same so synonymous mutation
                                # add 1 to count and append all relevant data to smaller list
                                syn = syn + 1
                                record.append(chrom)
                                record.append(entry.POS)
                                record.append(entry.REF)
                                record.append(i)
                                record.append('Synonymous')
                                record.append(parent.id)
                                record.append(aapos)
                                record.append(refaa)
                                record.append('NA')
                            else:
                                # codons have changed so non-synonymous mutation
                                # add 1 to count and append all relevant data to smaller list
                                nonsyn = nonsyn + 1
                                altaa = muttrans[aapos - 1]
                                record.append(chrom)
                                record.append(entry.POS)
                                record.append(entry.REF)
                                record.append(i)
                                record.append('Non-Synonymous')
                                record.append(parent.id)
                                record.append(aapos)
                                record.append(refaa)
                                record.append(altaa)

                            # append smaller list into bigger list so that it becomes a list of lists of all mutations
                            # and their relevant data
                            allrecords.append(record)
                    # another loop to iterate through the pseudogenes since even
                    # though they have a CDS, they are non coding
                    record3 = []
                    for pseudos in db.parents(snps.id, featuretype='pseudogenic_transcript'):
                        noncode += 1
                        for p in entry.ALT:
                            record3.append(chrom)
                            record3.append(entry.POS)
                            record3.append(entry.REF)
                            record3.append(p)
                            record3.append('Non-Coding')
                            record3.append('NA')
                            record3.append('NA')
                            record3.append('NA')
                            record3.append('NA')
                        allrecords.append(record3)

            else:
                # entry does not have CDS features so is non-coding
                # add 1 to count, then for each alt base iteration, append data to list
                noncode += 1
                record2 = []
                for j in entry.ALT:
                    record2.append(chrom)
                    record2.append(entry.POS)
                    record2.append(entry.REF)
                    record2.append(j)
                    record2.append('Non-Coding')
                    record2.append('NA')
                    record2.append('NA')
                    record2.append('NA')
                    record2.append('NA')
                # append lists to same bigger list as one used for coding variants
                allrecords.append(record2)
        else:
            # variants that did not pass quality filter, add 1 to count
            failfilter += 1
    except TypeError:
        logger.error('Quality not a float\n')
        raise SystemExit(1)

logger.info('Count of variants with quality less than {}: {}'.format(args.quality,failfilter))

# open tsv file then write each small list in bigger list as a separate variant
with open('{}.tsv'.format(args.output), 'w') as out:
    filewriter = csv.writer(out, delimiter='\t')
    # include headers
    filewriter.writerow(['CHROM', 'POS', 'REF', 'ALT', 'Type', 'ProteinLocation', 'RefAA', 'AltAA'])
    for lst in allrecords:
        filewriter.writerow(lst)

# creating lists to make the bar plot
types = ['Non-Coding', 'Synonymous', 'Non-Synonymous']
variants = [noncode, syn, nonsyn]

plt.bar(types, variants, color=['red', 'green', 'blue'])

# add appropriate labels for each observation and the name for the plot
plt.xlabel("Types")
plt.ylabel("Variants")
plt.title("Frequency of Variants in each Type")
# save the graph
plt.savefig('{}.png'.format(args.output))

# log info for output files generated
logger.info('The output files generated: {}.log\t{}.tsv\t{}.png\n'.format(args.output, args.output, args.output))
