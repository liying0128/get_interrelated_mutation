from Bio import SeqIO
from Bio.Seq import Seq
import swalign
from Bio.SeqRecord import SeqRecord
def filter (file,grabscore):
    name=file.replace('.fasta','')
    allseq = [seq_record.seq for seq_record in SeqIO.parse(file, "fasta")]
    template=allseq[0]
    allid= [seq_record.id for seq_record in SeqIO.parse(file, "fasta")]
    template.id=allid[0]
    match=2
    mismatch=-1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    fastafile=tuple()
    fastafile=(SeqRecord(template, id=template.id),)
    for i in range (len(allseq)):
        score=sw.align(template,allseq[i]).score
        if score>=grabscore:
            current_gene=(SeqRecord(allseq[i], id=allid[i]),)
            fastafile=fastafile+current_gene
    SeqIO.write(fastafile, name+str(-grabscore)+'.fasta', "fasta")
    
def combine (fastalist):
    allseq = [seq_record.seq for seq_record in SeqIO.parse(fastalist[0], "fasta")]
    template=allseq[0]
    allid= [seq_record.id for seq_record in SeqIO.parse(fastalist[0], "fasta")]
    template.id=allid[0]
    fastafile=tuple()
    fastafile=(SeqRecord(template, id=template.id),)
    for i in range(1,len(fastalist[0])):
        current_gene=(SeqRecord(allseq[i], id=allid[i]),)
        fastafile=fastafile+current_gene
    for i in range(1,len(fastalist)):
        for seq_record in SeqIO.parse(fastalist[i], "fasta"):
            if seq_record.id!='template':
                current_gene=(SeqRecord(seq_record.seq,id=seq_record.id),)
                fastafile=fastafile+current_gene
    SeqIO.write(fastafile, 'combined.fasta', "fasta")