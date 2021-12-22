from Bio import SeqIO
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import swalign
from Bio.SeqRecord import SeqRecord
def paralogfinder (proj_name,totalpart,localpart):    
    path1='/home/ly/paralog-finder/mutation/'
    path2='/home/ly/paralog-finder/kpgenome/'
    path3='/home/ly/paralog-finder/mutation/'+proj_name+'/'
    #search
    match=2
    mismatch=-1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    
    os.chdir(path3)
    for tem in SeqIO.parse(path3+proj_name+'.fasta','fasta'):
        template=tem.seq
        templateid=tem.id
        fastafile=tuple()
    fastafile=(SeqRecord(template, id=templateid),)
    os.chdir(path2)
    filelist=os.listdir()
    number=int(len(filelist)/totalpart)
    for i in range(math.floor(localpart*number),math.ceil((localpart+1)*number)):
        name=filelist[i]
        tempscore=[]
        allseq = [seq_record.seq for seq_record in SeqIO.parse(name, "fasta")]
        allid= [seq_record.id for seq_record in SeqIO.parse(name, "fasta")]
        for k in range(len(allseq)):
            score=sw.align(template,allseq[k]).score
            tempscore.append(score)
            if score>len(template)*0.1:
                current_gene=(SeqRecord(allseq[k], id=allid[k]),)
                fastafile=fastafile+current_gene
        fig,ax=plt.subplots()
        ax.plot(tempscore,color='g',linewidth=4)
        ax.set_ylim(10,int(len(template)*2.1))
        plt.savefig(path3+proj_name+str(name)+'.png',dpi=500)
 
    SeqIO.write(fastafile, path1+proj_name+str(localpart)+'.fasta', "fasta")