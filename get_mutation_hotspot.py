from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
def getspot(file):
    mutation=[]
    similarity=[]
    Protein_length=[]
    paralog_num=0
    for i in SeqIO.parse(file,'fasta'):
        paralog_num=paralog_num+1
        Protein_length.append(len(i))
    protein_length=max(Protein_length)
    for i in range(protein_length):
        a=[]
        for seq_record in SeqIO.parse(file,'fasta'):
            seq=str(seq_record.seq)
            seq=seq.replace('*','')
            seq=Seq(seq)
            if len(seq)<i+1:
                a.append(0)
            else:
                a.append(seq[i])
        result=pd.value_counts(a)
        similar=(1-result[0]/paralog_num)*100
        similarity.append(similar)
        if similar > 2.5:
            localmutation=[]
            localmutation.append(i)
            for k in range(len(result)):
                localmutation.append(result.index[k])
                localmutation.append(result.iloc[k])
            mutation.append(localmutation)
    plt.plot(range(protein_length),similarity,color='red')
    plt.ylabel('Mutation rate (%)',fontsize=12)
    plt.xlabel('Position',fontsize=12)
    plt.savefig('C:\\Users\\Administrator\\Desktop\\result.png',dpi=1000)
    return(mutation)
    print(mutation)    
    
def combine(a,b):
    n=len(a)-1
    templist=[]
    l2=[]
    for i in range(n):
        templist.append(a[i])
        templist.append(b[i])
    for i in templist:
        if i not in l2:
            l2.append(i)
    return (l2)
    
def getpair (imputlist,minsup,file):
    houxuan=[]
    for i in range(len(imputlist)):
        for k in range(4,len(imputlist[i]),2):
            if imputlist[i][k]>=minsup:
                templist=[]
                templist.append(imputlist[i][0])
                templist.append(imputlist[i][(k-1)])
                templist.append(imputlist[i][(k)])
                houxuan.append(templist)            
    
    double=[]
    for i in range(len(houxuan)):
        for k in range(i+1,len(houxuan)):
            m=0
            for seq_record in SeqIO.parse(file,'fasta'):
                seq=str(seq_record.seq)
                seq=seq.replace('*','')
                seq=Seq(seq)
                if seq[houxuan[i][0]]==houxuan[i][1] and seq[houxuan[k][0]]==houxuan[k][1]:
                    m=m+1
            if m>=minsup:
                double.append([houxuan[i],houxuan[k],m])
    return([houxuan,double])

def getover (double,minsup,file):               
   fif=[]
   final=[]
   for i in range(len(double)):
       for k in range(i+1,len(double)):
           new=combine(double[i],double[k])
           m=0
           if len(new)==len(double[0]):
               for seq_record in SeqIO.parse(file,'fasta'):
                   seq=str(seq_record.seq)
                   seq=seq.replace('*','')
                   seq=Seq(seq)
                   p=0
                   for q in range(len(new)):
                       if seq[new[q][0]]==new[q][1]:
                           p=p+1
                   if p==len(new):
                       m=m+1
           if m>=minsup:
               new.append(m)
               fif.append(new)
   for i in fif:
       if i not in final:
           final.append(i)
   return(final)
