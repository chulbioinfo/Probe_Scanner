# Probe Scanner
# version: 20230123
# Developed by Chul Lee (clee03@rockefeller.edu)
print("# Start - Prove Scanner (v20230123)")

# library
import sys
import os
import glob

# set global variable
print("# Set global variables")
try:
    fNAME_genome = sys.arg[1]
    fNAME_transcript = sys.arg[2]
    iMer = int(sys.arg[3])
except:
    fNAME_genome = "ref/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna"
    fNAME_transcript = "ref/GCF_003957565.2_bTaeGut1.4.pri_rna.fna"
    iMer = 16
try:
    os.system("mkdir scaff")
    os.system("mkdir transcript")
    os.system("mkdir kmer_genome_"+str(iMer))
    os.system("mkdir kmer_transcriptome_"+str(iMer))
except:
    pass

oPATH_scaff = "./scaff/"
oPATH_transcript = "./transcript/"
oPATH_kmer_genome = "./kmer_genome_"+str(iMer)+"/"
oPATH_kmer_transcriptome = "./kmer_transcriptome_"+str(iMer)+"/"
oNAME_kmer_genome = "./kmer_genome_"+str(iMer)+".txt"
oNAME_kmer_transcriptome = "./kmer_transcriptome_"+str(iMer)+".txt"
oNAME_gMer_wo_tMer = "ProbScan.out."+str(iMer)+".txt"
oNAME_gMer_wo_tMer_filtered = "ProbScan.out."+str(iMer)+".filtered.txt"

# functions
def split_fasta(fNAME_fasta,oPATH_fasta):
    print("# split fasta: ",fNAME_fasta)
    fpin = open(fNAME_fasta,'r')
    for line in fpin:
        if line[0]==">":
            eachFasta = line[1:].split(" ")[0]
            fpout = open(oPATH_fasta + eachFasta + ".fa",'w')
            fpout.write(line)
        else:
            fpout.write(line.strip())
    fpin.close()
    fpout.close()

def scan_kmer(oPATH_fasta,oPATH_kmer):
    print("# scan k-mer: ",oPATH_fasta)
    flist = glob.glob(oPATH_fasta+"*.fa")
    iCNT_fasta = len(flist)
    iCNT = 0
    iPercent = 0
    for fNAME_fasta in flist:
        iCNT += 1
        if iCNT in range(int(iCNT_fasta/10),iCNT_fasta,int(iCNT_fasta/10)):
            iPercent += 10
            print("# scan k-mer - "+str(iPercent)+'%')
        nFasta = os.path.basename(fNAME_fasta).split(".fa")[0]
        fpin = open(fNAME_fasta,'r')
        fpin.readline()
        nSeq = fpin.readline()
        fpin.close()
        nMer_iCNT_dic = {}
        for i in range(0,len(nSeq)-iMer):
            nMer = nSeq[i:i+iMer].upper()
            nMer_iCNT_dic.setdefault(nMer,0)
            nMer_iCNT_dic[nMer]+=1
        fpout = open(oPATH_kmer+nFasta,'w')
        for nMer in nMer_iCNT_dic.keys():
            tmpline = nMer +'\t'+ str(nMer_iCNT_dic[nMer]) +'\n'
            fpout.write(tmpline)
        fpout.close()

def sum_kmer(oPATH_kmer, oNAME_kmer):
    print("# sum_kmer: ",oPATH_kmer)
    flist = glob.glob(oPATH_kmer+"*")
    nMer_iCNT_dic = {}
    for fNAME in flist:
        fpin = open(fNAME,'r')
        for line in fpin:
            part = line.strip().split("\t")
            nMer = part[0]
            iCNT = int(part[1])
            nMer_iCNT_dic.setdefault(nMer,0)
            nMer_iCNT_dic[nMer]+=iCNT
        fpin.close()
    fpout = open(oNAME_kmer,'w')
    for nMer in nMer_iCNT_dic.keys():
        tmpline = nMer+'\t'+str(nMer_iCNT_dic[nMer])+'\n'
        fpout.write(tmpline)
    fpout.close()


def comp_kmer(oNAME_kmer_transcriptome, oNAME_kmer_genome, oNAME_gMer_wo_tMer):
    print("# compare kmer in genome with kmer in transcriptome")
    tMer_dic = {}
    fpin = open(oNAME_kmer_transcriptome,'r')
    for line in fpin:
        part = line.strip().split("\t")
        nMer = part[0]
        tMer_dic.setdefault(nMer,'')
    fpin.close()

    # compare k-mers in genome with k-mers in transcriptome
    fpout = open(oNAME_gMer_wo_tMer,'w')
    fpin = open(oNAME_kmer_genome,'r')
    for line in fpin:
        part = line.strip().split("\t")
        gMer = part[0]
        iCNT = int(part[1])
        if not gMer in tMer_dic.keys():
            if iCNT >= 300:
                if iCNT <= 600:
                   fpout.write(gMer+'\t'+str(iCNT)+'\n')
    fpin.close()
    fpout.close()


def check_self_complementary(oNAME_gMer_wo_tMer, oNAME_gMer_wo_tMer_filtered):
    print("# check self complementary")
    nBP_cpBP_dic = {"A":"T","G":"C","T":"A","C":"G"}
    fpin = open(oNAME_gMer_wo_tMer,'r')
    fpout = open(oNAME_gMer_wo_tMer_filtered,'w')
    for line in fpin:
        part = line.strip().split("\t")
        nSeq_probe = part[0]
        iFlag = False
        if not 'N' in nSeq_probe:
            iFlag = True
            for i in range(2,9):
                for j in range(0,len(nSeq_probe)-i):
                    nMer = nSeq_probe[j:j+i]
                    cpMer = ''
                    for nBP in nMer[::-1]:
                        cpBP = nBP_cpBP_dic[nBP]
                        cpMer+=cpBP
                    if cpMer in nSeq_probe:
                        iFlag = False
        if iFlag == True:
            fpout.write(line)
    fpin.close()
    fpout.close()


def main():
    split_fasta(fNAME_genome,oPATH_scaff)
    split_fasta(fNAME_transcript,oPATH_transcript)
    scan_kmer(oPATH_scaff,oPATH_kmer_genome)
    scan_kmer(oPATH_transcript,oPATH_kmer_transcriptome)
    sum_kmer(oPATH_kmer_genome, oNAME_kmer_genome)
    sum_kmer(oPATH_kmer_transcriptome, oNAME_kmer_transcriptome)
    comp_kmer(oNAME_kmer_transcriptome, oNAME_kmer_genome, oNAME_gMer_wo_tMer)
    check_self_complementary(oNAME_gMer_wo_tMer, oNAME_gMer_wo_tMer_filtered)


# main
if __name__ == "__main__":
    main()
    
print("# Finish: Probe Scanner (v20230123")
