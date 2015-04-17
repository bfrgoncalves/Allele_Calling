from Bio import SeqIO
from BCBio import GFF
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import sys
import os
import re
import HTSeq
import subprocess
from Bio.Seq import Seq
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from CommonFastaFunctions import runBlastParserTAB
from os import listdir
from os.path import isfile, join
from datetime import datetime
from Bio import SearchIO


def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:

        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]

def translateSeq(DNASeq):
    seq=DNASeq
    try:
        myseq= Seq(seq)
        #print myseq
        protseq=Seq.translate(myseq, table=11,cds=True)
    except:
        try:
            seq=reverseComplement(seq)
            myseq= Seq(seq)
            #print myseq
            protseq=Seq.translate(myseq, table=11,cds=True)
                        
        except:
            try:
                seq=seq[::-1]
                myseq= Seq(seq)
                #print myseq
                protseq=Seq.translate(myseq, table=11,cds=True)
            except:
                try:
                    seq=seq[::-1]                           
                    seq=reverseComplement(seq)
                    myseq= Seq(seq)
                    #print myseq
                    protseq=Seq.translate(myseq, table=11,cds=True)
                except:
                    raise
    return protseq


def Create_Blastdb( questionDB, overwrite, dbtypeProt, Path ):
    isProt=dbtypeProt

    if not os.path.isfile(Path + ".nin") and not os.path.isfile(Path + ".nhr") and not os.path.isfile(Path + ".nsq"):
        
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + Path + " -dbtype nucl -logfile " + Path + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + Path + " -dbtype prot -logfile " + Path + "_blast.log" )

    elif overwrite:
        if not isProt:
            os.system( "makeblastdb -in " + questionDB + " -out " + Path + " -dbtype nucl -logfile " + Path + "_blast.log" )
        else:
            os.system( "makeblastdb -in " + questionDB + " -out " + Path + " -dbtype prot -logfile " + Path + "_blast.log" )

    else:
        print "BLAST DB files found. Using existing DBs.."  
    return( Path )


def runProdigalContig(FASTAfile):

        # ------------ #
        # RUN PRODIGAL #
        # ------------ #  #/scratch/NGStools/prodigal-2.50/prodigal -i /scratch/spneumoniae_600/ERR067978/velvet_assembly/contigs.fa -c -m -g 11 -p single -f sco -q > test_all.txt
    path = "cd Prodigal/Prodigal-2.60 && ./prodigal -i ../../"+FASTAfile+" -c -m -g 11 -p single -f sco -q > ../../Presults.txt";
    os.system(path)
    cdsDict = {}
    tempList = []
    line = ' '
    x = open("Presults.txt","r")
    countContigs=0
    while line != '':
        
        # when it finds a contig tag
        if "seqhdr" in line:
        # add contig to cdsDict and start new entry
            
            if len(tempList) > 0:
                
                cdsDict[contigTag.replace("\r","")] = tempList
                tempList = []
            
            countContigs+=1
            contigTag = str(countContigs)
            
                # when it finds a line with cds indexes
        elif line[0] == '>':
            
                        # parsing
            cdsL = line.split('_')

                        # --- each element of this list is a pair of indices - the start and the end of a CDS --- #

            tempList.append([ int(cdsL[1]) - 1 , int(cdsL[2]) ])        # start index correction needed because prodigal indexes start in 1 instead of 0
        line = x.readline()
    if len(tempList) > 0:
         # --- brute force parsing of the contig tag - better solution is advisable --- #
        i=0
        for l in contigTag:
            if l == ' ':
                break
            i+=1
        contigTag=contigTag[:i]
        cdsDict[contigTag.replace("\r","")] =tempList

    os.remove('Presults.txt')

    return cdsDict


def getFASTAarray(FASTAfile, genomeArray):
    g_fp = HTSeq.FastaReader(FASTAfile)
    countContigs=0
    for contig in g_fp:
        countContigs+=1
        genomeArray[str(countContigs)]=contig.seq
    return genomeArray

def CreateNewAlleleFile(FASTAfile, ToNewAllele):

    fG = open( FASTAfile, 'a' )
    for i in ToNewAllele:
        fG.write(i)
    fG.close()

def CreateProteomeContig(cdsDict, referenceCDS, referenceGenomeArray):
    alleleProt=''
    proteome=""
    countP=0
    countCDS=0
    for j in cdsDict:
        for i in cdsDict[j]:
            countP+=1
            countCDS+=1
            AlignBegin=i[0]
            AlignEnd=i[1]
            IdCDS=str(j)+"--"+str(countCDS)
            referenceCDS[IdCDS] = [int(AlignBegin),int(AlignEnd)]
                    #print contigTag,value

                        
                        #print contigTag, protein[0], protein[1]
                        #print currentGenomeDict[ contigTag ]
            seq= referenceGenomeArray[j][ int(AlignBegin):int(AlignEnd) ].upper()
            try:
                protseq=translateSeq(seq)
            except:
                continue
            proteome+=">"+IdCDS+"\n"+str(protseq)+"\n"
                #proteome+=str(protseq)
    dirName = 'Referenceproteome.fasta'
    if os.path.isfile(dirName):
        os.remove(dirName)

    with open(dirName, "wb") as v:
        v.write(proteome)
    return dirName, referenceCDS


def CreateReferenceDatabase(referenceFASTA, databasePath):
    referenceGenomeArray = {}
    referenceCDS = {}
    referenceGenomeArray = getFASTAarray(referenceFASTA, referenceGenomeArray)
    CDSdict = runProdigalContig(referenceFASTA)
    refProteomePath, referenceCDS = CreateProteomeContig(CDSdict, referenceCDS, referenceGenomeArray)
    referenceDatabase = os.path.join(databasePath,'reference_db')
    Create_Blastdb(refProteomePath,1,True, referenceDatabase)

    return referenceDatabase, referenceCDS, referenceGenomeArray


def CreateQueryDatabase(FASTAfile, databasePath,queryProteomeName):
    gene_fp = HTSeq.FastaReader(FASTAfile)
    names=""
    alleleProt=''
    proteome=""
    isEmpty = True
    countAlleles = 0
    for allele in gene_fp: #new db for each allele to blast it against himself
        try:
            x = str(translateSeq(allele.seq))
            countAlleles+=1
            isEmpty = False
        except:
            print 'Could not translate'
            if countAlleles==0:
                isEmpty = True
            continue
        alleleProt+=">"+str(allele.name)+"\n"+x+"\n"
        proteome+=">"+str(allele.name)+"\n"+x+"\n"
    # with open(pathRef+'allAllelesAA.fasta', "wb") as f:
    #     f.write(alleleProt)
    databasePath = os.path.join(databasePath,queryProteomeName)
    databasePath = databasePath.split('.')[0]
    databasePath = databasePath+'_db'
    with open(queryProteomeName, "wb") as v:
        v.write(proteome)
    Gene_Blast_DB_name = Create_Blastdb(queryProteomeName,1,True, databasePath)
    
    return databasePath, isEmpty



def getOwnBlastScore(FASTAfile, databasePath, queryProteomeName, isXML):

    databasePath, isEmpty = CreateQueryDatabase(FASTAfile, databasePath,queryProteomeName)
        # --- get BLAST score ratio --- #

    allelescores = {}
    queryDef = ''

    if isEmpty:
        return allelescores, queryDef, isEmpty

    if isXML == 'True':

        blast_out_file = 'BLASTresults.xml'

        #cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=5, num_alignments=7000, num_descriptions=7000)

        print 'OwnBLAST:'
            #print cline
        allelescore=0

        #blast_records = runBlastParser(cline,blast_out_file, False)

        startTime = datetime.now()

        allelescores = {}

        #allelescores, queryDef = parseOwnBLASTRecordsXML(blast_records)

        print 'CheckResults:' + str(datetime.now() - startTime)

    else:
        blast_out_file = 'BLASTresults.tab'

        cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=6, num_alignments=7000, num_descriptions=7000)

        blast_records = runBlastParserTAB(cline,blast_out_file, False)
        
        startTime = datetime.now()
        
        allelescores = parseOwnBLASTRecordsTAB(blast_records)

        print 'CheckResults:' + str(datetime.now() - startTime)

    os.remove(databasePath+ ".pin")
    os.remove(databasePath+ ".phr")
    os.remove(databasePath+ ".psq")

    return allelescores, queryDef, isEmpty


def getBlastScoreRatios(FASTAfile, allelescore, queryDef, databasePath, queryProteomeName, referenceGenomeArray, referenceCDS, isXML):
    
    alleleProt=''
    proteome=""
    countP=0
    countCDS=0

    if isXML == 'True':
        blast_out_file = 'BLASTresults.xml'

        cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=5, num_alignments=7000, num_descriptions=7000)

        print 'BSR:'
            #print cline
        blast_records = runBlastParser(cline,blast_out_file, False)

        startTime = datetime.now()
        ToNewAllele = parseBLASTRecordsXML(blast_records, allelescore, queryDef, referenceGenomeArray, referenceCDS)
        print 'CheckResults:' + str(datetime.now() - startTime)

    else:
        blast_out_file = 'BLASTresults.tab'

        cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=6, num_alignments=7000, num_descriptions=7000)

        print 'BSR:'
        blast_records = runBlastParserTAB(cline,blast_out_file, False)

        startTime = datetime.now()
        ToNewAllele = parseBLASTRecordsTAB(blast_records, allelescore, queryDef, referenceGenomeArray, referenceCDS)
        print 'CheckResults:' + str(datetime.now() - startTime)

    os.remove(queryProteomeName)

    return ToNewAllele


def parseOwnBLASTRecordsXML(blast_records):

    allelescores = {}
    queryDef = ''

    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                for match in alignment.hsps:
                    try:
                        if allelescores[str(alignment.hit_def)] < match.score:
                            allelescores[str(alignment.hit_def)] = int(match.score)
                            queryDef = str(alignment.hit_def)
                            #break
                    except KeyError:
                        allelescores[str(alignment.hit_def)] = int(match.score)
                        queryDef = str(alignment.hit_def)
                        #break
            else:
                break

    return allelescores, queryDef

def parseOwnBLASTRecordsTAB(blast_records):

    allelescores = {}
    queryDef = ''

    for qresult in blast_records:
        for match in qresult:
            continue

    return allelescores


def parseBLASTRecordsXML(blast_records, allelescore, queryDef, referenceGenomeArray, referenceCDS):

    newAlleles = []
    ToNewAllele = []
    queryName = ''
    for blast_record in blast_records:
        queryName = blast_record.query
        for alignment in blast_record.alignments:
            isThere = False
            for match in alignment.hsps:

                blastScoreRatio = float(match.score) / float(allelescore[queryDef])
                
                if blastScoreRatio == 1:
                    #print 'Perfect match'
                    continue
                    #break
                
                elif blastScoreRatio>=0.6:
                    #print 'New Allele'
                    if isThere == False:
                        newAlleles.append(str(alignment.hit_def))
                        isThere = True
                    #print match.query
                    #print match.sbjct
                    #print allelescores
                    #break

    prevAllele = queryName

    for i in newAlleles:
        
        try:
            nextAllele=int(prevAllele.split("--")[1])+1
        except IndexError:
            nextAllele = 1
        newName = prevAllele.split("--")[0] + '--' + str(nextAllele)
        contig = i.split('--')[0]
        ToNewAllele.append(">" + str(newName) +"\n"+ referenceGenomeArray[contig][ int(referenceCDS[i][0])-1:int(referenceCDS[i][1]) ].upper()+ "\n")

    return ToNewAllele


def parseBLASTRecordsTAB(blast_records, allelescore, queryDef, referenceGenomeArray, referenceCDS):

    ToNewAllele = []

    for qresult in blast_records:
        for match in qresult:
            continue

    return ToNewAllele


def concatAllQueries(queryPath, maxBP, maxalleles):


    try:
        fg = open(os.path.join(queryPath,'Allalleles.fasta'),'w')
        fg.close()
        queryFilesOnDir = [ f for f in listdir(queryPath) if isfile(join(queryPath,f)) ]
        countAlleles = 0
        ToWrite = []
        for queryFile in queryFilesOnDir:
            if queryFile == 'Allalleles.fasta':
                ToWrite = []
                continue
            ToWrite = []
            g_fp = HTSeq.FastaReader(os.path.join(queryPath,queryFile))

            count = 0
        
            for allele in g_fp:
                if maxalleles != None and int(maxalleles) == count:
                    break
                else:
                    if maxBP != None:
                        if len(str(allele.seq)) > maxBP:
                            continue
                    countAlleles += 1
                    ToWrite.append(">" + str(countAlleles) +"\n"+ str(allele.seq).upper() + "\n")

                count+=1
            CreateNewAlleleFile(os.path.join(queryPath,'Allalleles.fasta'), ToWrite)
    except Exception:
        print 'An error occurred'
        return False
    
    return True


def writeFilteredAlleleFile(FASTAfile, newFile, maxalleles):

    Done = True
    try:
        g_fp = HTSeq.FastaReader(FASTAfile)
        newFile = open(newFile,'w')
        countAlleles=0
        for alleles in g_fp:
            if maxalleles != None and maxalleles == countAlleles:
                break
            countAlleles+=1
            newFile.write(">" + str(alleles.name) +"\n"+ str(alleles.seq).upper()+ "\n")
    except Exception:
        Done = False

    return Done