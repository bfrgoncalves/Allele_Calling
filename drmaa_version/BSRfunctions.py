#!/usr/bin/python
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import sys

import os
import re
import HTSeq
import subprocess
from Bio.Seq import Seq
#from CommonFastaFunctions import LoadAlelleFasta
#from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from CommonFastaFunctions import runBlastParserTAB
from os import listdir
from os.path import isfile, join
from datetime import datetime
from Bio import SearchIO
import csv



def SelfBLAST(args):
    alleleNumbers, sameAlleles, countNumberOfLocus, prevAlleleNames, isEmpty = func_Self_BLAST(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])

    return alleleNumbers, sameAlleles, countNumberOfLocus, prevAlleleNames, isEmpty

def BSR(args):
     resultsList, referenceCDSsequences, addNewAlleles = func_BSR(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])
     return resultsList, referenceCDSsequences, args[1], addNewAlleles

def check_ToWrite(args):
    ToWrite = TowriteToAlleleFiles(args[0], args[1])
    return ToWrite

def CreateNewAlleleFiles(args):
    NewAlleleFiles(args[0], args[1])

def get_correctAlleles(args):
    newArrayOfResults, genomeID = correctAlleles(args[0], args[1], args[2])

    return newArrayOfResults, genomeID



def func_Self_BLAST(queryFile, countNumberOfLocus, alleleFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleScoreFile, LocusToUse):
    queryFile = queryFile
    queryProtName = str(countNumberOfLocus)+'.fasta'



    alleleScores, isEmpty, proteinSeqToConcat, alleleNumbers, sameAlleles, prevAlleleNames = getOwnBlastScore(os.path.join(alleleFolder,queryFile), databaseP, queryProtName, str(countNumberOfLocus), blastResultsPath, LocusToUse, queryFile) #Self-BLAST to get Own Score
    writeAlleleScoresToFile(alleleScores,os.path.join(resultsFolder,alleleScoreFile))
    if proteinSeqToConcat != False:
        writeQuerySeqToFile(proteinSeqToConcat,os.path.join(resultsFolder,queryFileWithAll))

    os.remove(queryProtName)

    return alleleNumbers, sameAlleles, countNumberOfLocus, prevAlleleNames, isEmpty


def func_BSR(referenceGenome, countNumberOfGenomes, databaseFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleFolder, alleleScores, LocusToUse):
    isDone, alleleList, bestmatches = gatherAllQueries(alleleFolder) #MUDAR
    if isDone:
        queryFileWithAll = 'Allalleles.fasta'
    else:
        print 'There was an error concatenating the files'
        sys.exit()

    databaseName = str(countNumberOfGenomes)+'_db'

    startTime = datetime.now()

    TimeFile = os.path.join(resultsFolder,'timeResults.txt')
    print "--------------------"
    referenceDatabasePath, referenceCDS, referenceGenomeArray, referenceCDSsequences = CreateReferenceDatabase(os.path.join(databaseFolder,referenceGenome), databaseP, databaseName, str(countNumberOfGenomes)) #Create Database with one genome
    print "_______________"
    resultsList, addNewAlleles = getBlastScoreRatios(alleleScores, alleleList, referenceDatabasePath, os.path.join(resultsFolder,queryFileWithAll), referenceGenomeArray, referenceCDS, bestmatches, referenceCDSsequences, referenceGenome, str(countNumberOfGenomes), blastResultsPath, LocusToUse) #Get BLAST Score Ratio and call alleles
    print resultsList

    INF, EXM, LNF = checkNumberOfResults(resultsList)

    textToFile = 'BLAST of '+ queryFileWithAll + ' against ' + referenceGenome + ': '
    
    WriteTimeResults(textToFile, str(datetime.now() - startTime), TimeFile)

    textToFile = 'INF: '+ str(INF) + '\nEXM: '+ str(EXM) + '\nLNF: '+ str(LNF) + '\n'

    WriteTimeResults(textToFile, False, TimeFile)



    return resultsList, referenceCDSsequences, addNewAlleles


def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:

        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]

def translateSeq(DNASeq):
    seq=DNASeq
    reversedSeq=False
    try:
        myseq= Seq(seq)
        #print myseq
        protseq=Seq.translate(myseq, table=11,cds=True)
    except:
        reversedSeq=True
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
                reversedSeq=False
                try:
                    seq=seq[::-1]                           
                    seq=reverseComplement(seq)
                    myseq= Seq(seq)
                    #print myseq
                    protseq=Seq.translate(myseq, table=11,cds=True)
                except:
                    raise
    return protseq, reversedSeq


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


def runProdigalContig(FASTAfile, numberOfGenomes):

        # ------------ #
        # RUN PRODIGAL #
        # ------------ #  #/scratch/NGStools/prodigal-2.50/prodigal -i /scratch/spneumoniae_600/ERR067978/velvet_assembly/contigs.fa -c -m -g 11 -p single -f sco -q > test_all.txt
    pathResults = str(numberOfGenomes) + "_Presults.txt"
    path = "prodigal -i "+FASTAfile+" -c -m -g 11 -p single -f sco -q > "+ pathResults;
    os.system(path)
    cdsDict = {}
    tempList = []
    line = ' '
    x = open(pathResults,"r")
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

    os.remove(pathResults)

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

def CreateProteomeContig(cdsDict, referenceCDS, referenceGenomeArray, nameDB):
    alleleProt=''
    proteome=""
    countP=0
    countCDS=0
    referenceCDSsequences = {}
    for j in cdsDict:
        countCDS=0
        for i in cdsDict[j]:
            countP+=1
            countCDS+=1
            AlignBegin=i[0]
            AlignEnd=i[1]
            IdCDS = str(j) + '--' + str(countCDS)
            referenceCDS[IdCDS] = [int(AlignBegin),int(AlignEnd)]
                    #print contigTag,value

                        
                        #print contigTag, protein[0], protein[1]
                        #print currentGenomeDict[ contigTag ]
            seq= referenceGenomeArray[j][ int(AlignBegin):int(AlignEnd) ].upper()

            try:
                protseq=translateSeq(seq)
            except:
                continue

            proteome+=">"+IdCDS+"\n"+str(protseq[0])+"\n"

            referenceCDSsequences[IdCDS] = seq
                #proteome+=str(protseq)
    dirName = nameDB.split('_')[0] + '_GenomeProteome.fasta'
    if os.path.isfile(dirName):
        os.remove(dirName)

    with open(dirName, "wb") as v:
        v.write(proteome)
    
    return dirName, referenceCDS, referenceCDSsequences


def CreateReferenceDatabase(referenceFASTA, databasePath, nameDB, numberRefGenomes):
    referenceGenomeArray = {}
    referenceCDS = {}
    referenceGenomeArray = getFASTAarray(referenceFASTA, referenceGenomeArray)
    CDSdict = runProdigalContig(referenceFASTA, numberRefGenomes)
    refProteomePath, referenceCDS, referenceCDSsequences = CreateProteomeContig(CDSdict, referenceCDS, referenceGenomeArray, nameDB)
    referenceDatabase = os.path.join(databasePath,nameDB)
    Create_Blastdb(refProteomePath,1,True, referenceDatabase)

    os.remove(refProteomePath)

    return referenceDatabase, referenceCDS, referenceGenomeArray, referenceCDSsequences


def CreateQueryDatabase(FASTAfile, databasePath,queryProteomeName):
	
	
    gene_fp = HTSeq.FastaReader(FASTAfile)
    names=""
    alleleProt=''
    proteome=""
    isEmpty = True
    countAlleles = 0
    geneNumber = queryProteomeName.split('.')[0]
    proteinsToQueryFile = {} 
    queryAlleleList = {}
    prevAlleleName = {}
    for allele in gene_fp:
        try:
            x = translateSeq(allele.seq)
            countAlleles+=1
            isEmpty = False
        except:
            #print 'Could not translate'
            if countAlleles==0:
                isEmpty = True
            continue
        alleleProt+=">"+ str(geneNumber) + '--' + str(countAlleles)+"\n"+str(x[0])+"\n"
        proteome+=">"+ str(geneNumber) + '--' + str(countAlleles)+"\n"+str(x[0])+"\n"
        proteinsToQueryFile[str(geneNumber) + '--' + str(countAlleles)] = str(x[0])
        queryAlleleList[str(geneNumber) + '--' + str(countAlleles)] = str(allele.seq)
        prevAlleleName[allele.name] = str(geneNumber) + '--' + str(countAlleles)
    # with open(pathRef+'allAllelesAA.fasta', "wb") as f:
    #     f.write(alleleProt)
    databasePath = os.path.join(databasePath,queryProteomeName)
    databasePath = databasePath.split('.')[0]
    databasePath = databasePath+'_db'
    with open(queryProteomeName, "wb") as v:
        v.write(proteome)
    Gene_Blast_DB_name = Create_Blastdb(queryProteomeName,1,True, databasePath)


    if len(proteinsToQueryFile) == 0:
        proteinsToQueryFile = False

    return databasePath, isEmpty, proteinsToQueryFile, queryAlleleList, prevAlleleName



# def getOwnBlastScore(FASTAfile, databasePath, queryProteomeName):

#     databasePath, isEmpty, proteinsToQueryFile, queryAlleleList = CreateQueryDatabase(FASTAfile, databasePath,queryProteomeName)


#     allelescores = []

#     if isEmpty:
#         return allelescores, isEmpty, proteinsToQueryFile

#     blast_out_file = 'BLASTresults.xml'

#     cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=5, num_alignments=7000, num_descriptions=7000)

#     print 'OwnBLAST:'
#         #print cline
#     allelescore=0

#     blast_records = runBlastParser(cline,blast_out_file, False)

#     startTime = datetime.now()

#     allelescores = parseOwnBLASTRecordsXML(blast_records)

#     print 'CheckResults:' + str(datetime.now() - startTime)


#     os.remove(databasePath+ ".pin")
#     os.remove(databasePath+ ".phr")
#     os.remove(databasePath+ ".psq")
#     os.remove(databasePath+ "_blast.log")

#     return allelescores, isEmpty, proteinsToQueryFile


def getOwnBlastScore(FASTAfile, databasePath, queryProteomeName, numberOfLocus, blastResultsPath, LocusToUse, queryFile):
	
	
    print FASTAfile
    allelescores = []
    alleleNumbers = {}
    sameAlleles = {}
    prevAlleleName = {}

    databasePath, isEmpty, proteinsToQueryFile, queryAlleleList, prevAlleleName = CreateQueryDatabase(FASTAfile, databasePath,queryProteomeName)


    if isEmpty:
        return allelescores, isEmpty, proteinsToQueryFile, alleleNumbers, sameAlleles, prevAlleleName
	
    print blastResultsPath
    print numberOfLocus
	
    blast_out_file = blastResultsPath + numberOfLocus + '_BLASTresults.xml'

    cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=5, num_alignments=7000, num_descriptions=7000)

        #print cline
    allelescore=0

    blast_records = runBlastParser(cline,blast_out_file, False)


    allelescores, alleleList, alleleNumbers, sameAlleles = parseOwnBLASTRecordsAndDuplicates(blast_records, FASTAfile, queryAlleleList)

    proteinsToQueryFile = translateAlleleList(alleleList, queryProteomeName, LocusToUse, queryFile)


    os.remove(databasePath+ ".pin")
    os.remove(databasePath+ ".phr")
    os.remove(databasePath+ ".psq")
    os.remove(databasePath+ "_blast.log")
    os.remove(blast_out_file)

    return allelescores, isEmpty, proteinsToQueryFile, alleleNumbers, sameAlleles, prevAlleleName


def translateAlleleList(alleleList, queryProteomeName, LocusToUse, FASTAfile):
    
    if LocusToUse[FASTAfile] == False:
        return False
    
    proteinsToQueryFile = {}
    geneNumber = queryProteomeName.split('.')[0]
    countAlleles = 0
    for i in alleleList:
        x = translateSeq(alleleList[i])
        countAlleles+=1

        proteinsToQueryFile[str(geneNumber) + '--' + str(countAlleles)] = str(x[0])

    return proteinsToQueryFile



def getBlastScoreRatios(allelescore, alleleList, databasePath, queryProteomeName, referenceGenomeArray, referenceCDS, bestmatches, referenceCDSsequences,referenceFileName, countNumberOfGenomes, blastResultsPath, LocusToUse):
    
    alleleProt=''
    proteome=""
    countP=0
    countCDS=0

    blast_out_file = blastResultsPath + countNumberOfGenomes + '_BLASTresults.xml'


    cline = NcbiblastpCommandline(query=queryProteomeName, db=databasePath, out=blast_out_file, outfmt=5, num_alignments=7000, num_descriptions=7000)

        #print cline
    blast_records = runBlastParser(cline,blast_out_file, False)

    resultsList, addNewAlleles = parseBLASTRecordsXML(blast_records, allelescore, alleleList, referenceGenomeArray, referenceCDS, bestmatches, referenceCDSsequences, referenceFileName, LocusToUse)

    os.remove(blast_out_file)

    return resultsList, addNewAlleles


def parseOwnBLASTRecordsXML(blast_records):

    allelescores = {}

    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                #print blast_record.query, alignment.hit_def
                for match in alignment.hsps:

                    if alignment.hit_def == blast_record.query:
                        allelescores[alignment.hit_def] = int(match.score)

            else:
                break

    return allelescores

def parseOwnBLASTRecordsAndDuplicates(blast_records, FASTAfile, queryAlleleList):

    allelescores = {}
    alleleList = {}


    for blast_record in blast_records:
        found=False 
        for alignment in blast_record.alignments:
            if found is False:
                #print blast_record.query, alignment.hit_def
                for match in alignment.hsps:

                    if alignment.hit_def == blast_record.query:
                        allelescores[alignment.hit_def] = int(match.score)
                        alleleList[alignment.hit_def] = queryAlleleList[alignment.hit_def]

            else:
                break

    allelescores, alleleList, sameAlleles = checkForDuplicates(allelescores, alleleList)

    Done, alleleNumbers = writeAlleleFileWithoutDuplicates(alleleList, FASTAfile)



    return allelescores, alleleList, alleleNumbers, sameAlleles


def checkForDuplicates(allelescores, alleleList):

    newAlleleScores = {}
    newAlleleList = {}
    inserted = {}
    sameAlleles = {}

    for allele in allelescores:
        alleleList, inserted, newAlleleScores, newAlleleList, sameAlleles = compareLists(alleleList,allelescores, inserted, newAlleleScores, newAlleleList, allele, sameAlleles)

    return newAlleleScores, newAlleleList, sameAlleles


def compareLists(alleleList,allelescores, inserted, newAlleleScores, newAlleleList, allele, sameAlleles):

    toDelete = []
    for toCompare in alleleList:
            try:
                if allelescores[allele] == allelescores[toCompare] and alleleList[allele] == alleleList[toCompare]:
                    try:
                        if inserted[allele]:
                            toDelete.append(toCompare)
                            sameAlleles[allele].append(toCompare)
                    except KeyError:
                        inserted[allele] = True
                        sameAlleles[allele] = []
                        newAlleleScores[allele] = allelescores[allele]
                        newAlleleList[allele] = alleleList[allele]
            except KeyError:
                continue

    for j in toDelete:
        del alleleList[j]

    return alleleList, inserted, newAlleleScores, newAlleleList, sameAlleles


def parseOwnBLASTRecordsTAB(blast_records):

    allelescores = {}
    queryDef = ''

    for qresult in blast_records:
        for match in qresult:
            continue

    return allelescores


def parseBLASTRecordsXML(blast_records, allelescores, alleleList, referenceGenomeArray, referenceCDS, bestmatches, referenceCDSsequences, referenceFileName, LocusToUse):

    resultList={}
    addNewAlleles = False

    count = 0
    countInferred = 0

    for blast_record in blast_records:
                
            
            for alignment in blast_record.alignments:
                
                for match in alignment.hsps:

                    cdsStrName=str(blast_record.query)
                    ToSearchOnMatches = cdsStrName.split('--')[0] #Locus name

                    scoreRatio=float(match.score)/float(allelescores[cdsStrName])
                    queryDNAstr=alleleList[cdsStrName]
                    hitDef = str(alignment.hit_def)
                    referenceDNAstr=referenceCDSsequences[hitDef]
                    compare=False
                    if queryDNAstr==referenceDNAstr is False:
                        try:
                            queryDNAstr=reverseComplement(queryDNAstr)
                            if queryDNAstr==referenceDNAstr is False:
                                pass
                            else:
                                compare=True
                        except:
                            pass
                    else:
                        compare=True
                    
                    if "N" in queryDNAstr or "K" in queryDNAstr or "R" in queryDNAstr:
                        pass
                    
                    elif(scoreRatio == 1 and bestmatches[ToSearchOnMatches][2] is False and compare is True):
                        #print 'EXISTS'
                        queryFile = bestmatches[ToSearchOnMatches][9]
                        AllelesOnLocus = bestmatches[ToSearchOnMatches][8]
                        bestmatches[ToSearchOnMatches]=[match.score,scoreRatio,True,cdsStrName,hitDef,match,len(queryDNAstr),blast_record.query_letters, AllelesOnLocus, queryFile]
                        #print alignment
                        #print match

                    elif(scoreRatio == 1 and match.score>bestmatches[ToSearchOnMatches][0] and compare is True):
                        queryFile = bestmatches[ToSearchOnMatches][9]
                        AllelesOnLocus = bestmatches[ToSearchOnMatches][8]
                        bestmatches[ToSearchOnMatches]=[match.score,scoreRatio,True,cdsStrName,hitDef,match,len(queryDNAstr),blast_record.query_letters, AllelesOnLocus, queryFile]
                        #print match

                    elif(scoreRatio == 1 and bestmatches[ToSearchOnMatches][2] is False and compare is False):
                        queryFile = bestmatches[ToSearchOnMatches][9]
                        AllelesOnLocus = bestmatches[ToSearchOnMatches][8]
                        bestmatches[ToSearchOnMatches]=[match.score,scoreRatio,False,cdsStrName,hitDef,match,len(queryDNAstr),blast_record.query_letters, AllelesOnLocus, queryFile]
                        #print alignment
                        #print match
                    elif(scoreRatio == 1 and match.score>bestmatches[ToSearchOnMatches][0] and compare is False):
                        queryFile = bestmatches[ToSearchOnMatches][9]
                        AllelesOnLocus = bestmatches[ToSearchOnMatches][8]
                        bestmatches[ToSearchOnMatches]=[match.score,scoreRatio,False,cdsStrName,hitDef,match,len(queryDNAstr),blast_record.query_letters, AllelesOnLocus, queryFile]
                        #print match
                    elif(match.score>bestmatches[ToSearchOnMatches][0] and scoreRatio>0.6 and scoreRatio>bestmatches[ToSearchOnMatches][1] and bestmatches[ToSearchOnMatches][2] is False):
                        
                        #print 'AQUI'
                        queryFile = bestmatches[ToSearchOnMatches][9]
                        AllelesOnLocus = bestmatches[ToSearchOnMatches][8]
                        count += 1
                        #print match.query
                        #print match.sbjct
                        #print allelescores
                        bestmatches[ToSearchOnMatches]=[match.score,scoreRatio,False,cdsStrName,hitDef,match,len(queryDNAstr),blast_record.query_letters, AllelesOnLocus, queryFile]
                        #print match
    
    #print bestmatches
    for i in bestmatches:
        
        if LocusToUse[bestmatches[i][9]] == False:
            resultList[i] = ['LNIU',False,False,False,False,False,False] #Locus not in Use
            #print "Locus not in use"
        
        #print ("Classifying the match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))        
        elif bestmatches[i][0]==0 or "N" in referenceDNAstr or "K" in referenceDNAstr or "R" in referenceDNAstr :
                     #if no best match was found
                    
                     ###################
                     # LOCUS NOT FOUND #
                     ###################
            if  bestmatches[i][0]==0:        
                resultList[i] = ['LNF',False,False,False,False,False,False]
                #printinfo(genomeFile,geneFile)
                #print "Locus " + str(bestmatches[i][3].split('--')[0]) + "not found, no matches \n"
            else:
                resultList[i] = ['LNF',False,False,False,False,False,False]
                #print "Locus has strange base (N, K or R) \n"
        
        
        elif bestmatches[i][2] is True:
            referenceCDSName=bestmatches[i][4]
            queryName = bestmatches[i][3]
            
            matchLocation=referenceCDS[referenceCDSName]
            #matchLocation=matchLocation.split("-")
            
            referenceStr = referenceCDSsequences[referenceCDSName]
            protSeq, Reversed = translateSeq(referenceStr)
                    #if a perfect match was found
                    
                        ################################################
                        # EXACT MATCH --- MATCH == GENE --- GENE FOUND #
                        ################################################
                    
            if not Reversed:
                resultList[i] = ['EXM',queryName,str(matchLocation),"+", referenceCDSName,False]
            else:
                resultList[i] = ['EXM',queryName,str(matchLocation),"-", referenceCDSName,False]

            #print 'Exact match for ' + queryName
        
        else:
            
            match=bestmatches[i][5]
            #print match
            QueryalleleLen=bestmatches[i][6]
            referenceCDSName=bestmatches[i][4]

            contigName = referenceCDSName.split('--')[0]

            alleleName=bestmatches[i][3] 
            
            matchLocation=[match.sbjct_start, match.sbjct_end]
            
            seq = referenceGenomeArray[contigName]

            bestMatchRefLen=len(seq)
            
            alleleStr = alleleList[bestmatches[i][3]]
            protSeq, Reversed = translateSeq(alleleStr)
            
            
    #         print match
            # print bestmatches[i][1]
            # print match.match
            #print matchLocation
            # print bestMatchAlleleLen
            
            # get extra space to the right and left between the contig and match  
            leftmatchContig = referenceCDS[referenceCDSName][0] + (int(matchLocation[0])*3)
            rightmatchContig = bestMatchRefLen + referenceCDS[referenceCDSName][0] + (int(matchLocation[1])*3) 
            
            if Reversed:
                aux=rightmatchContig
                rightmatchContig=leftmatchContig
                leftmatchContig=aux
    #         """else:
    #             rightmatchContig=bestMatchContigLen-int(matchLocation[0])   
    #             leftmatchContig=int(matchLocation[1])"""
            
            #print leftmatchContig, rightmatchContig
            
            
            # get extra space to the right and left between the allele and match
            leftmatchAllele = (int(match.query_start)*3-2)
            rightmatchAllele = QueryalleleLen - (int(match.query_end)*3)   
            
            if Reversed: 
                aux=rightmatchAllele
                rightmatchAllele=leftmatchAllele
                leftmatchAllele=aux
                
    #         """else:
    #             rightmatchAllele=geneLen-(int(match.sbjct_start )*3)
    #             leftmatchAllele=(int(match.sbjct_end)*3)"""
            
            #if a best match was found but it's not an exact match  

                     ###########################
                     # LOCUS ON THE CONTIG TIP #
                     ###########################
            
            #print leftmatchAllele, rightmatchAllele
    #         print geneLen
            
            
            
    #         #if bestMatchContigLen <= geneLen:
            if leftmatchContig < leftmatchAllele and  rightmatchContig < rightmatchAllele:
            
                resultList[i] = ['LOTSC',False,False,False, False,False]

                #print "Locus is bigger than the contig \n"
            
    #         #if match.query_start ==1 and len(match.query) < geneLen:       
            elif leftmatchContig<leftmatchAllele:
                
                resultList[i] = ['LOT3',False,False,False, False,False]
           
                #print "Locus is on the 3' tip of the contig \n"
            
            
    #        #elif match.query_end == bestMatchContigLen and len(match.query) < geneLen:
            elif rightmatchContig < rightmatchAllele:
                
                resultList[i] = ['LOT5',False,False,False, False,False]

                #print "Locus is on the 5' tip of the contig \n"
            
                    
            else:
                         #######################
                         # ADD INFERRED ALLELE #     # a new allele 
                         #######################
                        
                                                
                     #print "infered allele has location : "+(CDSType)
                     #printinfo(genomeFile,geneFile) 
                addNewAlleles = True

                tagAux='INF'
                
                alleleName = 'allele_from_' + str(referenceFileName)

                #print bestmatches[i][3]
                
                referenceCDSName = bestmatches[i][4]
                locationToPass = str(referenceCDS[referenceCDSName][0]) + ',' + str(referenceCDS[referenceCDSName][1])

                if not Reversed:
                    resultList[i] = ['INF',alleleName,str(locationToPass),"+",referenceCDSName, bestmatches[i][9], False]
                else:
                    resultList[i] = ['INF',alleleName,str(locationToPass),"-",referenceCDSName, bestmatches[i][9], False]

                             
                #print "New allele! Adding allele "+ alleleName +" to the database\n"
                                                                                    

    return resultList , addNewAlleles          
                



def parseBLASTRecordsTAB(blast_records, allelescore, queryDef, referenceGenomeArray, referenceCDS):

    ToNewAllele = []

    for qresult in blast_records:
        for match in qresult:
            continue

    return ToNewAllele


def gatherAllQueries(queryPath):

    alleleList = {}
    bestmatches = {}

    try:
        queryFilesOnDir = [ f for f in listdir(queryPath) if isfile(join(queryPath,f)) ]
        countFiles = 0
        for queryFile in queryFilesOnDir:
            AllqueryFile = os.path.join(queryPath,queryFile)
            if queryFile == 'Allalleles.fasta':
                continue

            countFiles += 1
            g_fp = HTSeq.FastaReader(os.path.join(queryPath,queryFile))
            countAlleles = 0
            for allele in g_fp:
                countAlleles += 1
                #ToWrite.append(">" + str(countFiles) + '--' + str(countAlleles) +"\n"+ str(allele.seq).upper() + "\n")
                alleleList[str(countFiles) + '--' + str(countAlleles)] = str(allele.seq).upper()

            bestmatches[str(countFiles)] = [0,0,False,'','','',0,'', str(countAlleles), AllqueryFile] #To be used when searching for new alleles. On instance for each locus
                                            #Score, ScoreRatio, Found, queryName, HitName, MatchObject, lengthReference, lengthQuery, numberOfExistingAllelesForThatLocus

        #CreateNewAlleleFile(os.path.join(queryPath,'Allalleles.fasta'), ToWrite)
    except Exception:
        print 'An error occurred'
        return False , alleleList, bestmatches
    
    return True , alleleList, bestmatches


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


def writeAlleleScoresToFile(alleleScores,fileToWrite):
    fg = open(fileToWrite, 'a')
    countAlleles = 0
    for allele in alleleScores:
        countAlleles += 1
        locusNumber = allele.split('--')[0]
        fg.write(str(locusNumber)+'--'+str(countAlleles) + ',' + str(alleleScores[allele]) + '\n')
    fg.close()


def readCSValleleScores(CSVpath):

    allelescores = {}

    with open(CSVpath, mode='r') as infile:
        reader = csv.reader(infile)
        for rows in reader:
            allelescores[rows[0]] = rows[1]
    
    return allelescores


def writeQuerySeqToFile(proteinToAdd, fileToWrite):
    fg = open(fileToWrite, 'a')
    for protein in proteinToAdd:
        fg.write(">" + protein + '\n' + str(proteinToAdd[protein]) + '\n')
    fg.close()



def TowriteToAlleleFiles(dictOfResults, dictOfCDS):
    writeToFile = {}
    for i in dictOfResults[0]:
        if dictOfResults[0][i][0] == 'INF':

            referenceCDSName = dictOfResults[0][i][4]
            try:
                writeToFile[i].append(['>' + dictOfResults[0][i][1] + '\n' + dictOfCDS[0][referenceCDSName] + '\n', dictOfResults[0][i][5], i])
            except KeyError:
                writeToFile[i] = []
                writeToFile[i].append(['>' + dictOfResults[0][i][1] + '\n' + dictOfCDS[0][referenceCDSName] + '\n', dictOfResults[0][i][5], i])

    return writeToFile



def NewAlleleFiles(writeToFile, queryPath):

    for i in writeToFile:
        fG = open(os.path.join(queryPath, i[1]),'a')
        break

    for i in writeToFile:
        fG.write(i[0])

    fG.close()

    return True



def writeAlleleFileWithoutDuplicates(alleleList, FASTAfile):
    
    fG = open(FASTAfile,'w')
    countAlleles = 0

    alleleNumbers = {}

    #it = iter(sorted(alleleList.iteritems()))
    for i in iter(sorted(alleleList.iteritems())):
        countAlleles += 1
        alleleNumbers[i[0]] = i[0].split('--')[0] + '--' + str(countAlleles)
        #locusName = i.split('--')[0]
        fG.write('>'+ i[0].split('--')[0] + '--' + str(countAlleles) + '\n' + i[1] + '\n')

    fG.close()
    return True, alleleNumbers


def WriteTimeResults(whatToWrite, time, filePath):
    fG = open(filePath,'a')
    if time == False:
        fG.write(whatToWrite)
    else:
        fG.write(whatToWrite + time + '\n')

    fG.close()


def checkNumberOfResults(resultsList):
    INF = 0
    EXM = 0
    LNF = 0

    for i in resultsList:
        if resultsList[i][0] == 'INF':
            INF += 1
        elif resultsList[i][0] == 'EXM':
            EXM += 1
        elif resultsList[i][0] == 'LNF':
            LNF += 1

    return INF, EXM, LNF


def createPhylovizFile(arrayOfResults, GenomesID, LocusID, resultsFolder):
    fG = open(os.path.join(resultsFolder, 'inputPhyloviz.tab'), 'w')

    lineToWrite = 'FILE\t'
    profileOrder = []

    for i in arrayOfResults[1][0]:
        profileOrder.append(i)
        lineToWrite += LocusID[i] + '\t'

    lineToWrite = lineToWrite.strip('\t')
    fG.write(lineToWrite+ '\n')


    for i in GenomesID:
        lineToWrite = GenomesID[i] + '\t'
        for j in profileOrder:
            if arrayOfResults[i][0][j][0] == 'EXM' or arrayOfResults[i][0][j][0] == 'INF':
                lineToWrite += arrayOfResults[i][0][j][1].split('--')[1] + '\t'

            elif arrayOfResults[i][0][j][0] == 'LNF':
                lineToWrite += 'LNF' + '\t'

        lineToWrite = lineToWrite.strip('\t')

        fG.write(lineToWrite+ '\n')

    fG.close()


def checkLocusToUse(arrayOfResults, LocusID, countIterations):
    LocusToUse = {}

    for k in LocusID:
        LocusToUse[LocusID[k]] = False
    
    for i in arrayOfResults:
        for j in arrayOfResults[i][0]:
            if countIterations > 1:
                if arrayOfResults[i][0][j][0] == 'INF':
                    LocusToUse[LocusID[j]] = True
            else:
                if arrayOfResults[i][0][j][0] == 'LNF':
                    LocusToUse[LocusID[j]] = True

    return LocusToUse


def correctAlleles(arrayOfResults, dictOfChanges, genomeID):


    newArrayOfResults = {}
    results = {}


    for i in arrayOfResults[0]:
        newArrayOfResults[i] = arrayOfResults[0][i]



    for j in dictOfChanges:
        
        if dictOfChanges[j] == False:  #In case of errors during translation
            continue
        if newArrayOfResults[j][0] == 'LNF' or newArrayOfResults[j][0] == 'LNIU':
            continue

        alleleName = dictOfChanges[j][2][newArrayOfResults[j][1]]


        dictOfSameAlleles = dictOfChanges[j][1] 
        for allele in dictOfSameAlleles:
            for sameAllele in dictOfSameAlleles[allele]:
                if sameAllele == alleleName:
                    newArrayOfResults[j][1] = allele



    for j in dictOfChanges:
        if dictOfChanges[j] == False:  #In case of errors during translation
            continue
        if newArrayOfResults[j][0] == 'LNF':
            continue
        try:
            alleleName = dictOfChanges[j][2][newArrayOfResults[j][1]]
        except KeyError:
            alleleName = newArrayOfResults[j][1]
        
        CorrectAlleleNumbers = dictOfChanges[j][0] 
        for allele in CorrectAlleleNumbers:
            if allele == alleleName:
                #print 'PrevName: ' + alleleName
                newArrayOfResults[j][1] = CorrectAlleleNumbers[allele]
                #print 'New Name: ' + newArrayOfResults[j][1]



    return newArrayOfResults, genomeID





