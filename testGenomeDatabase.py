from BSRfunctions import CreateReferenceDatabase
from BSRfunctions import getOwnBlastScore
from BSRfunctions import getBlastScoreRatios
from BSRfunctions import CreateNewAlleleFile
from BSRfunctions import concatAllQueries
from BSRfunctions import writeFilteredAlleleFile
import os
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys


def main():


	parser = argparse.ArgumentParser(description="This program calls for alleles in a genome database")
	parser.add_argument('-a', nargs='?', type=str, help="folder with single allele files", required=True)
	parser.add_argument('-d', nargs='?', type=str, help='genome to be used as database', required=True)
	parser.add_argument('-p', nargs='?', type=str, help="Database Path", required=True)
	#parser.add_argument('-c', nargs='?', type=int, help='Number of BSR and nucleotide similarity calculations.', required=True)
	parser.add_argument('-r', nargs='?', type=str, help="Results folder", required=True)
	parser.add_argument('-l', nargs='?', type=int, help="Max BP length", required=False)
	parser.add_argument('-f', nargs='?', type=int, help="Max alleles per gene", required=False)
	parser.add_argument('-x', nargs='?', type=str, help="Is XML output format (If True)", required=True)
	parser.add_argument('-c', nargs='?', type=str, help="Concat queries into one file (if True)", required=True)

	args = parser.parse_args()


	alleleFolder=args.a
	databaseGenome = args.d
	resultsFolder=args.r
	databaseP = args.p
	maxBP = args.l
	maxalleles = args.f
	ToConcat = args.c
	isXML = args.x

	if not os.path.isdir(databaseP):
		os.makedirs(databaseP)

	if not os.path.isdir(resultsFolder):
		os.makedirs(resultsFolder)

	queryProtName = 'Queryproteome.fasta'

	if ToConcat == 'True':
		isDone = concatAllQueries(alleleFolder, maxBP, maxalleles)
		if isDone:
			queryFilesOnDir = ['Allalleles.fasta']
		else:
			sys.exit()
	else:
		queryFilesOnDir = [ f for f in listdir(alleleFolder) if isfile(join(alleleFolder,f)) ]

	startTime = datetime.now()

	referenceDatabasePath, referenceCDS, referenceGenomeArray = CreateReferenceDatabase(databaseGenome, databaseP)


	for queryFile in queryFilesOnDir:
		startTimeBLAST = datetime.now()
		if queryFile == 'Allalleles.fasta' and ToConcat != 'True':
			continue

		if maxalleles != None and ToConcat != 'True':
			queryFileProv = 'QueryProv.fasta'
			writeFilteredAlleleFile(os.path.join(alleleFolder,queryFile), os.path.join(alleleFolder,queryFileProv), maxalleles)
			alleleScores, hitDef, isEmpty = getOwnBlastScore(os.path.join(alleleFolder,queryFileProv), databaseP, queryProtName, isXML)

			if isEmpty:
				print 'There was an error translating the file'
			else:
				ToNewAllele = getBlastScoreRatios(os.path.join(alleleFolder,queryFileProv), alleleScores, hitDef, referenceDatabasePath, queryProtName, referenceGenomeArray, referenceCDS, isXML)
				if len(ToNewAllele)>0:
					print len(ToNewAllele)
					#CreateNewAlleleFile(os.path.join(alleleFolder,queryFile), ToNewAllele)
			print 'End of BSR for query File: ' + str(datetime.now() - startTimeBLAST)
			os.remove(os.path.join(alleleFolder,queryFileProv))

		else:
			alleleScores, hitDef, isEmpty = getOwnBlastScore(os.path.join(alleleFolder,queryFile), databaseP, queryProtName, isXML)
			if isEmpty:
				print 'There was an error translating the file'
			else:
				ToNewAllele = getBlastScoreRatios(os.path.join(alleleFolder,queryFile), alleleScores, hitDef, referenceDatabasePath, queryProtName, referenceGenomeArray, referenceCDS, isXML)
				if len(ToNewAllele)>0:
					print len(ToNewAllele)
					#CreateNewAlleleFile(os.path.join(alleleFolder,queryFile), ToNewAllele)
			print 'End of BSR for query File: ' + str(datetime.now() - startTimeBLAST)



	print 'Total:' + str(datetime.now() - startTime)

if __name__ == "__main__":
    main()