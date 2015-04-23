
from BSRfunctions import SelfBLAST
from BSRfunctions import BSR
from BSRfunctions import check_ToWrite
from BSRfunctions import CreateNewAlleleFiles
from BSRfunctions import readCSValleleScores
from BSRfunctions import WriteTimeResults
from BSRfunctions import createPhylovizFile
from BSRfunctions import checkLocusToUse
from BSRfunctions import get_correctAlleles
import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys
from multiprocessing import Pool


def main():

	startTime = datetime.now()


	parser = argparse.ArgumentParser(description="This program calls for alleles in a genome database")
	parser.add_argument('-a', nargs='?', type=str, help="folder with single allele files", required=True)
	parser.add_argument('-d', nargs='?', type=str, help='folder with genomes to be used as reference database', required=True)
	parser.add_argument('-p', nargs='?', type=str, help="Database Path", required=True)
	#parser.add_argument('-c', nargs='?', type=int, help='Number of BSR and nucleotide similarity calculations.', required=True)
	parser.add_argument('-r', nargs='?', type=str, help="Results folder", required=True)
	parser.add_argument('-x', nargs='?', type=str, help="Is XML output format (If True)", required=True)
	parser.add_argument('-v', nargs='?', type=str, help="Folder for BLAST results", required=True)


	args = parser.parse_args()


	alleleFolder=args.a
	databaseFolder = args.d
	resultsFolder=args.r
	databaseP = args.p
	isXML = args.x
	blastResultsPath = args.v

	TimeFile = os.path.join(resultsFolder,'timeResults.txt')

	fg = open(TimeFile, 'w')
	fg.close()

	countIterations = 0
	LocusToUse = 'All'
	#getProfile = True
	finalResults = {}
	arrayOfResults = {}


	def runIteration(alleleFolder, databaseFolder, resultsFolder, databaseP, isXML, runFirstBLAST, blastResultsPath, countIterations, LocusToUse, finalResults, arrayOfResults):


		countIterations += 1

		print str(countIterations) + ' Iteration'

		if not os.path.isdir(databaseP):
			os.makedirs(databaseP)

		if not os.path.isdir(resultsFolder):
			os.makedirs(resultsFolder)

		if not os.path.isdir(blastResultsPath):
			os.makedirs(blastResultsPath)

		alleleScoreFile = 'alleleScores.txt'
		queryFileWithAll = 'Allalleles.fasta'

		TimeFile = os.path.join(resultsFolder,'timeResults.txt')

		fg = open(os.path.join(resultsFolder,alleleScoreFile), 'w')
		fg.close()

		fg = open(os.path.join(resultsFolder,queryFileWithAll), 'w')
		fg.close()


		if runFirstBLAST == True:
		
			####Self-BLAST (To be parallelized by locus file)

			queryFilesOnDir = [ f for f in listdir(alleleFolder) if isfile(join(alleleFolder,f)) ]
			countQuery = 0
			try:
				indexToRemove = queryFilesOnDir.index('Allalleles.fasta')
				del queryFilesOnDir[indexToRemove]
			except ValueError:
				continue

			print "Self-BLAST and check for duplicates"

			poolJobs = Pool()
			startTimeBLAST = datetime.now()

			dictOfChanges = {}


			job_args = []
			countNumberOfLocus = 0
			LocusID = {}
			Tocheck = False
			if LocusToUse == 'All':
				LocusToUse = {}
				Tocheck = True

			for queryFile in queryFilesOnDir:
				if Tocheck == True:
					LocusToUse[queryFile] = True
				countNumberOfLocus += 1
				LocusID[str(countNumberOfLocus)] = queryFile
				job_args.append((queryFile, countNumberOfLocus, alleleFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleScoreFile, LocusToUse))

			it = poolJobs.imap(SelfBLAST, job_args)

			while 1:
				try:
					x = it.next()
					if x[4] == False:	
						dictOfChanges[str(x[2])] = [x[0], x[1], x[3]] #correctAlleleNumber, sameAlleles, prevAlleleNames  VER O QUE ESTA MAL
					else:
						dictOfChanges[str(x[2])] = False
				except StopIteration:
					break

			poolJobs.close()
			poolJobs.join()

			print 'End of Self-BLAST: ' + str(datetime.now() - startTimeBLAST)


		job_args = []

		###############################################VER###############################

		if countIterations >= 2:

			poolJobs = Pool()
			
			for i in arrayOfResults:
					job_args.append((arrayOfResults[i], dictOfChanges, i))

			it = poolJobs.imap(get_correctAlleles, job_args)

			while 1:
				try:
					res = it.next()
					if countIterations == 2:
						finalResults[res[1]] = []
						finalResults[res[1]].append(res[0])
					else:
						for locus in finalResults[res[1]][0]:
							if finalResults[res[1]][0][locus][0] == 'LNF' and (res[0][locus][0] == 'INF' or res[0][locus][0] == 'EXM'):
								finalResults[res[1]][0][locus] = res[0][locus]

				# except KeyError:
				# 	res = it.next()
				# 	print res[0]
				except StopIteration:
					break

			poolJobs.close()
			poolJobs.join()


		alleleScores = {}
		alleleList = {}
		arrayOfResults = {} 
		arrayOfGenomes = {}
		addNewAllelesArray = []
		writeToFile = {}


		####Allele-Calling (To be parallelized by reference genome)

		alleleScores = readCSValleleScores(os.path.join(resultsFolder,alleleScoreFile))

		referenceGenomesOnDir = [ f for f in listdir(databaseFolder) if isfile(join(databaseFolder,f)) ] #Check for all genomes in reference Database folder

		print "BSR"

		startTimeBLAST = datetime.now()
		
		countReference = 0

		job_args = []
		countNumberOfGenomes = 0
		GenomesID = {}
		for referenceGenome in referenceGenomesOnDir:
			countNumberOfGenomes += 1
			GenomesID[countNumberOfGenomes] = referenceGenome
			job_args.append((referenceGenome, countNumberOfGenomes, databaseFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleFolder, alleleScores, LocusToUse))

		poolJobs = Pool()

		startTimeBLAST = datetime.now()
		textToFile = '\nBegin of BSR at iteration ' + str(countIterations) + '\n'
		textToFile += '----------------------------------------\n'

		WriteTimeResults(textToFile, False, TimeFile)

		it = poolJobs.imap(BSR, job_args)

		while 1:
			try:
				x = it.next() # Results from BSR - [resultList, GenomeCDS, GenomeNumber]
				arrayOfResults[x[2]] = []
				arrayOfResults[x[2]].append(x[0])
				arrayOfGenomes[x[2]] = []
				arrayOfGenomes[x[2]].append(x[1])
				addNewAllelesArray.append(x[3])		
			except StopIteration:
				break

		textToFile = 'End of BSR at: '
		WriteTimeResults(textToFile, str(datetime.now() - startTime), TimeFile)
		
		poolJobs.close()
		poolJobs.join()

		addNewAlleles = False
		for i in addNewAllelesArray:
			if i == True:
				addNewAlleles = True
				break

		print 'End of BSR: ' + str(datetime.now() - startTimeBLAST)

	    ###(To be parallelized by instance on arrayOfResults)
		if addNewAlleles == True:
			job_args = []
			ToWriteArray = []
			for i in arrayOfResults:
				job_args.append((arrayOfResults[i], arrayOfGenomes[i]))

			poolJobs = Pool()

			it = poolJobs.imap(check_ToWrite, job_args)

			ToWriteDict = {}

			while 1:
				try:
					x = it.next() # Gather what to write
					for i in x:
						try: 
							ToWriteDict[x[i][0][2]].append([x[i][0][0], x[i][0][1]])
						except KeyError:
							ToWriteDict[x[i][0][2]] = []
							ToWriteDict[x[i][0][2]].append([x[i][0][0], x[i][0][1]])
				except StopIteration:
					break

			poolJobs.close()
			poolJobs.join()

				

		###(To be parallelized by instance on writeToFileArray)
		if addNewAlleles == True:
			job_args = []
			for locus in ToWriteDict:
				job_args.append((ToWriteDict[locus], alleleFolder))

			poolJobs = Pool()

			it = poolJobs.imap(CreateNewAlleleFiles, job_args)

			while 1:
				try:
					it.next() # Write new allele files
				except StopIteration:
					break

			poolJobs.close()
			poolJobs.join()


		if addNewAlleles == True:
			LocusToUse = checkLocusToUse(arrayOfResults, LocusID, countIterations)
			runIteration(alleleFolder, databaseFolder, resultsFolder, databaseP, isXML, True, blastResultsPath, countIterations, LocusToUse, finalResults, arrayOfResults)
		
		# elif getProfile == True:
		# 	getProfile = False
		# 	LocusToUse = 'All'
		# 	runIteration(alleleFolder, databaseFolder, resultsFolder, databaseP, isXML, True, blastResultsPath, countIterations, LocusToUse, getProfile)
		else:
			print 'Creating PhyloviZ input file...'

			createPhylovizFile(finalResults, GenomesID, LocusID, resultsFolder)

			print 'DONE!'
			print 'Total:' + str(datetime.now() - startTime)
			sys.exit()

	runIteration(alleleFolder, databaseFolder, resultsFolder, databaseP, isXML, True, blastResultsPath, countIterations, LocusToUse, finalResults, arrayOfResults)

if __name__ == "__main__":
	countNumberOfLocus = 0
	alleleScoreFile = 'alleleScores.txt'
	queryFileWithAll = 'Allalleles.fasta'
	main()

