
from BSRfunctions import readCSValleleScores
from BSRfunctions import WriteTimeResults
from BSRfunctions import createPhylovizFile
from BSRfunctions import checkLocusToUse
import os
import shutil
from os import listdir
from os.path import isfile, join
import argparse
from datetime import datetime
import sys


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

			print "Self-BLAST and check for duplicates"

			startTimeBLAST = datetime.now()

			dictOfChanges = {}


			job_args = []
			countNumberOfLocus = 0
			LocusID = {}
			Tocheck = False
			if LocusToUse == 'All':
				LocusToUse = {}
				Tocheck = True


			genepath=''
			basepath=''
			allQueryBasePaths = []
			action = "Locus"


			for queryFile in queryFilesOnDir:
				countNumberOfLocus += 1
				if Tocheck == True:
					LocusToUse[queryPath] = True
				LocusID[str(numberOfQueries)] = queryPath

				listOfArgs = (queryFile, countNumberOfLocus, alleleFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleScoreFile, LocusToUse)
				job_args, allQueryBasePaths = create_pickle(listOfArgs, queryPath, countNumberOfLocus, job_args, allQueryBasePaths, action)


			joblist =[]

			create_Jobs(job_args, 'Self_BLAST.py', allQueryBasePaths)

			countResults = 0
			for i in allQueryBasePaths:
				countResults += 1
				filepath=os.path.join( i,str(countResults)+"_"+ action + "_results.txt")

				with open(filepath,'rb') as f:
					x = pickle.load(f)

				if x[4] == False:	
						dictOfChanges[str(x[2])] = [x[0], x[1], x[3]] #correctAlleleNumber, sameAlleles, prevAlleleNames  VER O QUE ESTA MAL
					else:
						dictOfChanges[str(x[2])] = False


			print 'End of Self-BLAST: ' + str(datetime.now() - startTimeBLAST)


		job_args = []


		if countIterations >= 2:

			job_args = []
			allQueryBasePaths = []
			countNumberOfResults = 0
			action = "Correct"
			
			for i in arrayOfResults:

				listOfArgs = (arrayOfResults[i], dictOfChanges, i)

				queryPath = os.path.join(os.cwd(), 'correctAlleles')

				countNumberOfResults += 1

				job_args, allQueryBasePaths = create_pickle(listOfArgs, queryPath, countNumberOfResults, job_args, allQueryBasePaths, action)

			create_Jobs(job_args, 'correct_Alleles.py', allQueryBasePaths)


			countResults = 0
			for i in allQueryBasePaths:
				countResults += 1
				filepath=os.path.join( i,str(countResults)+"_"+ action + "_results.txt")

				with open(filepath,'rb') as f:
					x = pickle.load(f)

				if countIterations == 2:
						finalResults[res[1]] = []
						finalResults[res[1]].append(res[0])
					else:
						for locus in finalResults[res[1]][0]:
							if finalResults[res[1]][0][locus][0] == 'LNF' and (res[0][locus][0] == 'INF' or res[0][locus][0] == 'EXM'):
								finalResults[res[1]][0][locus] = res[0][locus]


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
		allQueryBasePaths = []
		action = "BSR"
		
		for referenceGenome in referenceGenomesOnDir:
			countNumberOfGenomes += 1
			GenomesID[countNumberOfGenomes] = referenceGenome
			
			listOfArgs = (referenceGenome, countNumberOfGenomes, databaseFolder, resultsFolder, queryFileWithAll, databaseP, blastResultsPath, alleleFolder, alleleScores, LocusToUse)

			queryPath = os.path.join(os.cwd(), 'BSR')
			job_args, allQueryBasePaths = create_pickle(listOfArgs, queryPath, countNumberOfGenomes, job_args, allQueryBasePaths, action)


		create_Jobs(job_args, 'BSR.py', allQueryBasePaths)

		startTimeBLAST = datetime.now()
		textToFile = '\nBegin of BSR at iteration ' + str(countIterations) + '\n'
		textToFile += '----------------------------------------\n'

		WriteTimeResults(textToFile, False, TimeFile)

		countResults = 0
		for i in allQueryBasePaths:
			countResults += 1
			filepath=os.path.join( i,str(countResults)+"_"+ action + "_results.txt")

			with open(filepath,'rb') as f:
				x = pickle.load(f)

			arrayOfResults[x[2]] = []
			arrayOfResults[x[2]].append(x[0])
			arrayOfGenomes[x[2]] = []
			arrayOfGenomes[x[2]].append(x[1])
			addNewAllelesArray.append(x[3])



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
			allQueryBasePaths = []
			countNumberOfResults = 0
			action = "Check"

			for i in arrayOfResults:
				countNumberOfResults += 1
				listOfArgs = (arrayOfResults[i], arrayOfGenomes[i])
				queryPath = os.path.join(os.cwd(), 'Check_ToWrite')
				job_args, allQueryBasePaths = create_pickle(listOfArgs, queryPath, countNumberOfResults, job_args, allQueryBasePaths, action)

			create_Jobs(job_args, 'filter_Results.py', allQueryBasePaths)


			ToWriteDict = {}

			countResults = 0
			for i in allQueryBasePaths:
				countResults += 1
				filepath=os.path.join( i,str(countResults)+"_"+ action + "_results.txt")

				with open(filepath,'rb') as f:
					x = pickle.load(f)

				for i in x:
					try: 
						ToWriteDict[x[i][0][2]].append([x[i][0][0], x[i][0][1]])
					except KeyError:
						ToWriteDict[x[i][0][2]] = []
						ToWriteDict[x[i][0][2]].append([x[i][0][0], x[i][0][1]])


				

		###(To be parallelized by instance on writeToFileArray)
		if addNewAlleles == True:
			job_args = []
			allQueryBasePaths = []
			countNumberOfResults = 0
			action = "Create"

			for locus in ToWriteDict:
				countNumberOfResults += 1
				listOfArgs = (ToWriteDict[locus], alleleFolder)
				queryPath = os.path.join(os.cwd(), 'Create_NewAlleleFiles')
				job_args, allQueryBasePaths = create_pickle(listOfArgs, queryPath, countNumberOfResults, job_args, allQueryBasePaths, action)

			create_Jobs(job_args, 'create_Allele_Files.py', allQueryBasePaths)


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

