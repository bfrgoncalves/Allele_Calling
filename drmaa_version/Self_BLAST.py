#!/usr/bin/python
from BSRfunctions import func_Self_BLAST


def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)


	def SelfBLAST(args):
	    alleleNumbers, sameAlleles, countNumberOfLocus, prevAlleleNames, isEmpty = func_Self_BLAST(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])

	    final =	(alleleNumbers, sameAlleles, countNumberOfLocus, prevAlleleNames, isEmpty)

	    filepath=os.path.join(temppath , args[1]+"_Locus_result.txt")

	    with open(filepath, 'wb') as f:
			pickle.dump(final, f)

	    return True


	SelfBLAST(argumentList)

if __name__ == "__main__":
    main()