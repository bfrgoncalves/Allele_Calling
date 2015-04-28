#!/usr/bin/python
from BSRfunctions import func_BSR


def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)

	def BSR(args):
		resultsList, referenceCDSsequences, addNewAlleles = func_BSR(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])
     	final =	(resultsList, referenceCDSsequences, args[1], addNewAlleles)

	    filepath=os.path.join(temppath , args[1]+"_BSR_result.txt")

	    with open(filepath, 'wb') as f:
			pickle.dump(final, f)

	    return True

	BSR(argumentList)

if __name__ == "__main__":
    main()