#!/usr/bin/python
from BSRfunctions import func_BSR
import sys
import pickle
import os

def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)
	BSR(argumentList,temppath)

def BSR(args,temppath):
    resultsList, referenceCDSsequences, addNewAlleles = func_BSR(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9])
    print resultsList
    
    final =	(resultsList, referenceCDSsequences, args[1], addNewAlleles)

    filepath=os.path.join(temppath , str(args[1])+"_BSR_result.txt")
    print filepath
    with open(filepath, 'wb') as f:
     	pickle.dump(final, f)

    return True

	

if __name__ == "__main__":
    main()
