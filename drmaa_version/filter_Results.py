#!/usr/bin/python
from BSRfunctions import TowriteToAlleleFiles
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
	check_ToWrite(argumentList,temppath)

def check_ToWrite(args,temppath):
	    ToWrite = TowriteToAlleleFiles(args[0], args[1])

	    final =	(ToWrite)

	    filepath=os.path.join(temppath , str(args[2])+"_Check_result.txt")

	    with open(filepath, 'wb') as f:
			pickle.dump(final, f)

	    return True

	
if __name__ == "__main__":
    main()
