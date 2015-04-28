#!/usr/bin/python
from BSRfunctions import NewAlleleFiles
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
	CreateNewAlleleFiles(argumentList)
	
def CreateNewAlleleFiles(args):
    NewAlleleFiles(args[0], args[1])

    return True

	

if __name__ == "__main__":
    main()
