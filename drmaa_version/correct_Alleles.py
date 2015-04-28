#!/usr/bin/python
from BSRfunctions import correctAlleles

def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)

	def get_correctAlleles(args):
	    newArrayOfResults, genomeID = correctAlleles(args[0], args[1], args[2])

	    final =	(newArrayOfResults, genomeID)

	    filepath=os.path.join(temppath , args[1]+"_Correct_result.txt")

	    with open(filepath, 'wb') as f:
			pickle.dump(final, f)

	    return True


	get_correctAlleles(argumentList)

if __name__ == "__main__":
    main()