from BSRfunctions import CreateNewAlleleFiles

def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)

	def CreateNewAlleleFiles(args):
    	NewAlleleFiles(args[0], args[1])

    	return True

	CreateNewAlleleFiles(argumentList)

if __name__ == "__main__":
    main()