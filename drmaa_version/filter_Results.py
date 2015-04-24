from BSRfunctions import check_ToWrite

def main():

	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)

	def check_ToWrite(args):
	    ToWrite = TowriteToAlleleFiles(args[0], args[1])

	    final =	(ToWrite)

	    filepath=os.path.join(temppath , args[1]+"_Check_result.txt")

	    with open(filepath, 'wb') as f:
			pickle.dump(final, f)

	    return True

	check_ToWrite()

if __name__ == "__main__":
    main()