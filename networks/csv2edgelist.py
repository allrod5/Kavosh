import sys
import getopt

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print('Error: incorrect parameters.\n Usage: csv2edgelist.py -i <inputfile> -o <outputfile>\n')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('Usage: test.py -i <inputfile> -o <outputfile>\n')
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg

	edgelist = open(outputfile, "w")

	with open(inputfile, "r") as csv:

		beginning = csv.tell()

		line = csv.readline()

		if line[0] != '0':
			mod = 1
		else:
			mod = 0

		print("Listing nodes...")

		nodelist = []
		csv.seek(beginning)
		counter = 0
		for line in csv:
			counter += 1
			nodes = line.split(';')
			if len(nodes) > 1:
				for node in nodes:
					if (int(node)-mod) not in nodelist:
						nodelist.append(int(node)-mod)
			if counter % 10000 == 0:
				print("\tline "+str(counter))

		print("\tdone!\n\nMarking gaps...")

		markslist = []
		nodelist.sort()
		
		for i in range(1,len(nodelist)):
			if nodelist[i-1]+1 != nodelist[i]:
				markslist.append([nodelist[i],nodelist[i]-nodelist[i-1]-1])

		print("\tdone!\n\nConverting...")
		csv.seek(beginning)
		counter = 0
		for line in csv:
			counter += 1
			nodes = line.split(';')
			source = str(int(nodes[0])-mod)
			
			gap = 0
			for mark in markslist:
				if mark[0] <= int(source):
					gap += mark[1]
			source = str(int(source)-gap)
			
			for j in range(1,len(nodes)):
				target = str(int(nodes[j])-mod)
				
				gap = 0
				for mark in markslist:
					if mark[0] <= int(target):
						gap += mark[1]
				target = str(int(target)-gap)
				
				print(''.join(source+' '+target), file = edgelist)

			if counter%10000 == 0:
				print("\tline "+str(counter))





if __name__ == "__main__":
	main(sys.argv[1:])

