import argparse
import sys
import TreeCode, Renamer, NameNavigate


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    description="A program to help label trees based on genbank taxonomy.")
    parser.add_argument("-d", "--description", action="store_true",
                        help="Gives a longer description of the program and different modes.")
    parser.add_argument("--tree",type = str,
                        help="input")
    parser.add_argument("--name_file",type = str,
                        help="genbank name file, typically called names.dmp")
    parser.add_argument("--node_file",type = str,
                        help="genbank node file, typically called node.dmp")
    parser.add_argument("--outfile", type=str, default="out",
                        help="A prefix for all of the outfiles produced by the program in whichever mode you run it in. ")

    if len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()

    args = parser.parse_args()
    if args.description:
        print("This program identifies labels tree nodes based on the most inclusive \n\
name available in genbank assuming the group is monophyletic \n")
        sys.exit(0)

    input_tree = args.tree
    outfile = args.outfile
    input_names = args.name_file
    input_nodes = args.node_file
    
    #Procedure for getting the code associated with nodes
    if input_tree and input_names and input_nodes:
    	tree_file = open(input_tree,"r")
    	NewickInput = TreeCode.NewickString()
    	NewickInput.string = tree_file.readline()
    	root = TreeCode.Node()
    	root.recurse_tree(NewickInput)
    	name_array = []
    	
    	#get the tips of the tree
    	root.get_nms(name_array)
    	
    	#print(name_array)
    	#sys.exit()
    	#get the tree bipartitions
    	bip_array = []
    	root.post_order_bips(bip_array)

    	#associate the tips with genbank data
    	taxa_array = [] #contains an array of arrays [TreeName,GenbankName,GenbankCode]
    	CodeNameHash = {} #contains the code and associated name for all taxa
    	CodeNameHash, taxa_array = Renamer.make_names(name_array,input_names) #send the name_array and the name file loc
    	
    	#get the full list of taxa affiliation
    	SpInfoHASH = NameNavigate.find_codes(taxa_array,input_nodes,CodeNameHash)
    	
    	#associate taxa info with bipartition
    	bipart_and_name = []
    	bipart_and_name = NameNavigate.BipartitionToCode(bip_array,name_array,SpInfoHASH,CodeNameHash)
    	
    	Choose = input("Choose the name (True) or have random selection (False): ")
    	root.associate_name(bipart_and_name,Choose)
    	print(root.get_newick_repr(True)+";")

    else:
    	print("Missing arguments")
    	

