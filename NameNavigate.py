import sys


#associate taxa with full code
def find_codes(taxa_array,node_file,CodeNameHASH):
	
	nodes = open(node_file,"r")
	
	HASH = {} #contains the code [key] and the parent code [value]
	RankHASH = {} #contains the code [key] and the rank [value]
	counter = 0
	
	#Get the information from the node file
	for i in nodes:	
		
		counter += 1
		i = "".join(i.split())
		array = i.split("|")
		HASH[array[0]] = array[1]
		RankHASH[array[0]] = array[2]
		sys.stderr.write("Running at: " + str(counter) + "\r")
	
	nodes.close()
	
	InfoHASH = {} #code associated with a series of array
	#navigate the hash for info associated with each taxa
	for i in taxa_array:
	
		rank_line = []
		code_line = []
		names_line = []
		cur_node = i[2]
		cur_rank = RankHASH[i[2]]
		cur_name = CodeNameHASH[i[2]]
		while cur_node != "1":
			rank_line.append(cur_rank)
			code_line.append(cur_node)
			names_line.append(cur_name)
			cur_node = HASH[cur_node]
			cur_rank = RankHASH[cur_node]
			cur_name = CodeNameHASH[cur_node]
		
		InfoHASH[i[0]] = [i,code_line,rank_line,names_line]

	return InfoHASH

#find where the earliest matching point in the arrays is
def get_meeting_point(bipartition):

	meet = ""
	for i in range(len(bipartition[0])):
		counter = 0
		for j in bipartition[1:]:
			if bipartition[0][i] in j:
				counter += 1
				if counter == len(bipartition[1:]):
					meet = bipartition[0][i]
					break
		if meet != "":
			break
	return(meet)


#find info for a bipartition (need to have it factor in other side)
def BipartitionToCode(bip_array,nms_array,SpInfoHASH,CodeNameHash):

	bipart_and_name = []
	#for each bipart find the lowest place in taxonomy that they match
	for i in bip_array:
		
		b_a_n = []
		right_side = list(set(nms_array)-set(i)) #get the other side of the bipartition
		#print(str(i) + "=========" + str(right_side))
		#find the earliest number they all share and is not shared by anyone else
		taxa_list = []
		meeting_point = ""
		for j in i:
			taxa_list.append(SpInfoHASH[j][1])
		meeting_point = get_meeting_point(taxa_list) #get the code of where the ancestor is
		
		#check to make sure the meeting point isn't in anything else
		for j in right_side:
			if meeting_point in SpInfoHASH[j][1]:
				meeting_point = ""
		
		if meeting_point != "":
			b_a_n = [i,CodeNameHash[meeting_point]]
		else:
			b_a_n = [i,[""]]
		#print("MATCH")
		#print(CodeNameHash[meeting_point])
		bipart_and_name.append(b_a_n)
	return bipart_and_name
				
		
		
		
		
		
		
		
		