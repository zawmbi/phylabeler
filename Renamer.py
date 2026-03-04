import sys

#Function designed to try and find the names in a name file
#Associated with names on genbank
def make_names(list_of_names,nf):
	
	#first put the names into a hash, key = name, value = code
	names_file = open(nf,"r")
	array = []
	HASH = {}
	HASH2 = {} #value will be an array since you can have multiple names
	counter = 0
	for i in names_file:
		counter += 1
		sys.stderr.write("Hashing: " + str(counter) + "\r")
		array = i.split("|")
		array[0] = array[0].strip()
		array[1] = array[1].strip()
		#print(array[1]+"\t"+array[3]+"\n") #This will contain the type of name it is
		sys.stderr.write("Hashing: " + str(counter) + "\r")
		HASH[array[1]] = array[0]
		if array[0] in HASH2:
			HASH2[array[0]].append(array[1])
		else:
			HASH2[array[0]] = []
			HASH2[array[0]].append(array[1])
	
	names_file.close()
	array_of_taxa = []
	#Try to identify the correct names
	print(list_of_names)
	for i in list_of_names:
		
		name = i
		b = False
		counter = 0
		while b == False:
			if name in HASH:
				array = [i,name,HASH[name]]
				#print("Found: " + i + " Name in genbank is: " + name + " Code is: " + HASH[name])
				array_of_taxa.append(array)
				b = True
			else:
				name = name.replace("_"," ") #test if problem is simple
				if 1 < counter:
					name = input(str(name) + " not found try another name: ")
				counter += 1
	return HASH2, array_of_taxa