import sys
import re

class NewickString:
    def __init__(self):
        self.string = ""
        self.pos = 0

#A class to create and manipulate hierarchichal structures
class Node:
	def __init__(self):
		self.label = ""
		self.length = ""
		self.parent = None
		self.children = []
		self.istip = False
		self.counter = 0
		self.newick = ""
		self.bipart = []
		
	
	def add_child(self,child):
		
		assert child not in self.children
		self.children.append(child)
		child.parent = self
	
	def remove_child(self, child):

		assert child in self.children
		self.children.remove(child)
		child.parent = None
		
	#gets the clade names
	def get_nms(self,nms):
		if self.istip == True:
			nms.append(self.label)
		for child in self.children:
			child.get_nms(nms)		
	
	
	def post_order_bips(self,bips):
		if self.istip != True:
			array = []
			self.get_nms(array)
			bips.append(array)
		for child in self.children:
			child.post_order_bips(bips)
			
	def associate_name(self,bipart_and_name,Choose):
		if self.istip != True:
			array = []
			self.get_nms(array)
			for x in bipart_and_name[1:]:
				if sorted(array) == sorted(x[0]):
					if len(x[1]) == 1:
						x[1][0] = x[1][0].replace(" ","")
						x[1][0] = x[1][0].replace(".","")
						x[1][0] = x[1][0].replace(",","")
						x[1][0] = x[1][0].replace("'","")
						x[1][0] = x[1][0].replace(")","")
						x[1][0] = x[1][0].replace("(","")
						self.label = x[1][0]
					elif Choose == "False":
						x[1][1] = x[1][1].replace(" ","")
						x[1][1] = x[1][1].replace(".","")
						x[1][1] = x[1][1].replace(",","")
						x[1][1] = x[1][1].replace("'","")
						x[1][1] = x[1][1].replace(")","")
						x[1][1] = x[1][1].replace("(","")
						self.label = x[1][1]
					else:
						a = input(str(x[1]))
						a = a.replace(" ","")
						a = a.replace(".","")
						a = a.replace(",","")
						a = a.replace("'","")
						a = a.replace(")","")
						a = a.replace("(","")
						self.label = a
		for child in self.children:
			child.associate_name(bipart_and_name,Choose)
	
	def get_newick_repr(self, showbl=False):
		ret = ""
		for i in range(len(self.children)):
			if i == 0:
				ret += "("
			ret += self.children[i].get_newick_repr(showbl)
			if i == len(self.children)-1:
				ret += ")"
			else:
				ret += ","
		if self.label != None:
			ret += str(self.label)
		if showbl == True:
			ret += ":" + str(self.length)
		return ret
	
	#Part of recurse tree
	def child_props(self, newick):
		
		while re.match(r"^[A-Za-z0-9\.@_-]*$", newick.string[newick.pos]):
			self.label += newick.string[newick.pos]
			newick.pos += 1
		if newick.string[newick.pos] == ":":
			newick.pos += 1
			while re.match(r"^[0-9e\.-]*$", newick.string[newick.pos]):
				self.length += str(newick.string[newick.pos])
				newick.pos += 1
		if newick.string[newick.pos] == ")":
			newick.pos += 1
		return newick.pos
	
	#Recursively make a tree structure
	def recurse_tree(self, newick):

		#check opening bracket and if this is true then that
		#means an opening
		if newick.string[newick.pos] == "(":
			newick.pos += 1
			if newick.string[newick.pos] == "(":
				nd = Node()
				nd.recurse_tree(newick)
				self.add_child(nd)
			#Its a tip node
			if re.match(r"^[A-Za-z0-9\.:@_-]*$", newick.string[newick.pos]):
				nd = Node()
				nd.istip = True
				newick.pos = nd.child_props(newick)
				self.add_child(nd)
		#This means it is either a sister or this is a situation where there are no
		#branch lengths or support/clade labels
		if newick.string[newick.pos] == ",":
			newick.pos += 1
			if re.match(r"^[A-Za-z0-9\.:@_-]*$", newick.string[newick.pos]):
				nd = Node()
				nd.istip = True
				newick.pos = nd.child_props(newick)
				self.add_child(nd)
			if newick.string[newick.pos] == "(":
				nd = Node()
				nd.recurse_tree(newick)
				self.add_child(nd)
			if newick.string[newick.pos] == ",":
				self.recurse_tree(newick)
			
		if re.match(r"^[A-Za-z0-9:\.@_-]*$", newick.string[newick.pos]):
			newick.pos = self.child_props(newick)
		
		if newick.string[newick.pos] == ")":
			newick.pos += 1
		
	#Adds the newick string at every node
	def get_newick(self):

		self.newick = ""
		for child in range(0,len(self.children)):
			if child == 0:
				self.newick += "("
			self.children[child].get_newick()
			self.newick += self.children[child].newick
			if child == (len(self.children)-1):
				self.newick += ")"
			else:
				self.newick += ","
		if self.label != "":
			self.newick += self.label
		if self.length != "":
			self.newick += ":" + self.length

if __name__ == "__main__":

    
	#if len(sys.argv) != 3:
	#	print("To run your own: python3 TreeCode.py Names FormattedTable")
	#	sys.exit()   
    
	tree_string = "(Spirodela_polyrhiza:0.0551939487,(Butomus_umbellatus:0.0483004826,Zostera_marina:0.0965913638)100:0.1212575939,((((((((((((Anthoceros_agrestis:0.0007600629,Anthoceros_punctatus:0.0010565928)100:0.0256554783,Anthoceros_angustus:0.0028050002)100:0.1153257700,Nothoceros_aenigmaticus:0.1314108925)100:0.0386647168,Leiosporoceros_dussii:0.0799416598)100:0.1492412636,(((((((Bartramia_pomiformis:0.0188044986,(Pohlia_nutans:0.0245015901,(((((Sanionia_uncinata:0.0100047285,Climacium_dendroides:0.0080624310)99:0.0010062270,(Callicladium_imponens:0.0066195059,Pseudanomodon_attenuatus:0.0079610314)99:0.0009364054)100:0.0026904614,Myurella_julacea:0.0095040092)100:0.0116586865,Ptychomnion_cygnisetum:0.0629123140)100:0.0189721506,(((Orthotrichum_obtusifolium:0.0022698224,Stoneobryum_bunyaense:0.0027515810)96:0.0002714822,Orthotrichum_stellatum:0.0021561149)100:0.0006879414,Ulota_hutchinsiae:0.0023782770)100:0.0183651456)100:0.0074008266)98:0.0024511502)100:0.0167403543,Syntrichia_filaris:0.0404411737)96:0.0021114580,(Funaria_hygrometrica:0.0056631840,Physcomitrium_patens:0.0093548591)100:0.0115392293)100:0.0308856165,Buxbaumia_aphylla:0.1125083079)100:0.0116766741,((Atrichum_angustatum:0.0144330181,Polytrichum_commune:0.0031060655)100:0.0448518143,Tetraphis_pellucida:0.0603425331)98:0.0039920220)100:0.0417167218,Sphagnum_palustre:0.1046800551)100:0.0623417493,((((Nowellia_curvifolia:0.0492713080,(Diplophyllum_taxifolium:0.0084862519,(Douinia_plicata:0.0066028911,Scapania_ampliata:0.0054419238)100:0.0008540633)100:0.0215315625)100:0.0119790203,Gymnomitrion_concinnatum:0.0283045601)100:0.0122100225,Aneura_pinguis:0.0730168041)100:0.0456545468,((Dumortiera_hirsuta:0.0153586176,(Riccia_fluitans:0.0076491037,Wiesnerella_denudata:0.0093523538)100:0.0012703094)100:0.0015482889,(Marchantia_paleacea:0.0105744748,Marchantia_polymorpha_subsp._ruderalis:0.0024331568)100:0.0055587032)98:0.0042950064)100:0.1842733808)98:0.0199562158)100:0.1553640741,(Ophioglossum_californicum:0.2971337661,Psilotum_nudum:0.2194976067)100:0.2166405068)100:0.4614425837,(Cycas_taitungensis:0.0638305711,(Ginkgo_biloba:0.0607298447,(Pinus_taeda:0.0598713034,Welwitschia_mirabilis:0.5825322776)100:0.0589609878)100:0.0336779760)100:0.0585096669)100:0.1242436430,(Nymphaea_colorata:0.0100883581,Nymphaea_hybrid_cultivar:0.0132497069)100:0.0783223805)100:0.0116842686,Schisandra_sphenanthera:0.0629326589)100:0.0148312810,(Liriodendron_tulipifera:0.0050617456,(Magnolia_biondii:0.0075154154,Magnolia_officinalis:0.0010697989)100:0.0022407865)100:0.0274891634)87:0.0037392560,((((((((((((((Apium_graveolens:0.0148940768,Saposhnikovia_divaricata:0.0127446897)100:0.0084888016,Coriandrum_sativum:0.0062260521)100:0.0410143520,Triosteum_pinnatifidum:0.0247541387)91:0.0022626871,((Arctium_lappa:0.0087335408,((((((((Bidens_bipinnata:0.0000784020,Bidens_biternata:0.0000024791)100:0.0001565092,Bidens_pilosa:0.0011405230)100:0.0010913656,Bidens_parviflora:0.0016341499)100:0.0020710270,Bidens_tripartita:0.0008889193)100:0.0091315655,(Helianthus_annuus:0.0036929880,(Helianthus_grosseserratus:0.0014571655,(Helianthus_strumosus:0.0001981572,Helianthus_tuberosus:0.0000028313)100:0.0002222249)100:0.0002122352)100:0.0094610074)100:0.0105001556,Diplostephium_hartwegii:0.0076969142)93:0.0008230371,Lactuca_sativa:0.0066728272)100:0.0086420000,Chrysanthemum_boreale:0.0103912992)92:0.0030885579)100:0.0486817936,Platycodon_grandiflorus:0.1517094131)100:0.0260232339)100:0.0055839117,Ilex_pubescens:0.0243516507)98:0.0033927037,(((((Asclepias_syriaca:0.0319278472,Cynanchum_auriculatum:0.0031113415)100:0.0289654000,Rhazya_stricta:0.0072539515)100:0.0151492103,Scyphiphora_hydrophyllacea:0.0294838690)100:0.0133652057,((((Calystegia_soldanella:0.0088009406,Ipomoea_nil:0.0068447906)78:0.0026364111,Evolvulus_alsinoides:0.0212712053)78:0.0023745364,Cuscuta_japonica:0.0471515897)100:0.0331286562,(((Capsicum_annuum:0.0102808758,((Solanum_aethiopicum:0.0011612407,Solanum_melongena:0.0008862905)100:0.0196898595,((Solanum_lycopersicum:0.0010246883,Solanum_pennellii:0.0006900082)100:0.0018780558,Solanum_tuberosum:0.0032179237)100:0.0076192114)92:0.0026512712)100:0.0159829959,(Hyoscyamus_niger:0.0154208123,Physochlaina_orientalis:0.0014064811)100:0.0045981486)100:0.0043938226,(Nicotiana_attenuata:0.0026099433,(Nicotiana_sylvestris:0.0000024791,Nicotiana_tabacum:0.0000024791)100:0.0014798648)100:0.0026852825)100:0.0311112521)100:0.0138567371)89:0.0025962860,((Dorcoceras_hygrometricum:0.0688199472,(Ajuga_reptans:0.2008429902,Castilleja_paramensis:0.0178408174)100:0.0111336453)100:0.0218444500,(Hesperelaea_palmeri:0.0059694283,Osmanthus_fragrans:0.0041252362)100:0.0128761315)100:0.0227927395)100:0.0117205247)100:0.0036691889,(((Rhododendron_simsii:0.0161608180,Vaccinium_macrocarpon:0.0313177514)100:0.0333355101,Aegiceras_corniculatum:0.0560940366)100:0.0047173466,Camellia_sinensis:0.0119333626)100:0.0129574697)94:0.0088887311,((((Agrostemma_githago:0.0094150108,Silene_latifolia:0.0177705910)100:0.0390871508,(((Beta_macrocarpa:0.0000910881,Beta_vulgaris_subsp._vulgaris:0.0006590551)100:0.0000010737,Beta_vulgaris_subsp._maritima:0.0000981301)100:0.0105616464,((Chenopodium_quinoa:0.0150349307,Spinacia_oleracea:0.0186163381)100:0.0084917680,Suaeda_glauca:0.0338384127)100:0.0064048806)100:0.0107617615)100:0.0103903458,(Bougainvillea_spectabilis:0.0044683441,(Mirabilis_himalaica:0.0152399048,Mirabilis_jalapa:0.0173755827)100:0.0098314913)100:0.0161811474)100:0.0288137377,Nepenthes_ventricosa_x_Nepenthes_alata:0.0212673716)100:0.0408846747)90:0.0032796893,((((((((((((((Bupleurum_falcatum:0.0015642195,Glycyrrhiza_uralensis:0.0007855101)100:0.0090647058,((Medicago_truncatula:0.0118168773,((Trifolium_aureum:0.0030187901,Trifolium_grandiflorum:0.0038382188)100:0.0024256645,(Trifolium_meduseum:0.0038989101,Trifolium_pratense:0.0048998855)100:0.0033213065)100:0.0026516031)100:0.0014790239,(Pisum_abyssinicum:0.0019755097,Pisum_fulvum:0.0026521202)100:0.0049558298)100:0.0147872093)100:0.0024453350,Lotus_japonicus:0.0115612368)100:0.0037878586,(((Glycine_max:0.0000028371,Glycine_soja:0.0016365218)100:0.0051550311,(Phaseolus_vulgaris:0.0065524042,(Vigna_angularis:0.0036696678,Vigna_radiata:0.0026305717)100:0.0059810465)100:0.0047800638)100:0.0064752780,Pongamia_pinnata:0.0098336655)100:0.0045040951)100:0.0050585967,Ammopiptanthus_nanus:0.0156257534)100:0.0125890732,Styphnolobium_japonicum:0.0150669007)100:0.0185243981,(((Acacia_ligulata:0.0064705688,Leucaena_trichandra:0.0123176611)100:0.0169588152,(Haematoxylum_brasiletto:0.0079146152,Libidibia_coriaria:0.0018661582)100:0.0045830728)100:0.0030687913,((Senna_occidentalis:0.0036163170,Senna_tora:0.0143149670)100:0.0070058530,Gleditsia_sinensis:0.0098229012)100:0.0063507720)100:0.0100624813)79:0.0020420360,Tamarindus_indica:0.0210562261)100:0.0176002337,Epirixanthes_elongata:0.1870581562)100:0.0082736540,(Geranium_maderense:0.2449695924,((Mangifera_longipes:0.0136588343,(Spondias_mombin:0.0186177750,Spondias_tuberosa:0.0021104623)100:0.0185719150)100:0.0389191679,((Citrus_maxima:0.0070968777,Citrus_sinensis:0.0278179024)100:0.0471974938,(Acer_yangbiense:0.0141503579,Sapindus_mukorossi:0.0139356440)100:0.0117614282)100:0.0036736096)100:0.0201979438)30:0.0043248443)34:0.0135612953,Fagus_sylvatica:0.0284206728)29:0.0145027909,((Citrullus_lanatus:0.0059960544,(Cucumis_sativus:0.0594603560,Cucurbita_pepo:0.0383390294)100:0.0074846865)100:0.0354051737,(((Cannabis_sativa:0.0192421479,Morus_notabilis:0.0169352091)100:0.0209316751,Ziziphus_jujuba:0.0225637363)100:0.0046575884,((((((Fragaria_gracilis:0.0000024791,Fragaria_tibetica:0.0001578890)100:0.0001594047,Fragaria_nubicola:0.0008607146)100:0.0012019798,(Fragaria_iinumae:0.0027550473,Fragaria_iturupensis:0.0016383585)95:0.0006443114)98:0.0004004776,((Fragaria_moschata:0.0000024791,Fragaria_orientalis:0.0001938190)100:0.0006349987,Fragaria_viridis:0.0003304020)98:0.0002186014)99:0.0009642168,Fragaria_nilgerrensis:0.0016451691)100:0.0391329298,(((Malus_domestica:0.0009569525,(Rhaphiolepis_bibas:0.0040786139,Sorbus_aucuparia:0.0012506894)100:0.0014095596)60:0.0006452131,Torminalis_glaberrima:0.0019152372)100:0.0175076754,(Prunus_avium:0.0075889464,Prunus_salicina_x_Prunus_armeniaca:0.0048518056)100:0.0141058178)100:0.0131293319)100:0.0162463125)100:0.0073302302)85:0.0011832961)29:0.0041399454,(((((((Arabidopsis_thaliana:0.0088490249,Arabis_alpina:0.0149291074)96:0.0110054829,Capsella_rubella:0.0062479500)89:0.0000908997,Boechera_stricta:0.0040045028)96:0.0032635639,((((Brassica_carinata:0.0000028944,Brassica_nigra:0.0112915500)100:0.0013111321,Sinapis_arvensis:0.0006361962)100:0.0008664685,((Brassica_juncea:0.0000800184,(Brassica_napus:0.0003151783,Brassica_oleracea:0.0000775380)94:0.0000024846)61:0.0000020490,Brassica_rapa:0.0107392668)100:0.0010471623)95:0.0004073285,Raphanus_sativus:0.0031758239)100:0.0061263426)100:0.0631471151,Carica_papaya:0.0237271980)100:0.0087596353,(((Bombax_ceiba:0.0202449289,(((Gossypium_arboreum:0.0001324196,((Gossypium_davidsonii:0.0002046367,((Gossypium_harknessii:0.0000795326,Gossypium_hirsutum:0.0000021271)100:0.0002275717,Gossypium_raimondii:0.0043745357)98:0.0000677419)100:0.0001302510,(Gossypium_thurberi:0.0000025639,Gossypium_trilobum:0.0000666909)100:0.0005326083)100:0.0004027253)97:0.0002131577,Gossypium_barbadense:0.0000564631)100:0.0055804725,Hibiscus_cannabinus:0.0150072668)100:0.0032733082)100:0.0159746289,(Corchorus_capsularis:0.0258380945,Corchorus_olitorius:0.0402260574)100:0.0724398165)100:0.0236309347,Aquilaria_sinensis:0.0531205933)100:0.0082066377)100:0.0090795981,((Manihot_esculenta:0.0193779449,Ricinus_communis:0.0156933462)100:0.0098173577,(Passiflora_edulis:0.0612861017,((Populus_alba:0.0003134438,(Populus_davidiana:0.0002236203,(Populus_tremula_x_Populus_alba:0.0000024791,Populus_tremula:0.0000747160)100:0.0002231279)100:0.0003655432)99:0.0008604971,((Salix_brachista:0.0009352362,(Salix_purpurea:0.0018124946,Salix_suchowensis:0.0006828492)99:0.0006150336)100:0.0020008115,Salix_dunnii:0.0021179042)100:0.0038160814)100:0.0647394613)100:0.0138409627)100:0.0140161099)95:0.0068736219)64:0.0025748041,(Lagerstroemia_indica:0.0402255769,Eucalyptus_grandis:0.0422219315)100:0.0184704430)98:0.0072947597)91:0.0054610239,((Tolypanthus_maclurei:0.0562256900,Viscum_album:0.8674572535)94:0.0138504572,Malania_oleifera:0.0296183561)94:0.0082450428)100:0.0073307129,Vitis_vinifera:0.0437443236)100:0.0256621383,Macadamia_integrifolia:0.0205287686)95:0.0033408902,Nelumbo_nucifera:0.0123396436)99:0.0024488370,(Aconitum_kusnezoffii:0.0324098043,Hepatica_maxima:0.0365369587)100:0.0252387043)100:0.0182122312)100:0.0350351242,(((Cocos_nucifera:0.0271425126,Phoenix_dactylifera:0.0049465734)100:0.0142178877,(Cyperus_esculentus:0.2199792432,((((Chrysopogon_zizanioides:0.0137752884,((Saccharum_officinarum:0.0004903787,Sorghum_bicolor:0.0038397413)100:0.0074279252,(Tripsacum_dactyloides:0.0025566146,(Zea_luxurians:0.0008739480,Zea_perennis:0.0010336432)100:0.0070916257)100:0.0230707743)99:0.0091621163)100:0.0172776173,Eleusine_indica:0.0232681990)100:0.0078855558,(((Oryza_rufipogon:0.0000313579,Oryza_sativa_Japonica_Group:0.0088694019)100:0.0143151198,Oryza_sativa_Indica_Group:0.0003665615)100:0.0021273062,Oryza_minuta:0.0016678483)100:0.0184406219)99:0.0046477678,(Triticum_aestivum:0.0006988807,Triticum_timopheevii:0.0038494298)100:0.0317741488)100:0.0384299999)100:0.0311540477)100:0.0080732476,(Allium_cepa:0.0707807752,Asparagus_officinalis:0.0154904111)100:0.0279091720)100:0.0247782510)100:0.0061973708);"

	root = Node()
	newick = NewickString()
	newick.string = tree_string
	root.recurse_tree(newick)
	full_array = []
	root.get_nms(full_array) #contains all taxa
    
	'''
    #create a hash to try and look up names
	names_file = open(sys.argv[1],"r")
	array = []
	HASH = {}
	counter = 0
	for i in names_file:
		counter += 1
		array = i.split("|")
		array[0] = array[0].strip()
		array[1] = array[1].strip()
		sys.stderr.write("Hashing: " + str(counter) + "\r")
		HASH[array[1]] = array[0]
	
	CodeHash = {}
	for i in full_array:
		
		name = i
		b = False
		counter = 0
		while b == False:
			if name in HASH:
				print("Found: " + name + " Name in dataset is: " + i + " Code is: " + HASH[name])
				CodeHash[HASH[name]] = i
				
				b = True
			else:
				name = name.replace("_"," ") #test if problem is simple
				if counter == 1:
					name = input(str(i) + "not found try another name: ")
				counter += 1
	test_file = open(sys.argv[2],"r")
	counter = 0
	for i in test_file:
		counter += 1
		test = i.split("\t")[0]
		if test in CodeHash:
			print(i.strip("\n"))
		
		sys.stderr.write("Running at line: " + str(counter) + "\r")
			
	'''		
			
			    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

