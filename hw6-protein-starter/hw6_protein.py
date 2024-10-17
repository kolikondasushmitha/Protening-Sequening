# """
# 15-110 Hw6 - Protein Sequencing Project
# Name: KOLIKONDA SUSHMITHA
# AndrewID: 2023501078
# """

# import hw6_protein_tests as test

# project = "Protein" # don't edit this

# ### WEEK 1 ###

# '''
# readFile(filename)
# #1 [Check6-1]
# Parameters: str
# Returns: str
# '''
# def readFile(filename):
#     with open(filename, 'r') as file:
#         #removing the new lines from the text
#         data = file.read().replace('\n', '')
#     return data


# '''
# dnaToRna(dna, startIndex)
# #2 [Check6-1]
# Parameters: str ; int
# Returns: list of strs
# '''
# def dnaToRna(dna, startIndex):
#     dna=dna.replace("T","U")        # we are replacing every "T" in dna with "U" , inorder to convert it into rna.
#     rna=[]
#     rna_codons = []                 # list to store rna codons.

#     for i in range(startIndex,len(dna),3):
#         rna.append(dna[i:i+3])      # creating rna list by adding strings of size 3.
#     for  rna_String in rna:
#         if rna_String in ["UAA","UAG","UGA"]: # if rna_string is a end codon , we append
#             rna_codons.append(rna_String)
#             break                   # once we find a end codon , we break the loop.
#         else:
#             rna_codons.append(rna_String)

#     return rna_codons


# '''
# makeCodonDictionary(filename)
# #3 [Check6-1]
# Parameters: str
# Returns: dict mapping strs to strs
# '''
# def makeCodonDictionary(filename):
#     import json
#     file=open(filename,"r")        # opening the file 
#     json_file=json.load(file)  
#     # print(json_file)    # loading the json file json.load() , here json file is a Dictionary.
#     codon_to_aminoacid={}          # dictionary to store codons as keys and amino acids as values. 
#     for amino_acid in json_file:
#         for codon in json_file[amino_acid]:
#             # print(json_file[amino_acid])
#             codon=codon.replace("T","U")    # replacing T with U.
#             if codon not in codon_to_aminoacid: # if codon not in dictionary, we add codon as key and amino acid a value.
#                 codon_to_aminoacid[codon]=amino_acid 

#     return codon_to_aminoacid


# '''
# generateProtein(codons, codonD)
# #4 [Check6-1]
# Parameters: list of strs ; dict mapping strs to strs
# Returns: list of strs
# '''
# def generateProtein(codons, codonD):
#     print(codons)
#        # this function is to extract protein from the RNA sequence.
#     protein=[]                          # to store the chain of amino acids.
#     for i in range(len(codons)):
#         if i==0 and codons[i]=="AUG" :  # checking if the first element of the codons is "AUG" , then we append start to the Protein list.
#             protein.append("Start")
#         elif codonD[codons[i]] == "Stop": # if the key in codonD is Stop , then we append the amino acids and break the loop , since we reached the STOP condition. 
#             protein.append(codonD[codons[i]])
#             break
#         else:
#             protein.append(codonD[codons[i]])   # append the amino acid in the protein.

#     return  protein


# '''
# synthesizeProteins(dnaFilename, codonFilename)
# #5 [Check6-1]
# Parameters: str ; str
# Returns: 2D list of strs
# '''
# def synthesizeProteins(dnaFilename, codonFilename):
#     # Read the DNA sequence from the file
#     dna_sequence = readFile(dnaFilename)
    
#     # Create the codon to amino acid dictionary from the codon file
#     codon_dict = makeCodonDictionary(codonFilename)
    
#     # List to store all proteins synthesized
#     proteins = []
    
#     # Initialize index to search for "ATG" start codons
#     start = dna_sequence.find("ATG")
    
#     # Loop through the DNA sequence and find all protein sequences
#     while start != -1:
#         # Convert DNA sequence starting from "ATG" to RNA sequence
#         rna = dnaToRna(dna_sequence, start)
        
#         # Generate the protein from the RNA sequence using the codon dictionary
#         protein = generateProtein(rna, codon_dict)
        
#         # Add the synthesized protein to the list
#         proteins.append(protein)
        
#         # Look for the next "ATG" start codon after the current one
#         start = dna_sequence.find("ATG", start + 3)
    
#     # Return the list of synthesized proteins
#     return proteins


# def runWeek1():
#     print("Human DNA")
#     humanProteins = synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/human_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")
#     print("Elephant DNA")
#     elephantProteins = synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/elephant_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")


# ### WEEK 2 ###

# '''
# commonProteins(proteinList1, proteinList2)
# #1 [Check6-2]
# Parameters: 2D list of strs ; 2D list of strs
# Returns: 2D list of strs
# '''
# def commonProteins(proteinList1, proteinList2):
#     return


# '''
# combineProteins(proteinList)
# #2 [Check6-2]
# Parameters: 2D list of strs
# Returns: list of strs
# '''
# def combineProteins(proteinList):
#     return


# '''
# aminoAcidDictionary(aaList)
# #3 [Check6-2]
# Parameters: list of strs
# Returns: dict mapping strs to ints
# '''
# def aminoAcidDictionary(aaList):
#     return


# '''
# findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
# #4 [Check6-2]
# Parameters: 2D list of strs ; 2D list of strs ; float
# Returns: 2D list of values
# '''
# def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
#     return


# '''
# displayTextResults(commonalities, differences)
# #5 [Check6-2]
# Parameters: 2D list of strs ; 2D list of values
# Returns: None
# '''
# def displayTextResults(commonalities, differences):
#     return


# def runWeek2():
#     humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
#     elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

#     commonalities = commonProteins(humanProteins, elephantProteins)
#     differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
#     displayTextResults(commonalities, differences)


# ### WEEK 3 ###

# '''
# makeAminoAcidLabels(proteinList1, proteinList2)
# #2 [Hw6]
# Parameters: 2D list of strs ; 2D list of strs
# Returns: list of strs
# '''
# def makeAminoAcidLabels(proteinList1, proteinList2):
#     return


# '''
# setupChartData(labels, proteinList)
# #3 [Hw6]
# Parameters: list of strs ; 2D list of strs
# Returns: list of floats
# '''
# def setupChartData(labels, proteinList):
#     return


# '''
# createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
# #4 [Hw6] & #5 [Hw6]
# Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
# Returns: None
# '''
# def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
#     import matplotlib.pyplot as plt
#     return


# '''
# makeEdgeList(labels, biggestDiffs)
# #5 [Hw6]
# Parameters: list of strs ; 2D list of values
# Returns: list of strs
# '''
# def makeEdgeList(labels, biggestDiffs):
#     return


# '''
# runFullProgram()
# #6 [Hw6]
# Parameters: no parameters
# Returns: None
# '''
# def runFullProgram():
#     return


# ### RUN CODE ###

# # This code runs the test cases to check your work
# if __name__ == "__main__":
#     print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
#     test.week1Tests()
#     print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
#     runWeek1()

#     ## Uncomment these for Week 2 ##
#     """
#     print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
#     test.week2Tests()
#     print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
#     runWeek2()
#     """

#     ## Uncomment these for Week 3 ##
#     """
#     print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
#     test.week3Tests()
#     print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
#     runFullProgram()
#     """
"""
15-110 Hw6 - Protein Sequencing Project
Name: kolikonda Sushmitha
ROLL NUMBER: 2023501078
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):    
    with open(filename, 'r') as file:
        #removing the new lines from the text
        data = file.read().replace('\n', '')
    return data


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):      # this function is used to convert dna to rna.
    dna=dna.replace("T","U")        # we are replacing every "T" in dna with "U" , inorder to convert it into rna.
    rna=[]
    rna_codons = []                 # list to store rna codons.
    for i in range(startIndex,len(dna),3):
        rna.append(dna[i:i+3])      # creating rna list by adding strings of size 3.
    for  rna_String in rna:
        if rna_String in ["UAA","UAG","UGA"]: # if rna_string is a end codon , we append
            rna_codons.append(rna_String)
            break                   
        else:
            rna_codons.append(rna_String)

    return rna_codons


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):  # this function returns a dictionary which maps codons to amino acids 
    import json
    file=open(filename,"r")        # opening the file 
    json_file=json.load(file)      # loading the json file json.load() , here json file is a Dictionary.
    codon_to_aminoacid={}          # dictionary to store codons as keys and amino acids as values. 
    for amino_acid in json_file:
        for codon in json_file[amino_acid]:
            codon=codon.replace("T","U")    # replacing T with U.
            if codon not in codon_to_aminoacid: # if codon not in dictionary, we add codon as key and amino acid a value.
                codon_to_aminoacid[codon]=amino_acid 

    return codon_to_aminoacid


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):    # this function is to extract protein from the RNA sequence.
    protein=[]                          # to store the chain of amino acids.
    for i in range(len(codons)):
        if i==0 and codons[i]=="AUG" :  # checking if the first element of the codons is "AUG" , then we append start to the Protein list.
            protein.append("Start")
        elif codonD[codons[i]] == "Stop": # if the key in codonD is Stop , then we append the amino acids and break the loop , since we reached the STOP condition. 
            protein.append(codonD[codons[i]])
            break
        else:
            protein.append(codonD[codons[i]])   # append the amino acid in the protein.

    return  protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna=readFile(dnaFilename)
    codon_dict=makeCodonDictionary(codonFilename)
    # print(dna,codon_dict)
    proteins=[]
    count=0 # unused base
    index=0
    while index<len(dna):
        if dna[index:index+3]=="ATG":       # if the string is ATG i.e start of the sequence.
            rna=dnaToRna(dna,index)         # converting dna to rna using the function we have written earlier.
            protein=generateProtein(rna,codon_dict)     # convertion of rna codons to amino acids. combination of amino acids is a Protein.
            proteins.append(protein)
            index+= 3*len(rna)      # incrementing index by 3 times length of rna.
        else:
            index+=1
            count+=1  # incrementing unused bases count 
            
    
    total_bases=len(dna)        
    total_proteins=len(proteins)

    print("Total number of Bases: ",total_bases)
    print("Number of un-used bases :",count)
    print("Total number of Proteins :",total_proteins)

    return proteins
    


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/human_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")
    
    print()
    print("************************************")
    print()
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/elephant_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")
    

### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):     # to get common elements from two lists.
    common_proteins=[]
    for protein in proteinList1:        
        if protein in proteinList2:                 # if the element is present in proteinlist2 , we append to common_proteins.
            common_proteins.append(protein)

    return common_proteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''                                         

# purpose is to find the differences between two proteins by dumping all amino acids into a single list.
def combineProteins(proteinList):           # this function takes a 2D list as input and merges all the sublists into 1D list.
    combined_proteins=[]                    # list to store the combined proteins.
    for protein in proteinList:
        combined_proteins.extend(protein)   # extend method is used to add specified list elements to the end of the existing list.  

    return combined_proteins    


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):            # this function takes a List and creates a dictionary to count the frequency of occurence. 
    AminoAcidDict={}                        # dictionary to hold amino acids. 
    for i in aaList:
        if i in AminoAcidDict:
            AminoAcidDict[i]+=1
        else:
            AminoAcidDict[i]=1

    return AminoAcidDict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):       # this function creates a list with amino acids in two protein lists which have frequencies greater than cutoff.
    # list format [amino acid, freq in list1, freq in list2]
    list1=combineProteins(proteinList1) 
    list2=combineProteins(proteinList2)
    dict1=aminoAcidDictionary(list1)
    dict2=aminoAcidDictionary(list2)

    len1=len(list1)       # total number of amino acids in the gene
    len2=len(list2)       # total number of amino acids in the gene

    differences=[]              # list to be returned. consists of amino acid name , frequencies in both lists.

    for amino_acid in dict1:    
        if amino_acid not in dict2:     # for every amino acid in dict 1 , if it is not present in dict 2, we are assigning it a count 0.
            dict2[amino_acid]=0     

    for amino_acid in dict2:
        if amino_acid not in dict1:
            dict1[amino_acid]=0
    
    for amino_acid in dict1.keys():
        if amino_acid not in ["Start" , "Stop"] and (abs(frequency(dict1[amino_acid],len1) - frequency(dict2[amino_acid],len2))>cutoff):
            differences.append([amino_acid,frequency(dict1[amino_acid],len1),frequency(dict2[amino_acid],len2)])
    
    return differences

def frequency(a,b):  # calculates frequency
    return a/b

'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The Proteins which occured in both DNA Sequences are: ")
    l=[]
    for list in commonalities:
        for j in list[1:-1]:
            if j not in l:
                l.append(j)
    for i in l:
        print(i,end=" ")
    
    print("============================================================================================")
    print("The Amino Acids which occured at most Differently rates in both DNA Sequences are: ")
    for list in differences:
        print(f"Amino Acid  {list[0]} : {round(list[1]*100,2)} % in Sequence 1, {round(list[2]*100,2)} % in Sequence 2 ")

    return


def runWeek2():
    humanProteins = synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/human_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "Protening-Sequening/hw6-protein-starter/data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):        # this function is used to generate a list of amino acids which are present in both the protein lists.
    a=combineProteins(proteinList1)
    b=combineProteins(proteinList2)
    dict1=aminoAcidDictionary(a)
    dict2=aminoAcidDictionary(b)

    a= list(dict1.keys()) + list(dict2.keys())
    a=set(a)                    # to deal with the duplicates.

    return sorted(list(a))


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):    # this function is used to create a frequency list which consists of frequencies of each gene. 
    protein=combineProteins(proteinList)
    dict=aminoAcidDictionary(protein)
    frequency_list=[]
    for amino_acid in labels:
        if amino_acid in protein:        
            frequency_list.append(dict[amino_acid]/len(protein))     # appending frequency 
        else:
            frequency_list.append(0.0)
    
    return frequency_list


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):  # function to plot Bar chart.
    import matplotlib.pyplot as plt
    w=0.40
    plt.bar(xLabels,freqList1,width=-w,align='edge',label=label1,edgecolor=edgeList)
    plt.bar(xLabels,freqList2,width=w,align='edge',label=label2,edgecolor=edgeList)
    plt.legend()
    plt.show()
    return

 
'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):         # this function returns a 1D list where each element of the list is "black" if the corr esponding amino acid is in the biggestDiffs list and "white" otherwise.
    edge_colors=[]         # list with Black and white colors for index.
    amino_acid=[item[0] for item in biggestDiffs]
    for i in labels:
        if i in amino_acid:     # if amino_acid is present in biggestDiffs, black color edgelist is added.
            edge_colors.append("black")
        else:
            edge_colors.append("white")     # else , white color edgelist is added.

    return edge_colors


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():   # Main function to generate graphs,  flow of process to get the analysis
    human_protein=synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/human_p53.txt","Protening-Sequening/hw6-protein-starter/data/codon_table.json")
    elephant_protein=synthesizeProteins("Protening-Sequening/hw6-protein-starter/data/elephant_p53.txt","Protening-Sequening/hw6-protein-starter/data/codon_table.json")
    similarities=commonProteins(human_protein,elephant_protein)
    differences=findAminoAcidDifferences(human_protein,elephant_protein,0.005) # cutoff given as 0.5% which is 0.005 
    displayTextResults(similarities,differences)
    labels_on_x=makeAminoAcidLabels(human_protein,elephant_protein)
    frequency1=setupChartData(labels_on_x,human_protein)
    frequency2=setupChartData(labels_on_x,elephant_protein)
    edgelist=makeEdgeList(labels_on_x,differences)
    createChart(labels_on_x,frequency1,"Human Genes",frequency2,"Elephant Genes",edgelist)

    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    # ## Uncomment these for Week 2 ##
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    # print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    
    # print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    