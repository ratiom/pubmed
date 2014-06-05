'''
Outline
journals - Define journals to search in
pubsearch - Download records from Pubmed (journals)
split - split records into fields (pubsearch)
gender - assign male of female and sort (split)
tidy - decide what fields to have in final (gender)
'''

from Bio import Entrez
from Bio import Medline
#is Biopython installed where this python executable is pointing to?
#this was set up to run under Anaconda (C:/Users/Conor/Anaconda)

import requests, json

def journals():
    journal_list = {"Angewandte":"0370543", "JACS":"7503056", "J Med Chem":"9716531"}

    jids = []
    for value in journal_list.values():
        jids.append(value)
    
    return jids


def getGenders(names):
    url = ""
    cnt = 0
    for name in names:
        if url == "":
            url = "name[0]=" + name
        else:
            url = url + "&name[" + str(cnt) + "]=" + name
        cnt += 1

    req = requests.get("http://api.genderize.io?" + url)
    results = json.loads(req.text)
    
    retrn = []
    for result in results:
        if result["gender"] != None:
            retrn.append((result["gender"], result["probability"], result["count"]))
        else:
            retrn.append((u'None',u'0.0',0.0))
        return retrn


    
    
def pubsearch():
    Entrez.email = "skippydoesntknow@gmail.com"
    #always let Entrez know who is calling    
    
    pubterm = ""
    for i in journals():
        pubterm += i+"[JID] or "
    
    IDhandle = Entrez.esearch(db="pubmed", term="peptide AND ("+pubterm+" ", mindate="2011", maxdate="2014", retstart=2500, retmax=100)
    #for documentation on esearch, see
    #http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
    #max number for retmax is 100k. Use retstart to get more than this.
    #Date range used to limit a search result by the date specified by datetype. These two parameters (mindate, maxdate) must be used together to specify an arbitrary date range. The general date format is YYYY/MM/DD, and these variants are also allowed: YYYY, YYYY/MM.
    
    record = Entrez.read(IDhandle)
    # record is returned as a dictionary. Lists search terms, all ID numbners etc
    
    idlist = record["IdList"]
    #return a list of ID numbers from the record dictionary
    
    recordHandle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    #search pubmed for records with idlist as input
    
    records = Medline.parse(recordHandle)
    #create dictionary from recordHandle
    
    return records
    
   

#comp: fields of interest get bundled into ta tuple, and then each tuple is appended to a list ("the list")
def pubmedRecords(records):
    pubmed_records = []
    for i in records:
        JournalID = i.get("JID", "?")
        authors = i.get("FAU", "?")
        abstracts = i.get("AB", "?")
        journal = i.get("JT", "?")
        tup = (JournalID, authors, journal)
        pubmed_records.append(tup)
    return pubmed_records
    #returns a list
   

    
def recordSort(records):
    pubmed_records = comp(records)
    # a new list called pubmed_records
    
    pubmed_auths = [l[1] for l in pubmed_records if l[0] in jids]
    #take the authors and fill up a new list
    
    lastFirstauthor_list = []
    #create list to hold the lastname, firstname group
    
    for auth in pubmed_auths:
        lastFirstauthor_list.append(auth[-1])
    
    # Go through records. If the words chemistry etc appear in the Journal Title (JT), grab the list of Author Full names (FAU) and append the last element of that list (hopefully the corresponding author)
        
    author_list_clean=[]
    for k in lastFirstauthor_list:
        if "," in k:
            last, first = k.split(',')
            author_list_clean.append(first[1:] + " " + last)
        
    first_names = []
    for i1 in author_list_clean:
        i2 = i1.split(" ")
        tup2 = (i1, (i2[0])) # gives name, firstname tuple. is this necessary?
        first_names.append(tup2) # new tuple (name, firstname)
    return first_names
    
    
pubsearch()

pubmedRecords(records)

first_names = recordSort(records)
print first_names

names = ["conor", "peter"]
for i in names:
    print getGenders(i)


#if __name__ == '__main__':
'''
l1=[] 
for elem in first_names:
    l1.append(list(elem))
    #new list containing names, firstnames

l2 = []
for elem in getGenders("conor", "piera"):
    l2.append(list(elem))
# new list containing "male, female"
    
z=zip(l1,l2)
#zip these into a three partite list of lists

znew = []
for i in z:
    znew.append(i[0] + i[1])
# just taking the names and first names   

wimmin = []
for i in znew:
    if "female" in i[2]:
        wimmin.append(i[0])
        
print wimmin
'''   