#We want to create a script that will take as input a uniprot ID
#And will give as output a list of all publications from pubmed.

import Bio
from Bio import Entrez
import numpy as np
import pandas as pd
import scholarly


def searchPubMed(query):
#This gets a list of IDs of articles associated with the protein.
    Entrez.email = "adamo.mancino@ucsf.edu"  #Q14524  #"Nav1.5"
    handle = Entrez.esearch(db='pubmed', sort='Most Recent', retmax = 100000, retmode='xml', term=query) #, field = "Title") #maximum allowed is 100,000 by Entrez documentation
    #This will return a dictionary, get the "ID List" value
    searchResults = Entrez.read(handle)
    id_list = searchResults["IdList"]
#    n_articles = len(id_list)
    #print(n_articles)
    #quit()
    if(len(id_list) == 0):
        noneDF = pd.DataFrame()
        return noneDF

    listPMID = []
    listYear = []
    listTitle = []
#The efetch function takes as input a comma-separated string
#eFetch has a difficult time handling large requests, take 1000 at a time...
    ids = ",".join(id_list)
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    fetchResults = Entrez.read(handle)
#    print(fetchResults["PubmedArticle"][0])
#    quit()


#It returns a dictionary of lits.
#Use the key "PubmedArticle" to get the list of articles.
#Each entry in the list is a dictionary.

    #Initialize lists of PMIDs, years, and titles

    #print(fetchResults["PubmedArticle"][0]["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"])   #["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"]))
    #quit()
    n_articles = (len(fetchResults["PubmedArticle"]))
#    quit()
#    print(fetchResults["PubmedArticle"][920:927])
#    quit()

    #For each entry...
        #print(fetchResults)
    #This is to extract PMID, which is redundant but I don't know if the order changed at all...
        #listPMID.append(
    for i in range(0, n_articles):
        try:
            listPMID.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["PMID"])
        except:
            listPMID.append(None)

    #This is to get the YYYY-MM-DD of the article's publication.
    #First try using the journal issue date.
    #If that fails, then try the publication date.
    #If that fails, then try the completed date.
    #Eventually, one of the correct dates will come up...
        try:
            listYear.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]["Year"])
        except:
            try:
                listYear.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"])
            except:
                try:
                    listYear.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["DateRevised"]["Year"])
                except:
                    print("I cri. The problematic index is: " + str(i))
                    listYear.append(None)
            #print(fetchResults["PubmedArticle"][i])
            #quit()

    #This is to get the article title.
        try:
            listTitle.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["ArticleTitle"])
        except:
            listTitle.append(None)

    outputDF = pd.DataFrame(index = listPMID)
    outputDF["Year of Pub"] = listYear
    outputDF["Title"] = listTitle
    return outputDF

#Pass the list of IDs to efetch'

def searchGoogleScholar(query):

    #print(scholarly.search_author("Marty Banks"))

    search_query = scholarly.search_pubs_query(query)
    nCount = 0
    for x in search_query:
        nCount = nCount + 1
        print(nCount)
    print("Final count is : " + str(nCount))

#    nCount = 0
#    for entry in search_query:
#        print(next(entry))
#        nCount =+ 1
#    print(next(search_query))
#    quit()

def multiJournalSearch(query, jStr):
    jList = jStr.replace(" ","").split(",")
    DF = pd.DataFrame()
    initial = 0
    for j in jList:
        searchTerm = query + " " + j + "[Journal]"
        intDF = searchPubMed(searchTerm)
        if initial == 0 and not intDF.empty:
            DF = intDF
            initial = 1
        else:
            DF = DF.append(intDF)
    DF = DF.sort_values(by = ["Year of Pub"], ascending = False)
    return DF



if __name__ == "__main__":
    #Lets user choose prompt.
    #prompt = "Please choose your protein of interest. "
    #querry = input(prompt)

    #Assume that the search querry is set for now...
    searchTerm = "Nav1.1"
    #Looking for simple search parameters.
    inTitle = 0
    if (inTitle == 1):
        searchTerm = searchTerm + "[Title]"
    inAbstract = 0
    if (inAbstract == 1):
        searchTerm = searchTerm + "[Title/Abstract]"
    inclJournal = 0
    if (inclJournal == 1):
        journal = "Science"
        searchTerm = searchTerm + " " + journal + "[Journal]"
    #querryDF = searchGoogleScholar(searchTerm)


    multiJournal = 0
    if multiJournal == 1 and inclJournal != 1:
        strList = "Nature, Science, Cell"
        querryDF = multiJournalSearch(searchTerm, strList)
    else:
        querryDF = searchPubMed(searchTerm)

    #output
    if not querryDF.empty:
        print(querryDF)
        print(querryDF.shape)
        querryDF.to_csv("Test.csv")
    else:
        print("Sorry, this search yielded no results.")
