#We want to create a script that will take as input a uniprot ID
#And will give as output a list of all publications from pubmed.

import Bio
from Bio import Entrez
import numpy as np
import pandas as pd
import scholarly
from habanero import counts


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
    listJournal = []
    listLeadAuthor = []
    listSeniorAuthor = []
    listWebLink = []
    listDOI = []
#The efetch function takes as input a comma-separated string
#eFetch has a difficult time handling large requests, take 1000 at a time...
    ids = ",".join(id_list)
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    fetchResults = Entrez.read(handle)
#    print(fetchResults["PubmedArticle"][0]["PubmedData"]["ArticleIdList"][1])
#    quit()
#    s = (fetchResults["PubmedArticle"][0]["PubmedData"]["ArticleIdList"])
#    print(s[2].attributes["IdType"])
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


    #quit()

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
    #This is to get the journal name
        try:
            listJournal.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["Journal"]["Title"])
        except:
            listJournal.append(None)
    #This is to get the article title.
        try:
            listTitle.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["ArticleTitle"])
        except:
            listTitle.append(None)
    #This is to get the lead author's name.
        try:
            listLeadAuthor.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"] + " " + fetchResults["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"][0]["Initials"])
        except:
            listLeadAuthor.append(None)
    #This is to get the senior author's name.
        try:
            listSeniorAuthor.append(fetchResults["PubmedArticle"][i]["MedlineCitation"]["Article"]["AuthorList"][-1]["LastName"] + " " + fetchResults["PubmedArticle"][0]["MedlineCitation"]["Article"]["AuthorList"][-1]["Initials"])
        except:
            listSeniorAuthor.append(None)
    #This will get the DOI and thus the web link.
        try:
            refArray = fetchResults["PubmedArticle"][i]["PubmedData"]["ArticleIdList"]
            doiValue = ""
            for entry in refArray:
                if entry.attributes["IdType"] == "doi":
                    doiValue = entry
            if doiValue != "":
                listDOI.append(doiValue)
                listWebLink.append("https://doi.org/" + doiValue)
            else:
                listDOI.append(None)
                listWebLink.append(None)
        except:
            listDOI.append(None)
            listWebLink.append(None)
    #This will get the number of citations associated with a DOI


    outputDF = pd.DataFrame(index = listPMID)
    outputDF["Year of Pub"] = listYear
    outputDF["Title"] = listTitle
    outputDF["Lead Author"] = listLeadAuthor
    outputDF["Senior Author"] = listSeniorAuthor
    outputDF["Journal"] = listJournal
    outputDF["DOI"] = listDOI
    outputDF["Web Link"] = listWebLink

    print("PubMed Search Complete!")
    return outputDF

#Pass the list of IDs to efetch'
"""
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
"""

def getCitationCount(inputDF):
    doiList = inputDF["DOI"].tolist()
    listCitationsCount = []
    count = 1
    for entry in doiList:
        print(count)
        try:
            if entry != "":
                nCit = counts.citation_count(doi = entry)
                listCitationsCount.append(nCit)
            else:
                listCitationsCount.append(None)
        except:
            listCitationsCount.append(None)
        count = count + 1

    inputDF["Citation Count"] = listCitationsCount
    print("CrossRef Search Complete!")
    return inputDF


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
    #queryDF = searchGoogleScholar(searchTerm)


    multiJournal = 0
    if multiJournal == 1 and inclJournal != 1:
        strList = "Nature, Science, Cell"
        queryDF = multiJournalSearch(searchTerm, strList)
    else:
        queryDF = searchPubMed(searchTerm)

    queryDF = queryDF.sort_values(by = ["Year of Pub"], ascending = False)

    getCitCount = 1
    if getCitCount == 1:
        queryDF = getCitationCount(queryDF)

    #Print and saved output
    if not queryDF.empty:
        print(queryDF)
        print(queryDF.shape)
        queryDF.to_csv("Test.csv")
    else:
        print("Sorry, this search yielded no results.")
