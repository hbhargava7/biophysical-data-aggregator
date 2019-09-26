#We want to create a script that will take as input a uniprot ID
#And will give as output a list of all publications from pubmed.

#To add
#API documentation for my functions, should be able to do help.
#Communication to the issue tracker, in a single Git Hub.
#Set retmax to 10


import Bio
from Bio import Entrez
import numpy as np
import pandas as pd
from habanero import counts
from altmetric import Altmetric


def searchPubMed(query):
    """This function takes the user input, the protein/term to search, and
    as output returns a dataframe with the results of the PubMed search.

    Parameters:
        query(str): protein or term to be searched in PubMed

    Returnds:
        Pandas dataframe, with PubMed IDs as indices and the following search
        parameters as columns:
            - Year of publication
            - Lead author
            - Senior author
            - Journal
            - DOI
            - WebLink based on DOI

    Example: queryDF = searchPubMed("PKC")"""

    #This gets a list of IDs of articles associated with the protein.
    Entrez.email = "adamo.mancino@ucsf.edu"  #Q14524  #"Nav1.5"
    handle = Entrez.esearch(db='pubmed', sort='Most Relevant', retmax = 100000, retmode='xml', term=query) #, field = "Title") #maximum allowed is 100,000 by Entrez documentation
    #This will return a dictionary, get the "ID List" value
    searchResults = Entrez.read(handle)
    id_list = searchResults["IdList"]

    if(len(id_list) == 0):
        noneDF = pd.DataFrame()
        return noneDF

    #Initializing some important lists...
    listPMID = []
    listYear = []
    listTitle = []
    listJournal = []
    listLeadAuthor = []
    listSeniorAuthor = []
    listWebLink = []
    listDOI = []

    #The efetch function takes as input a comma-separated string
    ids = ",".join(id_list)
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    fetchResults = Entrez.read(handle)

    n_articles = (len(fetchResults["PubmedArticle"]))

    #From the efetch, now let us take out some important parameters.
    for i in range(0, n_articles):
    #This is to extract the PubMed ID
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

    #This is to get the journal name.
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

    #Now take all those lists and combine them into a dataframe.
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


def getCitsNumber(inputDF):
    """This function reads in an already existing dataframe, typically the one
    obtained from the searchPubMed function. It returns said dataframe with an
    additional column, giving the number of citations for a given article DOI in the
    table according to CrossRef.

    Parameters:
        inputDF(Pandas.DataFrame): this is a Pandas dataframe with PubMed ID as indices,
        and a column labelled "DOI" with all the DOIs for the associated PubMed ID.

    Returns:
        Pandas.DataFrame, identical to input with addition of column "Citation Count"
        giving the number of citations registered in CrossRef for a given article.

    Example:
        queryDF = getCitsNumber(queryDF)"""

    #First gather the DOI for each article in the dataframe.
    doiList = inputDF["DOI"].tolist()
    listCitationsCount = []

    count = 1
    #Now for every DOI, use CrossRef to find the number of associated citations.
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

    #Now add a column of citation numbers to the dataframe.
    inputDF["Citation Count"] = listCitationsCount
    print("CrossRef Search Complete!")

    return inputDF



def getAltScore(inputDF):
    """This function reads in an already existing dataframe, typically the one
    obtained from the searchPubMed function. It returns said dataframe with an
    additional column, giving the Altmetric score for a given article DOI in the
    table according to Altmetrics.

    Parameters:
        inputDF(Pandas.DataFrame): this is a Pandas dataframe with PubMed ID as indices,
        and a column labelled "DOI" with all the DOIs for the associated PubMed ID.

    Returns:
        Pandas.DataFrame, identical to input with addition of column "Altmetric Score,"
        giving the Altmetric score, calculated by Altmetrics, for a given article.

    Example:
        queryDF = getAltScore(queryDF)"""


    #First gather the DOI for each article in the dataframe.
    doiList = inputDF["DOI"].tolist()
    listAltScore = []

    count = 1
    print("Getting Altmetrics Score")

    #Now for every DOI, use CrossRef to find the associated Altmetrics score.
    for entry in doiList:
        print(count)
        try:
            if entry != "":

                a = Altmetric()
                thisDOI = a.doi(entry)
                extractScore = thisDOI["score"]
                listAltScore.append(extractScore)

            else:
                listAltScore.append(None)
        except:
            listAltScore.append(None)

        count = count + 1

    #Now add a column of Altmetrics scores to the dataframe.
    inputDF["Altmetric Score"] = listAltScore
    print("Altmetrics Search Complete!")

    return inputDF



def multiJournalSearch(query, jStr):
    """This function allows one to iteratively run the function searchPubMed
    in the case when there are multiple journals to check. It takes as input the
    search term and a string with a comma-separated list of journals to check, and returns as
    output a dataframe similar to what is listed in searchPubMed documentation.

    Parameters:
        query(string) = protein or term to be searched in PubMed
        jStr(string) = comma-separated list of journal names

    Returns:
        Pandas dataframe, with PubMed IDs as indices and the following search
        parameters as columns:
            - Year of publication
            - Lead author
            - Senior author
            - Journal
            - DOI
            - WebLink based on DOI
        This differs from searchPubMed because the dataframe will be a combination
        of several dataframes, one for each journal in the input string jStr.

    Example:
        querryDF = multiJournalSearch("PKC", "Cell, Nature, Science")"""

    #First separate the comma-separated journal list into a list of journals.
    jList = jStr.replace(" ","").split(",")
    #Initialize a dataframe.
    DF = pd.DataFrame()
    initial = 0
    #For each journal in the list, now run searchPubMed and extract a dataframe for the search
    #term in that given journal. Append it to the existing dataframe.
    for j in jList:
        searchTerm = query + " " + j + "[Journal]"
        intDF = searchPubMed(searchTerm)
        if initial == 0 and not intDF.empty:
            DF = intDF
            initial = 1
        else:
            DF = DF.append(intDF)

    #Finally, return that dataframe.
    return DF



def searchCoordinator(searchTerm):
    """This function coordinates the search through PubMed, CrossRed, and Altmetrics.
    It takes as input a string containing the intended search query, and returns as output
    a dataframe. There are options to perform advanced searches, e.g. to search for a term
    in titles only, in titles & abstracts only, in one specific journal, in multiple specific
    journals, etc. But right now, they are set to the default: search for the term loosely
    throughout the article, restrict the journals to Cell, Nature, and Science.

    Parameters:
        searchTerm(str) = protein or term to be searched in PubMed

    Output:
        Pandas dataframe, with PubMed IDs as indices and the following search
        parameters as columns:
            - Year of publication *sorted in descending order*
            - Lead author
            - Senior author
            - Journal
            - DOI
            - WebLink based on DOI
            - CrossRef citations count
            - Altmetrics score

    Example:
        queryDF = searchCoordinator("Nav1.1")
    """


    #This lets user choose prompt.
    #prompt = "Please choose your protein of interest. "
    #querry = input(prompt)

    #Assume that the search querry is set for now...
    #searchTerm = "Nav1.1"
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

    #Designing for searches over several journals.
    multiJournal = 1
    if multiJournal == 1 and inclJournal != 1:
        print("Implementing CNS Search.")
        strList = "Nature, Science, Cell"
        queryDF = multiJournalSearch(searchTerm, strList)
    else:
        queryDF = searchPubMed(searchTerm)

    #This will organize the output dataframe by year of publication.
    queryDF = queryDF.sort_values(by = ["Year of Pub"], ascending = False)

    #This will get both the number of citations and the Altmetrics score, if requested.
    getCitCount = 1
    if getCitCount == 1:
        queryDF = getCitsNumber(queryDF)
        queryDF = getAltScore(queryDF)

    #Print and saved output
    if not queryDF.empty:
        print(queryDF)
        print(queryDF.shape)
        queryDF.to_csv("Test.csv")
    else:
        print("Sorry, this search yielded no results.")
    return queryDF



if __name__ == "__main__":
    searchCoordinator("Hsp90")
