"""This script will enable the user to extract all structural information
avaiable for a given protein search."""

import urllib.request as urllib
import requests
import numpy as np
import pandas as pd


def getPDBIDs(input):
    """This function takes in a protein name, finds all related PDB IDs, and exports that as a list.
    Parameters:
        input(str): protein name or molecule name to be searched.
    Returns:
        list of unique protein PDB IDs associated to that particular name.
    Example:
        listPDBs = getPDBIDs(searchTerm)"""



    url = 'http://www.rcsb.org/pdb/rest/search'

    queryText = """<?xml version="1.0" encoding="UTF-8"?><orgPdbQuery><queryType>org.pdb.query.simple.MoleculeNameQuery</queryType><description>Molecule Name Search : Molecule Name=""" + input + """</description><macromoleculeName>""" + input + """</macromoleculeName></orgPdbQuery>"""
    queryBytes = str.encode(queryText)

    req = urllib.Request(url, data=queryBytes)

    f = urllib.urlopen(req)
    resultBytes = f.read()
    resultsText = str(resultBytes, 'utf-8')
    print(type(resultsText))

    resultList = resultsText.split("\n")
    if (resultList[-1] == ""):
        del(resultList[-1])
    updatedPDBIDs = []
    for name in resultList:
        updatedPDBIDs.append(name[0:4])

    #Get rid of all redundant IDs
    setPDBIDs = set(updatedPDBIDs)
    updatedPDBIDs = list(setPDBIDs)
    return updatedPDBIDs


def infoByID(inputIDs):
    """This function takes in a list of PDB IDs and finds all relevant PDB structures, tabulating all their information.
    Parameters:
        input(list): list of PDB IDs.
    Returns:
        Pandas dataframe, with the information for every PDB ID, organized in columns...
            - PDB ID
            - Structure Title
            - Resolution
            - Date of deposit
            - Method Used for structural biology
            - DOI
            - Weblink to protein in RCSB
    Example:
        listPDBs = getPDBIDs(searchTerm)"""


    listPDBID = []
    listStructureTitle = []
    listResolution = []
    listDepositDate = []
    listMethod = []
    listLigand = []
    listDOI = []
    listWebLink = []
    listAuthors = []

    for name in inputIDs:
        #Go in, and find another way to list, not by commas

        url = "http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=" + name + "&CustomReportColumns=structureTitle,resolution,depositionDate,experimentalTechnique,pdbDoi,structureAuthor&format=csv"
        rString = requests.get(url) #, allow_redirects=True)
        rList = str(rString.content).split(" />")

        parameterList = rList[0]
        del(rList[0])
        del(rList[-1])

        for entry in rList:
            #print(entry)
            #quit()
            entry = entry.replace("<br","")
            valuesList = entry.split('","')


            listPDBID.append(valuesList[0].replace('"',''))
            listStructureTitle.append(valuesList[1])
            listResolution.append(valuesList[2])
            listDepositDate.append(valuesList[3])
            listMethod.append(valuesList[4])
            listDOI.append(valuesList[5].replace('"',''))
            url = "https://www.rcsb.org/structure/" + valuesList[0]
            listWebLink.append(url)
            listAuthors.append(valuesList[6])

        #print(listPDBID[0] + ", " + listChains[0])

    outputDF = pd.DataFrame()
    outputDF["PDB ID"] = listPDBID
    outputDF["Structure ID"] = listStructureTitle
    outputDF["Resolution"] = listResolution
    outputDF["Date of Deposit"] = listDepositDate
    outputDF["Method Used"] = listMethod
    outputDF["DOI"] = listDOI
    outputDF["Web Link"] = listWebLink
    outputDF["Author"] = listAuthors

    return outputDF



def search(searchTerm):
    """This function takes in a protein name, and organizes a search for its information in the RCBI.
    It calls two functions, one that will convert the protein name to a list of PDB IDs (getPDBIDs)
    another will convert the list of PDB IDs to a dataframe containing the relevant structural information
    for each of those PDB IDs (infoByID). This function ultimately returns a dataframe.
    Parameters:
        input(str): protein name or molecule name to be searched.
    Returns:
        This is the same output as the function infoByID.
        It returns a Pandas dataframe, with the information for every PDB ID, organized in columns...
            - PDB ID
            - Structure Title
            - Resolution
            - Date of deposit
            - Method Used for structural biology
            - DOI
            - Weblink to protein in RCSB
    Example:
        listPDBs = getPDBIDs(searchTerm)"""


    listPDBs = getPDBIDs(searchTerm)
    print(listPDBs)
    if (len(listPDBs) != 0):
        DF = infoByID(listPDBs)
    else:
        DF = pd.DataFrame()
    return(DF)
