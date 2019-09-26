import urllib.request as urllib
import requests
import numpy as np
import pandas as pd


def getPDBIDs(input):
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
    listPDBID = []
    listStructureTitle = []
    listResolution = []
    listDepositDate = []
    listMethod = []
    listLigand = []
    listDOI = []
    listWebLink = []

    for name in inputIDs:
        #Go in, and find another way to list, not by commas

        url = "http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=" + name + "&CustomReportColumns=structureTitle,resolution,depositionDate,experimentalTechnique,pdbDoi&format=csv"
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

        #print(listPDBID[0] + ", " + listChains[0])

    outputDF = pd.DataFrame()
    outputDF["PDB ID"] = listPDBID
    outputDF["Structure ID"] = listStructureTitle
    outputDF["Resolution"] = listResolution
    outputDF["Date of Deposit"] = listDepositDate
    outputDF["Method Used"] = listMethod
    outputDF["DOI"] = listDOI
    outputDF["Web Link"] = listWebLink

    return outputDF



if __name__ == "__main__":
    searchTerm = "AMPA"
    listPDBs = getPDBIDs(searchTerm)
    print(listPDBs)
    DF = infoByID(listPDBs)
    DF = DF.sort_values(by = ["Date of Deposit"], ascending = False)
    print(DF)
    DF.to_csv("PDB_Search.csv")
