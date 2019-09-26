import urllib.request as urllib


def getPDBIDs(input):
    url = 'http://www.rcsb.org/pdb/rest/search'

    input = "Nav1.1"

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

    return updatedPDBIDs


if __name__ == "__main__":
    searchTerm = "Hsp90"
    listPDBs = getPDBIDs(searchTerm)
    print(listPDBs)

    
