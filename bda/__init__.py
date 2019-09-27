from flask import Flask, render_template, redirect, request

import PubMedSearch as pms
import PDB as pdb

app = Flask(__name__)

@app.route('/')
def main():
    return render_template('index.html')

@app.route('/startQuery', methods = ['POST'])
def query():
    _query = request.form['inputQuery']
    print('pubmed query started')
    _result = pms.searchCoordinator(_query)
    _resultDict = _result.T.to_dict().values()
    print('pubmed query completed with ' + str(len(_resultDict)) + ' results.')

    print('PDB query started')

    _pdbResult = pdb.search(_query)
    _pdbDict = _pdbResult.T.to_dict().values()

    print('PDB query completed with ' + str(len(_pdbDict)) + ' results.')

    _summaryDict = {}
    _summaryDict["numStructures"] = str(len(_pdbDict))
    _summaryDict["numPapers"] = str(len(_resultDict))
    
    _pdbResult.dropna()
    _highResStructure = _pdbResult.loc[_pdbResult['Resolution'].replace("","1000").astype(float).idxmin()]
    
    _summaryDict["highResID"] = _highResStructure["PDB ID"]
    _summaryDict["highResolution"] = _highResStructure["Resolution"]
    
    _summaryDict["bestAuthor"] = _result["Senior Author"].mode()[0]
    
    return render_template('results.html', query=_query, result=_resultDict, pdb=_pdbDict,summary=_summaryDict)


if __name__ == '__main__':
    app.run()   
