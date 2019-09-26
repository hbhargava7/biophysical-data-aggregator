from flask import Flask, render_template, redirect, request

import PubMedSearch as pms

app = Flask(__name__)

@app.route('/')
def main():
    return render_template('index.html')

@app.route('/startQuery', methods = ['POST'])
def query():
    _query = request.form['inputQuery']
    print('pubmed query started')
    _result = pms.searchPubMed(_query)
    _resultDict = _result.T.to_dict().values()
    print('pubmed query completed with ' + str(len(_resultDict)) + ' results.')
    return render_template('results.html', query=_query, result=_resultDict)

if __name__ == '__main__':
    app.run()   
