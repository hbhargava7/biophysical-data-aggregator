import urllib2


url = 'http://www.rcsb.org/pdb/rest/search'

search = "Nav1.1"

queryText = """<?xml version="1.0" encoding="UTF-8"?><orgPdbQuery><queryType>org.pdb.query.simple.MoleculeNameQuery</queryType><description>Molecule Name Search : Molecule Name=""" + search + """</description><macromoleculeName>""" + search + """</macromoleculeName></orgPdbQuery>"""


print "query:\n", queryText

print "querying PDB...\n"

req = urllib2.Request(url, data=queryText)

f = urllib2.urlopen(req)

result = f.read()
print(result)
result = result.split("\n")
#del(result[-1])
print(len(result))
for name in result:
    print(name[0:4])
quit()

if result:

    print "Found number of PDB entries:", result.count('\n')

else:

    print "Failed to retrieve results"
