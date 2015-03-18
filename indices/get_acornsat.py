import urllib

stationfile = open('ACORN_STN_list.txt','r')
stations = []
for n in range(112):
    stations.append(stationfile.readline()[:6])

for idn in stations:
    filename = 'acorn.sat.maxT.%s.daily.txt' %(idn)
    url = 'http://www.bom.gov.au/climate/change/acorn/sat/data/%s' %(filename)
    urllib.urlretrieve(url, filename)
    filename = 'acorn.sat.minT.%s.daily.txt' %(idn)
    url = 'http://www.bom.gov.au/climate/change/acorn/sat/data/%s' %(filename)
    urllib.urlretrieve(url, filename)

