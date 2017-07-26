#code to return the number of individuals covered at 20x for the genomic_coordinate submitted from the Exac API
#variables in square brackets should be altered as required
#useage:

#python exac_coverage_array.py [Chromosome] [genomic coordinate]

import sys, requests

def exac_coverage_array(chrom, coord):
    server = "http://exac.hms.harvard.edu/"
    headers = {"Content-Type" : "application/json"}
    coord_plus = int(coord) + 1
    ext = "rest/region/coverage_array/" + chrom + "-" + coord + "-" + str(coord_plus)
    print(server + ext)
    r = requests.get(server + ext, headers=headers)
    r_json = r.json()
    for item in r_json:
        if str(item['pos']) == (coord):
            if item['has_coverage'] == False:
                print(0)
                return '0'
            else:
                x = int(float(item['20'])*60706)
                print(x)
                return str(x)
        else:
            pass
            
x = exac_coverage_array(sys.argv[1], sys.argv[2])
