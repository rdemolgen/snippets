import sys, subprocess, requests
from datetime import datetime

class AlamutVariant():

    def __init__(self, **entries):
        self.__dict__.update(entries)

class MpcVariant():
	
	def __init__(self, **entries):
	    self.__dict__.update(entries)

class AddMpcData():

    MPC_headers = ["chrom", "pos", "ref", "alt", "gene_name", "ENST", "ENSG", "CCDS" "Consequence",\
    	           "HGVSc", "HGVSp", "Amino_acid", "context", "SIFT", "PolyPhen", "obs_exp", "mis_badness",\
    	           "fitted_score", "MPC"]
    alamut_headers = ""
    alamut_instance_array = []

    def __init__(self, alamut_file, mpc_file):
    	self.MPC_file = mpc_file
    	self.alamut_file = alamut_file
        self.extract_to_instance(alamut_file)
        self.get_MPC_for_instances()
        self.write_instances_to_file()
        
    def extract_to_instance(self, input_file):
    	switch = 0
        with open(input_file, 'r') as inf:
        	for line in inf:
        		spl = line.strip('\n').split('\t')
        		if switch == 0:
        		    self.alamut_headers = spl
        		    switch += 1
        		else:
        		    self.alamut_instance_array.append(AlamutVariant(**{k:v for k,v in zip(self.alamut_headers, spl)}))
        print(str(len(self.alamut_instance_array)) + ' alamut variants were read into memory')
    
    def get_MPC_for_instances(self):
    	updated_instance_array = []
        for i in self.alamut_instance_array:      
            if i.codingEffect == 'missense':
                location = i.chrom + ':' + i.inputPos + '-' + i.inputPos 
                cmd = ['tabix', self.MPC_file, location]
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                tabix_match = process.communicate()[0].split('\n')
                tabix_mpc = self.filter_for_variant(tabix_match, i)
                if tabix_mpc:
                    i.__dict__['MPC'] = tabix_mpc.MPC
                    updated_instance_array.append(i)
                else:
                    i.__dict__['MPC'] = ''
                    updated_instance_array.append(i)
            else:
                i.__dict__['MPC'] = ''
                updated_instance_array.append(i)
        self.alamut_instance_array = updated_instance_array
        self.alamut_headers.append("MPC")

    def filter_for_variant(self, tabix_match, alamut_variant):
    	tabix_match = [i for i in tabix_match if i]
        for i in tabix_match:
            spl = i.split('\t')
            mpc = MpcVariant(**{k:v for k,v in zip(self.MPC_headers, spl)})
            if mpc.alt == alamut_variant.gNomen[-1]:
                if self.return_NM(mpc.ENST, alamut_variant.transcript):
                    return mpc
            else:
                continue

    def return_NM(self, ENST, required_match):
        server = "http://grch37.rest.ensembl.org/"
        ext = "/xrefs/id/" + ENST + "?object_type=transcript"
        r = requests.get(server + ext, headers={ "Content-Type" : "application/json"})
        json_r = r.json()
        for item in json_r:
            for k,v in item.iteritems():
            	if 'NM_' in str(v) or 'NM_' in str(k):
            	    if required_match in str(v):
            	        return True
    
    def write_instances_to_file(self):
        output = 'mpc_' + self.alamut_file
        with open(output, 'w') as ouf:
            ouf.write('{0}'.format('\t'.join(self.alamut_headers)) + '\n')
            for i in self.alamut_instance_array:
                for header in self.alamut_headers:
                    ouf.write(i.__dict__[header] + '\t')
                ouf.write('\n')    
        
if __name__ == "__main__":
    Add_mpc = AddMpcData(sys.argv[1], sys.argv[2])
    
