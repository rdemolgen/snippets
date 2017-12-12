import requests, json, sys, re, time, mechanicalsoup
from bs4 import BeautifulSoup

class Uniprot_api():

    def __init__(self):
        self.upload_ext = 'http://www.uniprot.org/uploadlists/'
        self.entry_ext = 'http://www.uniprot.org/uniprot/'

    #given an ensembl_id returns the uniprot responses as a list of dictionaries
    def get_pid_from_gene(self, ensembl_id):
        params = {'from':'ENSEMBL_ID','to':'ACC','format':'tab','query': ensembl_id }
        r = requests.get(self.upload_ext, params=params)
        dict_entry_list = []
        split_r = r.text.split('\n')
        headers = split_r[0].split('\t')
        for entry in split_r:
            if entry:
                entry_list = entry.split('\t')
                if entry_list != headers:
                    dict_entry_list.append(dict(zip(headers, entry_list)))
        return dict_entry_list

    def get_entry_gff(self, pid):
        r = requests.get(self.entry_ext + str(pid) + '.gff')
        return r

    def parse_gff(self, gff_response, required_list): 
        gff_objects = []
        gff_list = gff_response.split('\n')
        for annotation in gff_list:
            spl_anno = annotation.split('\t')
            if spl_anno != [''] and '#' not in spl_anno[0]:
                if spl_anno[2] in required_list:
                    gff_objects.append(Gff_object(spl_anno))
        return gff_objects  

class Gff_object():
  
    def __init__(self, gff_line):
        self.ID = gff_line[0]
        self.db_annotation = gff_line[1]
        self.anno_type = gff_line[2]
        self.start = gff_line[3]
        self.stop = gff_line[4]

class Ensembl_api():

    def __init__(self):
        self.server = "http://grch37.rest.ensembl.org/"

    #returns ensembl id from HGNC gene symbol
    def query_HGNC(self, gene_symbol):
        ext = "/xrefs/symbol/homo_sapiens/" + gene_symbol + "?external_db=HGNC"
        r = requests.get(self.server + ext, headers={ "Content-Type" : "application/json"})    
        json_r = r.json()
        for transcript in json_r:
            if "ENST" in transcript["id"]:              
                print("Ensembl Transcript ID: "+ transcript["id"])
        for reference in json_r:
            if "ENSG" in reference["id"]:              
                print("Ensembl Gene ID for " + gene_symbol + ": " + reference["id"] + "\n")
                return reference["id"]
                
    def tester(self, string):
        print(string)
                
class Exac_api():

    def __init__(self):
        self.server = "http://exac.hms.harvard.edu/"
        self.json_headers = {"Content-Type" : "application/json"}

    def variants_in_gene(self, ensembl_id):
        ext = "/rest/gene/variants_in_gene/" + ensembl_id
        r = requests.get(self.server + ext, headers=self.json_headers)
        json_r = r.json()
        return json_r
    
    def canonical_transcript(self, ensembl_id):
        ext = "/rest/gene/" + ensembl_id
        r = requests.get(self.server + ext, headers=self.json_headers)
        json_r = r.json()
        return json_r

    def variants_in_region(self, chr, start, stop):
        ext = "/rest/region/variants_in_region/" + chr + '-' + start + '-' + stop
        r = requests.get(self.server + ext, headers=self.json_headers)
        json_r = r.json()
        return json_r

   #filters variants in gene with a given key and value. keep or remove option
    def filter_variants(self, variant_list, key, value, remove=False):
        pass_vars = []
        for variant in variant_list:
            if remove == False:
                if variant[key] == value:
                    pass_vars.append(variant)
            #only keeps those not meeting the criteria
            else:
                if variant[key] == value:
                    pass
                else:
                    pass_vars.append(variant)
        return pass_vars

    #provide the starting list and dict {key : value} for each to be filtered.
    def filter_by_dict(self, starting_list, filter_dict):
        iteration_list = []
        count = 0
        for k,v in filter_dict.items():
            for variant_type in v:
                if count == 0:
                    this_filter = self.filter_variants(starting_list, k, variant_type, remove=True)
                    iteration_list.append(this_filter)
                    count += 1
                elif count > 0:
                    this_filter = self.filter_variants(iteration_list[-1], k, variant_type, remove=True)
                    iteration_list.append(this_filter)
        return iteration_list[-1]

    def position_frequency(self, var_list, homo=False, hemi=False):
        if homo==True:
            var_freq_pos_dict = self.dict_extractor('homo_freq', 'homo_pos', var_list)
        elif hemi==True:
            var_freq_pos_dict = self.dict_extractor('hemi_freq', 'hemi_pos', var_list)
        else:
            var_freq_pos_dict = self.dict_extractor('het_freq', 'het_pos', var_list)
        return var_freq_pos_dict

    def dict_extractor(self, freq_key, pos_key, var_list):
        freq_pos = { freq_key : [],
                     pos_key : []}
        for variant_dict in var_list:
            freq_pos[freq_key].append(variant_dict["allele_freq"])
            freq_pos[pos_key].append(self.extract_protein_position(variant_dict["HGVSp"]))
        return freq_pos


    #takes protein position from HGVS protein nomenclature for syn/miss
    def extract_protein_position(self, HGVSp):
        m = re.search(r'(p\.)([^0-9]*)(\d{1,})(.*)', HGVSp)
        #if matches returns the format
        if m:
            return m.group(3)
        else:
            return HGVSp

#creates class based on key, value **kwargs
class HGMD_variant():
    
    def __init__(self, **entries):
        self.__dict__.update(entries)

class HGMD_pro():

    #could make the object the login and methods what to do after
    def __init__(self, gene_name):
        self.gene = gene_name

    def scrape_HGMD_all_mutations(self,hgmd_username,hgmd_password):
        browser = mechanicalsoup.Browser()
        login_page = browser.get("http://portal.biobase-international.com/cgi-bin/portal/login.cgi")
        time.sleep(1)
        login_form = mechanicalsoup.Form(login_page.soup.select_one('#login_form'))
        time.sleep(2)
        # login username and user_password required as strings
        #login_form.input({"login" : "username", "password" : "user_password"})
        login_form.input({"login" : hgmd_username, "password" : hgmd_password})
        time.sleep(2)
        r = browser.submit(login_form, login_page.url)
        ##print(r.text)
        time.sleep(1)
        try:
            soup = self.form_finder(browser, self.gene)
        except:
            print('HGMD exception executed')
            browser = mechanicalsoup.Browser()
            login_page = browser.get("http://portal.biobase-international.com/cgi-bin/portal/login.cgi")
            time.sleep(1)
            login_form = mechanicalsoup.Form(login_page.soup.select_one('#login_form'))
            time.sleep(2)
            # login username and user_password required as strings
            #login_form.input({"login" : "username", "password" : "user_password"})
            login_form.input({"login" : hgmd_username, "password" : hgmd_password})
            time.sleep(2)
            r = browser.submit(login_form, login_page.url)
            soup = self.form_finder(browser, sys.argv[3])
            print(soup)
        return soup
      
    def form_finder(self, browser, gene):
        gene_search = browser.get("https://portal.biobase-international.com/hgmd/pro/gene.php?gene=" + gene)
        time.sleep(0.5)
        soup = BeautifulSoup(gene_search.content)
        form = soup.find("form", attrs={ "action" : "all.php" })
        gene_id_element = form.find("input", attrs={"name" : "gene_id"})
        gene_id_value = gene_id_element['value']
        trans_element = form.find("input", attrs={"name" : "refcore"})
        HGMD_transcript_id = trans_element['value']
        print("HGMD transcript ID: " + HGMD_transcript_id)     
        params = {"gene" : self.gene, "inclsnp" : "N", "base" : "Z", "refcore" : HGMD_transcript_id, "gene_id" : gene_id_value, "database" : "Get all mutations"}
        url = "https://portal.biobase-international.com/hgmd/pro/all.php"
        response = browser.post(url, data=params)
        soup = BeautifulSoup(response.content)
        #print(soup)
        return soup
        
        ###################################remember to un comment the above######################################################################
        
    def extract_missense_nonsense(self, all_mutations_soup):
    	#may need to add HGMD accession to headers when it's live
        HGMD_headers = ["HGMD_codon_change", "amino_acid_change", "HGVS_nuc", "HGVS_prot", "var_class", "phenotype","reference", "additional"]
        HGMD_headers_legacy = ["HGMD_codon_change", "amino_acid_change", "legacy_change", "HGVS_nuc", "HGVS_prot", "var_class", "phenotype","reference", "additional"]
        HGMD_var_objs = []
        soup = all_mutations_soup
        #soup = BeautifulSoup(open("html_to_parse"), "html.parser")
        table = soup.find("table", attrs={ "class" : "gene" })
        rows = table.find_all('tr')
        num_rows = 0
        for row in rows:
            cols = row.findAll("td")
            cols = [ele.text.strip() for ele in cols]
            #remove empty strings
            cols = [ele for ele in cols if ele]
            if len(cols) == 8:
                variant_dict = {key:value for key,value in zip(HGMD_headers, cols)}
                var_instance = HGMD_variant(**variant_dict)
                HGMD_var_objs.append(var_instance)
            elif len(cols) == 9: # if a legacy naming column exists, for example gene INSR
                variant_dict = {key:value for key,value in zip(HGMD_headers_legacy, cols)}
                var_instance = HGMD_variant(**variant_dict)
                HGMD_var_objs.append(var_instance)
            elif len(cols) == 7:
                variant_dict = {key:value for key,value in zip(HGMD_headers, cols)}
                if 'additional' not in variant_dict:
                	variant_dict['additional'] = 'None'
                elif additional in variant_dict:
                	continue
                var_instance = HGMD_variant(**variant_dict)
                HGMD_var_objs.append(var_instance)
            else:
                if num_rows != 0:
                    print("Columns in HGMD table = " + str(len(cols)) + ". Cannot handle this many columns")
            num_rows =+1
        return HGMD_var_objs

    def write_DM_file(self, variant_instance_list):
        with open("temp_hgmd_file", 'w') as ouf:
            for item in variant_instance_list:
                try:
                    m = re.search(r'(p\.)([^0-9]*)(\d{1,})(.*)', item.hgvs_prot)
                    if m:
                        if item.variant_class == "DM":
                            string = m.group(3) + '\t' + '1\n'
                        elif item.variant_class == "DM?":
                            string = m.group(3) + '\t' + '\t' + '2\n'
                        else:
                            string = m.group(3) + '\t' + '\t' + '\t' + '3\n'
                        ouf.write(string)
                    else:
                    	print("not_m " + item.hgvs_prot)
                except:
                    pass
        return "temp_hgmd_file"


if __name__ == "__main__":
    Ex = Exac_api()
    variants_in_gene_json = Ex.variants_in_gene("ENSG00000077279")
    print(variants_in_gene_json)

