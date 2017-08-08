from api import *
import re, sys, subprocess, json, os
import numpy as np
import pandas as pd
from collections import defaultdict

#the purpose of this class is to collect the data required of a plot and then write to a file for retention and subsequent plotting
#there is also a metadata file containing logs, gene name provided, ensembl_id used 
#protein_residue (1-length of residue,  

class Graph_object():
	
    DM_likely_objs = {}
    DM_objs = {}
    domain_count = 0
    HGMD_DM_track_count = 0
    HGMD_DMq_track_count = 0
    total_phen_count = 0
    phen_index = {}
    uniprot_columns = []
    slice_start = 0
    slice_end = 0
	
    def __init__(self, gene_name, user_pos):
        self.gene = gene_name
        self.user_pos = user_pos
    	#Ensembl#############################
        self.Ens = Ensembl_api()
        print('Getting Ensembl geneID for ' + gene_name)
        print(' ')
        self.ensembl_id = self.get_ensemble_id(gene_name)
        #Uniprot##########################
        self.Up = Uniprot_api()
        #composite file_name to add data to as it is parsed from different sources
        self.plotting_file = gene_name + 'composite.data'
        self.zoomed_plot = gene_name + 'zoomed_composite.data'
        #dict of uniprot entries for this ensembl transcript
        print('Getting Uniprot annotation file') 
        self.all_uniprot_entries = self.get_uniprot_entries(self.ensembl_id)
        #returns human reviews uniprot annotations
        self.reviewed_uniprot_entries = self.reviewed_human_entries(self.all_uniprot_entries)
        ##print(self.reviewed_uniprot_entries)
        self.length = self.reviewed_uniprot_entries[0]['Length']
        #returns text of all gff annotations
        self.all_gff_annotation = self.get_gff_domains()
        #Gff objects where values are ranges of 'Domain', 'Region', 'DNA binding'
        self.requried_gff_annotations = self.specific_gff_annotations()
        #generates dict of arrays to be written to master gnuplot file
        self.plottable_domains = self.generate_plottable_domains()
        
        #Exac##############################
        print('Gathering exac data using ensemble geneID')
        print(' ')
        self.Ex = Exac_api()
        self.all_exac_variants = self.get_exac_data(self.ensembl_id)
        
        #Consurf###########################
        print('Gathering Consurf data stored on the server')
        print(' ')
        self.consurf_file = self.find_consurf_file(gene_name)
        #dict of consurfscore and position
        self.consurf_data = self.parse_consurf_grades(self.consurf_file)
        #datafram of all data
        self.write_consurf_data = self.write_consurf_grades(self.consurf_data)
        
        #HGMD##############################
        print('Gathering HGMD data from website')
        self.hgmd_data = self.get_HGMD_data(gene_name)             
        self.write_DM_data = self.write_HGMD_data(self.DM_objs)
        self.write_DM_data2 = self.write_HGMD_data(self.DM_likely_objs, DM=False)
        
        
        #Gnomad_data######################
        #self.get_gnomad(
        
        #GNU plotter#######################
        print('Plotting all data')
        self.execute_gnuplot(gene_name)
        #self.create_smaller_graph_file()
        #self.execute_zoomed_gnuplot(gene_name)
        

    ### Ensembl_id   ####################################################
    def get_ensemble_id(self, HGNC_gene):
        id = self.Ens.query_HGNC(HGNC_gene)
        return id

    ### Uniprot_data #####################################################
    def get_uniprot_entries(self, ensembl_id):
        Up = Uniprot_api()
        up_entry_list = Up.get_pid_from_gene(ensembl_id)
        return up_entry_list
    
    def reviewed_human_entries(self, entry_list):
        human_match_list = []
        for entry_dict in entry_list:
            try:
                if entry_dict['Organism'] == 'Homo sapiens (Human)' and entry_dict['Status'] == 'reviewed':
                    human_match_list.append(entry_dict)
            except:
               print('check reviewd trascript is identified. see interpretation_assist.py reviewed_human_entries /n returned list :')
               print(human_match_list)
        return human_match_list
        
    def get_gff_domains(self):
        if len(self.reviewed_uniprot_entries) == 1:
            gff_response = self.Up.get_entry_gff(self.reviewed_uniprot_entries[0]["Entry"])
            return gff_response.text
        elif len(human_match_list) == 0:
            print('no_up_id_matches -- see interpretation_assist.py gff_domains')
        elif len(human_match_list) > 1:
            print('more than 1 reviewd human entry interpretation_assist.py gff_domains')

    def specific_gff_annotations(self):
    	required_list = ["Domain", "Region", "DNA binding", "Zinc finger", "Motif", "Transmembrane", "Intramembrane", "Compositional bias", "Topological domain"]
    	gff_objects = self.Up.parse_gff(self.all_gff_annotation, required_list)
    	#print(self.all_gff_annotation)
    	return gff_objects
    	
    def generate_plottable_domains(self):
        length = self.reviewed_uniprot_entries[0]['Length']
        master_dict = {}
        domain_array = {}
        annotation_index = 0
        arrays_to_save = []
        #generate a dictionary from the required gff objects 
        for item in self.requried_gff_annotations:
            if item.__dict__['anno_type'] in master_dict.keys():
                master_dict[item.__dict__['anno_type']].append(item.__dict__)
            else:
                master_dict[item.__dict__['anno_type']] = [item.__dict__]
        for k in master_dict.keys():
            self.uniprot_columns.append(k)
            annotation_index += 1
            #function to return the array on NA or integer
            this_type_column = self.array_creator(master_dict, k, annotation_index, length)
            np_array = np.array(this_type_column, dtype=np.float)
            arrays_to_save.append(np_array)    
        result_array = np.empty((0, int(length)))
        for item in arrays_to_save:
            result_array = np.append(result_array, [item], axis=0)
        with open(self.plotting_file, 'wb') as f:
            headers = '{0}'.format('\t'.join(self.uniprot_columns))
            f.write(bytes(headers, 'utf-8') + bytes('\n', 'utf-8'))
            np.savetxt(f, np.transpose([result_array]), delimiter='\t')
        self.domain_count = len(master_dict.keys())
        return arrays_to_save
    
    def array_creator(self, master_dict, key, annotation_index, length):
        ##print('array creator initiated for ' + key)
        this_annotation_array = []
        domain_count = 0
        position = 0
        ##print('iterating through annotations for ' + key)
        master_dict_iter = iter(master_dict[key])
        for v in master_dict[key]:
            #print(v)
            if domain_count == 0:
                interval_length = int(v['start']) - 1  
                ##print("interval " + str(interval_length))
                domain_range = int(v['stop']) - int(v['start']) + 1
                ##print("domain " + str(domain_range))
                domain_count += 1
                for i in range(interval_length):
                    position += 1
                    this_annotation_array.append(None)
                for i in range(domain_range):
                    position += 1
                    if domain_range == 1:
                        this_annotation_array.append(annotation_index)
                        this_annotation_array.append(annotation_index)
                    else:
                        this_annotation_array.append(annotation_index)
                        #next(master_dict_iter, None)
            elif domain_count > 0:
            	#to calculate the gap from last domain to the next
                interval_length = int(v['start']) - position
                ##print("interval " + str(interval_length))
                domain_range = int(v['stop']) - int(v['start']) + 1
                ##print("domain " + str(domain_range))
                domain_count += 1
                for i in range(interval_length):
                    position += 1
                    this_annotation_array.append(None)
                for i in range(domain_range):
                    position += 1
                    if domain_range == 1:
                        this_annotation_array.append(annotation_index)
                        this_annotation_array.append(annotation_index)
                    else:
                        this_annotation_array.append(annotation_index)
        final_interval = int(length) - position
        for i in range(final_interval):
            position += 1
            this_annotation_array.append(None)
        if len(this_annotation_array) > position:
            slice_amount = int((len(this_annotation_array) - position) / 2)
            #print(slice_amount)
            start = 0 + slice_amount
            #print(start)
            mod_annotation_array = this_annotation_array[start:-slice_amount]
            return mod_annotation_array
        return this_annotation_array

    ### Exac data  ##################################################################
    #gets exac variants in the ensembl gene
    #filters to retain exac pass filter, missense, initiator_codon_variant, stop_gained, frameshift_variant, inframe_deletion, stop_lost, canonical
    def get_exac_data(self, ensembl_id):
        filter_dict = {"major_consequence" : ["synonymous_variant",
                                              "splice_region_variant",
                                              "intron_variant",
                                              "splice_acceptor_variant",
                                              "5_prime_UTR_variant",
                                              "3_prime_UTR_variant",
                                              "non_coding_transcript_exon_variant",
                                              "splice_donor_variant",
                                              "inframe_deletion"] }
        all_variants = self.Ex.variants_in_gene(ensembl_id)
        pass_variants = self.Ex.filter_variants(all_variants, "filter", "PASS")
        canonical_transcript = self.Ex.filter_variants(pass_variants, "CANONICAL", "YES")
        missense_only = self.Ex.filter_by_dict(canonical_transcript, filter_dict)
        #list of homo entries
        homozygotes = self.Ex.filter_variants(missense_only, "hom_count", 0, remove=True)
        #list of hetero entries
        heterozygotes = self.Ex.filter_variants(missense_only, "hom_count", 0)
        #dicts of homo and het frequency and position
        homo_freq_pos = self.Ex.position_frequency(homozygotes, homo=True)
        het_freq_pos = self.Ex.position_frequency(heterozygotes)
        exac_to_composite = self.add_exac_to_composite(het_freq_pos, indexed=False)
        exac_to_composite = self.add_exac_to_composite(homo_freq_pos)
        return all_variants
    
    #indexed to indicate whether there is an index for pandas in column 0
    def add_exac_to_composite(self, freq_pos_data, indexed=True):
        if indexed:
            df = pd.read_csv(self.plotting_file, delimiter='\t', index_col=0)
        else:
            df = pd.read_csv(self.plotting_file, delimiter='\t')
        pos = [k for k,v in freq_pos_data.items() if 'pos' in k.lower()]
        pos_column = pd.DataFrame({pos[0] : freq_pos_data[pos[0]]})
        df = pd.concat([df, pos_column], axis=1)
        freq = [k for k,v in freq_pos_data.items() if 'freq' in k.lower()]
        freq_column = pd.DataFrame({freq[0] : freq_pos_data[freq[0]]})
        df = pd.concat([df, freq_column], axis=1)
        df.to_csv(self.plotting_file, sep='\t', na_rep='?')       

    ### Consurf data #################################################################
    def find_consurf_file(self, gene_name):
        #print(gene_name)
        base = "/home/cannons/interpretation_graphs/consurf_scores/"
        for dirs, subdirs, files in os.walk(base):
            if dirs.endswith(gene_name):
                print(os.path.join(dirs, files[0]))
                return os.path.join(dirs, files[0])

    def parse_consurf_grades(self, consurf_file):
        cons_pos = { "cons" : [],
    	             "pos" : []}
        if consurf_file:
            with open(consurf_file, 'r') as inf:
                count = 0
                for line in inf:
                    spl = line.strip(' ').split('\t')
                    if len(spl) > 10 and 'COLOR' not in spl[4]:
                        count += 1
                        cons_pos["pos"].append(count)
                        if '*' in spl[4]:
                            cons_pos["cons"].append(0)
                        else:
                            cons_pos["cons"].append(spl[4].strip(' '))
        return cons_pos                        	      

    #opens data, creates a new df with conservation, running mean and cumulative mean. Then adds position to be used with gnuplot
    def write_consurf_grades(self, cons_pos_dict, indexed=True):
        df = pd.read_csv(self.plotting_file, delimiter='\t', index_col=0)
        cons_df = pd.DataFrame({'cons' : cons_pos_dict['cons']})
        cons_df['cons_MA'] = cons_df.rolling(window=3).mean()
        cons_df['cons_cumsum'] = cons_df['cons'].expanding().mean()
        df = pd.concat([df, cons_df], axis=1)
        pos_df = pd.DataFrame({'pos' : cons_pos_dict['pos']})
        df = pd.concat([df, pos_df], axis=1)
        df.to_csv(self.plotting_file, sep='\t', na_rep='?')
        return df    
    
### HGMD data  ###################################################################
    def get_HGMD_data(self, gene_name):
        HGMD = HGMD_pro(gene_name)
        all_mutations_soup = HGMD.scrape_HGMD_all_mutations()
        #list of objects
        variant_instances = HGMD.extract_missense_nonsense(all_mutations_soup)
        #returns dict of mutation class and list of objects in that class
        separate_var_class = self.var_class_separator(variant_instances)
        for var_class,objs in separate_var_class.items():
            #separate_each_list_by returning dictionaries
            if var_class == 'DM':
                self.DM_objs = self.phenotype_lists(var_class, objs)
            elif var_class == 'DM?':
                self.DM_likely_objs = self.phenotype_lists(var_class, objs)
        return variant_instances
        #write_hgmd = HGMD.write_DM_file(variant_instances)

    def write_HGMD_data(self, obj_d, DM=True):
        #for k,v in obj_d.items():
    	#    print(k,v)
        df = pd.read_csv(self.plotting_file, delimiter='\t', index_col=0)
        #phen_index where { phen : count }
        track = []
        phenotype = []
        position = []
        if DM:
            count = 0
            for phen, var_list in obj_d.items():
                var_pos_list = self.var_pos_extractor(var_list)
                if var_pos_list != []:
                    #generate columns
                    count += 1
                    count_index = [count] * len(var_pos_list)
                    track = track + count_index
                    phen_length = [phen] * len(var_pos_list)
                    phenotype = phenotype + phen_length
                    position = position + var_pos_list
                    #create index for use with DMq
                    self.phen_index[phen] = count    
            #create dataframe and write to file
            track_df = pd.DataFrame({'DM_track': track})
            df = pd.concat([df, track_df], axis=1)
            phen_df = pd.DataFrame({'DM_phenotype' : phenotype})  
            df = pd.concat([df, phen_df], axis=1)
            pos_df = pd.DataFrame({'DM_Residue' : position}) 
            df = pd.concat([df, pos_df], axis=1)
            #df.to_csv(self.plotting_file, sep='\t', na_rep="?")
            ##print(df)
            df.to_csv(self.plotting_file, sep='\t', na_rep=list(self.phen_index.keys())[list(self.phen_index.values()).index(1)])
            	#(list(mydict.keys())[list(mydict.values()).index(16)])
            self.HGMD_DM_track_count = count
            self.total_phen_count = count
        elif not DM:
            count = 0
            for phen, var_list in obj_d.items():
                var_pos_list = self.var_pos_extractor(var_list)
                if var_pos_list != []:
                    count += 1
                    #see if the phenotype has been plotted and use that index
                    try:
                        count_index = [self.phen_index[phen]] * len(var_pos_list)
                    #otherwise make index total number of phenotypes
                    except:
                        self.total_phen_count += 1
                        count_index = [self.total_phen_count] * len(var_pos_list)
                        self.phen_index[phen] = count_index[0]
                    track = track + count_index
                    phen_length = [phen] * len(var_pos_list)
                    phenotype = phenotype + phen_length
                    position = position + var_pos_list
            track_df = pd.DataFrame({'DMq_track': track})
            df = pd.concat([df, track_df], axis=1)
            phen_df = pd.DataFrame({'DMq_phenotype' : phenotype})  
            df = pd.concat([df, phen_df], axis=1)
            pos_df = pd.DataFrame({'DM_qResidue' : position}) 
            df = pd.concat([df, pos_df], axis=1)
            
            #sorted_keys = sorted(self.phen_index, key=self.phen_index.__getitem__)
            #phen_uniq_df = pd.DataFrame({'phenotype_labels' : sorted_keys})
            #df = pd.concat([df, phen_uniq_df], axis=1)
            df.to_csv(self.plotting_file, sep='\t', na_rep=list(self.phen_index.keys())[list(self.phen_index.values()).index(1)])
            self.HGMD_DMq_track_count = count
           

    def var_pos_extractor(self, var_list):
        ##print(var_list)
        pos_list = []
        if type(var_list) == list:
            for var in var_list:
                try:
                    m = re.search(r'(p\.)([^0-9]*)(\d{1,})(.*)', var.HGVS_prot)
                    if m:
                    	pos_list.append(m.group(3))
                except:
                    print( )
                    print(var.__dict__)
                    print(var.HGVS_prot)
        return pos_list 

    def var_class_separator(self, obj_list):
        class_groups = defaultdict(list)
        for obj in obj_list:
            class_groups[obj.var_class].append(obj)
        return class_groups

    def phenotype_lists(self, var_class, obj_list, DM=False):
        phen_groups = defaultdict(list)
        for obj in obj_list:
            phen_groups[obj.phenotype].append(obj)
        phen_groups['phen_count'] = len(phen_groups.keys())
        return phen_groups
             

### plot all data ###############################################################

    #constructs command line argument to feed variables into gnuplot script
    #gene name as title, length of x axis
    #returns output_variable=variable
    def construct_gnuplot_command(self, var_name, var):
        variable = "'" + var + "'"
        cmd = var_name + '=' + variable
        #print(cmd)
        return(cmd)
        
    def execute_zoomed_gnuplot(self, gene_name):
        length = int(self.slice_end) - int(self.slice_start)
        #print(length)
        svg_name_string = gene_name + '_zoomed_composite.svg'
        svg_name = self.construct_gnuplot_command("svg_name", svg_name_string)
        print(svg_name_string)
        x_length = self.construct_gnuplot_command("x_length" , str(length))
        gene_name = self.construct_gnuplot_command("gene_name", gene_name)
        data = self.construct_gnuplot_command("data", self.zoomed_plot)
        #print(data)
        if self.domain_count == 0:
            domain_count = self.coonstruct_gnuplot_command("domain_count", "1")
        else:
            domain_count = self.construct_gnuplot_command("domain_count", str(self.domain_count))
            domain_gnu = self.construct_gnuplot_command("domain_gnu", str(self.domain_count + 1))
        #self.HGMD_track_count = count
        #MAY NEED SOME LOGIC TO PLOT THE LARGEST AXIS FOR PHENOTYPES CURRENTLY PLOTTED BASED ON DM COUNT 
        DM_phen_count = self.construct_gnuplot_command("DM_phen_count", str(self.HGMD_DM_track_count))
        #print(DM_phen_count)
        DMq_phen_count = self.construct_gnuplot_command("DMq_phen_count", str(self.HGMD_DMq_track_count))
        #print(DMq_phen_count)
        total_phen_count = self.construct_gnuplot_command("total_phen_count", str(self.total_phen_count))
        user_pos = self.construct_gnuplot_command("user_pos", str(self.user_pos))
        gnuplot_command = ['gnuplot',
                           '-e', svg_name,
                           '-e', gene_name,
                           '-e', user_pos,
                           '-e', x_length,
                           '-e', data,
                           '-e', domain_count,
                           '-e', domain_gnu,
                           '-e', DM_phen_count,
                           '-e', DMq_phen_count,
                           '-e', total_phen_count,
                           "multiplot_zoomed"]
        #print(gnuplot_command)
        subprocess.call(gnuplot_command)

    #call gnuplot script
    def execute_gnuplot(self, gene_name):
        ##print(len(self.plottable_domains))
        adjusted_length = int(self.length)
        svg_name_string = gene_name + '_composite.svg'
        svg_name = self.construct_gnuplot_command("svg_name", svg_name_string)
        x_length = self.construct_gnuplot_command("x_length" , str(adjusted_length))
        gene_name = self.construct_gnuplot_command("gene_name", gene_name)
        data = self.construct_gnuplot_command("data", self.plotting_file)
        if self.domain_count == 0:
            domain_count = self.construct_gnuplot_command("domain_count", "1")
        else:
            domain_count = self.construct_gnuplot_command("domain_count", str(self.domain_count))
            domain_gnu = self.construct_gnuplot_command("domain_gnu", str(self.domain_count + 1))
        #self.HGMD_track_count = count
        #MAY NEED SOME LOGIC TO PLOT THE LARGEST AXIS FOR PHENOTYPES CURRENTLY PLOTTED BASED ON DM COUNT 
        DM_phen_count = self.construct_gnuplot_command("DM_phen_count", str(self.HGMD_DM_track_count))
        #print(DM_phen_count)
        DMq_phen_count = self.construct_gnuplot_command("DMq_phen_count", str(self.HGMD_DMq_track_count))
        #print(DMq_phen_count)
        total_phen_count = self.construct_gnuplot_command("total_phen_count", str(self.total_phen_count))
        user_pos = self.construct_gnuplot_command("user_pos", str(self.user_pos))
        #domains = construct_gnuplot_command("domains", "temp_domain_file")
        #domain_ints = construct_gnuplot_command("domain_ints", "temp_dom_ints_file")
        #exac = construct_gnuplot_command("exac", "temp_freq_file")
        #consurf = construct_gnuplot_command("consurf", "temp_cons_file")
        #hgmd = construct_gnuplot_command("hgmd", "temp_hgmd_file")
        #gnuplot_command = ['gnuplot', '-e', svg_name, '-e', gene_name, '-e', x_length, '-e', domains, '-e', domain_ints, '-e', exac, '-e', consurf, '-e', hgmd, "plot"]
        gnuplot_command = ['gnuplot',
                           '-e', svg_name,
                           '-e', gene_name,
                           '-e', user_pos,
                           '-e', x_length,
                           '-e', data,
                           '-e', domain_count,
                           '-e', domain_gnu,
                           '-e', DM_phen_count,
                           '-e', DMq_phen_count,
                           '-e', total_phen_count,
                           "multiplot_final"]
        #print(gnuplot_command)
        subprocess.call(gnuplot_command)
                
    def create_smaller_graph_file(self):
        #read in dataframe from file.data
        df = pd.read_csv(self.plotting_file, delimiter='\t', index_col=0)
        #extract columns into new dataframe; domains based on n columns from start, "cons", "cons_MA", "cons_cumsum"
        columns = self.uniprot_columns + ["cons", "cons_MA", "cons_cumsum", "pos"]
        df_mini = df.filter(columns, axis=1)
        #slice rows +/- 300 residues around user_pos
        self.slice_start = int(self.user_pos) - 200
        if self.slice_start < 0:
            self.slice_start = 0
        self.slice_end = int(self.user_pos) + 200
        if self.slice_end > int(self.length):
            self.slice_end = int(self.length)
        df_slice = df_mini[self.slice_start:self.slice_end]
        self.columns_in_range(df_slice, df, columns)
        
    	
    def columns_in_range(self, df_sliced, df_whole, columns):    
    	#drop unwanted_columns
    	df_other = df_whole.drop(columns, axis=1)
    	#return columns of interest
    	exac_het_df = self.df_filter_columns(df_other, ["het_pos", "het_freq"], "het_pos")
    	df_sliced = pd.concat([df_sliced, exac_het_df], axis=1) 	
    	exac_hom_df = self.df_filter_columns(df_other, ["homo_pos", "homo_freq"], "homo_pos")
    	df_sliced = pd.concat([df_sliced, exac_hom_df], axis=1)
    	DM_df = self.df_filter_columns(df_other, ["DM_track", "DM_phenotype", "DM_Residue"], "DM_Residue")
    	df_sliced = pd.concat([df_sliced, DM_df], axis=1)
    	DMq_df = self.df_filter_columns(df_other, ["DMq_track", "DMq_phenotype", "DM_qResidue"], "DM_qResidue")
    	df_sliced = pd.concat([df_sliced, DMq_df], axis=1)
    	df_sliced.to_csv(self.zoomed_plot, sep='\t', na_rep='?')
    	
    def df_filter_columns(self, df, include_list, column_name):
        filtered_df = df.filter(include_list, axis=1)
        this_range = []
        for i in range(self.slice_start, self.slice_end):
            this_range.append(str(i))
        sliced_filtered_df = filtered_df.loc[df[column_name].isin(this_range)]
        return sliced_filtered_df   
        return filtered_df
        
if __name__ == "__main__":
    Graph_object(sys.argv[1], sys.argv[2])
