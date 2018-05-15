from api import *
import re, sys, subprocess, json, os, getpass
import numpy as np
import pandas as pd
from collections import defaultdict

#the purpose of this class is to model and collect the data required of a plot and then write to a file for retention and subsequent plotting

hgmd_username = input("\nEnter HGMD Pro licence username: ")
hgmd_password = getpass.getpass(prompt="Enter HGMD Pro licence password (hidden): ")

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
        self.hgmd_username = hgmd_username
        self.hgmd_password = hgmd_password

   	#Ensembl#############################
        self.Ens = Ensembl_api()
        print('\nGetting Ensembl gene ID for ' + gene_name + '...')
        self.ensembl_id = self.get_ensembl_id(gene_name)
        #Uniprot##########################
        self.Up = Uniprot_api()
        #composite file_name to add data to as it is parsed from different sources
        self.plotting_file = gene_name + '_composite_' + user_pos + '.data'
        #self.plotting_file = gene_name + '_composite.data'
        self.zoomed_plot = gene_name + 'zoomed_composite.data'
        #dict of uniprot entries for this ensembl transcript
        print('Getting Uniprot annotation file...') 
        self.all_uniprot_entries = self.get_uniprot_entries(self.ensembl_id)
        #returns human reviews uniprot annotations
        self.reviewed_uniprot_entries = self.reviewed_human_entries(self.all_uniprot_entries) 
        ## N.B. Failure of the script at this point may be due to lack of connection to Uniprot
        ##Traceback (most recent call last):  File "PM1_plotter.py", line 41, in __init__ self.length = self.reviewed_uniprot_entries[0]['Length'] IndexError: list index out of range
        ## Try running script again
        ##print(self.reviewed_uniprot_entries)
        self.length = self.reviewed_uniprot_entries[0]['Length']
        #returns text of all gff annotations
        self.all_gff_annotation = self.get_gff_domains()
        #Gff objects where values are ranges of 'Domain', 'Region', 'DNA binding'
        self.required_gff_annotations = self.specific_gff_annotations()
        #generates dict of arrays to be written to master gnuplot file
        self.plottable_domains = self.generate_plottable_domains(self.length)
        
        #Exac##############################
        print("Gathering ExAC data for " + self.ensembl_id + " ...")
        self.Ex = Exac_api()
        self.all_exac_variants,self.chrom = self.get_exac_data(self.ensembl_id)
        
        #Consurf###########################
        print('\nGathering Consurf data stored on the server...')
        try:
            self.consurf_file = self.find_consurf_file(gene_name)
        except:
            self.consurf_file = "no_file"
        self.consurf_data = self.parse_consurf_grades(self.consurf_file, self.length)
        self.write_consurf_data = self.write_consurf_grades(self.consurf_data)

        #HGMD##############################
        print('\nGathering HGMD data from website...')
        self.hgmd_data = self.get_HGMD_data(gene_name)             
        if self.DM_objs == {} and self.DM_likely_objs == {}:
            print("No variants found in HGMD")
        elif self.DM_objs == {}:
            self.write_DM_data2 = self.write_HGMD_data(self.DM_likely_objs, DM=False)
        else:
            self.write_DM_data = self.write_HGMD_data(self.DM_objs)
            self.write_DM_data2 = self.write_HGMD_data(self.DM_likely_objs, DM=False)

        #Gnomad_data######################
        
        #GNU plotter#######################
        print('\nPlotting all data...\n')
        if self.chrom=='X':
            self.execute_gnuplot(gene_name, user_pos, self.chrom, hemi=True)
        elif self.chrom=='Y':
            self.execute_gnuplot(gene_name, user_pos, self.chrom, chrY=True)
        else:
            self.execute_gnuplot(gene_name, user_pos, self.chrom)
        print("Data plotted.\n")
        #self.create_smaller_graph_file()
        #self.execute_zoomed_gnuplot(gene_name)
        
    ### Ensembl_id   ####################################################
    #required for ExAC query
    def get_ensembl_id(self, HGNC_gene):
        id = self.Ens.query_HGNC(HGNC_gene)
        return id

    ### Uniprot_data #####################################################
    def get_uniprot_entries(self, ensembl_id):
        Up = Uniprot_api()
        up_entry_list = Up.get_pid_from_gene(ensembl_id)
        return up_entry_list
    
    #returns reviewed human entries 
    def reviewed_human_entries(self, entry_list):
        human_match_list = []
        for entry_dict in entry_list:
            try:
                if entry_dict['Organism'] == 'Homo sapiens (Human)' and entry_dict['Status'] == 'reviewed':
                    human_match_list.append(entry_dict)
                    #uniprot_revd_entry = entry_dict['Entry']
                    #print(uniprot_revd_entry)
                    #return uniprot_revd_entry
                    #print(human_match_list)
            except:
                print('Check reviewed transcript is identified:')
                print(human_match_list)
                print('Check gene name at www.uniprot.org and/or try re-running script, as Uniprot connection possibly failed')
        return human_match_list
 
    #returns the gff annotation data for gene of interest
    def get_gff_domains(self):
        if len(self.reviewed_uniprot_entries) == 1:
            gff_response = self.Up.get_entry_gff(self.reviewed_uniprot_entries[0]["Entry"])
            #print(gff_response)
            print('One reviewed Uniprot entry found, connection to Uniprot successful')
            return gff_response.text
        elif len(human_match_list) == 0:
            print('no_up_id_matches -- see PM1_plotter.py gff_domains')
        elif len(human_match_list) > 1:
            print('More than 1 reviewed human Uniprot transcript entry PM1_plotter.py gff_domains')

    #filters gff data EDIT the list here to change which annotations are included
    # N.B. adding too many areas of interest may cause the plot to become too crowded
    def specific_gff_annotations(self):
        # Possible options to add the required_list include:
        #required_list = ["Beta strand", "Helix", "Turn", "Repeat","Motif", "Domain", "Region", "Transmembrane",  "DNA binding", "Zinc finger", "Disulfide bond", "Nucleotide binding", "Coiled coil"]
        # Adding "Topological domain" causes the script to fail
        try:
            required_list = ["Region","Repeat", "Domain", "Transmembrane", "DNA binding", "Zinc finger", "Disulfide bond", "Nucleotide binding", "Coiled coil"]
            gff_objects = self.Up.parse_gff(self.all_gff_annotation, required_list)
            if len(gff_objects) <1:
                sys.exit()
            else:
                return gff_objects
        except:
            print("\nNo specified domains in required list found in GFF file. \n\nCheck required_list and Uniprot gff file. Some areas of interest may need to be removed\n")
            sys.exit()

    def generate_plottable_domains(self,length):
        #length = self.reviewed_uniprot_entries[0]['Length']
        print("Uniprot canonical transcript protein length: " + str(length) + " amino acids\n")
        master_dict = {}
        domain_array = {}
        annotation_index = 0
        arrays_to_save = []
        #generate a dictionary from the required gff objects 
        for item in self.required_gff_annotations:
            if item.__dict__['anno_type'] in master_dict.keys():
                master_dict[item.__dict__['anno_type']].append(item.__dict__)
            else:
                master_dict[item.__dict__['anno_type']] = [item.__dict__]
        for k in master_dict.keys(): # k represents each of the plottable domains found e.g. transmembrane, helix etc. from required_list
            self.uniprot_columns.append(k)
            annotation_index += 1
            #function to return the array on NA or integer
            this_type_column = self.array_creator(master_dict, k, annotation_index, length)
            if len(this_type_column) == int(length)+1:
                this_type_column = this_type_column[:-1]
            elif len(this_type_column) > int(length)+1:
                print(k + " feature exceeds length of protein by more than one amino acid or multiple " + k + " features overlap. Program will now close")
                sys.exit()
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
    
    #generates arrays as specified above, 
    #note the additional logic to deal with the first and last differently 
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
    #gets ExAC variants in the ensembl gene
    #filters to retain ExAC pass filter, missense, initiator_codon_variant, stop_gained, frameshift_variant, inframe_deletion, stop_lost, canonical
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
        # identify canonical transcript
        canon_trans_info = self.Ex.canonical_transcript(ensembl_id)        
        self.exac_canon_transcript_id = canon_trans_info['gene']['canonical_transcript']
        print("ExAC canonical transcript ID: " + self.exac_canon_transcript_id)
        self.chrom = canon_trans_info['gene']['chrom']
        updated_missense_only = self.Ex.update_variant(missense_only)        

        # list of heterozygous entries
        heterozygotes = self.Ex.filter_variants(updated_missense_only, "het_count", 0, remove=True)
        # dict of het frequency and position
        het_freq_pos = self.Ex.position_frequency(heterozygotes)
        # save het data to cmposite data file
        exac_to_composite = self.add_exac_to_composite(het_freq_pos, indexed=False)

        # list of homozygous entries (if they exist - not available on Y chromosome records)
        try:
            homozygotes = self.Ex.filter_variants(updated_missense_only, "hom_count", 0, remove=True)
            # dict of hom frequency and position
            hom_freq_pos = self.Ex.position_frequency(homozygotes, hom=True)
            #save homozygous data to composite data file
            exac_to_composite = self.add_exac_to_composite(hom_freq_pos)
        except:
            pass

        # list of hemizygotes entries (if they exist)
        try:
            hemizygotes = self.Ex.filter_variants(updated_missense_only, "hemi_count", 0, remove=True)
            # dict of hemi frequency and position
            hemi_freq_pos = self.Ex.position_frequency(hemizygotes, hemi=True)
            # save hemi data to composite data file
            exac_to_composite = self.add_exac_to_composite(hemi_freq_pos)
        except: 
            pass
        return (all_variants,self.chrom)
    
    #indexed to indicate whether there is an index for pandas in column 0
    #uses pandas and enables column lengths to differ
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
        base = "consurf_scores/"
        for dirs, subdirs, files in os.walk(base):
            if dirs.endswith(gene_name):
                return os.path.join(dirs, files[0])

    #excludes consurf scores marked with * to indicate low confidence
    def parse_consurf_grades(self, consurf_file, length):
        cons_pos = { "cons" : [],
    	             "pos" : []}
    	# consurf.grades file for gene does not exist in consurf_outputs folder, create mock file
        if consurf_file == "no_file":
            cons_pos["pos"] = list(range(int(length)+1))
            del(cons_pos["pos"][0])
            cons_pos["cons"] = [0]
            cons_pos["cons"] = cons_pos["cons"]*int(length)
        elif consurf_file:
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

    #opens data, creates a new df with conservation, running mean and cumulative mean.
    #Then adds position to be used with gnuplot
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
        all_mutations_soup = HGMD.scrape_HGMD_all_mutations(hgmd_username,hgmd_password)
        #list of objects
        variant_instances = HGMD.extract_missense_nonsense(all_mutations_soup)
        #for i in variant_instances:
           #print(i.__dict__) #returns dict of mutation class and list of objects in that class
        separate_var_class = self.var_class_separator(variant_instances)
        for var_class,objs in separate_var_class.items():
            #separate_each_list_by returning dictionaries
            if var_class == 'DM':
                self.DM_objs = self.phenotype_lists(var_class, objs)
            elif var_class == 'DM?':	
                self.DM_likely_objs = self.phenotype_lists(var_class, objs)
        return variant_instances
        #write_hgmd = HGMD.write_DM_file(variant_instances)

    #Writes HGMD data into columns by phentype, separates by DM and DM? annotation
    #assigns an index number to each phentotype so they can be plotted separately without many columns in the datafile
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
            df.to_csv(self.plotting_file, sep='\t', na_rep="?")
            #print(df)
            #print(list(self.phen_index.keys())) # print list of phenotypes associated with a gene
            #print(list(self.phen_index.keys())[list(self.phen_index.values()).index(1)]) # print the first phenotype in the phenotype list
            #print([list(self.phen_index.values())])
            df.to_csv(self.plotting_file, sep='\t', na_rep=list(self.phen_index.keys())[list(self.phen_index.values()).index(1)])
            #print(df)
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

    #extracts the amino acid position of each variant 
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

    #retain the dictionary order
    def var_class_separator(self, obj_list):
        class_groups = defaultdict(list)
        for obj in obj_list:
            class_groups[obj.var_class].append(obj)
        return class_groups

    #returns a list of phenotypes based on phenotype
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
    def execute_gnuplot(self, gene_name, user_pos, chrom, hemi=False, chrY=False):
        adjusted_length = int(self.length)
        svg_name_string = gene_name + '_composite_' + user_pos + '.svg'
        svg_name = self.construct_gnuplot_command("svg_name", svg_name_string)
        x_length = self.construct_gnuplot_command("x_length" , str(adjusted_length))
        gene_name = self.construct_gnuplot_command("gene_name", gene_name)
        data = self.construct_gnuplot_command("data", self.plotting_file)
        chrom = self.construct_gnuplot_command("chrom", chrom)
        if self.domain_count == 0:
            print("No Uniprot domains to plot")
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

        if hemi==True or chrY==True:
            gnuplot_command = ['gnuplot',
                               '-e', svg_name,
                               '-e', gene_name,
                               '-e', user_pos,
                               '-e', x_length,
                               '-e', data,
                               '-e', chrom,
                               '-e', domain_count,
                               '-e', domain_gnu,
                               '-e', DM_phen_count,
                               '-e', DMq_phen_count,
                               '-e', total_phen_count,
                               "multiplot_final_hemi"]
        else:
            gnuplot_command = ['gnuplot',
                               '-e', svg_name,
                               '-e', gene_name,
                               '-e', user_pos,
                               '-e', x_length,
                               '-e', data,
                               '-e', chrom,
                               '-e', domain_count,
                               '-e', domain_gnu,
                               '-e', DM_phen_count,
                               '-e', DMq_phen_count,
                               '-e', total_phen_count,
                               "multiplot_final"]
        #print(gnuplot_command)
        subprocess.call(gnuplot_command)

    #inital method to create a zoomed in data file for plotting with an additional gnuplot script
    def create_smaller_graph_file(self):
        print("create_smaller_graph function in use")
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
        print("columns_in_range function in use")
    	#drop unwanted_columns
        df_other = df_whole.drop(columns, axis=1)
    	#return columns of interest
        exac_het_df = self.df_filter_columns(df_other, ["het_pos", "het_freq"], "het_pos")
        df_sliced = pd.concat([df_sliced, exac_het_df], axis=1) 	
        exac_hom_df = self.df_filter_columns(df_other, ["hom_pos", "hom_freq"], "hom_pos")
        df_sliced = pd.concat([df_sliced, exac_hom_df], axis=1)
        #exac_hemi_df = self.df_filter_columns(df_other, ["hemi_pos", "hemi_freq"], "hemi_pos")
  	#df_sliced = pd.concat([df_sliced, exac_hemi_df], axis=1)       
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
