#!/usr/bin/env python3

import os
import pandas as pd

from pathlib import Path
from PopColors import PopColor

class SubjectPops:
    """
    Read ethnicities or ancestry group IDs from file and set colors to
    population groups
    
    """

    def __init__(self, args):
        self.sbj_anc_file = args.input_file  # GrafAnc result file
        self.sbj_pop_file = args.input_file  # optional sbj race file
        self.ga_sbj_col = "Sample"           # sbj column in the GrafAnc result file
        self.ga_anc_col = "AncGroupID"       # population column in GrafAnc result file
        self.ga_snp_col = "SNPs"             # num SNPs column in GrafAnc result file
        self.show_ancs = []                  # Anc IDs to be included in the plot
        self.num_show_ancs = 0               # number of the above Anc IDs
        self.sbj_col = ""                    # sbj column in sbj pop file
        self.pop_col = ""                    # population column in sbj pop file
        self.pop_type = "SuperAnc"           # type can be Anc, SuperAnc, Pop

        if args.sub_anc:
            self.pop_type = "SubAnc"

        if self.pop_type != "Pop":
            self.sbj_col = self.ga_sbj_col
            self.pop_col = self.ga_anc_col

        # Subjects are colored by:
        # 1. Races/ethnicities if sbj race file is specified
        # 2. Super AncGroup IDs if pop_type = SuperAnc 
        # 3. AncGroup IDs if pop_type = SubAnc 
        if args.sbj_pop_file:
            self.sbj_pop_file = args.sbj_pop_file
            self.pop_type = "Pop"
            self.sbj_col = args.sbj_col
            self.pop_col = args.pop_col

        # Subjects and Anc IDs are read from the GrafAnc result file
        self.sbj_anc_ids = {}        # subject Anc ID dict
        self.sbj_sup_anc_ids = {}    # subject Super Anc ID dict
        self.all_anc_ids = {}        # all Anc IDs
        self.all_sup_anc_ids = {}    # all Super Anc IDs
        self.tot_sbjs = 0            # total num sbjs in GrafAnc file
        self.num_anc_sbjs = 0        # num sbjs with selected Anc IDs
        self.num_anc_ids = 0         # num selected Anc IDs
        self.num_sup_anc_ids = 0     # num selected Super Anc IDs
        self.sbj_mean_snps = 0       # mean genotyped SNPs per sbj
        self.sbj_min_snps = 0        # min genotyped SNPs per sbj
        self.sbj_max_snps = 0        # max genotyped SNPs per sbj
        
        # Subject races (or pops) are read from the user-specified sbj race file.
        # Only the subjects included in the GrafAnc results are read
        self.sbj_pops = {}           # sbj race dict
        self.all_pops = {}           # list of populations
        self.pop_cnts = {}           # population, #sbjs dict
        self.sorted_pops = []        # sorted (by count) population list
        self.selected_pops = []      # selected populations to be colored
        self.num_sbjs = 0            # number of subjects
        self.num_pops = 0            # number of populations
        self.anc_disp_names = {}     # Ancestry display names
        self.sup_anc_disp_names = {} # Super ancestry display names

        self.pop_colors = {}         # population color dict
        self.color_list = []         # list of colors for populations

        self.lgd_pops = []           # list of populations shown on legend
        self.lgd_colors = []         # list of colors for above populations

        pop_color = PopColor(args)
        if args.color_file:
            pop_color.SetColorListFromFile(args.color_file)

        self.color_list = pop_color.color_list
        self.background_color = pop_color.background_color        

        self.GetAncestryDisplayNames()
        # When '--show_ancs' is specified, read only the specified Anc IDs from file
        if args.show_ancs:
            anc_list = args.show_ancs.split (",")
            for anc in anc_list:
                anc_id = int(anc)
                if anc_id in self.sup_anc_disp_names or anc_id in self.anc_disp_names:
                    self.show_ancs.append(anc_id)
                else:
                    print(f'\nERROR: {anc} is an invalid AncGroupID.\n')
                    exit()

        self.num_show_ancs = len(self.show_ancs)
        if self.num_show_ancs > 0:
            print(f'The following Anc IDs will be displayed: {self.show_ancs}')
        
        # Read subject Anc IDs from GrafAnc result, keep only some Anc IDs if specified 
        self.ReadSubjectAncs(args)

        # Read subject races from race file, ignoring sbjs not in the above list
        if self.pop_type == "Pop":
            self.ReadSubjectPopsFromRaceFile(args)
        else:
            # No self-reported races. Pops are Anc IDs from the GrafAnc file
            self.GetSubjectPopsFromGrafAncFile()

        self.SetPopulationColors(args)            
        self.ShowSummary()

    def ReadSubjectAncs(self, args):
        """
        Read subjects and AncGroupIDs from GrafAnc result file
        """
        if not os.path.isfile(self.sbj_anc_file):
            print(f"ERROR: didn't find file {self.sbj_anc_file}\n")
            exit()

        print(f'\nReading GrafAnc results from {self.sbj_anc_file}')

        pf = None
        try:
            pf = pd.read_csv(self.sbj_anc_file, sep="\t")
        except Exception as error:
            print(f'\nERROR: Cannot open {self.sbj_anc_file}: {error}\n')
            exit()
            
        self.tot_sbjs = len(pf[self.ga_sbj_col])
        if self.tot_sbjs < 1:
            print(f'\nERROR: No results found in {self.sbj_anc_file}')
            exit()
            
        print(f'\tTotal {self.tot_sbjs} subjects in column {self.ga_sbj_col}')

        tot_snps = 0
        min_snps, max_snps = 500000, 0
        
        for i in range(self.tot_sbjs):
            sbj = pf[self.ga_sbj_col][i]
            anc = pf[self.ga_anc_col][i]
            num_snps = pf[self.ga_snp_col][i]
            sbj = str(sbj)
            anc_id = int(anc)
            sup_anc_id = int(anc_id/100) * 100
            tot_snps += num_snps
            if num_snps < min_snps:
                min_snps = num_snps
            if num_snps > max_snps:
                max_snps = num_snps
                
            # If show_ancs is specified, only read subjects with these Ancs 
            flag = True
            if self.num_show_ancs > 0:
                if anc_id not in self.show_ancs and sup_anc_id not in self.show_ancs:
                    flag = False

            if flag:
                self.sbj_anc_ids[sbj] = anc_id
                self.sbj_sup_anc_ids[sbj] = sup_anc_id

                if self.all_anc_ids.get(anc_id) is None:
                    self.all_anc_ids[anc_id] = 1
                else:
                    self.all_anc_ids[anc_id] += 1

                if self.all_sup_anc_ids.get(sup_anc_id) is None:
                    self.all_sup_anc_ids[sup_anc_id] = 1
                else:
                    self.all_sup_anc_ids[sup_anc_id] += 1

        mean_snps = tot_snps * 1.0 / self.tot_sbjs
        print(f'\tSNPs/subject: {min_snps} - {max_snps}, mean = {mean_snps:.1f}')

        self.sbj_min_snps = min_snps
        self.sbj_max_snps = max_snps
        self.sbj_mean_snps = mean_snps
        
        self.num_anc_sbjs = len(self.sbj_anc_ids)
        self.num_anc_ids = len(self.all_anc_ids)
        self.num_sup_anc_ids = len(self.all_sup_anc_ids)

        sel_ancs = []
        for anc in self.all_anc_ids.keys():
            cnt = self.all_anc_ids[anc]
            sel_ancs.append(anc)
        sel_ancs.sort()

        sel_sup_ancs = []
        for anc in self.all_sup_anc_ids.keys():
            cnt = self.all_sup_anc_ids[anc]
            sel_sup_ancs.append(anc)
        sel_sup_ancs.sort()

        if self.num_show_ancs > 0:
            print(f'\tSelected subjects: {self.num_anc_sbjs}')
            print(f'\tSelected AncIDs: {self.num_anc_ids}: {sel_ancs}')
            print(f'\tSelected Super AncIDs: {self.num_sup_anc_ids}: {sel_sup_ancs}')

        print("")

    def ReadSubjectPopsFromRaceFile(self, args):
        """
        Read subjects and self-reported races/ethnicities from subject race file.
        Only those included in those from ReadSubjectAncs will be read.
        """
        if not os.path.isfile(self.sbj_pop_file):
            print(f"ERROR: didn't find file {self.sbj_pop_file}\n")
            exit()

        num_file_sbjs = 0
        pf = None
        has_header = True
        sbj_col, pop_col = self.sbj_col, self.pop_col
        
        if self.sbj_col and self.pop_col:            
            try:
                pf = pd.read_csv(self.sbj_pop_file, sep="\t")
            except Exception as error:
                print(f'\nERROR: Cannot open {self.sbj_pop_file}: {error}\n')
                exit()

            cols = pf.columns
            if self.sbj_col not in cols:
                print(f'\nERROR: did not find column {self.sbj_col} in {self.sbj_pop_file}\n')
                exit()
            if self.pop_col not in cols:
                print(f'\nERROR: did not find column {self.pop_col} in {self.sbj_pop_file}\n')
                exit()
                
            num_file_sbjs = len(pf[self.sbj_col])
            print(f'Reading sbj races/pops from {self.sbj_pop_file} columns {self.sbj_col} and {self.pop_col}')  
        else:
            # If sbj and race columns are not specified, assume first two columns are sbj and race cols
            sbj_col, pop_col = "sbj", "pop"
            try:
                pf = pd.read_csv(self.sbj_pop_file, sep="\t", usecols=[0,1], names=[sbj_col, pop_col])
            except Exception as error:
                print(f'\nERROR: Cannot open {self.sbj_pop_file}: {error}\n')
                exit()

            num_file_sbjs = len(pf)
            has_header = False
            print(f'Read races/pops of {num_file_sbjs} subjects from first two columns of {self.sbj_pop_file}')  

        print(f'\tTotal {num_file_sbjs} subjects in file {self.sbj_pop_file}')

        for i in range(num_file_sbjs):
            sbj, pop = "", "NO VALUE"
            sbj = pf[sbj_col][i]
            pop = pf[pop_col][i]
            sbj = str(sbj)

            anc = self.sbj_anc_ids.get(sbj)

            if self.sbj_anc_ids.get(sbj) is not None:
                self.sbj_pops[sbj] = pop
                if self.all_pops.get(pop) is None:
                    self.all_pops[pop] = 1
                else:
                    self.all_pops[pop] += 1

        self.num_sbjs = len(self.sbj_pops)    
        self.num_pops = len(self.all_pops)    

        print(f'\tRead {self.num_sbjs} subjects {self.num_pops} races from {self.sbj_pop_file}\n')
        if self.num_sbjs < 1:
            err = "\tNo race/ethnicity info found in race/ethnicity file for subjects in GrafAnc results.\n"
            err += "\tPlease make sure the file is tab-delimited and correct column headers are specified."
            
            print(f'\nERROR: {err}\n')
            exit()

    def GetSubjectPopsFromGrafAncFile(self):
        """
        By default, subjects are colored by Anc IDs in the GrafAnc result file.
        So, sbj_pops are the same as sbj_anc_ids.
        """
        for sbj in self.sbj_anc_ids:
            anc = self.sbj_anc_ids[sbj]
            sup_anc = self.sbj_sup_anc_ids[sbj]
            pop = sup_anc
            if self.pop_type == "SubAnc":
                pop = anc
                
            self.sbj_pops[sbj] = pop
            if self.all_pops.get(pop) is None:
                self.all_pops[pop] = 1
            else:
                self.all_pops[pop] += 1

        self.num_sbjs = len(self.sbj_pops)    
        self.num_pops = len(self.all_pops)    
        sorted_pop_cnts = sorted(self.all_pops.items(), key=lambda item: item[1], reverse=True)

        print(f'Read {self.num_sbjs} subjects {self.num_pops} races/pops from {self.sbj_pop_file}\n')
                    
    def SetPopulationColors(self, args):
        """
        By default, populations are sorted by subject counts, and the top ones
        with colors are selected to plot on front. Others are plotted to 
        the background with background color
        """
        num_colors = len(self.color_list)

        sorted_pop_cnts = sorted(self.all_pops.items(), key=lambda item: item[1], reverse=True)
        
        # Assign colors in the color list to pops in the sorted list.
        # Set all other pops to background color if pop list is longer.
        for i in range(self.num_pops):
            pop_cnt = sorted_pop_cnts[i]
            pop = pop_cnt[0]
            cnt = pop_cnt[1]
            
            self.sorted_pops.append(pop)
            self.pop_cnts[pop] = cnt
            color = self.background_color
            
            if i < num_colors:
                color = self.color_list[i]
                self.selected_pops.append(pop)
                disp_pop = pop
                if self.pop_type != "Pop":
                    disp_pop = self.anc_disp_names[pop]

                self.lgd_pops.append(disp_pop)
                self.lgd_colors.append(color)
                
            self.pop_colors[pop] = color

        other_pop = "Other ancestry groups"
        if self.pop_type == "Pop":
            other_pop = "Other populations"

        # When --pops option is set, the selected pops are specified by user
        # through this comma-delimited integer list, with each integer be
        # the pop num in the sorted pop list
        if args.pops:
            sel_pop_strs = args.pops.split(",")
            num_sel_pops = len(sel_pop_strs)
            sel_pop_ids = []
            for i in range(num_sel_pops):
                pop_str = sel_pop_strs[i]
                if pop_str != "":
                    pop_no = int(pop_str)
                if pop_no > 0:
                    sel_pop_ids.append(pop_no)

            self.selected_pops = []
            self.lgd_pops = []
            self.lgd_colors = []

            sel_pop_no = 0
            for i in range(self.num_pops):
                pop_cnt = sorted_pop_cnts[i]
                pop = pop_cnt[0]
                pop_no = i + 1

                color = self.background_color

                if pop_no in sel_pop_ids:
                    if sel_pop_no < num_colors:
                        color = self.color_list[sel_pop_no]
                        self.selected_pops.append(pop)
                        disp_pop = pop
                        if self.pop_type != "Pop":
                            disp_pop = self.anc_disp_names[pop]

                        self.lgd_pops.append(disp_pop)
                        self.lgd_colors.append(color)
                        sel_pop_no += 1
                    
                self.pop_colors[pop] = color

        self.lgd_pops.append(other_pop)
        self.lgd_colors.append(self.background_color)

    def GetAncestryDisplayNames(self):
        self.anc_disp_names[100] = "100: AFR"
        self.anc_disp_names[200] = "200: MENA"
        self.anc_disp_names[300] = "300: EUR"
        self.anc_disp_names[400] = "400: SAS"
        self.anc_disp_names[500] = "500: EAS"
        self.anc_disp_names[600] = "600: AMR"
        self.anc_disp_names[700] = "700: OCN"
        self.anc_disp_names[800] = "800: MIX"

        self.anc_disp_names[101] = "101: Nigeria"
        self.anc_disp_names[102] = "102: Western Africa"
        self.anc_disp_names[103] = "103: Central Africa"
        self.anc_disp_names[104] = "104: Kenya"
        self.anc_disp_names[105] = "105: Southern Africa"
        self.anc_disp_names[106] = "106: Northeastern Africa"
        self.anc_disp_names[107] = "107: African American"
        self.anc_disp_names[108] = "108: Other Africa"
        self.anc_disp_names[201] = "201: Northern Africa"
        self.anc_disp_names[202] = "202: Middle East 1"
        self.anc_disp_names[203] = "203: Middle East 2"
        self.anc_disp_names[301] = "301: Finland"
        self.anc_disp_names[302] = "302: Northern Europe"
        self.anc_disp_names[303] = "303: Western Europe"
        self.anc_disp_names[304] = "304: Southern Europe"
        self.anc_disp_names[305] = "305: Northeastern Europe"
        self.anc_disp_names[306] = "306: Southeastern Europe"
        self.anc_disp_names[307] = "307: Balkans"
        self.anc_disp_names[308] = "308: Other Europe"
        self.anc_disp_names[401] = "401: Asian India"
        self.anc_disp_names[402] = "402: Gujarati in India"
        self.anc_disp_names[403] = "403: Northern South Asia"
        self.anc_disp_names[404] = "404: Sri Lanka"
        self.anc_disp_names[405] = "405: Bangladesh"
        self.anc_disp_names[501] = "501: Japan Ryukyu"
        self.anc_disp_names[502] = "502: Japan Main Islands"
        self.anc_disp_names[503] = "503: Korea"
        self.anc_disp_names[504] = "504: Northern Asia"
        self.anc_disp_names[505] = "505: Northern China 1"
        self.anc_disp_names[506] = "506: Northern China 2"
        self.anc_disp_names[507] = "507: Southern China 1"
        self.anc_disp_names[508] = "508: Southern China 2"
        self.anc_disp_names[509] = "509: Southeast Asia"
        self.anc_disp_names[510] = "510: Thailand"
        self.anc_disp_names[511] = "511: Other East Asia"
        self.anc_disp_names[601] = "601: Latin American 1"
        self.anc_disp_names[602] = "602: Latin American 2"
        self.anc_disp_names[603] = "603: Native American"

    def ShowSummary(self):
        num_colors = len(self.color_list)
        num_pops = len(self.sorted_pops)
        print(f'Total {self.num_sbjs} subjects from {self.num_pops} populations')

        for i in range(self.num_pops):
            pop_no = i + 1
            pop = self.sorted_pops[i]
            cnt = self.pop_cnts[pop]
            color = self.pop_colors[pop]
            if pop == "":
                pop = "NO VALUE"
            print(f'\tPop No {pop_no}: {pop} (n={cnt} c={color})')
