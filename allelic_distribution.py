#!/usr/bin/env  python3

"""
This script can be used to assess the distribution of cg/wgMLST alleles across the groups of any metadata variable of interest

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import datetime as datetime
import pandas

version = "1.0.0"
last_updated = "2025-09-22"

def concat_df(df_alleles, df_groups, group_col):
	"""
	concatenate allele and metadata matrices
	"""
	
	df_groups = df_groups.set_index(df_groups.columns[0], drop = True)
	df_alleles = df_alleles.set_index(df_alleles.columns[0], drop = False)

	df = pandas.concat([df_groups[group_col], df_alleles], axis=1).dropna()
	df = df.reset_index(drop = True)

	return df

def set_variables(df, group_col):
    """
    initiate dictionaries for the output matrices 
    """

    groups = df[group_col].values.tolist()
    summary = {"locus": [], "allele": []}
    dist = {"locus": []}
    order_columns = ["locus", "allele"]
    order_columns2 = ["locus"]
    group_size = {}
    checked = {}
    for group in set(groups):
        summary[group] = []
        dist[group] = []
        group_size[group] = groups.count(group)
        if group not in order_columns:
            order_columns.append(group)
            order_columns2.append(group)
    
    return groups, summary, dist, order_columns, order_columns2, group_size, checked

def group_analysis(df, df_interest, summary, group_size, groups, checked):
	"""
    generate data for the summary analysis
    """

	for locus in df.columns[2:]:
		if locus not in checked.keys():
			checked[locus] = []
		alleles_group_interest = set(df_interest[locus].values.tolist())
		for allele in alleles_group_interest:
			if allele not in checked[locus]:
				checked[locus].append(allele)
				summary["locus"].append(locus)
				summary["allele"].append(allele)
				df_allele = df[df[locus] == allele]
				observed_group = df_allele[df_allele.columns[0]].values.tolist()
				for group in set(groups):
					if group not in observed_group:
						info = "0"
					else:
						count = observed_group.count(group)
						pct = count/group_size[group]
						val = str(pct)
						info = str(val)
					summary[group].append(info)
	
	return summary

def group_distribution(df_work, group_col, dist, groups):
	"""
    generate data for the distribution analysis
    """

	for locus in df_work.columns[2:]:
		dist["locus"].append(locus)
		for group in set(groups):
			counter = {}
			df_group = df_work[df_work[group_col] == group]
			alleles_group_interest = df_group[locus].values.tolist()
			for allele in set(alleles_group_interest):
				counter[allele] = alleles_group_interest.count(allele)
			info2report = []
			for v in sorted(counter, key=counter.get, reverse=True):
				rel_freq = float(counter[v]/len(alleles_group_interest))
				statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)"
				info2report.append(statistics)
			joint = ", ".join(info2report) + " (n = " + str(len(alleles_group_interest)) + ")"
			dist[group].append(joint)
	
	return dist

def main():
	parser = argparse.ArgumentParser(prog="allelic_distribution.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                           allelic_distribution.py                           #
									#                                                                             #
									###############################################################################
									                            
									This script assesses the distribution of each allele present in all the groups 
									of a given variable and also provides a summary table with information about 
									the proportion of isolates of each group that harbor a given allele.
									
									-------------------------------------------------------------------------------"""))
										
	parser.add_argument("-m", "--metadata", dest="metadata", required=True, type=str, help="[MANDATORY] Metadata file in .tsv format.")
	parser.add_argument("-a", "--allele", dest="alleles", required=True, type=str, help="[MANDATORY] Allele matrix in .tsv format.")
	parser.add_argument("-c", "--group-column", dest="group_column", required=True, type=str, help="[MANDATORY] Name of the metadata column with the groups of interest.")
	parser.add_argument("-g", "--group-interest", dest="group_interest", required=False, default="", type=str, help="[OPTIONAL] Comma-separated list of groups of interest to select the alleles under analysis. \
                        If nothing is indicated, all goups will be considered as of interest and all alleles will be analysed.")
	parser.add_argument("-o", "--output", dest="output", required=False, default="Allele_distribution", type=str, help="[OPTIONAL] Tag for output files. default: Allele_distribution")
	
	args = parser.parse_args()
    
    # starting logs	----------
    
	log_name = args.output + ".log"
	log = open(log_name, "w+")

	print("\n******************** running allelic_distribution.py ********************\n", log)
	print("version " + str(version) + " last updated on " + str(last_updated) + "\n")
	print(" ".join(sys.argv))
	
	start = datetime.datetime.now()
	print("start: " + str(start))
    
    #inputs
	df_alleles = pandas.read_table(args.alleles, dtype=str)
	df_groups = pandas.read_table(args.metadata, dtype=str)
	group_col = args.group_column
	
	df_work = concat_df(df_alleles, df_groups, group_col)
	groups, summary, dist, order_columns, order_columns2, group_size, checked = set_variables(df_work, group_col)
    
	if args.group_interest == "":
		group_interest = groups
	else:
		group_interest = args.group_interest.split(",")
        
	for group in group_interest:
		df_interest = df_work[df_work[group_col] == group]
		df_interest = df_interest.astype(str)
		summary = group_analysis(df_work, df_interest, summary, group_size, groups, checked)
	
	distribution = group_distribution(df_work, group_col, dist, groups)
    
	summary_df = pandas.DataFrame(data = summary, columns = order_columns)
	distribution_df = pandas.DataFrame(data = distribution, columns = order_columns2)
	summary_df.to_csv(args.output + "_" + group_col + "_summary.tsv", index = False, header=True, sep ="\t")
	distribution_df.to_csv(args.output + "_" + group_col + "_distribution.tsv", index = False, header=True, sep ="\t")

	end = datetime.datetime.now()
	
	print("end: " + str(end))
	print("\nDone!")
	log.close()
	print("************************************************************")

if __name__ == "__main__":
    main()