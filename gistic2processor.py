#!/usr/bin/python

import sys, os, json, argparse, ConfigParser, copy, tarfile
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_table_line(table_line):
    '''
    Get widepeak and maxpeak gene list from table lines.
    Return both in the form of a set, as well as the chromosome number.
    '''

    line = table_line.rstrip().split('\t')
    widepeak_genes = line[9][:-1]
    if widepeak_genes[0] == '[':
        widepeak_genes = widepeak_genes[1:]
    widepeak_genes = widepeak_genes.split(',')


    maxpeak_genes = line[11][:-1]
    if maxpeak_genes[0] == '[':
        maxpeak_genes = maxpeak_genes[1:]
    maxpeak_genes = maxpeak_genes.split(',')

    return widepeak_genes, maxpeak_genes, line[1]

def get_peak_targets_wrapper(target_genes, config):
    '''
    Wrapper allows underlying function to be used by other modules without
    the need for a configuration file.
    '''
    amp_table = os.path.join(config.get('GisticData','data'), config.get('GisticData','amp_file'))
    del_table = os.path.join(config.get('GisticData','data'), config.get('GisticData','del_file'))
    return get_peak_targets(target_genes, amp_table, del_table)

def get_peak_targets(target_genes, amp_table, del_table):
    '''
    Using gene targets and gistic2 amp/del output files, identify
    and return target genes found in peaks. If a peak occurs and
    there is only one gene, it will be returned regardless if it
    is a target gene or not.
    '''

    # Dictionary of target genes that are actually found in peaks, in addition
    # to single genes found in max peaks.
    found_targets = {'Amp':set(), 'Del':set(), 'amp_peaks':0,'del_peaks':0}

    files = {'Amp':amp_table, 'Del':del_table}

    # Creates a dict->dict->set that will automatically generate entries
    # on assignment. The order is CNA type -> chromosome -> gene
    genes_in_peaks = defaultdict(lambda: defaultdict(set))

    genes_to_peak = defaultdict(set)
    peak_names = defaultdict()


    # extract genes found in wide and max peaks
    for cna_type in files:
        with open(files[cna_type]) as table_file:
            for line in table_file:
                # skip column name row
                if line[0] == 'i':
                    continue

                if [cna_type] == 'Amp':
                    found_targets['amp_peaks'] += 1
                else:
                    found_targets['del_peaks'] += 1

                widepeak_genes, maxpeak_genes, chromosome = parse_table_line(line)
                widepeak_set = set(widepeak_genes)
                maxpeak_set = set(maxpeak_genes)


                # If there is only one gene in the widepeak, add it as a target
                # as it will also be in the max peak
                if len(widepeak_set) == 1:
                    found_targets[cna_type].update(widepeak_set)
                    short_name = widepeak_genes
                else:
                    present_gene_targets = set(target_genes[cna_type]).intersection(widepeak_set)
                    if len(present_gene_targets) >= 1:
                        # Gene targets found in wide peak
                        found_targets[cna_type].update(present_gene_targets)
                        short_name = list(present_gene_targets)
                    elif len(maxpeak_set) == 1:
                        # If one gene in the maxpeak and no intersection between
                        # targets and genes in widepeak, use maxpeak gene as target
                        found_targets[cna_type].update(maxpeak_set)
                        short_name = maxpeak_genes
                    else:
                        short_name = widepeak_genes

                ## Here we format and collect data used in the CoMEt output ##

                # map short names to full peak events
                short_name = ','.join(short_name)+ '('+cna_type[0]+')'
                peak_names[short_name] = ','.join(widepeak_genes)+ '('+cna_type[0]+')'

                # Associate single genes with full peak events
                for gene in widepeak_genes:
                    genes_to_peak[gene].add(short_name)


    return found_targets, genes_to_peak, peak_names

def load_target_genes(gene_target_file):
    '''
    Return a dictionary of two sets, amplified target genes
    and deleted target genes.
    '''

    targets = {'Amp':set(), 'Del':set()}

    with open(gene_target_file) as gene_file:
        for line in gene_file:
            gene, cna_type = line.rstrip().split('\t')
            if cna_type.lower() == 'amp':
                targets['Amp'].add(gene)
            elif cna_type.lower() == 'del':
                targets['Del'].add(gene)
            else:
                targets['Amp'].add(gene)
                targets['Del'].add(gene)

    return targets

def load_focal_matrix_data_wrapper(config, targets):
    '''
    Wrapper allows underlying function to be used by other modules without
    the need for a configuration file.
    '''
    amp_threshold = config.getfloat('GisticData', 'amplification_cutoff')
    del_threshold = config.getfloat('GisticData', 'deletion_cutoff')
    consistency_threshold = config.getfloat('GisticData', 'cna_consistency_threshold')
    focal_data_file = os.path.join(config.get('GisticData', 'data'), config.get('GisticData','focal_matrix'))
    process_comet = 'comet' in config.get('GisticData', 'type')
    return load_focal_matrix_data(targets, amp_threshold, del_threshold, consistency_threshold, focal_data_file, process_comet)


def load_focal_matrix_data(targets, amp_threshold, del_threshold, consistency_threshold, focal_data_file, process_comet):
    '''
    Read in focal matrix file and CNA thresholds, return two dictionaries,
    one mapping samples to genes, and another mapping genes to samples.
    NOTE: These will be filtered using the list of target genes generated
    by the user provided target gene file and processed by the peak finding function.
    '''

    gene_to_sample = defaultdict(lambda: defaultdict(set))

    sample_set = set()

    sample_list = []
    n_samples = 0

    with open(focal_data_file) as focal_data:
        for line in focal_data:
            line = line.rstrip().split('\t')

            # First line, use to fill out list of all samples
            if line[0] == 'Gene Symbol':
                sample_list = line[3:]
                n_samples = len(sample_list)
            else:
                if n_samples == 0 or sample_list == []:
                    continue
                gene = line[0]

                thresholds = line[3:]

                # If user asks for comet data, every gene must be considered, else
                # we can filter here for gene targets and increase performance.
                # Decided consistency is more important than a couple second performance 
                # increase. Will always process all genes now, filter later. Can uncomment
                # the below if statement to re-enable faster non-comet performance.
                # if process_comet or gene in targets['Amp'] or gene in targets['Del']:

                for sample_indice in range(0, n_samples):
                    sample = sample_list[sample_indice]

                    # take only first three segments of name if sample is from TCGA
                    if sample[:4] == 'TCGA':
                        sample = '-'.join((sample.split('-'))[:3])

                    sample_set.add(sample)

                    # Add samples if they meet magnitude threshold criteria
                    if float(thresholds[sample_indice]) >= amp_threshold:
                        gene_to_sample[gene]['Amp'].add(sample)

                    elif float(thresholds[sample_indice]) <= del_threshold:
                        gene_to_sample[gene]['Del'].add(sample)

    gene_to_sample_unfiltered = copy.deepcopy(gene_to_sample)
    gene_to_sample_filtered = filter_targets(copy.deepcopy(gene_to_sample), targets)
    gene_to_sample_filtered = filter_consistency(gene_to_sample_filtered, consistency_threshold)

    gene_to_sample = filter_consistency(gene_to_sample, consistency_threshold)
    sample_to_gene = create_sample_to_gene(gene_to_sample)

    sample_to_gene_filtered = create_sample_to_gene(filter_targets(gene_to_sample, targets))

    return sample_to_gene, sample_to_gene_filtered, gene_to_sample_filtered, gene_to_sample_unfiltered, sample_set

def create_sample_to_gene(gene_to_sample):
    '''
    Generate a sample to gene dictionary from gene to sample dictionary
    '''

    sample_to_gene = defaultdict(lambda: defaultdict(set))
    for gene in gene_to_sample:
        for cna_type in gene_to_sample[gene]:
            for sample in gene_to_sample[gene][cna_type]:
                sample_to_gene[sample][cna_type].add(gene)
    return sample_to_gene

def filter_consistency(gene_to_sample, consistency_threshold):
    '''
    Discard genes that are not consistently mutated in the same way at
    a threshold defined by the user (default is 0.75).
    '''

    gene_list = gene_to_sample.keys()
    for gene in gene_list:

        num_amp = float(len(gene_to_sample[gene]['Amp']))
        num_del = float(len(gene_to_sample[gene]['Del']))
        if num_amp == 0 or num_del == 0:
            continue
        if num_amp/(num_del+num_amp) >= float(consistency_threshold):
            del gene_to_sample[gene]['Del']
        elif num_del / (num_del+num_amp) >= float(consistency_threshold):
            del gene_to_sample[gene]['Amp']
        else:
            # Threshold not met by either amp or del, remove gene from consideration
            del gene_to_sample[gene]

    return gene_to_sample

def filter_targets(gene_to_sample_dict, targets):
    '''
    Discard genes that are not in the target list provided by the user.
    '''

    gene_list = gene_to_sample_dict.keys()
    for gene in gene_list:

        if gene not in targets['Amp']:
            if gene_to_sample_dict[gene]['Amp']:
                del gene_to_sample_dict[gene]['Amp']

        if gene not in targets['Del']:
            if gene_to_sample_dict[gene]['Del']:
                del gene_to_sample_dict[gene]['Del']

        if not gene_to_sample_dict[gene]['Del'] and not gene_to_sample_dict[gene]['Amp']:
            del gene_to_sample_dict[gene]

    return gene_to_sample_dict

def output_hotnet2(sample_to_gene, output_dir, prefix):
    """
    Writes to a file in the following tab separated format:
    Patient-ID Gene1(A or D) Gene2(A or D) ...
    """

    with open(os.path.join(output_dir,prefix+"_hotnet2.tsv"), 'w') as out_file:
        for sample in sorted(sample_to_gene):

            gene_list = [s + '(A)' for s in sample_to_gene[sample]['Amp']] + [s + '(D)' for s in sample_to_gene[sample]['Del']]

            out_file.write(sample+'\t'+'\t'.join(sorted(gene_list)) + '\n')

def output_samples(sample_set, output_dir, prefix):
    """
    Writes a single column files containing ALL samples (even those which were 
    filtered) from the processed GISTIC2 input.
    """

    with open(os.path.join(output_dir,prefix+"_sample_list.tsv"), 'w') as out_file:
        for sample in sample_set:
            out_file.write(sample+'\n')


def output_comet(genes_to_peak, peak_names, sample_to_gene, output_dir, prefix):
    """
    Writes two files, one with every gene found by gistic for each peak,
    another with the target genes found at each peak. Row location maps
    one to the other.
    """

    with open(os.path.join(output_dir,prefix+"_all_cna_comet.tsv"), 'w') as out_file:
        for sample in sorted(sample_to_gene):
            event_set = set()
            for gene in sample_to_gene[sample]['Del']:
                event_set.update(genes_to_peak[gene])
            for gene in sample_to_gene[sample]['Amp']:
                event_set.update(genes_to_peak[gene])

            if len(event_set) > 0:
                out_file.write(sample + '\t' + '\t'.join(sorted(list(event_set)))+'\n')

    with open(os.path.join(output_dir,prefix+"_name_map_comet.tsv"), 'w') as out_file:
        for short_name in peak_names:
            out_file.write(peak_names[short_name] + ' ' + short_name+'\n')


def output_magi_file(gene_to_sample, gene_location_dict, config, output_dir, prefix):
    '''
    Write cna data to MAGI file format.
    '''

    focal_segments_file = os.path.join(config.get('GisticData', 'data'), 
                                       config.get('GisticData', 'focal_segment'))

    amp_cutoff = config.getfloat('GisticData', 'amplification_cutoff')
    del_cutoff = config.getfloat('GisticData', 'deletion_cutoff')
    seg_slack = config.getint('GisticData', 'range_cna')



    magi_dict = defaultdict(list)
    segment_dict = load_focal_segment_file(focal_segments_file)

    for gene in gene_to_sample:
        # test if we have gene location data
        if gene in gene_location_dict:
            gene_start_location, gene_end_location, gene_chromosome = gene_location_dict[gene][0]

            for cna_type, samples in gene_to_sample[gene].iteritems():
                for sample in samples:
                    for segment_start, segment_end, segment_amplitude in segment_dict[sample][gene_chromosome]:

                        if (segment_amplitude > amp_cutoff and cna_type is "Amp"
                          and not (gene_start_location-seg_slack > segment_end
                          or gene_end_location+seg_slack < segment_start)):
                            magi_dict[gene].append((sample, cna_type, segment_start, segment_end))

                        if (segment_amplitude < del_cutoff and cna_type is "Del"
                          and not (gene_start_location-seg_slack > segment_end
                          or gene_end_location+seg_slack < segment_start)):
                            magi_dict[gene].append((sample, cna_type, segment_start, segment_end))

    header = '\t'.join(["Gene","Sample ID","CNA Type","Left","Right\n"])

    with open(os.path.join(output_dir,prefix+"_magi.tsv"), 'w') as out_file:
        out_file.write(header)
        for gene, data_list in sorted(magi_dict.items()):
            for data in data_list:
                out_file.write(gene+'\t'+'\t'.join([str(s) for s in data])+"\n")

def load_focal_segment_file(focal_segments_file):
    '''
    Parse and load focal segments file into a dictionary
    (sample id->chromosome->[focal segment information])
    '''
    segments_dict = defaultdict(lambda: defaultdict(list))

    with open(focal_segments_file) as focal_file:
        for line in focal_file:
            if line[:6] == 'Sample':
                continue
            line = line.rstrip().split('\t')
            sample = line[0]
            if sample[:4] == 'TCGA':
                sample = '-'.join((sample.split('-'))[:3])
            chromosome = line[1]

            # Add start bp location, end bp location, and CNA amplitude
            segments_dict[sample][chromosome].append((int(line[2]), int(line[3]), float(line[5])))

    return segments_dict

def summary_stats(config, target_genes, found_genes, sample_to_gene):

    # of each type of target gene

    # found genes from target list, # of non-target genes "found"

    targets_both = target_genes['Amp'].intersection(target_genes['Del'])
    targets_amp = target_genes['Amp'] - targets_both
    targets_del = target_genes['Del'] - targets_both

    non_target_peak_amp = found_genes['Amp'] - target_genes['Amp']
    non_target_peak_del = found_genes['Del'] - target_genes['Del']

    mutated_genes = 0
    for sample in sample_to_gene:
        mutated_genes += len(sample_to_gene[sample])

    output_list = []
    output_list.append("*** Summary statistics ***")
    output_list.append("*  Number of samples:\t" + str(len(sample_to_gene)))
    output_list.append("*  Total genes mutated:\t" + str(mutated_genes))
    output_list.append("*  Total gene amplification targets:\t" + str(len(targets_amp)))
    output_list.append("*  Total gene deletion targets:\t" + str(len(targets_del)))
    output_list.append("*  Total gene targets for amplification and deletion:\t" + str(len(targets_both)))
    output_list.append("*  New amplification targets found (total %s): " % str(len(non_target_peak_amp)))
    for gene in non_target_peak_amp:
        output_list.append("**\t" + gene)

    output_list.append("*  New deletion targets found (total %s): " % str(len(non_target_peak_del)))
    for gene in non_target_peak_del:
        output_list.append("**\t" + gene)




    out_dir = config.get('GisticData','output_directory')
    out_name = config.get('GisticData', 'prefix')+'_summary.txt'


    with open(os.path.join(out_dir, out_name), 'w') as write_file:
        for line in output_list:
            write_file.write(line+'\n')

def visualize_data(gene_to_sample, gene_to_sample_filtered, prefix, output_dir):

     ##################################################
    #      Unfiltered plot of top 20 mutated genes     # 
     ##################################################
    top_list = sorted(gene_to_sample, key=lambda k: len(gene_to_sample[k]['Amp'])+
                        len(gene_to_sample[k]['Del']), reverse=True)[:20]

    # Stored for later use
    most_mut_unfiltered = len(gene_to_sample[top_list[0]]['Amp'])+len(gene_to_sample[top_list[0]]['Del'])

    plot_items = []
    plot_valueA = []
    plot_valueD = []

    for gene in top_list:
        plot_items.append(gene)
        plot_valueD.append(len(gene_to_sample[gene]['Del']))
        plot_valueA.append(len(gene_to_sample[gene]['Amp']))

    ax = plt.subplot(111)
    width = 0.3

    y_pos = np.arange(len(plot_items))
    plt.yticks(y_pos, plot_items)
    plt.xlabel('Total mutations')
    plt.ylabel('Genes')
    plt.title('Unfiltered 20 Most Mutated Genes')

    ax.barh(y_pos, np.array(plot_valueA), width, align='center', alpha=0.8, color='r', label="Amp")
    ax.barh(y_pos+width, np.array(plot_valueD), width, align='center', alpha=0.8, color='b', label="Del")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_20MostMutatedUnfiltered.svg'))
    plt.close()

     ##################################################
    #       Filtered plot of top 20 mutated genes      # 
     ##################################################
 

    top_list = sorted(gene_to_sample_filtered, key=lambda k: len(gene_to_sample_filtered[k]['Amp'])+
                        len(gene_to_sample_filtered[k]['Del']), reverse=True)[:20]

    most_mut_filtered = len(gene_to_sample_filtered[top_list[0]]['Amp'])+len(gene_to_sample_filtered[top_list[0]]['Del'])

    plot_items = []
    plot_value = []

    for gene in top_list:
        plot_items.append(gene)
        plot_value.append(len(gene_to_sample_filtered[gene]['Del'])+len(gene_to_sample_filtered[gene]['Amp']))
    
    y_pos = np.arange(len(plot_items))
    plt.yticks(y_pos, plot_items)
    plt.xlabel('Total mutations')
    plt.ylabel('Genes')
    plt.title('Filtered 20 Most Mutated Genes')

    plt.barh(y_pos, np.array(plot_value), width, align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_20MostMutatedFiltered.svg'))
    plt.close()

     ##################################################
    #       Filtered plot of mutation count            # 
     ##################################################
    num_del = 0
    num_amp = 0

    for gene in gene_to_sample_filtered:
        num_del += len(gene_to_sample_filtered[gene]['Del'])
        num_amp += len(gene_to_sample_filtered[gene]['Amp'])

    plot_items = ['Deletions','Amplifications']
    plot_value = [num_del,num_amp]

    x_pos = np.arange(len(plot_items))
    plt.xticks(x_pos, plot_items)
    plt.xlabel('Mutation type')
    plt.ylabel('Number of mutations')
    plt.title('Filtered Mutation Count')

    plt.bar(x_pos, np.array(plot_value), align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_FilteredMutationCount.svg'))
    plt.close()

     ##################################################
    #       Unfiltered plot of mutation count          # 
     ##################################################

    num_del = 0
    num_amp = 0
    
    for gene in gene_to_sample:
        num_del += len(gene_to_sample[gene]['Del'])
        num_amp += len(gene_to_sample[gene]['Amp'])


    plot_items = ['Deletions','Amplifications']
    plot_value = [num_del,num_amp]

    x_pos = np.arange(len(plot_items))
    plt.xticks(x_pos, plot_items)
    plt.xlabel('Mutation type')
    plt.ylabel('Number of mutations')
    plt.title('Unfiltered Mutation Count')

    plt.bar(x_pos, np.array(plot_value), align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_UnfilteredMutationCount.svg'))
    plt.close()


     #####################################################
    #  Filtered number of mutations vs number of samples  # 
     #####################################################
    mutations_f = defaultdict(int, [(i, 0) for i in range(most_mut_filtered)])
     
    for gene in gene_to_sample_filtered:
        muts = 0
        muts += len(gene_to_sample_filtered[gene]['Del'])
        muts += len(gene_to_sample_filtered[gene]['Amp'])
        mutations_f[muts] += 1

    plt.figure(figsize=(12,6))

    x_pos = np.arange(len(mutations_f.keys()))
    plt.xticks(range(0,most_mut_filtered,4), range(0,most_mut_filtered,4),rotation=45)
    plt.xlabel('Number of mutations')
    plt.ylabel('Number of samples')
    plt.title('Filtered Mutation Count Distribution')


    plt.tick_params(axis='x', labelsize=6)
    plt.bar(x_pos, np.array(mutations_f.values()), align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_filteredNumMutationsDistribution.svg'))
    plt.close()


     #####################################################
    # Unfiltered number of mutations vs number of samples # 
     #####################################################

    mutations = defaultdict(int, [(i, 0) for i in range(most_mut_unfiltered)])
     
    for gene in gene_to_sample:
        muts = 0
        muts += len(gene_to_sample[gene]['Del'])
        muts += len(gene_to_sample[gene]['Amp'])
        mutations[muts] += 1

    plt.figure(figsize=(12,6))

    x_pos = np.arange(len(mutations.keys()))
    plt.xticks(range(0,most_mut_unfiltered,4), range(0,most_mut_unfiltered,4),rotation=45)
    plt.xlabel('Number of mutations')
    plt.ylabel('Number of samples')
    plt.title('Unfiltered Mutation Count Distribution')


    plt.tick_params(axis='x', labelsize=6)
    plt.bar(x_pos, np.array(mutations.values()), align='center', alpha=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, prefix+'_unfilteredNumMutationsDistribution.svg'))
    plt.close()


# read peaks from amp_genes or del_genes 
def read_peak_file(filein, peak_q_cutoff):

    peaks = dict()
    ind2bd = dict()
    bound = dict()
    for l in open(filein):
        if l.startswith("wide peak boundaries"):
            v = l.rstrip().split('\t')
            ind2bd = dict( (i, c) for i, c in enumerate(v[1:]))
            bound = dict((c, list()) for c in v[1:])

    for l in open(filein):
        if l.startswith("wide peak boundaries"):
            continue

        elif l.startswith("residual"):
            v = l.rstrip().split('\t')
            for i, q in enumerate(v[1:]):
                if float(q) > peak_q_cutoff:
                    del bound[ind2bd[i]]
        elif l.startswith("q value"):
            continue
        elif l.startswith("cytoband"):
            v = l.rstrip().split("\t")
            for i, b in enumerate(v[1:]):
                if ind2bd[i] in bound:
                    peaks[ind2bd[i]] = b

        elif l.startswith("genes in wide peak"):
            v = l.rstrip().split("\t")
            for i, g in enumerate(v[1:]):
                if g and g[-1] == "]":
                    g = g[1:-1]
                if g and ind2bd[i] in bound:
                    bound[ind2bd[i]].append(g)

        else: # genes
            v = l.rstrip().split("\t")
            for i, g in enumerate(v[1:]):
                if g and g[-1] == "]":
                    g = g[1:-1]
                if g and ind2bd[i] in bound:
                    bound[ind2bd[i]].append(g)

    return peaks, bound


# write peak in to table format ( each row as one peak )
def write_table(fout, genes, peak):
    out = list()
    header = "index\tchromosome\tregion_start\tregion_end\tpeak_start\tpeak_end\tenlarged_peak_start\tenlarged_peak_end\tn_genes_in_region\tgenes_in_region\tn_genes_in_peak\tgenes_in_peak\tn_genes_on_chip\tgenes_on_chip\ttop 3"
    out.append(header)
    for i, c in enumerate(genes.keys()):
        chrm, start, end = c.split(":")[0][3:], c.split(":")[1].split("-")[0], c.split(":")[1].split("-")[1]
        out.append("\t".join([str(i), chrm, start, end, start, end, start, end, str(len(genes[c])), ",".join(genes[c])+",", str(len(genes[c])), ",".join(genes[c])+",", str(len(genes[c])), ",".join(genes[c])+",", ""]))

    fout.write("\n".join(out))


def run(config):
    '''
    Main function.
    '''

    # Find and load in target genes and Gistic2 data for processing

    with open(config.get('GisticData', 'gene_dictionary')) as json_file:
        gene_location_dict = json.load(json_file)


    target_gene_list = load_target_genes(config.get('GisticData', 'target_genes'))

    if config.getboolean('GisticData', 'alternate'):

        ampPeaks = os.path.join(config.get('GisticData', 'data'), config.get('GisticData', 'alt_amp_file'))
        delPeaks = os.path.join(config.get('GisticData', 'data'), config.get('GisticData', 'alt_del_file'))

        amp_threshold = config.getfloat('GisticData', 'amplification_cutoff')
        del_threshold = config.getfloat('GisticData', 'deletion_cutoff')

        peak_amp, bound_amp = read_peak_file(ampPeaks, amp_threshold)
        peak_del, bound_del = read_peak_file(delPeaks, del_threshold)

        config.set('GisticData', 'amp_file', 'table_amp.conf_99.txt')
        config.set('GisticData', 'del_file', 'table_del.conf_99.txt')

        amp_loc =  os.path.join(config.get('GisticData', 'data'), 'table_amp.conf_99.txt')
        del_loc = os.path.join(config.get('GisticData', 'data'), 'table_del.conf_99.txt')

        with open(amp_loc, 'w') as fout:
            write_table(fout, bound_amp, peak_amp)
        with open(del_loc, 'w') as fout:
            write_table(fout, bound_del, peak_del)


    found_targets, genes_to_peak, peak_names = get_peak_targets_wrapper(target_gene_list, config)
    sample_to_gene, sample_to_gene_filtered, gene_to_sample_filtered, gene_to_sample, sample_set = load_focal_matrix_data_wrapper(config, found_targets)


    if 'hotnet2' in config.get('GisticData', 'type'):
        output_hotnet2(sample_to_gene_filtered, config.get('GisticData', 'output_directory'), config.get('GisticData', 'prefix'))
    if 'magi' in config.get('GisticData', 'type') and not config.getboolean('GisticData', 'alternate'):
        output_magi_file(gene_to_sample_filtered, gene_location_dict, config, 
                        config.get('GisticData', 'output_directory'), config.get('GisticData', 'prefix'))
    if 'comet' in config.get('GisticData', 'type'):
        output_comet(genes_to_peak, peak_names, sample_to_gene, config.get('GisticData', 'output_directory'), config.get('GisticData', 'prefix'))

    output_samples(sample_set, config.get('GisticData', 'output_directory'), config.get('GisticData', 'prefix'))


    # remove temp folder if it exists
    if config.getboolean('temp_folder_exists', 'value'):
        import shutil
        shutil.rmtree(config.get('GisticData', 'temp_folder'))

    if config.getboolean('GisticData','visualizations'):
        visualize_data(gene_to_sample, gene_to_sample_filtered, config.get('GisticData', 'prefix'), config.get('GisticData','output_directory'))

    if config.getboolean('GisticData','statistics'):
        summary_stats(config, target_gene_list, found_targets, sample_to_gene)

############################################################
#### Configuration, argument parsing, and main function ####
############################################################

def check_automatic_alternate(config):
    '''
    Attempt to detect if GISTIC data is using alternate format
    '''

    gistic_data = config.get('GisticData', 'data')

    alternate = False
    norm_list = []


    # If data is in a tar, search tar for alternate file names
    if gistic_data.endswith('.gz') or gistic_data.endswith('.tar'):
        tar = tarfile.open(gistic_data)
        file_list = tar.getmembers()
        for member in file_list:
            if member.name.endswith(config.get('GisticData','alt_amp_file')) or member.name.endswith(config.get('GisticData','alt_del_file')):
                alternate = True
            if member.name.endswith(config.get('GisticData','amp_file')) or member.name.endswith(config.get('GisticData','del_file')):
                norm_list.append(member.name)
                
    else:
        for _, _, file_list in os.walk(gistic_data):
            for member in file_list:
                if member.endswith(config.get('GisticData','alt_amp_file')) or member.endswith(config.get('GisticData','alt_del_file')):
                    alternate = True
                if member.endswith(config.get('GisticData','amp_file')) or member.endswith(config.get('GisticData','del_file')):
                    norm_list.append(member)

    if len(norm_list) == 2: # both required table files found
        return False
    elif alternate:
        return True
    else: 
        return None

def get_config(args):
    '''
    Find and parse the configuration file.
    '''

    # Get and set up config file
    config = ConfigParser.SafeConfigParser()

    cfg_loc = os.path.dirname(os.path.realpath(sys.argv[0]))
    found = config.read(os.path.join(cfg_loc, 'gistic.cfg'))

    if not found:
        raise IOError("Config file not found!")

    config.add_section('temp_folder_exists')
    config.set('temp_folder_exists', 'value','false')

    # config.set('GisticData', 'types', str(args.type))


    # Integrate command line arguments into configuration
    for arg in vars(args):
        attribute = getattr(args, arg)
        if attribute:
            config.set('GisticData', arg, str(attribute))

    if config.getboolean('GisticData', 'check_alternate'):
        if check_automatic_alternate(config):
            config.set('GisticData', 'alternate', 'True')
            print "GISTIC2 alt format found"
        config.set('GisticData', 'alternate', 'False')
    else:
        config.set('GisticData', 'alternate', 'False')


    # check if gistic data is given in a tar or gzip file
    if (config.get('GisticData', 'data').endswith('.gz')
            or config.get('GisticData', 'data').endswith('.tar')):
        config.set('temp_folder_exists', 'value', 'true')
        tar = tarfile.open(config.get('GisticData', 'data'))
        file_list = tar.getmembers()
        gistic_files = []
        for member in file_list:

            if not config.getboolean('GisticData', 'alternate'):
                if member.name.endswith(config.get('GisticData','amp_file')):
                    gistic_files.append(member)
                if member.name.endswith(config.get('GisticData','del_file')):
                    gistic_files.append(member)
            else:
                # These should only find a file if the GISTIC output is in the alternate format
                if member.name.endswith(config.get('GisticData','alt_amp_file')):
                    gistic_files.append(member)
                if member.name.endswith(config.get('GisticData','alt_del_file')):
                    gistic_files.append(member)
            if member.name.endswith(config.get('GisticData','focal_matrix')):
                gistic_files.append(member)
            if member.name.endswith(config.get('GisticData','focal_segment')):
                gistic_files.append(member)

        # flatten extraction directory
        for member in gistic_files:
            member.name = os.path.basename(member.name)
            tar.extract(member,config.get('GisticData', 'temp_folder'))

        config.set('GisticData', 'data', config.get('GisticData', 'temp_folder'))

    # Ensure required files exist
    file_loc = config.get('GisticData', 'data')

    if not file_loc.endswith('/'):
        file_loc = file_loc + '/'
        config.set('GisticData', 'data', file_loc)


    if not config.getboolean('GisticData', 'alternate'):
        if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','amp_file'))):
            raise IOError('Cannot find amplified peak file: ' + config.get('GisticData','amp_file'))

        if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','del_file'))):
            raise IOError('Cannot find deleted peak file: ' + config.get('GisticData','del_file'))
    else:
        if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','alt_amp_file'))):
            raise IOError('Cannot find amplified peak file: ' + config.get('GisticData','alt_amp_file'))

        if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','alt_del_file'))):
            raise IOError('Cannot find deleted peak file: ' + config.get('GisticData','alt_del_file'))


    if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','focal_matrix'))):
        raise IOError('Cannot find focal matrix file: ' + config.get('GisticData','focal_matrix'))

    if not config.getboolean('GisticData', 'alternate') and 'magi' in config.get('GisticData', 'type'):   
        if not os.path.isfile(os.path.join(file_loc, config.get('GisticData','focal_segment'))):
            raise IOError('Cannot find segment file: ' + config.get('GisticData','focal_segment'+'. Either check the segment file name is correct, '))

    if not os.path.isfile(config.get('GisticData', 'gene_dictionary')):
        raise IOError('Cannot find gene dictionary file, please check path')

    if not os.path.isfile(config.get('GisticData', 'target_genes')):
        raise IOError('Cannot find target genes file, please check path')

    if not os.path.exists(config.get('GisticData', 'output_directory')):
        os.makedirs(config.get('GisticData', 'output_directory'))

    if not config.has_option('GisticData','prefix') or config.get('GisticData','prefix') == '':
        config.set('GisticData','prefix', 'output')


    return config

def get_parser():
    '''
    Parse arguments.
    '''
    description = 'Creates CNA to TSV data for MAGI, HotNet2, and CoMEt.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--type', help="Output format options: magi,hotnet2,comet", type=str.lower,
                         choices=['magi', 'hotnet2', 'comet'],nargs='+')
    parser.add_argument('-o', '--output_directory', help='Output directory.')
    parser.add_argument('-g', '--gene_dictionary', help='Gene location file')
    parser.add_argument('-d', '--data', help='Directory or tar file containing gistic2 output data')
    parser.add_argument('-tf', '--temp_folder', help='Specify name and location of temporary \
                                folder to extract tar file to. This will be cleaned on script exit.')
    parser.add_argument('-tg', '--target_genes', help='Location of gene targets file. \
                                Two columns, tab separated, "gene_name Amp/Del/both"')
    parser.add_argument('-ap', '--amp_file', help='Name of Gistic2 output file (not path!) containing amplified gene peaks')
    parser.add_argument('-dp', '--del_file', help='Name of Gistic2 output file (not path!) containing deleted gene peaks')
    parser.add_argument('-altap', '--alt_amp_file', help='Name of Gistic2 output file (not path!) containing alternative amplified gene peaks')
    parser.add_argument('-altdp', '--alt_del_file', help='Name of Gistic2 output file (not path!) containing alternative deleted gene peaks')
    parser.add_argument('-fm', '--focal_matrix', help='Name of Gistic2 output file (not path!) containing focal matrix data')
    parser.add_argument('-fs', '--focal_segment', help='Name of Gistic2 output file (not path!) containing segments data')
    parser.add_argument('-r', '--range_cna', help='Tolerant range of CNA are included in the browser.')
    parser.add_argument('-ac', '--amplification_cutoff', help='Amplification changes to be considered.')
    parser.add_argument('-dc', '--deletion_cutoff', help='Deletion changes to be considered.')
    parser.add_argument('-cct', '--cna_consistency_threshold', help='CNA cna_consistency_threshold to be considered.')
    parser.add_argument('-p', '--prefix', help='String to prefix output files.')
    parser.add_argument('-v', '--visualizations', action='store_true')
    parser.add_argument('-s', '--statistics', action='store_true', help = "Output GISTIC2 processing statistics")
    parser.add_argument('-ca', '--check_alternate', action='store_true', help = "Flag to indicate if the script should check for aternate GISTIC format.")
    return parser

if __name__ == '__main__':
    run(get_config( get_parser().parse_args( sys.argv[1:])))
