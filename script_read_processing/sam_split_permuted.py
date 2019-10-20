import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def samfile_barcode_count(sam_file, barcode_file):
    
    #generate the barcode list and barcode dictionary
    barcodes = open(barcode_file)
    barcode_ls = []
    barcode_dic = {}
    for line in barcodes:
        barcode = line.strip()
        barcode_ls.append(barcode)
        barcode_dic[barcode] = 0
    barcodes.close()
    #read the sam file, and count the number per barcode
    sam = open(sam_file)
    for line in sam:
        if (line[0] == '@'):
            continue
        else:
            name = (((line.split('\t'))[0]).split(','))
            barcode = name[0]
            barcode_dic[barcode] += 1
    sam.close()
    return barcode_dic

def permute_samples(sam_file, barcode_count, barcode_list, output_folder):
    '''
    this function accept a sam file, a barcode count dictionary, a barcode list and a output folder.
    then for each barcode in the barcode list, it find the reads number n associated with the barcode, 
    and sample n reads from the samfile and then output the sampled reads to the 
    '''
    
    # Generate a list of output file
    file_name = (sam_file.split('/')[-1]).split('.')[0]
    output_files = {}
    for barcode in barcode_list:
        output_file = output_folder + '/' + file_name + '.' + barcode + '.permuted.sam'
        output_files[barcode] = open(output_file, 'w')
    
    # output the header into each file and count the line number in the header and total reads number
    header_number = 0
    input_file = open(sam_file)
    all_lines = input_file.readlines()
    for line in all_lines:
        if (line[0] == '@'):
            header_number += 1
            for barcode in barcode_list:
                output_files[barcode].write(line)
        else:
            break
    
    # for each barcode, generate the permuted line array
    # for each barcode, output the permuted lines to the output file
    permuted_lines = np.random.permutation(all_lines[header_number:])
    first_line = 0
    for barcode in barcode_list:
        end_line = first_line + barcode_count[barcode]
        output_lines = list(permuted_lines[first_line:end_line])
        for output_line in output_lines:
            output_files[barcode].write(output_line)
        first_line = end_line
    
    # close the output file and sam file
    for barcode in barcode_list:
        output_files[barcode].close()
    input_file.close()
    
def split_samfile(sam_file, barcode_file, output_folder, cutoff):
    '''
    this script accept a sam file, a barcode file, a output_file, a cutoff value,
    then it will call the samfile_barcode_count function and get the total read count per barcode,
    then it use the cutoff value to filter the barcode,
    and generate the output samfile for single cells, generate the sample_ID.txt in the output folder,
    generate the reads distribution in the output folder/read_distribution_barcode;
    '''
    
    # generate the count per barcode
    barcode_count = samfile_barcode_count(sam_file, barcode_file)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (sam_file.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('frequency')
    plt.xlabel('Number of unique reads')
    fig_output = output_folder + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    read_dist = open(output_folder + '/' + plot_name + '.txt', 'w')
    for barcode in barcode_count:
        line = barcode + ', %d\n' %(barcode_count[barcode])
        read_dist.write(line)
    read_dist.close()
    
    #filter the barcode based on the cutoff value
    barcode_filtered = []
    for barcode in barcode_count:
        if barcode_count[barcode] >= cutoff:
            barcode_filtered.append(barcode)
    #print barcode_filtered
     
    # for the barcode in the barcode filter list, generate permuted sample in the output folder
    print "Generatign permuted sequences..."
    permute_samples(sam_file, barcode_count, barcode_filtered, output_folder)
    
    #generate the output sam file and sample_list file
    sample_list_file = open(output_folder + '/' + plot_name + '.' + 'sample_list.txt', 'w')
    output_files = {}
    for barcode in barcode_filtered:
        output_file = output_folder + '/' + plot_name + '.' + barcode + '.sam'
        output_files[barcode] = open(output_file, 'w')
        sample_list_file.write(plot_name + '.' + barcode + '\n')
    
    # output the each read to the output sam file
    sam = open(sam_file)
    for line in sam:
        if (line[0] == '@'):
            for barcode in barcode_filtered:
                output_files[barcode].write(line)
        else:
            barcode = (((line.split('\t'))[0]).split(','))[0]
            if barcode in barcode_filtered:
                output_files[barcode].write(line)
    
    #close the files:
    sample_list_file.close()
    sam.close()
    for barcode in barcode_filtered:
        output_files[barcode].close()

if __name__ == '__main__':
    sam_file = sys.argv[1]
    barcode_file = sys.argv[2]
    output_folder = sys.argv[3]
    cutoff = int(sys.argv[4])
    split_samfile(sam_file, barcode_file, output_folder, cutoff)
