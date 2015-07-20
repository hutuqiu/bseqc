#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
import os
import sys
import logging
import numpy as np

# ------------------------------------
#logging object
# ------------------------------------

logging.basicConfig(level=20,
                    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w")
info = logging.info
error = logging.error
warning = logging.warning

try:
    '''
    The matplotlib, python module, is used for draw the Mbias plot.
    We recommend installing it.
    '''
    import matplotlib

    try:
        matplotlib.pyplot.figure()
        from matplotlib.lines import Line2D
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        import matplotlib.cm as cm
    except:
        matplotlib.use('Agg')
        from matplotlib.lines import Line2D
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        import matplotlib.cm as cm
    pdf_off = False
except:
    info('Please install the matplotlib module to draw the Mbias plot./'
         ' If not, you can only get the Mbias table')
    pdf_off = True


def Mbias_plot(strand_p, strand_t, strand, length, length_max, ymin, ymax, ax, color):
    '''
    Use Matplotlib to draw Mbias plot for every read length in every strand
    '''
    length_label = str(length) + 'bp'
    # M/(M+U)*100%
    Mcall = strand_p[strand][length][0,] / strand_p[strand][length][1,] * 100
    (trim_5_p, trim_5_M) = (strand_t[strand][length][0] + 1, Mcall[strand_t[strand][length][0]])
    (trim_3_p, trim_3_M) = (strand_t[strand][length][1] + 1, Mcall[strand_t[strand][length][1]])
    position = np.arange(1, length + 1)
    ax.plot(position, Mcall, label=length_label, color=color, lw=5.0)
    ax.add_line(Line2D([trim_5_p, trim_5_p], [ymin - 10, ymax + 10], c=color, ls='--', lw=4.0))
    ax.add_line(Line2D([trim_3_p, trim_3_p], [ymin - 10, ymax + 10], c=color, ls='--', lw=4.0))
    title = ax.set_title(strand + ' strand', fontsize=25)
    title.set_y(1.09)
    ax.set_ylim(ymin - 10, ymax + 10)
    ax.set_xlim(-2, length_max + 2)
    ax.set_ylabel('Methylation level (%)', fontsize=25)
    #ax.xaxis.set_minor_locator(MultipleLocator(5))
    #ax.xaxis.grid(True,'minor')
    ax.xaxis.grid(True, 'major', lw=2)
    box_p = ax.get_position()
    ax.set_position([box_p.x0, box_p.y0 + box_p.height * 0.05, box_p.width, box_p.height * 0.8])
    ax.legend(bbox_to_anchor=(1.01, 1.10))


def save_Mbias_table(strand_p, strand_t, name):
    '''
    export the Mcall biases table for every read length in every strand
    '''
    info("export the Mbias table...")
    table_d = name + '_Mbias_table'
    os.popen("mkdir " + table_d)
    rows = np.array(['Methylatin', 'All'], dtype='|S20')[:, np.newaxis]
    trim_file = open(table_d + '/' + name + '_trim_file.txt', 'w')
    for s in strand_p.keys():
        for l in strand_p[s].keys():
            t_file = table_d + '/' + name + '_' + s + '_' + str(l) + '.csv'
            np.savetxt(t_file, np.hstack((rows, strand_p[s][l])), delimiter=', ', fmt='%s')
            trim_file.write('%s\t%d\t%d\t%d\n' % (s, l, strand_t[s][l][0] + 1, strand_t[s][l][1] + 1))
    trim_file.close()


def mbias_generator(strand_p, strand_t, name):
    '''
    Using plots and tables to illustrate the Mcall biases
    '''
    info("draw Mbias plot and export the Mbias table...")
    strand_num = len(strand_p.keys())
    if strand_num == 2:
        strand_list = ['++', '-+']
    else:
        strand_list = ['++', '-+', '+-', '--']
    ymax = 0
    ymin = 100
    if pdf_off:
        save_Mbias_table(strand_p, strand_t, name)
        info("Can't import matplotlib package. Just save Mcall biases table!")
        return
    else:
        save_Mbias_table(strand_p, strand_t, name)
        #determine the max values of length and m% for the ranges of xaxis and yaxis respectively
    length_num = 0
    for s in strand_list:
        if length_num < len(strand_p[s].keys()):
            length_num = len(strand_p[s].keys())
        for l in strand_p[s].keys():
            #some read length only include a few reads, which are not enough to measure the Mbias
            if 0 in strand_p[s][l][1,]:
                zero_index = np.where(strand_p[s][l][1,] == 0)
                strand_p[s][l][1,][zero_index] = strand_p[s][l][1,][zero_index] + 1

            # M/(M+U)*100%
            Mcall = strand_p[s][l][0,] / strand_p[s][l][1,] * 100
            if ymax < Mcall.max():
                ymax = Mcall.max()
            if ymin > Mcall.min():
                ymin = Mcall.min()

    info("Draw Mbias plot...")
    # construct the rows and cols for the Mcall biases plot
    # single: lengthes*2: ; paired end: lengthes*4
    fig = plt.figure(figsize=(16 * strand_num, 8 * length_num))
    fig.subplots_adjust(hspace=0.1, wspace=0.15, left=0.1, right=0.9, top=0.95, bottom=0.05)
    # draw the plot
    for i in range(1, strand_num + 1):
        length_list = strand_p[strand_list[i - 1]].keys()
        length_max = max(length_list)
        row = 0
        for l in range(1, len(length_list) + 1):
            row += 1
            ax = fig.add_subplot(length_num, strand_num, i + (l - 1) * strand_num)
            color = cm.jet(1. * (l - 1) / length_num)
            Mbias_plot(strand_p, strand_t, strand_list[i - 1], length_list[l - 1], length_max, ymin, ymax, ax, color)
        ax.set_xlabel('Position in read (bp)', fontsize=25)
    plt.suptitle('Mbias Plot', fontsize=30)
    plt.savefig(name + '_Mbias_plot.pdf', format="PDF")
    info("Save Mbias plot and Mbias table!")
    return