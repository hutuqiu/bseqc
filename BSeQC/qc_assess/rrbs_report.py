#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
from collections import defaultdict
import sys
import logging


# ------------------------------------
#logging object
# ------------------------------------

logging.basicConfig(level=20,
                    format=' %(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
)
info = logging.info
error = logging.error
warning = logging.warning

try:
    '''
    The matplotlib, python module, is used for draw the Mbias plot.
    We recommend installing it.
    '''
    import matplotlib
    #matplotlib.use('Agg')
    try:
        matplotlib.pyplot.figure()
        from matplotlib.lines import Line2D
        import matplotlib.pyplot as plt
    except:
        matplotlib.use('Agg')
        from matplotlib.lines import Line2D
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        import matplotlib.cm as cm
    pdf_off = False
except:
    info('Please install the matplotlib module to draw the Mbias plot./'
         ' If not, you can not get the RRBS assessing report')
    pdf_off = True



def calculate_methlevel(dige_site_strand, strand):
    meth_level = []
    for site in dige_site_strand.keys():
        if dige_site_strand[site][1] == 0:
            if site == 's':
                error('Find zero restriction enzyme digestion site in the 5 prime of the %s strand read!!!' % strand)
                error('Please check the input file!!!')
            else:
                warning('Find zero restriction enzyme digestion site in the 5 prime of the %s strand read!!!' % strand)
                warning(
                    'Please check the input file!!! (It is normal for the short length read, but not for the long one)')

            meth_level.append(0.0)
        else:
            meth_level.append(dige_site_strand[site][0] / float(dige_site_strand[site][1]) * 100)
    return meth_level



def generator(dige_dict, single_on, name):
    '''
    generate the qc report of rrbs
    '''


    if single_on:
        strand_list = ['++', '-+']
        meth_level_list = [0] * 4
        meth_level_list_se = [0] * 2
        string_list = ["5'", "3'"]
    else:
        strand_list = ['++', '-+', '+-', '--']
        meth_level_list = [0] * 8
        meth_level_list_se = [0] * 4
        string_list = ["mate1 5'", "mate2 5'", "mate1 3'", "mate2 3'"]


    if pdf_off:
        info('Please install the matplotlib module to draw the Mbias plot./'
             ' If not, you can not get the rrbs QC report')
    else:
        #fig = plt.figure(16, 12)
        #plt = fig.add_subplot(1, 1, 1)
        max_meth_level = 0
        for i in range(len(strand_list)):
            meth_level = calculate_methlevel(dige_dict[strand_list[i]], strand_list[i])
            meth_level_list[i] = meth_level[0]
            meth_level_list[i + len(strand_list)] = meth_level[1]
            if max_meth_level < max(meth_level):
                max_meth_level = max(meth_level)
            #color = cm.jet(1. * i / len(strand_list))
            #plt.plot([1, 2], meth_level, label=strand_list[i], color=color, lw=1, linestyle='--', marker='o')
            #plt.plot([1, 2], meth_level, color='r', marker='o')
        for i in range(len(meth_level_list)/2):
            meth_level_list_se[i] = (meth_level_list[i*2] + meth_level_list[i*2+1])/2

        #x_location = range(1, len(meth_level_list) + 1)
        x_location = range(1, len(meth_level_list_se) + 1)
        #plt.scatter(x_location, meth_level_list, s=100, linewidths=2, c='lightgrey')
        #plt.plot(x_location, meth_level_list, linestyle='--', c='black')
        plt.scatter(x_location, meth_level_list_se, s=100, linewidths=2, c='lightgrey')
        plt.plot(x_location, meth_level_list_se, linestyle='--', c='black')
        plt.ylim(-10, max_meth_level + 20)
        #plt.xlim(0, len(meth_level_list) + 1)
        plt.xlim(0, len(meth_level_list_se) + 1)
        plt.ylabel('Methylation level (%)', fontsize=20)
        #plt.xticks([1, 2], ["5'", "3'"])
        plt.xlabel('The position of MspI site', fontsize=20)
        plt.title('The Mbias of the second C in MspI site', fontsize=20)
        #plt.legend(loc=1)
        plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        for i, txt in enumerate(string_list):
            plt.annotate(txt, (x_location[i], meth_level_list_se[i] + 2), size=15)
        plt.savefig(name + '_Mbias_MspI.pdf', format="PDF", figsize=(16, 12))
        info("Save Mbias MspI plot!")

    #calculate the frequency of the location of the restriction enzyme site
    #for key in dige_dict.keys():
    #   for s in dige_dict[key]:
    #        d = defaultdict(int)
    #        for item in dige_site[key][s][0]:
    #            d[item] += 1
    #            dige_site[key][s][0] = d.items()
    return
