#!/usr/bin/python

# ------------------------------------
#python package
# ------------------------------------
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
         ' If not, you can not get the nonuniform report')
    pdf_off = True


def nonuniform_generator(trim_position, name):
    '''
    generate a distribution of trimming position for nonuniform trimming
    '''

    if pdf_off:
        info('Please install the matplotlib module to draw the Mbias plot./'
             ' If not, you can not get the nonuniform QC report')
    else:
        #trim_p_array = np.array(trim_position)
        trim_p_percent = [i / float(sum(trim_position)) for i in trim_position]
        ind = np.arange(len(trim_position))
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(ind, trim_p_percent, color='red')
        ax.set_ylim(-0, max(trim_p_percent) + 0.1)
        ax.set_xlim(0, ind[-1])
        ax.set_xlabel('The number of trimmed nucleotides')
        ax.set_ylabel('The percentage of reads')
        plt.savefig(name + '_trimmed_nucleotides_dis.pdf', format='PDF')
        info('The distribution of the number of trimmed nucleotides in trimming has been drawn!')
    return