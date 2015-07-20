import sys
from distutils.core import setup
def main():
    if not float(sys.version[:3])>=2.6:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.6! python 2.7 is recommended!\n")
        sys.exit(1)
    setup(
        name='BSeQC',
        version='1.0.0',
        url='http://code.google.com/p/BSeQC/',
        license='Artistic License/GPL',
        package_dir={'BSeQC' : 'BSeQC'},
        packages = ['BSeQC', 'BSeQC.qc_assess', 'BSeQC.read', 'BSeQC.qc_filter'],
        platforms = ['Linux'],
        author='Xueqiu Lin',
        author_email='xueqiu.lin@gmail.com',
        description=' Quality Control of bisulfite sequencing experiments',
        scripts=['bin/bseqc'],
        )
if __name__ == '__main__':
    main()