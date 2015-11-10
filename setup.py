from setuptools import setup

setup(
    name = "boobook",
    version = "0.1",
    py_modules = ['boobook'],
    install_requires = [
        'click>=5.1',
        'biopython>=1.65',
        'HTSeq>=0.6.1'
    ],
    entry_points = '''
        [console_scripts]
        boobook = boobook:cli
        ''',
    author = "Anders Goncalves da Silva",
    author_email = "andersgs@gmail.com",
    description = "A tool to transform your RNAseq data into a format suitable for analysis with Degust.",
    license = "GPL3",
    keywords = "RNASeq Degust BWA HTSeq",
    url = "https://github.com/MDU-PHL/boobook"
)
