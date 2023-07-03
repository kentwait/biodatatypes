from setuptools import setup, find_packages

long_desc = open('README.md').read()
required = []

setup(
    name='biodatatypes',
    version='0.1.1',
    author='Kent Kawashima',
    author_email='kentkawashima@gmail.com',
    license='MIT',
    description='Nucleotide, amino acid, and codon datatypes for Python',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    url='https://github.com/kentwait/biodatatypes',
    project_urls={
        'Source': 'https://github.com/kentwait/biodatatypes',
        'Bug Tracker': 'https://github.com/kentwait/biodatatypes/issues',
    },
    keywords=['bioinformatics', 'biology', 'biological sequence', 'datatypes'],
    install_requires=required,
    packages=find_packages(),
    python_require='>=3.8',
)
