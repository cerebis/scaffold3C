import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='scaffold3C',
    description='Compute the order and orientation of assembly contigs using Hi-C sequencing data.',
    long_description=long_description,
    version='0.1a1',
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    platforms='Linux-86_x64',
    packages=setuptools.find_packages(),
    url='https://github.com/cerebis/scaffold3C',
    license='GNU Affero General Public License v3',

    install_requires=['proxigenomics_toolkit @ git+https://github.com/cerebis/proxigenomics_toolkit@master#egg=proxigenomics_toolkit'],

    classifiers=[
        'Programming Language :: Python :: 2.7',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 2 - Pre-Alpha'
    ]

)
