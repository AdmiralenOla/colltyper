from setuptools import setup
from __init__ import __version__ as cv

setup(name='colltyper',
    version=cv,
    description='Type TB strains using Coll (2015) scheme',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['TB','lineage','classification','typing','mycobacterium','tuberculosis'],
    url='https://github.com/AdmiralenOla/colltyper',
    author='Ola Brynildsrud',
    author_email='ola.brynildsrud@fhi.no',
    license='MIT',
    packages=['colltyper'],
    install_requires=[
        'argparse',
        'pyvcf'],
    entry_points={
        'console_scripts': ['colltyper=colltyper.colltyper:main']
    },
    package_data={'colltyper': ['LICENSE', 'setup.py','example/*','docs/*','data/*']},
    include_package_data=True
    )