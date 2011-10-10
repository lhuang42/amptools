from setuptools import setup

setup(
    name='amptools',
    version='0.1',
    description='Amplicon experiment tools',
    author='James Casbon',
    author_email='james.casbon@popgentech.com',
    packages=['amptools'],
    install_requires=[
        'pysam>0.6',
        'cmdln'
    ], 
    scripts=['amptools/amptools']

)