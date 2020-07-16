from setuptools import setup

setup(
    name='clumps',
    python_requires='>3.5.0',
    author='Shankara Anand - Broad Institute - Cancer Genome Computational Analysis',
    author_email='sanand@broadinstitute.org',
    description='CLUMPS driver gene discovery using 3D protein structure (Getz Lab).',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    packages = [
        'clumps',
        'clumps.samplers',
        'clumps.mapping',
        'clumps.db',
    ],
    install_requires=[
        "twobitreader>=3.1.7",
        "statsmodels>=0.9.0",
        "scipy>=1.3.0",
        "pyopenssl>=19.0.0",
        "prody>=1.10.10",
        "lxml>=4.3.4",
        "jpype1>=0.7.0",
        "canine",
        "biopython>=1.73",
        "tqdm>=4.32.1",
        "agutil"
    ],
    entry_points={
        'console_scripts':[
            'clumps-prep = clumps.prep:main',
            'clumps = clumps.clumps:main',
            'clumps-postprocess = clumps.postprocess:main'
        ]
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
