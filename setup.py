from setuptools import find_packages, setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(
    name="meta_transcriptomics_pipeline",
    version="0.2.2",
    license='MIT',
    packages=find_packages(include=['meta_transcriptomics_pipeline']),
    author="Yasir Kusay",
    author_email='yasirsaad1234@hotmail.com',
    description="To be added!.",
    ext_modules=cythonize(["meta_transcriptomics_pipeline/obtain_relevant_taxids.pyx"]),
    include_dirs=np.get_include(),
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
    entry_points = {
    "console_scripts": [
            "rapid_metatranscriptomics = meta_transcriptomics_pipeline.__main__:main",
        ]
    },
    python_requires='>=3.6.2',
    install_requires=[
        'pandas>=1.1.5',
        'numpy>=1.19.5',
        'matplotlib>=3.3.4',
        'ete3>=3.1.2',
        'cython>=0.26.0'
    ]
)

# The important parts with the setup call are the last three arguments. ext_modules specifies C extensions for the package, include_dirs needs to be specified so that the numpy-dependent C extensions can be compiled, and 
# install_requires specifies the packages thatclassification_library depends on.

# guide: https://levelup.gitconnected.com/how-to-deploy-a-cython-package-to-pypi-8217a6581f09