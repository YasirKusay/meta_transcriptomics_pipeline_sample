from setuptools import find_packages, setup
from Cython.Build import cythonize
import numpy as np

setup(
    name="meta_transcriptomics_pipeline_test",
    version="0.1.14",
    packages=find_packages(include=['meta_transcriptomics_pipeline']),
    author="Yasir Kusay",
    description="To be added!.",
    # long_description_content_type='text/markdown',
    # url="https://github.com/lol-cubes/classification-library",
    ext_modules=cythonize(["meta_transcriptomics_pipeline/__init__.pyx", "meta_transcriptomics_pipeline/obtain_relevant_taxids.pyx"]),
    # ext_modules=cythonize(["__init__.pyx", "obtain_relevant_taxids.pyx"]),
    # s
    # zip_safe=False, # from https://cython.readthedocs.io/en/latest/src/quickstart/build.html
    include_dirs=np.get_include(),
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
    #setup_requires=['pytest-runner'],
    #tests_require=['pytest==4.4.1'],
    #test_suite='tests',
)

# The important parts with the setup call are the last three arguments. ext_modules specifies C extensions for the package, include_dirs needs to be specified so that the numpy-dependent C extensions can be compiled, and 
# install_requires specifies the packages thatclassification_library depends on.

# guide: https://levelup.gitconnected.com/how-to-deploy-a-cython-package-to-pypi-8217a6581f09

# python3 setup.py sdist bdist_wheel
# twine upload --repository testpypi dist/*
# pip install python3 -m pip install --index-url https://test.pypi.org/simple/ meta_transcriptomics_pipeline


'''
setup(
    name='meta_transcriptomics_pipeline',
    packages=find_packages(include=['meta_transcriptomics_pipeline']),
    version='0.1.0',
    description='meta_transcriptomics_pipeline',
    author='Yasir Kusay',
    license='MIT',
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
'''