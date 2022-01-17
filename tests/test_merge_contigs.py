import pytest
import os
import sys

from meta_transcriptomics_pipeline_sample.meta_transcriptomics_pipeline.merge_contigs import merge_contigs

#sys.path.insert(1, "/Users/yasir/Desktop/honours2/meta_transcriptomics_pipeline_sample/meta_transcriptomics_pipeline")
#from merge_contigs import merge_contigs

nucl_file = "/sample/nucl.sam"
prot_file = "/sample/prot.sam"
nucl_out = "/sample/nucl_out.sam"
prot_out = "/sample/prot_out.sam"

nucl_out_test, prot_out_test = merge_contigs(nucl_file, prot_file, os.getcwd()) 