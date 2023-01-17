# list of Cython modules containing tests
'''
def test_join():
    cython_test_modules = ["test_join"]

    for mod in cython_test_modules:
        try:
            # For each callable in `mod` with name `test_*`,
            # set the result as an attribute of this module.
            mod = importlib.import_module(mod)
            assert(1 == 2)
            for name in dir(mod):
                item = getattr(mod, name)
                if callable(item) and name.startswith("test_"):
                    setattr(sys.modules[__name__], name, item)
        except ImportError:
            pass
'''

import os
import filecmp
from meta_transcriptomics_pipeline.obtain_relevant_taxids import obtain_relevant_taxids

def test_obtain_relevant_taxids(capsys):
    print("HE")
    path = "tests/samples/join/"
    #path = "samples/join/"
    #assert(os.path.exists(path + "to_map.txt"))
    #assert(os.path.exists(path + "to_map.txt"))
    #assert(os.path.exists(path + "mapping.txt"))
    #assert(os.path.exists(path + "actual1.txt"))
    #assert(os.path.exists(path + "to_map_2.txt"))
    #assert(os.path.exists(path + "mapping.txt"))
    #assert(os.path.exists(path + "actual1.txt") == False)
    #assert(os.path.exists(path + "actual2.txt") == False)
    #print()
    #obtain_relevant_taxids(path + "to_map.txt", path + "mapping.txt", path + "actual1.txt")
    #captured = capsys.readouterr()
    #print()

    if os.path.exists(path + "actual1.txt"):
        os.remove(path + "actual1.txt")
    if os.path.exists(path + "actual2.txt"):
        os.remove(path + "actual2.txt")

    obtain_relevant_taxids(path + "to_map.txt", path + "mapping.txt", path + "actual1.txt")
    obtain_relevant_taxids(path + "to_map_2.txt", path + "mapping.txt", path + "actual2.txt")
    assert(os.path.exists(path + "actual1.txt"))
    assert(os.path.exists(path + "actual2.txt"))

    assert(filecmp.cmp(path + "exp1.txt", path + "actual1.txt", shallow=False))
    assert(filecmp.cmp(path + "exp2.txt", path + "actual2.txt", shallow=False))
    #assert(os.path.exists(path + "to_map.txt"))
    print("TEST FONISHED")