#!/bin/env python
import xas_module
import os
import local_module

list1d_key = local_module.def_list1d_key()
dict_structure = local_module.def_dict_structure()

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    xas_module.def_atom_abworkflow(
        class_structure = class_structure
        )
