#!/bin/env python
import exp
import os
import exp_setup

list1d_key = exp_setup.def_list1d_key()
dict_structure = exp_setup.def_dict_structure()

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())
    
    exp.def_exp_info_json(
        class_structure = class_structure,
        str_jsonfile = 'exp_info.json',
        )
