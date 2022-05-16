#!/bin/env python
import os
from tf_xas_kit import mod_class_structure
from tf_xas_kit import class_paras_code2xas
from tf_xas_kit import pp

list1d_key = mod_class_structure.def_list1d_key()
dict_structure = mod_class_structure.def_dict_structure( )
class_paras = class_paras_code2xas.def_class_paras()

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())
    print(class_structure.dict_atom)

    pp.def_ave( 
        class_structure,
        class_paras,
        )
