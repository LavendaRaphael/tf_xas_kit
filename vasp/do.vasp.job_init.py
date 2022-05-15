import os
from tf_xas_kit import mod_class_structure
from tf_xas_kit import vasp

list1d_key = mod_class_structure.def_list1d_key()
dict_structure = mod_class_structure.def_dict_structure( )

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    vasp.def_vasp_jobinit( class_structure )
