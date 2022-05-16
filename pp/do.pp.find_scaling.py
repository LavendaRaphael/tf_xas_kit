import os
import def_dict_structure
from tf_xas_kit import pp

dict_structure = def_dict_structure.def_dict_structure_scaling_json()
list1d_key = def_dict_structure.def_list1d_key()

for str_key in list1d_key:
    class_structure = dict_structure[ str_key ]
    str_chdir = class_structure.str_chdir
    os.chdir(str_chdir)
    print(os.getcwd())

    pp.def_scaling_json(
        list1d_alignangle = [20, 90, 'trigonal'],
        list1d_scalingangle = [20, 90, 'trigonal'],
        class_structure = class_structure,
        )

