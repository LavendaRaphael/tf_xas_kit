def def_exp_info_json(
        class_structure,
        str_jsonfile,
        ):
#----------------------------------------------[]
    str_datfile = class_structure.str_datfile
    int_xcolumn = class_structure.list1d_column[0]
    int_ycolumn = class_structure.list1d_column[1]
    _, array2d_xdata_origin = def_extract(
        str_datfile = str_datfile,
        list1d_column = [int_xcolumn],
        )

    _, array2d_ydata_origin = def_extract(
        str_datfile=str_datfile,
        list1d_column = [int_ycolumn],
        )

    dict_scaling = def_findscaling_dict(
        array1d_xdata = array2d_xdata_origin,
        array1d_ydata = array2d_ydata_origin,
        class_structure = class_structure,
        )
    dict_peaks = def_findpeaks(
        array1d_xdata = array2d_xdata_origin,
        array1d_ydata = array2d_ydata_origin
        )
    dict_scaling.update({
        'float_onset': dict_peaks[ 'E(eV)' ][0]
         })
    #--------------------------------------------------[output]
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump(
            obj= dict_scaling,
            fp=obj_jsonfile,
            indent=4 )

    def_endfunc()
    return

