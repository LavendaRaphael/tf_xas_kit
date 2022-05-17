import pandas
import numpy
import json
from tf_xas_kit import io

def def_exp_scaling(
        class_structure,
        str_outfile,
        class_paras,
        ):
    str_datfile = class_structure.str_datfile
    int_xcolumn = class_structure.list1d_column[0]
    int_ycolumn = class_structure.list1d_column[1]
    df_xas = pandas.read_csv(
        str_datfile,
        header = 0,
        names= ['E(eV)', 'Intensity'],
        usecols = [int_xcolumn,int_ycolumn],
        comment='#',
        )
    df_xas = df_xas.dropna()
    print(df_xas)
    dict_scaling = data_tool.def_findscaling_dict(
        array1d_xdata = df_xas['E(eV)'].to_numpy(),
        array1d_ydata = df_xas['Intensity'].to_numpy(),
        class_structure = class_structure
        )
    float_scaling = dict_scaling[ class_paras.str_scalingmethod ]
    df_xas['Intensity'] *= float_scaling

    df_xas.to_csv( str_outfile )

def def_exp_xyzfit( list2d_alpha, str_outfile ):
#----------------------------------------------[]
# list2d_alpha = []
#----------------------------------------------[]
   io.def_startfunc( locals() )

    int_lenalpha = len(list2d_alpha)
    array1d_alpha = numpy.empty( shape=(int_lenalpha) )
    list2d_data = []
    for int_i in range(int_lenalpha):
        array1d_alpha[int_i] = list2d_alpha[ int_i ][ 0 ]
        str_datfile = list2d_alpha[ int_i ][ 1 ]
        list1d_xycolumn = list2d_alpha[ int_i ][ 2 ]
        list1d_xheader, array2d_xdata = io.def_extract( str_datfile=str_datfile, list1d_column=[ list1d_xycolumn[0] ] )
        list1d_yheader, array2d_ydata = io.def_extract( str_datfile=str_datfile, list1d_column=[ list1d_xycolumn[1] ] )
        list2d_data.append( [ array2d_xdata, array2d_ydata] )

    array1d_cosalpha2 = numpy.cos( numpy.radians( array1d_alpha ) ) **2
    def_print_paras( locals(), ['array1d_cosalpha2'] )

    array1d_xdata_interp, list1d_ydata_interp = data_tool.def_interp( list2d_data )

    int_len1dxdata = len( array1d_xdata_interp )
    def_print_paras( locals(), ['int_len1dxdata'] )

    array2d_ydata_fit = numpy.empty( shape=(int_len1dxdata,2) )
    array1d_temp = numpy.empty( shape=(int_lenalpha) )
    for int_i in range( int_len1dxdata ):
        for int_j in range( int_lenalpha ):
            array1d_temp[int_j] = list1d_ydata_interp[int_j][int_i]
        polyfit = numpy.polynomial.Polynomial.fit( array1d_cosalpha2, array1d_temp, 1 )
        b,k = polyfit.convert().coef
        array2d_ydata_fit[ int_i ][0] = b
        array2d_ydata_fit[ int_i ][1] = k+b

    io.def_writedata(
        list2d_header = [ ['E(eV)', 'fit_xy','fit_z'] ],
        list3d_data = [ array1d_xdata_interp, array2d_ydata_fit ],
        str_outfile=str_outfile
        )

def def_exp_info_json(
        class_structure,
        str_jsonfile,
        ):
#----------------------------------------------[]
    str_datfile = class_structure.str_datfile
    int_xcolumn = class_structure.list1d_column[0]
    int_ycolumn = class_structure.list1d_column[1]
    _, array2d_xdata_origin = io.def_extract(
        str_datfile = str_datfile,
        list1d_column = [int_xcolumn],
        )

    _, array2d_ydata_origin = io.def_extract(
        str_datfile=str_datfile,
        list1d_column = [int_ycolumn],
        )

    dict_scaling = data_tool.def_findscaling_dict(
        array1d_xdata = array2d_xdata_origin,
        array1d_ydata = array2d_ydata_origin,
        class_structure = class_structure,
        )
    dict_peaks = data_tool.def_findpeaks(
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

    io.def_endfunc()
    return

