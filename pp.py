import numpy
import math
import os
import json
from tf_xas_kit import io
from tf_xas_kit import vasp
from tf_xas_kit import feff
from tf_xas_kit import data_tool

def def_scaling_json(
        list1d_alignangle,
        list1d_scalingangle,
        instance_structure,
        str_datfile = 'xas.ave.csv'
        ):
#----------------------------------------------[]
# list_alignangle = [ alpha0, beta0 ]
# list_scalingangle = [ alpha1, beta1 ]
#----------------------------------------------[]
    io.def_startfunc( locals(), ['instance_structure'] )
    #--------------------------------------------------[extract]
    list1d_column = [0]
    _, array2d_xdata_origin = io.def_extract(
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    list1d_column = [1,2,3]
    _, array2d_ydata_origin = io.def_extract(
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    #--------------------------------------------------[align]
    _, array2d_ydata_alphabeta = def_alphabeta(
        list2d_angle = [ list1d_alignangle ],
        array2d_ydata = array2d_ydata_origin
        )

    dict_peaks = data_tool.def_findpeaks(
        array1d_xdata = array2d_xdata_origin,
        array1d_ydata = array2d_ydata_alphabeta)
    float_align = instance_structure.float_onset - dict_peaks[ 'E(eV)' ][0]
    array1d_xdata_align = array2d_xdata_origin + float_align
    #--------------------------------------------------[scaling]
    _, array2d_ydata_alphabeta = def_alphabeta(
        list2d_angle = [ list1d_scalingangle ],
        array2d_ydata = array2d_ydata_origin
        )

    dict_scaling = def_findscaling_dict(
        array1d_xdata = array1d_xdata_align,
        array1d_ydata = array2d_ydata_alphabeta,
        tuple_mainxrange = instance_structure.tuple_mainxrange,
        tuple_postxrange = instance_structure.tuple_postxrange,
        )

    #--------------------------------------------------[output]
    str_abname = def_abname( alpha=list1d_scalingangle[0], beta=list1d_scalingangle[1])
    str_jsonfile = 'xas.'+str_abname+'.scaling.json'
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump(
            obj=dict_scaling,
            fp=obj_jsonfile,
            indent=4 )

    io.def_endfunc()
    return

def def_weight (
        alpha,
        beta,
        str_symmetry,
        ):
#----------------------------------------------[]
#-------------------------[in]
# alpha: float, degree angle
# beta:  float, degree angle
#-------------------------[out]
# array1d_weight:  list, shape(3)
#               array1d_weight[1] = array1d_weight[0]
#               sigma = [sigma_x, sigma_y, sigma_z]*array1d_weight
#----------------------------------------------[]
    alpha = math.radians (alpha)
    beta = math.radians (beta)
    array1d_weight = numpy.empty( shape=(3) )
    array1d_weight[0]= math.cos(beta)**2
    array1d_weight[1]= math.sin(alpha)**2 * math.sin(beta)**2
    array1d_weight[2]= math.cos(alpha)**2 * math.sin(beta)**2

    if (str_symmetry=='orthorhombic'):
        pass
    elif (str_symmetry=='trigonal'):
        array1d_weight[0] = ( array1d_weight[0] + array1d_weight[1] ) / 2
        array1d_weight[1] = array1d_weight[0]
    else:
        raise ValueError('str_symmetry')

    return array1d_weight

def def_abname( alpha, beta ):
    str_abname = 'a'+str(alpha)+'_b'+str(beta)
    return str_abname

def def_alphabeta( list2d_angle, array2d_ydata ):
#----------------------------------------------[]
# list2d_angle = []
# list2d_angle.append( [alpha1, beta1] )
# list2d_angle.append( [alpha2, beta2] )
# array2d_ydata:
#   x,y,z
#   ...
#----------------------------------------------[]
    io.def_startfunc( locals(), ['array2d_ydata'] )

    int_shape2dydata0 = numpy.shape(array2d_ydata)[0]
    io.def_print_paras( locals(), ['int_shape2dydata0'] )
    array2d_ydata_alphabeta = numpy.empty( shape=(int_shape2dydata0, len(list2d_angle)) )
    list1d_yheader = []
    for int_i in range(len(list2d_angle)):
        alpha, beta, str_symmetry = list2d_angle[int_i]
        str_yheader = def_abname( alpha, beta )
        list1d_yheader.append( str_yheader )

        array1d_weight = def_weight (alpha,beta,str_symmetry)
        array1d_ydata_dot = numpy.dot(array2d_ydata, array1d_weight)

        array2d_ydata_alphabeta[:,int_i] = array1d_ydata_dot

    return list1d_yheader, array2d_ydata_alphabeta

def def_align(
        array1d_xdata,
        array2d_ydata,
        list1d_alignangle,
        float_onset,
        ):
#----------------------------------------------[]
# list_alignangle = [ alpha0, beta0 ]
# list_scalingangle = [ alpha1, beta1 ]
#----------------------------------------------[]
    io.def_startfunc( locals(), ['array1d_xdata', 'array2d_ydata'] )

    array1d_xdata_origin = array1d_xdata
    array2d_ydata_origin = array2d_ydata
    #--------------------------------------------------[align]
    _, array2d_ydata_alphabeta = def_alphabeta(
        list2d_angle = [ list1d_alignangle ],
        array2d_ydata = array2d_ydata_origin
        )

    dict_peaks = data_tool.def_findpeaks(
        array1d_xdata = array1d_xdata_origin,
        array1d_ydata = array2d_ydata_alphabeta)
    float_align = float_onset - dict_peaks[ 'E(eV)' ][0]
    array1d_xdata_align = array1d_xdata_origin + float_align

    str_abname = def_abname( alpha=list1d_alignangle[0], beta=list1d_alignangle[1])
    str_jsonfile = 'xas.'+str_abname+'.align.json'
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump( obj={'float_align': float_align}, fp=obj_jsonfile, indent=4 )

    io.def_endfunc()
    return array1d_xdata_align

def def_alphabeta_workflow(
        instance_structure,
        instance_paras,
        ):
#----------------------------------------------[]
# list_alignangle = [ alpha0, beta0 ]
# list_scalingangle = [ alpha1, beta1 ]
#----------------------------------------------[]
    io.def_startfunc( locals(), ['instance_structure', 'instance_paras'] )

    str_datfile = instance_paras.str_avefile
    str_outfile = instance_paras.str_alphafile
    list1d_alignangle = instance_paras.list1d_alignangle
    list2d_angle = instance_paras.list2d_angle
    #--------------------------------------------------[extract]
    list1d_column = [0]
    list1d_xheader, array2d_xdata_origin = io.def_extract(
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    list1d_column = [1,2,3]
    _, array2d_ydata_origin = io.def_extract(
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    #--------------------------------------------------[alignscaling]

    array1d_xdata_align = def_align(
        array1d_xdata = array2d_xdata_origin,
        array2d_ydata = array2d_ydata_origin,
        list1d_alignangle = list1d_alignangle,
        float_onset = instance_structure.float_onset,
        )

    array2d_ydata_scaling = array2d_ydata_origin * instance_structure.float_scaling

    list1d_yheader, array2d_ydata_alphabeta = def_alphabeta(
        list2d_angle = list2d_angle,
        array2d_ydata = array2d_ydata_scaling
        )

    io.def_writedata(
        list2d_header = [ list1d_xheader, list1d_yheader ],
        list3d_data = [ array1d_xdata_align, array2d_ydata_alphabeta ],
        str_outfile=str_outfile
        )

    io.def_endfunc()
    return

def def_sft( array1d_xdata, str_code):

    if ( str_code == 'vasp' ):
        float_sft = vasp.def_sft_vasp()
    elif ( str_code == 'feff' ):
        float_sft = 0
    else:
        raise

    array1d_xdata_sft = array1d_xdata + float_sft
    return array1d_xdata_sft

def def_ave(
        instance_structure,
        instance_paras,
        ):

    dict_atom = instance_structure.dict_atom
    str_code = instance_structure.str_code

    list2d_data = []
    for list1d_atom in dict_atom.values():

        str_chdir = list1d_atom[1]
        float_scaling = list1d_atom[0]
        os.chdir(str_chdir)
        print(os.getcwd())

        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = def_code2xas(
            str_code = str_code,
            str_outfile = instance_paras.str_xasfile,
            log_tm2xas = instance_paras.log_tm2xas,
            str_broadmethod = instance_paras.str_broadmethod,
            float_hwhm = instance_paras.float_hwhm,
            int_broadnbin = instance_paras.int_broadnbin,
            )

        array2d_xdata_sft = def_sft( array2d_xdata, str_code )

        list_ycolums = list(range(len(array2d_ydata[0])))
        list2d_data.append( [ array2d_xdata_sft, array2d_ydata, list_ycolums, float_scaling ] )

        os.chdir('..')

    array1d_xdata_mix, array2d_ydata_mix = data_tool.def_mix(list2d_data=list2d_data)
    io.def_writedata(
        list2d_header = [ list1d_xheader, list1d_yheader ],
        list3d_data = [ array1d_xdata_mix, array2d_ydata_mix ],
        str_outfile = instance_paras.str_avefile
        )

def def_code2xas(
        str_code = 'vasp',
        str_outfile = 'xas.csv',
        log_tm2xas = False,
        str_broadmethod = None,
        float_hwhm = None,
        int_broadnbin = None,
        ):
    if ( str_code == 'vasp' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = vasp.def_vasp2xas(
            log_tm2xas = log_tm2xas,
            str_broadmethod = str_broadmethod,
            float_hwhm = float_hwhm,
            int_broadnbin = int_broadnbin,
            )
    elif ( str_code == 'feff' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = feff.def_feff2xas()
    else:
        raise

    io.def_writedata(
        list2d_header = [list1d_xheader, list1d_yheader],
        list3d_data = [array2d_xdata, array2d_ydata],
        str_outfile = str_outfile)

    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

