#!/bin/env python
#===============================================[README]
# @FeifeiTian
# 2021.10.13
#===============================================<<
import csv
import json
from ase.calculators import vasp
import ase.io
import scipy.signal
import numpy
import copy
import inspect
import math
import os
import local_module
import pandas
from distutils.dir_util import copy_tree
import re
import group_module
import subprocess

def def_atom_findpeak( 
        list1d_angle,
        str_workdir,
        str_jsonfile
        ):
    def_startfunc( locals() )

    _, _, array2d_xdata, array2d_ydata = def_vasp_outcar2xas()

    float_eigcore = -array2d_xdata[0][0] / 0.9
    # float_eigcore = -514.703961144
    def_print_paras( locals(), ['float_eigcore'])

    _, array2d_ydata_alphabeta = def_alphabeta( 
        list2d_angle = [list1d_angle], 
        array2d_ydata = array2d_ydata 
        )
    dict_peaks = def_findpeaks(
        array1d_xdata=array2d_xdata, 
        array1d_ydata=array2d_ydata_alphabeta
        )
    
    str_cwddir = os.getcwd()
    os.chdir('..')
    with open(str_jsonfile) as obj_jsonfile:
        float_align = json.load( fp=obj_jsonfile )['float_align']
    dict_structure = local_module.def_dict_structure()
    str_chdir = dict_structure[ str_workdir ].dict_atom[1][1]
    os.chdir(str_chdir)
    float_finalenergy_1 = def_vasp_finalenergy()
    os.chdir( str_cwddir )
    float_finalenergy = def_vasp_finalenergy()
    float_sft = float_finalenergy-float_finalenergy_1
    float_sftplusalign = float_sft + float_align
    def_print_paras( locals(), ['float_align','float_sft','float_sftplusalign'])

    array1d_energy_origin = dict_peaks[ 'E(eV)' ]
    array1d_energy_pluscore = dict_peaks[ 'E(eV)' ] + float_eigcore
    array1d_energy_sftplusalign = dict_peaks[ 'E(eV)' ] + float_sftplusalign
    def_print_paras( locals(), ['array1d_energy_origin','array1d_energy_pluscore', 'array1d_energy_sftplusalign'])

    def_endfunc()

def def_chgrdf_workflow( 
        str_chgfile,
        str_outdir='',
        ):
    def_startfunc( locals() )

    if ( not def_has_numbers(str_chgfile) ):
        str_outfile = 'chgrdf.csv'
    elif ( str_chgfile[:6] == 'PARCHG' ):
        str_outfile = 'chgrdf.B' + str_chgfile[ 7:11] + '_K' + str_chgfile[12:16] + '.csv'
    elif ( str_chgfile[:3] == 'WFN'):
        str_outfile = 'chgrdf.B' + str_chgfile[13:17] + '_K' + str_chgfile[19:23] + '.csv'
    def_print_paras( locals(), ['str_outfile'] )

    array1d_r, array1d_rpd, array1d_rpi, array1d_rdf = def_chgrdf( str_chgfile=str_chgfile)

    def_writedata( 
            list2d_header = [ ['r(ang)'], ['chgrpd'], ['chgrpi'], ['chgrdf'] ],
            list3d_data = [ array1d_r, array1d_rpd, array1d_rpi, array1d_rdf ],
            str_outfile = str_outdir+str_outfile
            )

    def_endfunc()
    return

def def_has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

def def_chgrdf( 
        str_chgfile,    
        float_r0=3, 
        float_slice = 0.05
        ):
    def_startfunc( locals() )

    obj_chgcar = vasp.VaspChargeDensity(filename=str_chgfile)
    obj_atoms = obj_chgcar.atoms[0]
    array3d_chgdens = obj_chgcar.chg[0]

    array1d_cell_ngrid = numpy.array(numpy.shape(array3d_chgdens))
    int_ngrid = 1
    for int_i in array1d_cell_ngrid:
        int_ngrid *= int_i
    float_volume=obj_atoms.get_volume()
    if ( 'CHG' in str_chgfile ):
        float_chgsum = numpy.sum(array3d_chgdens) * float_volume / int_ngrid
    elif ( 'WFN' in str_chgfile ):
        float_chgsum = numpy.sum(array3d_chgdens) * float_volume * 2
    else:
        raise ValueError('str_chgfile not clear!')
    def_print_paras( locals(), ['array1d_cell_ngrid', 'float_volume', 'float_chgsum'] )

    array1d_atom1_pos = obj_atoms.get_positions()[0]
    array1d_cell_paras = obj_atoms.cell.cellpar()[0:3]
    array1d_grid_paras = array1d_cell_paras/array1d_cell_ngrid
    array1d_atom1_ngrid = array1d_atom1_pos/array1d_grid_paras
    def_print_paras( locals(), ['array1d_atom1_pos','array1d_grid_paras', 'array1d_cell_paras','array1d_atom1_ngrid'] )

    int_nslice = int(float_r0 // float_slice) + 1
    array1d_rpd = numpy.zeros( shape=(int_nslice) )
    float_r0_new = float_slice * int_nslice
    array1d_r0_ngrid = float_r0_new/array1d_grid_paras
    def_print_paras( locals(), ['int_nslice','float_r0_new','array1d_r0_ngrid'] )

    array1d_rangel = numpy.ceil(array1d_atom1_ngrid - array1d_r0_ngrid).astype(int)
    array1d_ranger = numpy.floor(array1d_atom1_ngrid + array1d_r0_ngrid).astype(int)
    def_print_paras( locals(), ['array1d_rangel','array1d_ranger'] )

    #array3d_chgdens_test = numpy.zeros( shape=array1d_cell_ngrid )

    for int_x in range( array1d_rangel[0], array1d_ranger[0]+1 ):
        print(int_x)
        for int_y in range( array1d_rangel[1], array1d_ranger[1]+1 ):
            for int_z in range( array1d_rangel[2], array1d_ranger[2]+1 ):
                array1d_ngrid = numpy.array( [int_x, int_y, int_z] )
                array1d_dist = ( array1d_ngrid - array1d_atom1_ngrid ) * array1d_grid_paras
                float_dist = numpy.sqrt( array1d_dist.dot( array1d_dist ) )
                if (float_dist >= float_r0_new): continue
                int_temp = int( float_dist // float_slice)
                array1d_ngrid = array1d_ngrid % array1d_cell_ngrid
                array1d_rpd[ int_temp ] += array3d_chgdens[ array1d_ngrid[0], array1d_ngrid[1], array1d_ngrid[2] ]
                #array3d_chgdens_test[ array1d_ngrid[0], array1d_ngrid[1], array1d_ngrid[2] ] = (
                #    array3d_chgdens[ array1d_ngrid[0], array1d_ngrid[1], array1d_ngrid[2] ]
                #    )

    del array3d_chgdens
    #obj_chgcar.chg[0] = array3d_chgdens_test
    #obj_chgcar.write( filename='CHG_test.vasp' )

    if ( 'CHG' in str_chgfile ):
        array1d_rpd *= float_volume / int_ngrid / float_slice
    elif ('WFN' in str_chgfile):
        array1d_rpd *= float_volume * 2 / float_slice
    array1d_rpd /= float_chgsum

    array1d_r = numpy.linspace( 0, float_r0_new, num=int_nslice, endpoint=False )
    array1d_r += float_slice/2

    array1d_rpi = numpy.empty( shape=(int_nslice) ) 
    array1d_rpi[ 0 ] = array1d_rpd[ 0 ]
    for int_i in range( 1, int_nslice ):
        array1d_rpi[ int_i ] = array1d_rpi[ int_i-1 ] + array1d_rpd[ int_i ]
    for int_i in range( int_nslice ):
        array1d_rpi[ int_i ] -= array1d_rpd[ int_i ]/2
    array1d_rpi *= float_slice

    array1d_rdf = numpy.empty( shape=(int_nslice) )
    for int_i in range( int_nslice ):
        array1d_rdf[ int_i ] = array1d_rpd[int_i] / ( 3*int_i**2 + 3*int_i + 1 )
    array1d_rdf /= 4/3 * math.pi * float_slice**2
    #float_rho0 = 1 / float_volume
    #array1d_rdf /= float_rho0

    def_endfunc()
    return array1d_r, array1d_rpd, array1d_rpi, array1d_rdf

def def_tm_findmax( 
        array1d_xdata, 
        array1d_ydata, 
        array1d_kb, 
        float_onset, 
        str_abname, 
        int_ntm=1, 
        float_xwidth=0.5 ):
    def_startfunc( locals(), ['array1d_xdata', 'array1d_ydata','array1d_kb'] )

    array1d_xdata = numpy.reshape( array1d_xdata, newshape=-1 )
    array1d_ydata = numpy.reshape( array1d_ydata, newshape=-1 )

    str_jsonfile = 'xas_tm.'+str_abname+'.findtm.json'
    list1d_index_xwidth = []
    for int_i in range( numpy.shape(array1d_xdata)[0] ):
        float_x = array1d_xdata[ int_i ]
        if ( float_x < float_onset-float_xwidth or float_x > float_onset+float_xwidth): continue
        list1d_index_xwidth.append( int_i )
    array1d_index_xwidth_topn = numpy.argpartition( array1d_ydata[ list1d_index_xwidth ], -int_ntm )[-int_ntm:]
    array1d_index_topn = [ list1d_index_xwidth[int_i] for int_i in array1d_index_xwidth_topn ]

    array1d_xdata_topn = array1d_xdata[ array1d_index_topn ]
    array1d_ydata_topn = array1d_ydata[ array1d_index_topn ]
    array1d_kb_topn = array1d_kb[ array1d_index_topn ]

    dict_jsonfile = {}
    dict_jsonfile[ 'array1d_index_topn' ] = array1d_index_topn
    dict_jsonfile[ 'array1d_xdata_topn' ] = array1d_xdata_topn
    dict_jsonfile[ 'array1d_ydata_topn' ] = array1d_ydata_topn
    dict_jsonfile[ 'array1d_kb_topn' ] = array1d_kb_topn
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump( obj=dict_jsonfile, fp=obj_jsonfile, indent=4, cls=NumpyEncoder )

    def_endfunc()
    return array1d_index_topn

def def_atom_abworkflow(  
        class_structure,
        ):

    class_paras = local_module.def_class_paras()

    dict_atom = class_structure.dict_atom

    list2d_angle = class_paras.list2d_angle
    float_tm_scaling = class_paras.float_tm_scaling

    str_alignfile = class_paras.str_alignfile
    with open(class_paras.str_alignfile) as obj_jsonfile:
        float_align = json.load( fp=obj_jsonfile )['float_align']

    str_chdir = dict_atom[ class_paras.int_atomkey ][1]
    os.chdir( str_chdir )

    int_l=str_alignfile.find('.')
    int_r=str_alignfile.rfind('.')
    str_abname = str_alignfile[(int_l+1):int_r]
    def_print_paras( locals(), ['str_abname'] )

    #----------------------------------------------[alignscaling]
    array2d_xdata_align, array2d_ydata_scaling, array2d_tm_xdata_align, array2d_tm_ydata_scaling, array2d_tm_kb = def_atom_alignscaling( 
        float_align=float_align, 
        float_scaling= class_structure.float_scaling, 
        )
    #----------------------------------------------[
    list1d_yheader, array2d_ydata_alphabeta = def_alphabeta( 
        list2d_angle = list2d_angle, 
        array2d_ydata = array2d_ydata_scaling )

    list2d_header = [['E(eV)'], list1d_yheader]
    list3d_data = [array2d_xdata_align, array2d_ydata_alphabeta]
    str_outfile = 'xas.'+str_abname+'.csv'
    def_writedata( list2d_header=list2d_header, list3d_data=list3d_data, str_outfile=str_outfile)
    #----------------------------------------------[
    _, array2d_tm_ydata_alphabeta = def_alphabeta( 
        list2d_angle=list2d_angle, 
        array2d_ydata=array2d_tm_ydata_scaling )

    list2d_header = [ ['Kpoint', 'Band'] ] + list2d_header
    list3d_data = [ array2d_tm_kb, array2d_tm_xdata_align, array2d_tm_ydata_alphabeta * float_tm_scaling]
    str_outfile = 'xas_tm.'+str_abname+'.csv'
    def_writedata( list2d_header=list2d_header, list3d_data=list3d_data, str_outfile=str_outfile)

    def_endfunc()

def def_atom_alignscaling(float_align, float_scaling):
    def_startfunc( locals() )

    class_paras = local_module.def_class_paras()
    #----------------------------------------------[extract]
    _, array2d_xdata = def_extract(
        str_datfile = class_paras.str_xasfile,
        list1d_column = [0]
        ) 
    _, array2d_ydata = def_extract(
        str_datfile = class_paras.str_xasfile,
        list1d_column = [1,2,3]
        )
    #----------------------------------------------[alignscaling]
    float_sft = float_align + def_vasp_finalenergy()
    array2d_xdata_align = array2d_xdata + float_sft
    array2d_ydata_scaling = array2d_ydata * float_scaling
    #----------------------------------------------[tm]
    array2d_tm_xdata, array2d_tm_ydata, array2d_tm_kb = def_tm_extract()
    array2d_tm_xdata_align = array2d_tm_xdata + float_sft
    array2d_tm_ydata_scaling = array2d_tm_ydata * float_scaling

    def_endfunc()
    return array2d_xdata_align, array2d_ydata_scaling, array2d_tm_xdata_align, array2d_tm_ydata_scaling, array2d_tm_kb

def def_tm_extract( str_datfile='MYCARXAS' ):
#------------------------------[]
#------------------------------[]
    def_startfunc( locals() )

    df_tm = pandas.read_csv( 
        str_datfile, 
        sep=' ', 
        skipinitialspace=True, 
        names= ['E(eV)','x','y','z','band'],
        usecols=[0,1,2,3,7],
        comment='#',
        )
    array1d_tm_xdata = df_tm['E(eV)'].to_numpy()
    array2d_tm_ydata = df_tm[['x','y','z']].to_numpy()
    array2d_tm_ydata *= def_vasp_volume()
    array1d_tm_band = df_tm['band'].to_numpy()

    int_lenline = numpy.shape( array1d_tm_band )[0]
    int_nb = numpy.amax( array1d_tm_band )
    int_nk = int_lenline//int_nb
    def_print_paras( locals(), ['int_lenline','int_nb','int_nk'] )
    array1d_tm_kpoint = numpy.empty( shape=(int_lenline), dtype= int )
    for int_k in range(int_nk):
        array1d_tm_kpoint[ int_nb*int_k:int_nb*(int_k+1) ].fill( int_k+1 )
    df_tm['kpoint'] = array1d_tm_kpoint

    array2d_tm_kb = df_tm[['kpoint','band']].to_numpy()
    def_endfunc()
    return array1d_tm_xdata, array2d_tm_ydata, array2d_tm_kb

