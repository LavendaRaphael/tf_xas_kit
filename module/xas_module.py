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

def def_vasp_tm2xas(
        str_broadmethod,
        float_hwhm,
        int_broadnbin
        ):
    array1d_tm_xdata, array2d_tm_ydata, _ = def_tm_extract()

    array2d_tm_ydata *= array1d_tm_xdata[ :,None ]

    array1d_xdata, array2d_ydata =  def_broad( 
        array1d_xdata = array1d_tm_xdata, 
        array2d_ydata = array2d_tm_ydata, 
        float_hwhm = float_hwhm,
        str_method = str_broadmethod,
        int_nbin = int_broadnbin
        )

    array2d_ydata /= array1d_xdata[:,None]

    list1d_xheader = ['E(eV)']
    list1d_yheader = ['x','y','z']

    return list1d_xheader, list1d_yheader, array1d_xdata, array2d_ydata
def def_broad(
        array1d_xdata,
        array2d_ydata,
        float_hwhm,
        str_method,
        int_nbin,
        ):
    float_xl = numpy.amin( array1d_xdata )
    float_xr = numpy.amax( array1d_xdata )
    array1d_xdata_new = numpy.linspace( 
        start = float_xl, 
        stop = float_xr, 
        num = int_nbin
        )

    array2d_delta =  array1d_xdata_new[:,None] - array1d_xdata[None,:]
    array2d_delta = def_lineshape(
        str_method = str_method,
        float_x = array2d_delta,
        float_hwhm = float_hwhm,
        )

    int_xdata_new_shape = array1d_xdata_new.shape[0]
    int_ydata_shape1 = array2d_ydata.shape[1]
    array2d_ydata_new = numpy.empty( shape=( int_xdata_new_shape, int_ydata_shape1 ) )
    for int_z in range( int_ydata_shape1 ):
        array1d_temp = array2d_delta * array2d_ydata[ :,int_z ][None,]
        array2d_ydata_new[ :,int_z ] = numpy.sum( array1d_temp, axis=1 )

    return array1d_xdata_new, array2d_ydata_new

def def_lineshape(
        str_method,
        float_x,
        float_hwhm,
        ):
    if ( str_method == 'gaussian' ):
        float_y = def_gaussian(
            float_x = float_x,
            float_hwhm = float_hwhm,
            )
    elif ( str_method == 'lorentzian' ):
        float_y = def_lorentzian(
            float_x = float_x,
            float_hwhm = float_hwhm,
            )
    else:
        raise ValueError()

    return float_y

def def_gaussian(
        float_x,
        float_hwhm,
        ):
    return numpy.sqrt(numpy.log(2) / numpy.pi) / float_hwhm\
                             * numpy.exp(-(float_x / float_hwhm)**2 * numpy.log(2))
def def_lorentzian(
        float_x,
        float_hwhm,
        ):
    return float_hwhm / numpy.pi / (float_x**2 + float_hwhm**2)

def def_exp_scaling( 
        class_structure,
        str_outfile 
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
    dict_scaling = def_findscaling_dict(
        array1d_xdata = df_xas['E(eV)'].to_numpy(),
        array1d_ydata = df_xas['Intensity'].to_numpy(),
        class_structure = class_structure
        )
    class_paras = local_module.def_class_paras()
    float_scaling = dict_scaling[ class_paras.str_scalingmethod ]
    df_xas['Intensity'] *= float_scaling
    
    df_xas.to_csv( str_outfile )

def def_code2xas(
        str_code = 'vasp',
        str_outfile = 'xas.csv',
        log_tm2xas = False,
        str_broadmethod = None,
        float_hwhm = None,
        int_broadnbin = None,
        ):
    if ( str_code == 'vasp' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = def_vasp2xas(
            log_tm2xas = log_tm2xas,
            str_broadmethod = str_broadmethod,
            float_hwhm = float_hwhm,
            int_broadnbin = int_broadnbin,
            )
    elif ( str_code == 'feff' ):
        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = def_feff2xas()
    else:
        raise

    def_writedata(
        list2d_header = [list1d_xheader, list1d_yheader],
        list3d_data = [array2d_xdata, array2d_ydata], 
        str_outfile = str_outfile)

    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

def def_vasp2xas(
        log_tm2xas = False,
        str_broadmethod = None,
        float_hwhm = None,
        int_broadnbin = None,
        ):
    if ( not log_tm2xas ):
        return def_vasp_outcar2xas()
    return def_vasp_tm2xas( 
        str_broadmethod = str_broadmethod,
        float_hwhm = float_hwhm,
        int_broadnbin = int_broadnbin
        )

def def_feff2xas():

    str_xheader = 'E(eV)'
    list1d_yheader = [ 'x','y','z' ]

    df_xas_x = pandas.read_csv( 
        'ellip_z/xmu.dat', 
        sep=' ', 
        skipinitialspace=True, 
        header=None, 
        names= [str_xheader, list1d_yheader[0]] ,
        usecols=[0,3],
        comment='#',
        )
    print(df_xas_x)

    df_xas_y = df_xas_x.rename( columns={ list1d_yheader[0]: list1d_yheader[1] } )

    df_xas_z = pandas.read_csv(       
        'polar_z/xmu.dat',
        sep=' ',
        skipinitialspace=True,
        header=None,
        names=[str_xheader, list1d_yheader[2]],
        usecols=[0,3],
        comment='#',
        )
    print(df_xas_z)

    df_xas = pandas.merge( df_xas_x, df_xas_y, on=str_xheader )
    df_xas = pandas.merge( df_xas, df_xas_z, on=str_xheader )
    print(df_xas)
    
    # df_xas.to_csv( str_outfile )

    list1d_xheader = [str_xheader]
    array2d_xdata = df_xas[list1d_xheader].to_numpy()
    array2d_ydata = df_xas[list1d_yheader].to_numpy()

    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

def def_vasp_sub( class_structure ):
    os.chdir('template/')
    with open( 'vasp_sub.py','r' ) as file_sub:
        str_subfile = file_sub.read()
    os.chdir('..')

    for list1d_atomi in class_structure.dict_atom.values():
        str_atomdir = list1d_atomi[1]
        os.chdir(str_atomdir)

        with open( 'vasp_sub.py','w' ) as file_sub:
            str_subfile_new = re.sub( 'xNUMx',str(int_atomi),str_subfile )
            file_sub.write( str_subfile_new )

        subprocess.run( ['python','vasp_sub.py'] )
        os.chdir('..')

def def_vasp_jobinit( class_structure ):
    os.chdir('template/')

    with open('POSCAR', 'r') as file_poscar:
        list1d_poscar = file_poscar.readlines()
        list1d_atomnum = [int(str_num) for str_num in list1d_poscar[6].split()]
        if (list1d_atomnum[0] != 1):
            list1d_poscar[5] = list1d_poscar[5].split()[0] +' '+ list1d_poscar[5]
            list1d_atomnum[0] -= 1
            list1d_poscar[6] = '1'
            for int_i in list1d_atomnum:
                list1d_poscar[6] += ' '+str(int_i)
            list1d_poscar[6] += '\n'

        if ( list1d_poscar[7].strip()[0] in 'CcKkDd' ):
            int_directline = 7
        elif ( list1d_poscar[8].strip()[0] in 'CcKkDd' ):
            int_directline = 8
        else:
            raise NameError ( 'POSCAR not format!' )
 
    int_nelect = 0
    list1d_atomname = list1d_poscar[5].split()
    list1d_atomnum = [int(str_num) for str_num in list1d_poscar[6].split()]
    dict_pot = group_module.def_dict_pot()
    for int_i in range(len(list1d_atomname)):
        int_val = dict_pot[ list1d_atomname[ int_i ] ][1]
        int_atomnum = list1d_atomnum[ int_i ]
        int_nelect += int_val * int_atomnum
    with open('INCAR','r+') as file_incar:
        list1d_incar = file_incar.readlines()
        for int_i in range(len(list1d_incar)):
            if (list1d_incar[ int_i ].strip()[0:6] == 'NBANDS'):
                list1d_incar[ int_i ] = 'NBANDS = '+str(int_nelect)
                break
        file_incar.seek(0)
        file_incar.truncate()
        file_incar.writelines( list1d_incar )

    os.chdir('..')

    for int_atomi, list1d_atomi in class_structure.dict_atom.items():
        str_atomdir = list1d_atomi[1]
        copy_tree( 'template/', str_atomdir )
        os.chdir(str_atomdir)
        list1d_poscar_atomi = list( list1d_poscar)
        list1d_poscar_atomi[ int_directline +1 ] = list1d_poscar[ int_directline +int_atomi ]
        list1d_poscar_atomi[ int_directline +int_atomi ] = list1d_poscar[ int_directline +1 ]
        with open('POSCAR', 'w') as file_poscar:
            file_poscar.writelines( list1d_poscar_atomi )

        os.chdir('..')

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

def def_exp_xyzfit( list2d_alpha, str_outfile ):
#----------------------------------------------[]
# list2d_alpha = []
#----------------------------------------------[]
    def_startfunc( locals() )

    int_lenalpha = len(list2d_alpha)
    array1d_alpha = numpy.empty( shape=(int_lenalpha) )
    list2d_data = []
    for int_i in range(int_lenalpha):
        array1d_alpha[int_i] = list2d_alpha[ int_i ][ 0 ]
        str_datfile = list2d_alpha[ int_i ][ 1 ]
        list1d_xycolumn = list2d_alpha[ int_i ][ 2 ]
        list1d_xheader, array2d_xdata = def_extract( str_datfile=str_datfile, list1d_column=[ list1d_xycolumn[0] ] )
        list1d_yheader, array2d_ydata = def_extract( str_datfile=str_datfile, list1d_column=[ list1d_xycolumn[1] ] )
        list2d_data.append( [ array2d_xdata, array2d_ydata] )

    array1d_cosalpha2 = numpy.cos( numpy.radians( array1d_alpha ) ) **2
    def_print_paras( locals(), ['array1d_cosalpha2'] )

    array1d_xdata_interp, list1d_ydata_interp = def_interp( list2d_data )

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

    def_writedata(
        list2d_header = [ ['E(eV)', 'fit_xy','fit_z'] ],
        list3d_data = [ array1d_xdata_interp, array2d_ydata_fit ],
        str_outfile=str_outfile
        )

class class_paras(object):

    @property
    def float_tm_scaling(self):
        return self._float_tm_scaling
    @float_tm_scaling.setter
    def float_tm_scaling(self, obj_temp):
        self._float_tm_scaling = obj_temp

    @property
    def str_alignfile(self):
        return self._str_alignfile
    @str_alignfile.setter
    def str_alignfile(self, obj_temp):
        self._str_alignfile = obj_temp

    @property
    def int_atomkey(self):
        return self._int_atomkey
    @int_atomkey.setter
    def int_atomkey(self, obj_temp):
        self._int_atomkey = obj_temp

    @property
    def str_scalingmethod(self):
        return self._str_scalingmethod
    @str_scalingmethod.setter
    def str_scalingmethod(self, str_temp):
        self._str_scalingmethod = str_temp

    @property
    def log_tm2xas(self):
        return self._log_tm2xas
    @log_tm2xas.setter
    def log_tm2xas(self, obj_temp):
        self._log_tm2xas = obj_temp

    @property
    def str_broadmethod(self):
        return self._str_broadmethod
    @str_broadmethod.setter
    def str_broadmethod(self, obj_temp):
        self._str_broadmethod = obj_temp

    @property
    def float_hwhm(self):
        return self._float_hwhm
    @float_hwhm.setter
    def float_hwhm(self, obj_temp):
        self._float_hwhm = obj_temp

    @property
    def str_xasfile(self):
        return self._str_xasfile
    @str_xasfile.setter
    def str_xasfile(self, obj_temp):
        self._str_xasfile = obj_temp

    @property
    def str_avefile(self):
        return self._str_avefile
    @str_avefile.setter
    def str_avefile(self, obj_temp):
        self._str_avefile = obj_temp

    @property
    def str_alphafile(self):
        return self._str_alphafile
    @str_alphafile.setter
    def str_alphafile(self, obj_temp):
        self._str_alphafile = obj_temp

    @property
    def int_broadnbin(self):
        return self._int_broadnbin
    @int_broadnbin.setter
    def int_broadnbin(self, obj_temp):
        self._int_broadnbin = obj_temp

    @property
    def list2d_anlge(self):
        return self._list2d_anlge
    @list2d_anlge.setter
    def list2d_anlge(self, obj_temp):
        self._list2d_anlge = obj_temp

    @property
    def list1d_alignangle(self):
        return self._list1d_alignangle
    @list1d_alignangle.setter
    def list1d_alignangle(self, obj_temp):
        self._list1d_alignangle = obj_temp

class class_structure(object):

    @property
    def list1d_column(self):
        return self._list1d_column
    @list1d_column.setter
    def list1d_column(self, list1d_temp):
        self._list1d_column = list1d_temp

    @property
    def str_datfile(self):
        return self._str_datfile
    @str_datfile.setter
    def str_datfile(self, str_temp):
        self._str_datfile = str_temp

    @property
    def tuple_mainxrange(self):
        return self._tuple_mainxrange
    @tuple_mainxrange.setter
    def tuple_mainxrange(self, tuple_temp):
        self._tuple_mainxrange = tuple_temp

    @property
    def tuple_postxrange(self):
        return self._tuple_postxrange
    @tuple_postxrange.setter
    def tuple_postxrange(self, tuple_temp):
        self._tuple_postxrange = tuple_temp
 
    @property
    def float_onset(self):
        return self._float_onset
    @float_onset.setter
    def float_onset(self, float_temp):
        self._float_onset = float_temp

    @property
    def str_chdir(self):
        return self._str_chdir
    @str_chdir.setter
    def str_chdir(self, str_temp):
        self._str_chdir = str_temp
  
    @property
    # dict_atom[1] = [
    #   [ 1.0, 'atom_1' ]
    #   ]
    def dict_atom(self):
        return self._dict_atom
    @dict_atom.setter
    def dict_atom(self, dict_temp):
        float_sum = 0
        for list_i in dict_temp.values(): float_sum += list_i[0]
        for list_i in dict_temp.values(): list_i[0] /= float_sum
        for int_i in dict_temp: dict_temp[ int_i ].append( 'atom_'+str(int_i) )
        self._dict_atom = dict_temp

    @property
    def list1d_bbox(self):
        return self._list1d_bbox
    @list1d_bbox.setter
    def list1d_bbox(self, list1d_temp):
        self._list1d_bbox = list1d_temp

    @property
    def str_cif(self):
        return self._str_cif
    @str_cif.setter
    def str_cif(self, str_temp):
        self._str_cif = str_temp

    @property
    def str_code(self):
        return self._str_code
    @str_code.setter
    def str_code(self, str_temp):
        self._str_code = str_temp

    @property
    def float_scaling(self):
        return self._float_scaling
    @float_scaling.setter
    def float_scaling(self, float_temp):
        self._float_scaling = float_temp

def def_ave( 
        class_structure, 
        ):
    def_startfunc( locals(), [ 'class_structure' ] )

    dict_atom = class_structure.dict_atom
    str_code = class_structure.str_code
    class_paras = local_module.def_class_paras()

    list2d_data = []
    for list1d_atom in dict_atom.values():

        str_chdir = list1d_atom[1]
        float_scaling = list1d_atom[0]
        os.chdir(str_chdir)
        print(os.getcwd())

        list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata = def_code2xas( 
            str_code = str_code,
            str_outfile = class_paras.str_xasfile,
            log_tm2xas = class_paras.log_tm2xas,
            str_broadmethod = class_paras.str_broadmethod,
            float_hwhm = class_paras.float_hwhm,
            int_broadnbin = class_paras.int_broadnbin,
            )
    
        array2d_xdata_sft = def_sft( array2d_xdata, str_code )

        list_ycolums = list(range(len(array2d_ydata[0])))
        list2d_data.append( [ array2d_xdata_sft, array2d_ydata, list_ycolums, float_scaling ] )

        os.chdir('..')

    array1d_xdata_mix, array2d_ydata_mix = def_mix(list2d_data=list2d_data)
    def_writedata( 
        list2d_header = [ list1d_xheader, list1d_yheader ],
        list3d_data = [ array1d_xdata_mix, array2d_ydata_mix ],
        str_outfile = class_paras.str_avefile
        )

def def_alphabeta_workflow( 
        class_structure,
        ):
#----------------------------------------------[]
# list_alignangle = [ alpha0, beta0 ]
# list_scalingangle = [ alpha1, beta1 ]
#----------------------------------------------[]
    def_startfunc( locals(), ['class_structure'] )

    class_paras = local_module.def_class_paras()
    str_datfile = class_paras.str_avefile
    str_outfile = class_paras.str_alphafile
    list1d_alignangle = class_paras.list1d_alignangle
    list2d_angle = class_paras.list2d_angle
    #--------------------------------------------------[extract]
    list1d_column = [0]
    list1d_xheader, array2d_xdata_origin = def_extract( 
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )
    
    list1d_column = [1,2,3]
    _, array2d_ydata_origin = def_extract( 
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    #--------------------------------------------------[alignscaling]

    array1d_xdata_align = def_align( 
        array1d_xdata = array2d_xdata_origin, 
        array2d_ydata = array2d_ydata_origin, 
        list1d_alignangle = list1d_alignangle, 
        float_onset = class_structure.float_onset, 
        )

    array2d_ydata_scaling = def_scaling(
        array2d_ydata = array2d_ydata_origin,
        class_structure = class_structure,
        )
    
    list1d_yheader, array2d_ydata_alphabeta = def_alphabeta( 
        list2d_angle = list2d_angle, 
        array2d_ydata = array2d_ydata_scaling 
        )

    def_writedata( 
        list2d_header = [ list1d_xheader, list1d_yheader ],
        list3d_data = [ array1d_xdata_align, array2d_ydata_alphabeta ],
        str_outfile=str_outfile
        )

    def_endfunc()
    return

def def_findscaling_dict(
        array1d_xdata,
        array1d_ydata,
        class_structure
        ): 
    float_mainscaling = def_findscaling(
        array1d_xdata = array1d_xdata, 
        array1d_ydata = array1d_ydata, 
        tuple_xrange = class_structure.tuple_mainxrange
        )

    float_postscaling = def_findscaling(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = class_structure.tuple_postxrange
        )
    dict_scaling = {
        'float_mainscaling' : float_mainscaling,
        'float_postscaling': float_postscaling
        }

    def_print_paras( locals(),['dict_scaling'] ) 
    return dict_scaling

def def_findscaling(
        array1d_xdata,
        array1d_ydata,
        tuple_xrange
        ):

    if ( numpy.amax( array1d_xdata ) < tuple_xrange[1] ):
        raise ValueError('tuple_xrange is too large!')
    float_area = def_findarea(
        array1d_xdata = array1d_xdata,
        array1d_ydata = array1d_ydata,
        tuple_xrange = tuple_xrange
        )
    float_scaling = ( tuple_xrange[1] - tuple_xrange[0] ) / float_area

    return float_scaling

def def_scaling( 
        array2d_ydata,
        class_structure,
        ):
    
    float_scaling = class_structure.float_scaling

    array2d_ydata_scaling = array2d_ydata * float_scaling

    return array2d_ydata_scaling

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
    def_startfunc( locals(), ['array1d_xdata', 'array2d_ydata'] )
    
    array1d_xdata_origin = array1d_xdata
    array2d_ydata_origin = array2d_ydata
    #--------------------------------------------------[align]
    _, array2d_ydata_alphabeta = def_alphabeta( 
        list2d_angle = [ list1d_alignangle ],
        array2d_ydata = array2d_ydata_origin 
        )
    
    dict_peaks = def_findpeaks( 
        array1d_xdata = array1d_xdata_origin, 
        array1d_ydata = array2d_ydata_alphabeta)
    float_align = float_onset - dict_peaks[ 'E(eV)' ][0]
    array1d_xdata_align = array1d_xdata_origin + float_align

    str_abname = def_abname( alpha=list1d_alignangle[0], beta=list1d_alignangle[1])
    str_jsonfile = 'xas.'+str_abname+'.align.json'
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump( obj={'float_align': float_align}, fp=obj_jsonfile, indent=4 )

    def_endfunc()
    return array1d_xdata_align

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

def def_scaling_json( 
        list1d_alignangle, 
        list1d_scalingangle, 
        class_structure, 
        str_datfile='xas.ave.csv'):
#----------------------------------------------[]
# list_alignangle = [ alpha0, beta0 ]
# list_scalingangle = [ alpha1, beta1 ]
#----------------------------------------------[]
    def_startfunc( locals(), ['class_structure'] )
    #--------------------------------------------------[extract]
    list1d_column = [0]
    _, array2d_xdata_origin = def_extract( 
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )
    
    list1d_column = [1,2,3]
    _, array2d_ydata_origin = def_extract( 
        str_datfile=str_datfile,
        list1d_column = list1d_column,
        )

    #--------------------------------------------------[align]
    _, array2d_ydata_alphabeta = def_alphabeta( 
        list2d_angle = [ list1d_alignangle ],
        array2d_ydata = array2d_ydata_origin 
        )
    
    dict_peaks = def_findpeaks( 
        array1d_xdata = array2d_xdata_origin, 
        array1d_ydata = array2d_ydata_alphabeta)
    float_align = class_structure.float_onset - dict_peaks[ 'E(eV)' ][0]
    array1d_xdata_align = array2d_xdata_origin + float_align
    #--------------------------------------------------[scaling]
    _, array2d_ydata_alphabeta = def_alphabeta(
        list2d_angle = [ list1d_scalingangle ], 
        array2d_ydata = array2d_ydata_origin 
        )
    
    dict_scaling = def_findscaling_dict(
        array1d_xdata = array1d_xdata_align, 
        array1d_ydata = array2d_ydata_alphabeta,
        class_structure = class_structure,
        )

    #--------------------------------------------------[output]
    str_abname = def_abname( alpha=list1d_scalingangle[0], beta=list1d_scalingangle[1])
    str_jsonfile = 'xas.'+str_abname+'.scaling.json'
    with open( str_jsonfile, 'w' ) as obj_jsonfile:
        json.dump( 
            obj=dict_scaling, 
            fp=obj_jsonfile, 
            indent=4 )

    def_endfunc()
    return

def def_alphabeta( list2d_angle, array2d_ydata ):
#----------------------------------------------[]
# list2d_angle = []
# list2d_angle.append( [alpha1, beta1] )
# list2d_angle.append( [alpha2, beta2] )
# array2d_ydata:
#   x,y,z
#   ...
#----------------------------------------------[]
    def_startfunc( locals(), ['array2d_ydata'] )

    int_shape2dydata0 = numpy.shape(array2d_ydata)[0]
    def_print_paras( locals(), ['int_shape2dydata0'] )
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

def def_vasp_outcar2xas():
    def_startfunc( locals() )

    str_tempfile = 'outcar2xas.tmp'
    with open( 'OUTCAR', 'r' ) as obj_datfile:
        for str_line in obj_datfile:
            if (str_line.strip() == 'frequency dependent IMAGINARY DIELECTRIC FUNCTION (independent particle, no local field effects) density-density'):
                break
        with open( str_tempfile, mode='w' ) as obj_tmp: 
            str_line = next(obj_datfile)
            obj_tmp.write( str_line )
            str_line = next(obj_datfile)
            for str_line in obj_datfile:
                if (not str_line.strip()): break
                obj_tmp.write( str_line )

    list1d_column = [0]
    str_datfile = str_tempfile
    list1d_xheader, array2d_xdata = def_extract( str_datfile=str_datfile, list1d_column=list1d_column )

    list1d_column = [1,2,3]
    list1d_yheader, array2d_ydata = def_extract( str_datfile=str_datfile, list1d_column=list1d_column )

    int_shapexdata = numpy.shape(array2d_xdata)[0]
    def_print_paras( locals(),['int_shapexdata'] )
    for int_i in range(int_shapexdata):
        array2d_ydata[int_i,:] *= array2d_xdata[int_i][0]

    array2d_ydata *= def_vasp_volume()

    os.remove( str_tempfile )

    #def_writedata(
    #    list2d_header = [ list1d_xheader, list1d_yheader ],
    #    list3d_data = [ array2d_xdata, array2d_ydata ],
    #    str_outfile = 'xas.csv'
    #)

    def_endfunc()
    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata

def def_vasp_volume():
    def_startfunc( locals() )

    float_volume = ase.io.read( filename='OUTCAR' ).get_volume()
    def_print_paras( locals(), ['float_volume'] )

    def_endfunc()
    return float_volume

def def_vasp_finalenergy():
    def_startfunc( locals() )

    float_finalenergy = ase.io.read( filename='OUTCAR' ).get_total_energy()
    def_print_paras( locals(), ['float_finalenergy'] )

    def_endfunc()
    return float_finalenergy

def def_sft( array1d_xdata, str_code='vasp' ):

    if ( str_code == 'vasp' ):
        float_sft = def_sft_vasp()
    elif ( str_code == 'feff' ):
        float_sft = 0
    else:
        raise

    array1d_xdata_sft = array1d_xdata + float_sft
    return array1d_xdata_sft

def def_sft_vasp( ): 
    float_finalenergy = def_vasp_finalenergy()
    float_sft = float_finalenergy

    return float_sft

def def_writedata( list2d_header, list3d_data, str_outfile):
    def_startfunc( locals(), ['list3d_data'] )

    with open( str_outfile, 'w', newline='' ) as obj_outfile:
        obj_outwriter = csv.writer( obj_outfile, delimiter=',' )
        # header
        list1d_header = []
        for list1d_temp in list2d_header: 
            list1d_header.extend( list1d_temp )
        obj_outwriter.writerow( list1d_header )
        
        int_lenline = numpy.shape( list3d_data[0] )[0]
        def_print_paras( locals(),['int_lenline'] )
        for int_i in range(len(list3d_data)):
            if ( list3d_data[int_i].ndim == 2 ): continue
            list3d_data[int_i] = list3d_data[int_i].reshape(int_lenline,1)
        for int_i in range(int_lenline):
            list1d_data = []
            for array2d_temp in list3d_data:
                list1d_data.extend( array2d_temp[int_i]  ) 
            obj_outwriter.writerow( list1d_data )

    def_endfunc()
    return

def def_findarea( array1d_xdata, array1d_ydata, tuple_xrange):
    def_startfunc( locals(), ['array1d_xdata', 'array1d_ydata'] )

    array1d_xdata = numpy.reshape( a=array1d_xdata, newshape=-1 )
    array1d_ydata = numpy.reshape( a=array1d_ydata, newshape=-1 )
 
    array_x_new=[]
    array_y_new=[]

    for int_i in range(len(array1d_xdata)):
        float_x = array1d_xdata[int_i]
        float_y = array1d_ydata[int_i]
        if ( float_x > tuple_xrange[0] and float_x < tuple_xrange[1] ):
            array_x_new.append( float_x )
            array_y_new.append( float_y )

    float_area = numpy.trapz( y=array_y_new, x=array_x_new )
    print(json.dumps( obj={'float_area':float_area}, indent=4))

    def_endfunc()
    return float_area

def def_findpeaks( 
        array1d_xdata, 
        array1d_ydata, 
        float_relheight = 0.4, 
        float_relprominence = 0.02
        ):
#------------------------------[]
#------------------------------[]
    def_startfunc( locals(), ['array1d_xdata', 'array1d_ydata'] )

    array1d_xdata = numpy.reshape( a=array1d_xdata, newshape=-1 )
    array1d_ydata = numpy.reshape( a=array1d_ydata, newshape=-1 )

    float_y_max = max(array1d_ydata)

    array1d_peak_indices, dict_properties = scipy.signal.find_peaks( 
        array1d_ydata, 
        height = float_relheight * float_y_max, 
        prominence = float_relprominence * float_y_max
        )
    
    dict_peaks = {}
    dict_peaks[ 'E(eV)' ] = array1d_xdata[ array1d_peak_indices ]
    dict_peaks[ 'relheight' ] = dict_properties[ 'peak_heights' ] / float_y_max
    dict_peaks[ 'relprominence' ] = dict_properties[ 'prominences' ] / float_y_max
    def_print_paras( locals(),['dict_peaks'])

    def_endfunc()
    return dict_peaks

def def_mix(list2d_data):
#------------------------------[]
# list2d_data = []
# list2d_data.append( [ array1d_xdata, array2d_ydata, [0,2], 0.7 ] )
#------------------------------[]
    dict_args = copy.deepcopy(locals()) 
    for list_temp in dict_args['list2d_data']:
        del list_temp[0:2]
    def_startfunc( dict_args )

    list2d_xydata = []
    for list1d_data in list2d_data:
        list2d_xydata.append( [ list1d_data[0], list1d_data[1] ] )
    array1d_xdata_interp, list1d_ydata_interp = def_interp( list2d_xydata )

    int_len1dxdata = len( array1d_xdata_interp )
    int_lenycolumn = len( list2d_data[0][2] )
    array2d_ydata_mix = numpy.zeros( shape=(int_len1dxdata, int_lenycolumn) )

    int_len2ddata = len(list2d_data)
    for int_i in range(int_len2ddata):
        list1d_ycolumn = list2d_data[int_i][2]
        array2d_ydata = list1d_ydata_interp[ int_i][:, list1d_ycolumn ]
        float_scaling = list2d_data[int_i][3]
        array2d_ydata_mix += array2d_ydata * float_scaling

    def_endfunc()
    return array1d_xdata_interp, array2d_ydata_mix

def def_interp(list2d_data):
#------------------------------[]
# list2d_data = []
# list_data.append( [ array1d_xdata, array2d_ydata ] )
#------------------------------[]
    def_startfunc()

    for list1d_data in list2d_data:
        list1d_data[0] = list1d_data[0].reshape( -1 )

    float_xl = float('-inf')
    float_xr = float('+inf')
    for list1d_data in list2d_data:
        array1d_xdata = list1d_data[0]
        float_xl = max( float_xl, numpy.amin( array1d_xdata ))
        float_xr = min( float_xr, numpy.amax( array1d_xdata ))

    dict_json = {}
    dict_json[ 'float_xl' ] = float_xl
    dict_json[ 'float_xr' ] = float_xr
    print( json.dumps( dict_json, indent=4 ) )

    int_shape1dxdata = numpy.shape(list2d_data[0][0])[0]
    dict_json = {}
    dict_json[ 'int_shape1dxdata' ] = int_shape1dxdata
    print( json.dumps( dict_json, indent=4 ) )

    array1d_xdata_interp = numpy.linspace( start=float_xl, stop=float_xr, num=int_shape1dxdata )

    list1d_ydata_interp = []
    for list1d_data in list2d_data:
        array1d_xdata = list1d_data[0]
        array2d_ydata = list1d_data[1]
        int_shape2dydata_1 = numpy.shape( array2d_ydata )[1]
        array2d_temp = numpy.empty( shape=(int_shape1dxdata, int_shape2dydata_1) )
        for int_i in range(int_shape2dydata_1):
            array2d_temp[:,int_i] = numpy.interp( x=array1d_xdata_interp, xp=array1d_xdata, fp=array2d_ydata[:,int_i] )

        list1d_ydata_interp.append( array2d_temp )

    def_endfunc()
    return array1d_xdata_interp, list1d_ydata_interp

def def_extract( str_datfile, list1d_column, log_head=True, dtype=float ):
#------------------------------[]
#------------------------------[]
    def_startfunc( locals() )

    with open( str_datfile, 'r', newline='' ) as obj_datfile:
        obj_datreader = csv.reader( filter( lambda row: (row.strip() and (row.strip()[0]!='#')), obj_datfile ), delimiter= ' ', skipinitialspace=True )
        if log_head:
            list1d_line = next(obj_datreader)
        list1d_line = next(obj_datreader)
        if (',' in list1d_line[0]):
            delimiter=','
        else:
            delimiter=' '
    def_print_paras( locals(),['delimiter'] )

    list1d_header = []
    list2d_data = []
    with open( str_datfile, 'r', newline='' ) as obj_datfile:
        obj_datreader = csv.reader( filter( lambda row: (row.strip() and (row.strip()[0]!='#')), obj_datfile ), delimiter=delimiter, skipinitialspace=True )
        if log_head:
            list1d_line = next(obj_datreader)
            list1d_header = [ list1d_line[i] for i in list1d_column ]
        else:
            list1d_header = [ '' for i in list1d_column ]
        def_print_paras( locals(),['list1d_header'] )

        for list1d_line in obj_datreader:
            if ( not list1d_line[ list1d_column[0] ] ): continue
            list_temp = []
            for int_i in list1d_column:
                list_temp.append( list1d_line[int_i] )
            list2d_data.append( list_temp )
    array2d_data = numpy.array( list2d_data, dtype=dtype )
    tup_shapedata = numpy.shape(array2d_data)
    def_print_paras( locals(),['tup_shapedata'] )
 
    def_endfunc()
    return list1d_header, array2d_data

def def_abname( alpha, beta ):
    str_abname = 'a'+str(alpha)+'_b'+str(beta)
    return str_abname

def def_print_paras( dict_localpara, list_paraname ):
    dict_paraprint = {}
    for str_paraname in list_paraname:
        dict_paraprint[ str_paraname ] = dict_localpara[ str_paraname ]
    print(json.dumps( dict_paraprint, indent=4, cls=NumpyEncoder ))

class NumpyEncoder(json.JSONEncoder):
#----------------------------------------------[]
# https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
#----------------------------------------------[]
    """ Special json encoder for numpy types """
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif isinstance(obj, type):
            return str(obj)
        return json.JSONEncoder.default(self, obj)

def def_startfunc( dict_args={}, list1d_del=[] ):
    this_function_name = inspect.currentframe().f_back.f_code.co_name
    print("#"+'-'*20+"["+this_function_name+"]\n")
    for str_del in list1d_del:
        del dict_args[ str_del ]
    print(json.dumps( obj=dict_args, indent=4, cls=NumpyEncoder ))

def def_endfunc():
    print("#"+'-'*20+"<<\n")
