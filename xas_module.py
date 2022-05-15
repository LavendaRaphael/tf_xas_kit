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

def def_abname( alpha, beta ):
    str_abname = 'a'+str(alpha)+'_b'+str(beta)
    return str_abname

