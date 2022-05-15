import os
import subprocess
import re
import numpy
from tf_xas_kit import misc
from tf_xas_kit import pp

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

def def_vasp_outcar2xas(): 
    misc.def_startfunc( locals() ) 
 
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
    list1d_xheader, array2d_xdata = misc.def_extract( str_datfile=str_datfile, list1d_column=list1d_column ) 
 
    list1d_column = [1,2,3] 
    list1d_yheader, array2d_ydata = misc.def_extract( str_datfile=str_datfile, list1d_column=list1d_column ) 
 
    int_shapexdata = numpy.shape(array2d_xdata)[0] 
    misc.def_print_paras( locals(),['int_shapexdata'] ) 
    for int_i in range(int_shapexdata): 
        array2d_ydata[int_i,:] *= array2d_xdata[int_i][0] 
 
    array2d_ydata *= def_vasp_volume() 
 
    os.remove( str_tempfile ) 
 
    #def_writedata( 
    #    list2d_header = [ list1d_xheader, list1d_yheader ], 
    #    list3d_data = [ array2d_xdata, array2d_ydata ], 
    #    str_outfile = 'xas.csv' 
    #) 
 
    misc.def_endfunc() 
    return list1d_xheader, list1d_yheader, array2d_xdata, array2d_ydata 

def def_vasp_volume():
    misc.def_startfunc( locals() )

    float_volume = ase.io.read( filename='OUTCAR' ).get_volume()
    misc.def_print_paras( locals(), ['float_volume'] )

    misc.def_endfunc()
    return float_volume

def def_vasp_tm2xas(
        str_broadmethod,
        float_hwhm,
        int_broadnbin
        ):
    array1d_tm_xdata, array2d_tm_ydata, _ = misc.def_tm_extract()

    array2d_tm_ydata *= array1d_tm_xdata[ :,None ]

    array1d_xdata, array2d_ydata =  pp.def_broad(
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


def def_vasp_sub( class_structure ):
    os.chdir('template/')
    with open( 'vasp_sub.py','r' ) as file_sub:
        str_subfile = file_sub.read()
    os.chdir('..')

    for int_atomi, list1d_atomi in class_structure.dict_atom.items():
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
        print( list1d_atomname[ int_i ], 'int_atomnum=', int_atomnum, 'int_val=', int_val )
        int_nelect += int_val * int_atomnum
    print('int_nelect=', int_nelect)
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

