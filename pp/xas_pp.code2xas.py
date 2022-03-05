#!/bin/env python
import xas_module
import os
import local_module

class_paras = local_module.def_class_paras()

xas_module.def_code2xas(
    str_code = 'vasp',
    str_outfile = class_paras.str_xasfile,
    log_tm2xas = class_paras.log_tm2xas,
    str_broadmethod = class_paras.str_broadmethod,
    float_hwhm = class_paras.float_hwhm,
    int_broadnbin = class_paras.int_broadnbin,
    )
 
