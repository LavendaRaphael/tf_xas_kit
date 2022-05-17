#!/bin/env python
import os
from tf_xas_kit import pp

class_paras = pp.class_paras_code2xas.def_class_paras()

pp.def_code2xas(
    str_code = 'vasp',
    str_outfile = class_paras.str_xasfile,
    log_tm2xas = class_paras.log_tm2xas,
    str_broadmethod = class_paras.str_broadmethod,
    float_hwhm = class_paras.float_hwhm,
    int_broadnbin = class_paras.int_broadnbin,
    )
 
