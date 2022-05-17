class class_structure():

    @property
    def dict_input(self):
        return self._dict_input
    @dict_input.setter
    def dict_input(self, temp):
        self._dict_input = temp

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
        for float_i in dict_temp.values(): float_sum += float_i
        self._dict_atom = {}
        for int_i in dict_temp: self._dict_atom[int_i] = [ dict_temp[int_i]/float_sum, 'atom_'+str(int_i) ]

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

