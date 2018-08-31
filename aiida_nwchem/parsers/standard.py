# -*- coding: utf-8 -*-
###########################################################################
# Copyright (c), The AiiDA team. All rights reserved.                     #
# This file is part of the AiiDA code.                                    #
#                                                                         #
# The code is hosted on GitHub at https://github.com/aiidateam/aiida_core #
# For further information on the license, see the LICENSE.txt file        #
# For further information please visit http://www.aiida.net               #
###########################################################################
from aiida_nwchem.parsers import BasenwcParser
from aiida_nwchem.calculations.standard import StandardCalculation
from aiida.orm.data.parameter import ParameterData

class StandardParser(BasenwcParser):
    """
    Parser for the output of nwchem.

    The goal for the standard parser is to parse all standard
    NWChem modules.

    Currently supported modules:
    - SCF
    - DFT
    - Geo-opt

    """
    def __init__(self,calc):
        """
        Initialize the instance of StandardParser
        """
        # check for valid input
        self._check_calc_compatibility(calc)
        super(StandardParser, self).__init__(calc)

    def _check_calc_compatibility(self,calc):
        from aiida.common.exceptions import ParsingError
        if not isinstance(calc,StandardCalculation):
            raise ParsingError("Input calc must be a StandardCalculation")

    def _get_output_nodes(self, output_path, error_path):
        """
        Extracts output nodes from the standard output and standard error
        files.
        """
        from aiida.orm.data.array.trajectory import TrajectoryData
        import re

        state = None
        step = None
        scale = None
        with open(output_path) as f:
            lines = [x.strip('\n') for x in f.readlines()]

        result_dict = dict()
        trajectory = None
        for line in lines:
            
            # 1. Determine which NWChem module we are using
            # SCF (HF)
            if state is None and re.match('^\s*NWChem SCF Module\s*$',line):
                state = 'nwchem-scf-module'
                continue
            # DFT
            if state is None and re.match('^\s*NWChem DFT Module\s*$',line):
                state = 'nwchem-dft-module'
                continue
            # Geometry Optimisation
            if state is None and re.match('^\s*NWChem Geometry Optimization\s*$',line):
                state = 'nwchem-geometry-optimisation'
                trajectory = TrajectoryData()
                continue
            
            # 2. Get the data from each module
            # SCF (HF)
            if state == 'nwchem-scf-module' and re.match('^\s*Final RHF \s*results\s*$',line):
                state = 'final-rhf-results'
                continue
            if re.match('^\s*\-*\s*$',line):
                continue
            if state == 'final-rhf-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    result_dict[key] = result.group(2)
                else:
                    state = 'nwchem-scf-module'
            # DFT 
            # Note the search for the Total DFT energy. NWChem doesn't otherwise
            # announce that the results are being printed. 
            if state == 'nwchem-dft-module' and re.match('^\s*Total DFT energy', line):
                state = 'final-dft-results'
            if state == 'final-dft-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    result_dict[key] = result.group(2)
                else: 
                    state = 'nwchem-dft-module'
            # Geometry Optimisation
            if state == 'nwchem-geometry-optimisation' and re.match('^\s*Step\s+\d+\s*$',line):
                result = re.match('^\s*Step\s+(\d+)\s*$',line)
                step = result.group(1)
                continue
            if state == 'nwchem-geometry-optimisation' and \
                re.match('^\s*Output coordinates in a.u.',line):
                state = 'nwchem-geometry-optimisation-coordinates'
                result = re.match('scale by \s(*[\-\d\.]+)',line)
                scale = result.group(1)
                continue
        return [('parameters', ParameterData(dict=result_dict))]
