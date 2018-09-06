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
from aiida.orm.data.structure import StructureData
import numpy as np
import ase 
import re

class StandardParser(BasenwcParser):
    """
    Parser for the output of nwchem.

    The goal for the standard parser is to parse all standard
    NWChem modules.

    Currently supported modules:
    - SCF
    - DFT
    - Geo-opt

    Multiple tasks are possible so we must parse each one.
    The output will be parameter data containing a list 
    of dictionaries in the order of the tasks. 
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

        state = None
        step = None
        scale = None
        with open(output_path) as f:
            lines = [x.strip('\n') for x in f.readlines()]
        
        # All of the parsed data 
        data = {}
        # A list to hold the dictionaries for each calculation task
        result_list = []
        # Dictionary to hold the data parsed from each task
        task_dict = {}
    
        # For geo opt, lists to hold atom symbols, positions and the cell
        symbols = []
        positions = []
        cell = []
    
        # List to hold all output nodes 
        output_nodes_list = []
        
        for line in lines:
            
            # 1. Determine which NWChem module we are using
            # SCF (HF)
            if re.match('^\s*NWChem SCF Module\s*$',line):
                if state is None:
                    state = 'nwchem-scf-module'
                    task_dict = {'module':'scf'}
                    continue
                elif state == 'nwchem-geoopt':
                    if not tast_dict:
                        task_dict = {'module':'geoopt_scf'}
                    continue
            # DFT
            if re.match('^\s*NWChem DFT Module\s*$',line):
                if state is None:
                    state = 'nwchem-dft-module'
                    task_dict = {'module':'dft'}
                    continue
                elif state == 'nwchem-geoopt':
                    if not task_dict:
                        task_dict = {'module':'geoopt_dft'}
                    continue
    
            # Geometry Optimisation
            if state is None and re.match('^\s*NWChem Geometry Optimization\s*$',line):
                state = 'nwchem-geoopt'
                task_dict = {}
                continue
    
            
            # 2. Get the data from each module
            # SCF (HF)
            if state == 'nwchem-scf-module' and re.match('^\s*Final [ROU]+HF\s*results\s*$',line):
                state = 'final-hf-results'
                print 'yes'
                continue
            if re.match('^\s*\-*\s*$',line):
                continue
            if state == 'final-hf-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict[key] = result.group(2)
                    
            # DFT 
            # Note the search for the Total DFT energy. NWChem doesn't otherwise
            # announce that the results are being printed. 
            if state == 'nwchem-dft-module' and re.match('^\s*Total DFT energy', line):
                state = 'final-dft-results'
            if state == 'final-dft-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict[key] = result.group(2)
    
            # Geometry Optimisation
            if state == 'nwchem-geoopt' and re.match('^\s*Optimization converged\s*$',line):
                state = 'final-geoopt-results'
                continue
            # - Parse step and energy
            if state == 'final-geoopt-results':
                result = re.match('^@\s*([\d]+)\s*([\-\d\.]+)',line)
                if result:
                    task_dict['final_step'] = result.group(1)
                    # Note, this energy may be unreliable - do another scf on the final structure
                    task_dict['final_opt_energy'] = result.group(2)
            # - Parse coords
            if state == 'final-geoopt-results' and re.match('^\s*Output coordinates in angstroms\s*',line):
                state = 'final-geoopt-coords'
                continue
            if state == 'final-geoopt-coords' and re.match('^[A-z\s\.-]*$', line):
                continue
            elif state == 'final-geoopt-coords':
                result = re.match('^\s*[\d]+\s*([a-zA-Z]+)\s*[\-\d]+\.[\d]+'\
                                  '\s*([\-\d\.]+)\s*([\-\d\.]+)\s*([\-\d\.]+)', line)
                if result:
                    symbols.append(result.group(1))
                    positions.append([result.group(2),result.group(3),result.group(4)])
                else:
                    state = 'final-geoopt-results'
            # - Parse cell 
            if state == 'final-geoopt-results' and re.match('^\s*lattice vectors in angstroms',line):
                state = 'nwchem-geoopt-cell'
                continue
            if state == 'nwchem-geoopt-cell':
                result = re.match('^\s*a[1-3]=<\s*([\d\.\d]+)\s*([\d\.\d]+)\s*([\d\.\d]+)', line) 
                if result:
                    cell.append([result.group(1),result.group(2),result.group(3)])
                else:
                    state = None
            
            # If see timings, we are exiting the module  
            if re.match('^ Task  times  cpu:\s*([\d\.\d]+)s\s*wall:\s*([\d\.\d]+)s', line):
                result = re.match('^ Task  times  cpu:\s*([\d\.\d]+)s\s*wall:\s*([\d\.\d]+)s', line)
                task_dict['cpu_time'] = result.group(1)
                task_dict['wall_time'] = result.group(2)
                result_list.append(task_dict.copy())
                state = None
                
        # 3. Create output nodes
        # - ParameterData
        output_nodes_list.append(('parameters', ParameterData(dict={'data':result_list})))
    
        # - StructureData
        if positions:
            postions = np.array(positions, np.float64)
            if not cell:
                # If the cell is specified, ASE defaults to 0,0,0, which throws an AiiDA
                # error for cell with volume of 0. Here we arbitrarily set the cell. 
                # This isn't really satisfactory.
                # TODO: Look into changing AiiDA test of cell volume.
                cell = (1.,1.,1.)
            else:
                cell = np.array(cell, np.float64)
            atoms = ase.Atoms(symbols=symbols, positions=positions, cell=cell)
            output_nodes_list.append(('output_structure', StructureData(ase=atoms)))
    
        return output_nodes_list

