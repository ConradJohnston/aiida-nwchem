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
    - Frequency analysis

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

        with open(output_path) as f:
            all_lines = [x.strip('\n') for x in f.readlines()]
        
        # A list to hold the dictionaries for each calculation task
        result_list = []
        # List to hold all output nodes 
        output_nodes_list = []

        # Supported task types - dict to map task names to functions
        parse_task = {
            'scf': parse_scf,
            'dft': parse_dft,
            'geoopt': parse_geoopt,
            'freq': parse_freq,
        }
                               
        # Cut the data in to lists
        task_types, task_lines_list = separate_tasks(all_lines)

        # Loop over tasks, parsing each one
        for task, task_lines in zip(task_types, task_lines_list):
            result_dict, nodes = parse_task[task](task_lines)
            result_list.append(result_dict)
            output_nodes_list = output_nodes_list + nodes
        
        # Create ParameterData node
        output_nodes_list.append(('parameters', ParameterData(dict={'data':result_list})))
        
        return output_nodes_list
 

    def separate_tasks(lines):
        """
        Parse the different task types and their associated lines

        args: lines: all the lines in the output file
        returns: task_types: a list of the tasks parsed
                 task_lines_list: a list of lists or lines from each task
        """

        # List holding all the lines from a given task
        task_lines = []
        # List to hold all of the lists of lines from the tasks
        task_lines_list = []
        # Task types
        task_types = []
        # State to track if we're in a task or not
        in_task = False
        
        for line in all_lines:

            # Determine which NWChem module we are using
            if not in_task:
                if re.match('^\s*NWChem SCF Module\s*$',line):
                    in_task = True
                    task_types.append('scf')
                elif re.match('^\s*NWChem DFT Module\s*$',line):
                    in_task = True
                    task_types.append('scf')
                elif re.match('^\s*NWChem Geometry Optimization\s*$',line):
                    in_task = True
                    task_types.append('geoopt')
                elif re.match('^\s*NWChem Nuclear Hessian and Frequency Analysis\s*$',line):
                    in_task = True
                    task_types.append('freq')
                else:
                    continue
            else:
                task_lines.append(line)
                if re.match('^ Task  times  cpu:', line):
                   in_task = False
                   task_lines_list.append(task_lines)
                   task_lines = []

        return task_types, task_lines_list

    def parse_scf(lines):
        """
        Parse a SCF(HF) task block
        
        args: lines: the lines to parse
        returns: task_dict : a dictionary of results for the task
                 nodes       : the nodes created by parsing
        """
         
        task_dict = {'task':'scf'}
        nodes = []
        state = None
         
        for line in lines:

            if re.match('^\s*Final [ROU]+HF\s*results\s*$',line):
                state = 'final-results'
            if state == 'final-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict[key] = result.group(2)
        
            # End of task
            if re.match('^ Task  times  cpu:', line):
                result = re.match('^ Task  times  cpu:\s*([\d\.\d]+)s\s*wall:\s*([\d\.\d]+)s', line)
                task_dict['cpu_time'] = result.group(1)
                task_dict['wall_time'] = result.group(2)
                break

        return task_dict, nodes


    def parse_dft(lines):
        """
        Parse a DFT task block
        
        args: lines: the lines to parse
        returns: task_dict : a dictionary of results for the task
                 nodes       : the nodes created by parsing
        """
         
        task_dict = {'task':'dft'}
        nodes = []
        state = None

        for line in lines:
            # Note the search for the Total DFT energy. NWChem doesn't otherwise
            # announce that the results are being printed. 
            if re.match('^\s*Total DFT energy', line):
                state = 'final-results'
            if state == 'final-results':
                result = re.match('^\s*([^=]+?)\s*=\s*([\-\d\.]+)$',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict[key] = result.group(2)
            
            # End of task
            if re.match('^ Task  times  cpu:', line):
                result = re.match('^ Task  times  cpu:\s*([\d\.\d]+)s\s*wall:\s*([\d\.\d]+)s', line)
                task_dict['cpu_time'] = result.group(1)
                task_dict['wall_time'] = result.group(2)
                break
         
        return task_dict, nodes

    
    def parse_geoopt(lines):
        """
        Parse a geometry optimisation task block
        
        args: lines: the lines to parse
        returns: task_dict : a dictionary of results for the task
                 nodes       : the nodes created by parsing
        """
          
        task_dict = {'task':'geo-opt'}
        nodes = []
        state = None
        symbols = []
        positions = []
        
        for line in lines:
        
            if re.match('^\s*Optimization converged\s*$',line):
                state = 'final-results'
                continue
            if state == 'final-results':
                # Parse step and energy
                result = re.match('^@\s*([\d]+)\s*([\-\d\.]+)',line)
                if result:
                    task_dict['final_step'] = result.group(1)
                    # Note, this energy may not be perfect - do another scf on the final structure
                    task_dict['final_opt_energy'] = result.group(2)
                    continue
                # Parse coords
                if re.match('^\s*Output coordinates in angstroms\s*',line):
                    state = 'final-coords'
                    continue
                # Parse cell
                if re.match('^\s*lattice vectors in angstroms',line):
                    state = 'final-cell'
                    continue
            if state == 'final-coords':
                result = re.match('^\s*[\d]+\s*([a-zA-Z]+)\s*[\-\d]+\.[\d]+'\
                                  '\s*([\-\d\.]+)\s*([\-\d\.]+)\s*([\-\d\.]+)', line)
                if result:
                    symbols.append(result.group(1))
                    positions.append([result.group(2),result.group(3),result.group(4)])
                else:
                    state = 'final-results'
                continue
            if state == 'final-cell':
                result = re.match('^\s*a[1-3]=<\s*([\d\.\d]+)\s*([\d\.\d]+)\s*([\d\.\d]+)', line) 
                if result:
                    cell.append([result.group(1),result.group(2),result.group(3)])
                else:
                    state = 'final-results'

            if re.match('^ Task  times  cpu:', line):
                result = re.match('^ Task  times  cpu:\s*([\d\.\d]+)s\s*wall:\s*([\d\.\d]+)s', line)
                task_dict['cpu_time'] = result.group(1)
                task_dict['wall_time'] = result.group(2)
                break
    
        # Create StructureData node 
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
        nodes.append(('output_structure', StructureData(ase=atoms)))
        
        return task_dict, nodes


    def parse_freq(lines):
        """
        Parse a frequency analysis task block
        
        args: lines: the lines to parse
        returns: task_dict : a dictionary of results for the task
                 nodes     : the nodes created by parsing
        """
          
        task_dict = {'task':'freq'}
        nodes = []
        state = None
        
        for line in lines:
             
            if re.match('^\s*Rotational Constants\s*$',line):
                state = 'final-results'
                continue
            if state == 'final-results':
                result = re.match('^\s([A-z\s\(\)]+)\s+=\s*([\d\.]+)',line)
                if result:
                    if result.group(1) == 'Total Entropy':
                        state = 'final-entropy'
                        task_dict['entropy'] = {}
                        key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                        task_dict['entropy'][key] = result.group(2)
                        continue
                    elif result.group(1) == 'Cv (constant volume heat capacity)':
                        state = 'final-cv'
                        task_dict['head_capacity'] = {}
                        task_dict['heat_capacity']['total_cv'] = result.group(2)
                        continue
                    else:
                        key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                        task_dict[key] = result.group(2)
                        continue
                # Derivative Dipole
                if re.search('Projected Derivative Dipole',line):
                    state = 'final-dipole'
                    dipole_list = []
                    frequencies = []
                    continue
                # Infrared 
                if re.search('Projected Infra Red',line):
                    state = 'final-freq-results-ir'
                    intensities = []
                    continue
            # Entropy
            if state == 'final-entropy':
                result = re.match('^\s*-\s([A-z\s\(\)]+)\s+=\s*([\d\.]+)',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict['entropy'][key] = result.group(2)
                else:
                    state = 'final-results'
                    continue
            # Heat capacity
            if state == 'final-cv':
                result = re.match('^\s*-\s([A-z\s\(\)]+)\s*=\s*([\d\.]+)',line)
                if result:
                    key = re.sub('[^a-zA-Z0-9]+', '_', result.group(1).lower())
                    task_dict['heat_capacity'][key] = result.group(2)
                else:
                    state = 'final-results'
                    continue
            # Parse dipole data
            if state == 'final-freq-results-dipole':
                result = re.match('^\s*[\d]\s*([\-\d\.]+)\s*\|\|'
                                  's*([-\d.]+)\s*([-\d.]+)\s*([-\d.]+)$',line)
                if result:
                    # Get vibrational eigenvalues (cm^-1)
                    frequencies.append(re.group(1))
                    # Get dipole moments (cartesian, debye/angs)
                    dipoles_list.append([re.group(2), re.group(3), re.group(4)])
                    continue
                if re.match('^\s-+$', line):
                    state = 'final-results' 
                    task_dict['frequencies'] = frequencies
                    task_dict['dipoles'] = np.array(dipoles_lisr, np.float64)
                    continue
            # Parse IR data 
            if state == 'final-freq-results-ir':
                result = re.match('^\s*[\d]\s*[\-\d\.]+\s*\|\|\s*([-\d.]+)'
                                  '\s*([-\d.]+)\s*([-\d.]+)\s*([-\d.]+)$',line)
                if result:
                    # Get intensity (arbitraty units)
                    intensities.append(re.group(4))
                    continue
                if re.match('^\s-+$', line):
                    state = 'final-results' 
                    task_dict['ir-intensities'] = intensities
                    continue
            # End of task
            if re.match('^ Task  times  cpu:', line):
                result = re.match('^ Task  times  cpu:\s*([\d\.\d]+)s'
                                  '\s*wall:\s*([\d\.\d]+)s', line)
                task_dict['cpu_time'] = result.group(1)
                task_dict['wall_time'] = result.group(2)
                break
        
        return task_dict, nodes
