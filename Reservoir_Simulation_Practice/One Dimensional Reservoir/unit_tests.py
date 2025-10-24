#!/usr/bin/env python

# Copyright 2025 Andrew Z. Miller
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import unittest
import nbconvert
import numpy as np
import scipy.sparse

with open("one_d_reservoir_simulator.ipynb") as f:
    exporter = nbconvert.PythonExporter()
    python_file, _ = exporter.from_file(f)


with open("one_d_reservoir_simulator.py", "w") as f:
    f.write(python_file)

from one_d_reservoir_simulator import OneDimReservoir

class TestSolution(unittest.TestCase):

    def setUp(self):
    
        self.inputs = {
            'conversion factor': 6.33e-3,
            'fluid': {
                'water': {
                    'compressibility': 1e-6, #psi^{-1}
                    'viscosity': 1, #cp
                    'formation volume factor': 1,
                },
            },
            'reservoir': {
                'permeability': 50, #mD
                'porosity': 0.2,
                'length': 10000, #ft
                'cross sectional area': 200000 #ft^2
            },
            'initial conditions': {
                'pressure': 1000 #psi
            },
            'boundary conditions': {
                'left': {
                    'type': 'prescribed pressure',
                    'value': 2000 #psi
                },
                'right': {
                    'type': 'prescribed flux',
                    'value': 0 #ft^3/day
                }
            },
            'numerical': {
                'solver': 'implicit',
                'number of grids': {'x': 4},
                'time step': 1, #day
                'number of time steps' : 3 
            },
            'plots': {
                'frequency': 1
            }
        }

        return

    def test_compute_transmissibility(self):

        problem = OneDimReservoir(self.inputs)

        np.testing.assert_allclose(problem.compute_transmissibility(0,0), 4000.0)

        return

    def test_compute_accumulation(self):

        problem = OneDimReservoir(self.inputs)

        np.testing.assert_allclose(problem.compute_accumulation(0), 100.0)

        return 

    def test_is_transmissiblity_matrix_sparse(self):

        problem = OneDimReservoir(self.inputs)

        assert scipy.sparse.issparse(problem.T)

        return

    def test_implicit_solve_one_step(self):

        implicit = OneDimReservoir(self.inputs)
        implicit.solve_one_step()
        np.testing.assert_allclose(implicit.get_solution(), 
                                   np.array([1295.1463, 1051.1036, 1008.8921, 1001.7998]), 
                                   atol=0.5)
        return

    def test_explicit_solve_one_step(self):

        self.inputs['numerical']['solver'] = 'explicit'

        explicit = OneDimReservoir(self.inputs)

        explicit.solve_one_step()

        np.testing.assert_allclose(explicit.get_solution(), 
                               np.array([ 1506., 1000.,  1000.,  1000.004]), 
                               atol=0.5)
        return 

    def test_mixed_method_solve_one_step_implicit(self):

        self.inputs['numerical']['solver'] = {'mixed method': {'theta': 0.0}}

        mixed_implicit = OneDimReservoir(self.inputs)

        mixed_implicit.solve_one_step()

        np.testing.assert_allclose(mixed_implicit.get_solution(),
                               np.array([1295.1463, 1051.1036, 1008.8921, 1001.7998]),
                               atol=0.5)
        return

    def test_mixed_method_solve_one_step_explicit(self):

        self.inputs['numerical']['solver'] = {'mixed method': {'theta': 1.0}}

        mixed_explicit = OneDimReservoir(self.inputs)

        mixed_explicit.solve_one_step()

        np.testing.assert_allclose(mixed_explicit.get_solution(),
                               np.array([ 1506., 1000.,  1000.,  1000.004]),
                               atol=0.5)
        return

    def test_mixed_method_solve_one_step_crank_nicolson(self):

        self.inputs['numerical']['solver'] = {'mixed method': {'theta': 0.5}}

        mixed = OneDimReservoir(self.inputs)

        mixed.solve_one_step()

        np.testing.assert_allclose(mixed.get_solution(),
                                   np.array([ 1370.4,  1037.8 ,  1003.8,  1000.4]),
                                   atol=0.5)
        return

    def test_implicit_solve(self):

        implicit = OneDimReservoir(self.inputs)
        implicit.solve()
        np.testing.assert_allclose(implicit.get_solution(),
                                   np.array([1582.9, 1184.8, 1051.5, 1015.9]),
                                   atol=0.5)
        return

    def test_implicit_solve_reverse_boundary_conditions(self):


        self.inputs['boundary conditions'] = {
                'right': {
                    'type': 'prescribed pressure',
                    'value': 2000 #psi
                },
                'left': {
                    'type': 'prescribed flux',
                    'value': 0 #ft^3/day
                }
            }

        implicit = OneDimReservoir(self.inputs)
        implicit.solve()
        np.testing.assert_allclose(implicit.get_solution(),
                                   np.array([1015.9, 1051.5, 1184.8, 1582.9]), 
                                   atol=0.5)
        return

    def test_explicit_solve(self):

        self.inputs['numerical']['solver'] = 'explicit'

        explicit = OneDimReservoir(self.inputs)

        explicit.solve()

        np.testing.assert_allclose(explicit.get_solution(),
                               np.array([1689.8, 1222.3, 1032.4, 1000.0]),
                               atol=0.5)
        return

    def test_mixed_method_solve_crank_nicolson(self):

        self.inputs['numerical']['solver'] = {'mixed method': {'theta': 0.5}}

        mixed = OneDimReservoir(self.inputs)

        mixed.solve()

        np.testing.assert_allclose(mixed.get_solution(),
                                   np.array([ 1642.0,  1196.5,  1043.8,  1009.1]),
                                   atol=0.5)
        return

    def test_implicit_heterogeneous_permeability_solve_one_step(self):

        self.inputs['reservoir']['permeability'] = [10., 100., 50., 20]

        implicit = OneDimReservoir(self.inputs)
        implicit.solve_one_step()
        np.testing.assert_allclose(implicit.get_solution(),
                                   np.array([1085.3,  1005.8,  1001.3,  1000.1]),
                                   atol=0.5)
        return

    def test_implicit_heterogeneous_permeability_and_grid_size_solve_one_step(self):

        self.inputs['reservoir']['permeability'] = [10., 100., 50., 20]
        self.inputs['numerical']['delta x'] = [2000., 3000., 1500., 3500]

        implicit = OneDimReservoir(self.inputs)
        implicit.solve_one_step()
        np.testing.assert_allclose(implicit.get_solution(),
                                   np.array([1123.0,  1008.5,  1003.1,  1000.2]),
                                   atol=0.5)
        return


    def test_implicit_heterogeneous_permeability_and_grid_size_solve(self):

        self.inputs['reservoir']['permeability'] = [10., 100., 50., 20]
        self.inputs['numerical']['delta x'] = [2000., 3000., 1500., 3500]

        implicit = OneDimReservoir(self.inputs)
        implicit.solve()
        np.testing.assert_allclose(implicit.get_solution(),
                                   np.array([1295.6,  1039.1,  1019.9,  1002.5]),
                                   atol=0.5)
        return

if __name__ == '__main__':
    unittest.main()