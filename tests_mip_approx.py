
import cdc_formulation as form
from tqdm import tqdm
import iss_checker
from random import random
import numpy as np
import math
from matplotlib import pyplot as plt

# This test ensures that a formulation for a relaxation of a constraint is indeed a relaxation,
# that is, that the points of a manifold are in the feasible set of its purported envelope.

# To do this, we generate random points p_i in the feasible set, and create an optimization problem
# where we try to minimize the distance between x and each p_i where x is in the relaxation.
# If the minimum distance is zero for all points p_i, the relaxation passes the test.

def test_circumference(mode, ntests, tolerance = 1e-8): # Test for the cone shell constraint over the circumference x²+y²=1
	exact_matches = 0
	max_error = tolerance
	for random_point in tqdm(range(ntests)):
		angle = random()*2*math.pi
		pointx, pointy = math.cos(angle), math.sin(angle)
		vars2 = {'var_x': {}, 'var_y': {}, 'var_z': {'lower': 1, 'upper': 1}, 'objective': {'lower': 0}}
		# cones2 = {'cone_obj': [{'objective': 1}, {var1: 1, 'constant_1': -pointx}, {var2: 1, 'constant_1': -pointy}, {}]}
		linear_ineq = {
			'obj1': {'objective': 1, 'var_x': -1, '>=': -pointx},
			'obj2': {'objective': 1, 'var_y': -1, '>=': -pointy},
			'obj3': {'objective': 1, 'var_x': 1, '>=': pointx},
			'obj4': {'objective': 1, 'var_y': 1, '>=': pointy},
		}
		shell = {'shell1': [{'var_z': 1}, {'var_x': 1}, {'var_y': 1}, {'z_bound': 1}]}
		formulation = form.QuadProblem('Cvx_hull', vars2, linear_ineq, {}, shell, {}, {'objective': 1})
		lp_results = formulation.solve(mode=mode, time_limit=None, helix_depth = 7)
		formulation.model.optimize()
		distance = abs(formulation.modelvars['var_x'].X - pointx) + abs(formulation.modelvars['var_y'].X - pointy)
		exact_matches += (distance < tolerance)
		max_error = max(max_error, distance)
	nerrors = ntests-exact_matches
	return nerrors, max_error

def test_helix(mode, ntests, tolerance=1e-8): # Test for the helix constraint over the graph of y=sin(x)
	exact_matches = 0
	max_error = tolerance
	for random_point in tqdm(range(ntests)):
		angle = random()*2*math.pi
		pointx, pointy = math.cos(angle), math.sin(angle)
		vars2 = {'var_x': {}, 'var_y': {}, 'var_z': {'lower': 1, 'upper': 1}, 'var_alpha': {'lower': 0, 'upper': 2*math.pi}, 'objective': {'lower': 0}}
		# cones2 = {'cone_obj': [{'objective': 1}, {var1: 1, 'constant_1': -pointx}, {var2: 1, 'constant_1': -pointy}, {}]}
		linear_ineq = {
			'obj1': {'objective': 1, 'var_x': -1, '>=': -pointx},
			'obj2': {'objective': 1, 'var_y': -1, '>=': -pointy},
			'obj3': {'objective': 1, 'var_x': 1, '>=': pointx},
			'obj4': {'objective': 1, 'var_y': 1, '>=': pointy},
			'obj5': {'objective': 1, 'var_alpha': 1, '>=': angle},
			'obj6': {'objective': 1, 'var_alpha': -1, '>=': -angle},
		}
		helix = {'helix1': [{'var_alpha': 1}, {'var_z': 1}, {'var_x': 1}, {'var_y': 1}, {'z_bound': 10}]}
		formulation = form.QuadProblem('Cvx_hull', vars2, linear_ineq, {}, {}, helix, {'objective': 1})
		lp_results = formulation.solve(mode='cdc', time_limit=None, helix_depth = 7)
		formulation.model.optimize()
		distance = abs(formulation.modelvars['var_x'].X - pointx) + abs(formulation.modelvars['var_y'].X - pointy) + abs(formulation.modelvars['var_alpha'].X - angle)
		exact_matches += (distance < tolerance)
		max_error = max(max_error, distance)
	errors = ntests-exact_matches
	return nerrors, max_error

if __name__ == '__main__':
	ntests = 1000
	nerrors, max_error = test_circumference('cdc', ntests)
	res1 = ('\033[91m{-} Test 1 failed' if nerrors else '\033[92m{+} Test 1 Passed', '(%d/%d) errors\033[0m'%(nerrors, ntests),
		'\n Maximum error bound: %.2E'%max_error)

	nerrors, max_error = test_helix('cdc', ntests)
	res2 = ('\033[91m{-} Test 2 failed' if nerrors else '\033[92m{+} Test 2 Passed', '(%d/%d) errors\033[0m'%(nerrors, ntests),
		'\n Maximum error bound: %.2E'%max_error)

	print(*res1)
	print(*res2)

