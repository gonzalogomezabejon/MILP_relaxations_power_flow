
import numpy as np
import cdc_formulation
import math
from matplotlib import pyplot as plt
from random import random
from gurobipy import GRB

# This script tests the different compact relaxations by generating
# a plot of their feasible regions.

# Get the closest point to randomly generated points in the feasible region
def generate_points(q, var1, var2, range1, range2, n_points, params, mode):
	points_x, points_y = [], []
	for ii in range(n_points):
		pointx, pointy = range1[0]+random()*(range1[1]-range1[0]), range2[0]+random()*(range2[1]-range2[0])
		vars2 = {**q.variables, 'objective': {'lower': 0}}
		cones2 = {**q.cone_ineq, 'cone_obj': [{'objective': 1}, {var1: 1, 'constant_1': -pointx}, {var2: 1, 'constant_1': -pointy}, {}]}
		formulation = cdc_formulation.QuadProblem('Cvx_hull', vars2, q.linear_ineq, cones2, q.cone_outeq, q.helix_cones, {'objective': 1}, bilinears=q.bilinears)
		lp_results = formulation.solve(mode=mode, **params)
		formulation.model.optimize()
		points_x.append(formulation.modelvars[var1].X)
		points_y.append(formulation.modelvars[var2].X)
	return points_x, points_y

# Get n points in the convex hull of a region by optimizing in all directions in the plane
def generate_convex_hull(q, var1, var2, n_directions, params, mode):
	points_x = []
	points_y = []
	for ii in range(n_directions+1):
		angle = 2*math.pi*ii/n_directions
		weights = {var1: math.sin(angle), var2: math.cos(angle)}
		formulation = cdc_formulation.QuadProblem('Cvx_hull', q.variables, q.linear_ineq, q.cone_ineq, q.cone_outeq, q.helix_cones, weights, bilinears=q.bilinears)
		lp_results = formulation.solve(mode=mode, **params)
		formulation.model.optimize()
		points_x.append(formulation.modelvars[var1].X)
		points_y.append(formulation.modelvars[var2].X)
	return points_x, points_y

def generate_contour(q, var1, var2, n_directions, params, binary_names):
	all_points = []
	for jj in range(2**len(binary_names)):
		print(jj)
		fixed_binaries = {binary_names[ii]: jj%(2**(ii+1))//(2**ii) for ii in range(len(binary_names))}
		points_x = []
		points_y = []
		for ii in range(n_directions+1):
			angle = 2*math.pi*ii/n_directions
			weights = {var1: math.sin(angle), var2: math.cos(angle)}
			formulation = cdc_formulation.QuadProblem('Contour_piece', q.variables, q.linear_ineq, q.cone_ineq, q.cone_outeq, q.helix_cones, weights, exponentials = q.exponentials)

			formulation.formulate_mixed(**params)
			for binary in fixed_binaries:
				formulation.model.addConstr(formulation.modelvars[binary] == fixed_binaries[binary])
			lp_results = formulation.solve(mode=None)
			# formulation.model.optimize()
			points_x.append(formulation.modelvars[var1].X)
			points_y.append(formulation.modelvars[var2].X)
		all_points.append((points_x, points_y))
	return all_points
	
def generate_branching_pattern(q, var1, var2, n_directions, params, binary_names):
	all_points = []
	piece = 10 #int(random()*2**len(binary_names))
	for jj in range(1+len(binary_names)):
		print(jj)
		fixed_binaries = {binary_names[ii]: piece%(2**(ii+1))//(2**ii) for ii in range(jj)}
		non_fixed_binaries = [binary_names[ii] for ii in range(jj, len(binary_names))]
		points_x = []
		points_y = []
		for ii in range(n_directions+1):
			angle = 2*math.pi*ii/n_directions
			weights = {var1: math.sin(angle), var2: math.cos(angle)}
			formulation = cdc_formulation.QuadProblem('Contour_piece', q.variables, q.linear_ineq, q.cone_ineq, q.cone_outeq, q.helix_cones, weights, exponentials = q.exponentials)

			formulation.formulate_mixed(**params)
			for binary in fixed_binaries:
				formulation.model.addConstr(formulation.modelvars[binary] == fixed_binaries[binary])
			for binary in non_fixed_binaries:
				formulation.modelvars[binary].setAttr(GRB.Attr.VType, GRB.CONTINUOUS)
				formulation.model.addConstr(formulation.modelvars[binary] <= 1)
				formulation.model.addConstr(formulation.modelvars[binary] >= 0)
			lp_results = formulation.solve(mode=None)
			# formulation.model.optimize()
			points_x.append(formulation.modelvars[var1].X)
			points_y.append(formulation.modelvars[var2].X)
		all_points.append((points_x, points_y))
	return all_points
	
def generate_error(q, n_points, params):
	max_error = 0.0
	for ii in range(n_points):
		angle_val = random()*2*math.pi
		weights = {'objective': 1 if random()>0.5 else -1}
		vars2 = {**q.variables, 'objective': {'lower': 0}}
		linear_ineq2 = {**q.linear_ineq, 'obj1': {'objective' : 1, 'my_y': 1, '==': math.sin(angle_val)}, 'fix': {'my_angle': 1, '==': angle_val}}
		formulation = cdc_formulation.QuadProblem('Contour_piece', vars2, linear_ineq2, q.cone_ineq, q.cone_outeq, q.helix_cones, weights)
		formulation.formulate_milp(**params)
		lp_results = formulation.solve(mode=None)
		# print(angle_val, lp_results)
		# formulation.model.optimize()
		obj_value = -formulation.model.objVal
		max_error = max(obj_value, max_error)
	return max_error
	

def test_helix(mode):
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_angle'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 2, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})

	for my_cone_depth in [2,3,4]:#[4,5,6]:
		points_x, points_y = generate_convex_hull(quad, target1, target2, 100, {"helix_depth": my_cone_depth}, mode)
		plt.plot(points_x, points_y)
		plt.plot([math.pi*ii/100. for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		# plt.savefig('test_helix_hull_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

		points_x, points_y = generate_points(quad, target1, target2, [-1, 7.3], [-1.4, 1.4], 2000, {"helix_depth": my_cone_depth}, mode)
		plt.scatter(points_x, points_y, s=2)
		plt.plot([math.pi*ii/100. for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], c='r')
		# plt.savefig('test_helix_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

def test_nmdt(mode):
	linear_ineq = {}
	variables = {'my_x': {'lower': -1, 'upper': 1.4}, 'my_y': {'lower': -1, 'upper': 2}, 'constant_1': {'lower': 1, 'upper': 1}}
	target1 = 'my_x'
	target2 = 'my_y'
	bilinears = {'parab1': [{'my_y': 1}, {'my_x': 1}, {'my_x': 1}, {}]}
	cones = {}
	quad = cdc_formulation.QuadProblem('Nmdt_parabola', variables, {}, {}, {}, {}, {'constant_1': 1}, bilinears = bilinears)

	for my_cone_depth in [2,3,4]:#[4,5,6]:
		points_x, points_y = generate_convex_hull(quad, target1, target2, 100, {"bilinear_depth": my_cone_depth}, mode)
		plt.plot(points_x, points_y)
		plt.plot([ii/100. for ii in range(201)], [ii*ii/10000. for ii in range(201)], c='r')
		# plt.plot([ii/200. for ii in range(201)], [ii*ii/40000. for ii in range(201)])
		# plt.savefig('test_helix_hull_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

		points_x, points_y = generate_points(quad, target1, target2, [-1.1, 1.6], [-1.1, 1.6], 2000, {"bilinear_depth": my_cone_depth}, mode)
		plt.scatter(points_x, points_y, s=2)
		plt.plot([ii/100. for ii in range(201)], [ii*ii/10000. for ii in range(201)], c='r')
		# plt.savefig('test_helix_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

def test_circle(mode):
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}}
	target1 = 'my_x'
	target2 = 'my_y'
	cones = {'cone1': [{'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'ang_max': math.pi/3, 'ang_min': -math.pi/3, 'z_bound': 10}]} # 
	quad = cdc_formulation.QuadProblem('Cone test', variables, {}, {}, cones, {}, {'constant_1': 1})

	for my_cone_depth in [2,3,4]:
		points_x, points_y = generate_convex_hull(quad, target1, target2, 100, {"cone_depth": my_cone_depth, "helix_depth": my_cone_depth}, mode)
		plt.plot(points_x, points_y)
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		# plt.savefig('test_helix_hull_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

		points_x, points_y = generate_points(quad, target1, target2, [-1.4, 1.4], [-1.4, 1.4], 2000, {"cone_depth": my_cone_depth, "helix_depth": my_cone_depth}, mode)
		plt.scatter(points_x, points_y, s=2)
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], c='r')
		# plt.savefig('test_helix_depth%d'%my_cone_depth)
		plt.show()
		plt.close()


def test_bilinear(mode):
	linear_ineq = {}
	variables = {'my_x': {'lower': 1, 'upper': 100}, 'my_y': {'lower': 1, 'upper': 100}, 'my_z': {}, 'constant_1': {'lower': 1, 'upper': 1}}
	target1 = 'my_x'
	target2 = 'my_z'
	lin_constraints = {'lin1': {'my_x': 1, 'my_y': -1, '==': 0}}
	# cones = {'cone1': [{'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'ang_max': math.pi/3, 'ang_min': -math.pi/3, 'z_bound': 10}]} # 
	bilinear = {'bil1': [{'my_z': 1}, {'my_y': 1}, {'my_x': 1}, {'lb2': 1,'ub2': 100,'lb3': 1,'ub3': 100}]}
	quad = cdc_formulation.QuadProblem('Cone test', variables, lin_constraints, {}, {}, {}, {'constant_1': 1}, bilinear)

	for my_cone_depth in [2,3,4]:
		points_x, points_y = generate_convex_hull(quad, target1, target2, 100, {"cone_depth": my_cone_depth, "helix_depth": my_cone_depth}, mode)
		plt.plot(points_x, points_y)
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		# plt.savefig('test_helix_hull_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

		points_x, points_y = generate_points(quad, target1, target2, [-1.4, 1.4], [-1.4, 1.4], 2000, {"cone_depth": my_cone_depth, "helix_depth": my_cone_depth}, mode)
		plt.scatter(points_x, points_y, s=2)
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], c='r')
		# plt.savefig('test_helix_depth%d'%my_cone_depth)
		plt.show()
		plt.close()

def test_sine():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_angle'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})
	plt.style.use("seaborn-v0_8")
	for helix_depth in [2,3,4]:
		all_points = generate_contour(quad, target1, target2, 40, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		for points_x, points_y in all_points:
			plt.fill(points_x, points_y, facecolor='#AAAAFF', edgecolor='b', linewidth=1)
		# plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		plt.plot([math.pi*ii/100. for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], 'k--', linewidth=1)
		plt.savefig('test_sine_contour_%d'%helix_depth)
		# plt.show()
		plt.close()

def test_sine_branching():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_angle'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})
	plt.style.use("seaborn-v0_8")
	for helix_depth in [4]:
		all_points = generate_branching_pattern(quad, target1, target2, 400, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		contour_points = generate_contour(quad, target1, target2, 40, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		colours = '158BDDDDD'
		cc = 0
		for points_x, points_y in all_points:
			plt.fill(points_x, points_y, facecolor='#'+colours[helix_depth-cc]*4+'FF', edgecolor='b', linewidth=0.5, alpha=0.5)
			cc += 1
		for points_x, points_y in contour_points:
			plt.fill(points_x, points_y, facecolor='#FFBBBB', edgecolor='r', linewidth=0.5, alpha=0.7)
		points_x, points_y = all_points[-1]
		plt.fill(points_x, points_y, facecolor='#BBFFBB', edgecolor='g', linewidth=0.5, alpha=0.8)
		# plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		plt.plot([math.pi*ii/100. for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], 'k--', linewidth=1)
		# plt.savefig('test_sine_contour_%d'%helix_depth)
		plt.show()
		plt.close()

def test_circle_contour():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_x'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})
	plt.style.use("seaborn-v0_8")
	for helix_depth in [2,3,4]:
		all_points = generate_contour(quad, target1, target2, 20, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		for points_x, points_y in all_points:
			plt.fill(points_x, points_y, facecolor='#AAAAFF', edgecolor='b', linewidth=1)
		# plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)])
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], 'k--', linewidth=1)
		plt.savefig('test_circle_contour_%d'%helix_depth)
		# plt.show()
		plt.close()

def test_circle_branching():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_x'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})
	plt.style.use("seaborn-v0_8")
	for helix_depth in [4]:
		all_points = generate_branching_pattern(quad, target1, target2, 400, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		contour_points = generate_contour(quad, target1, target2, 40, {"helix_depth": helix_depth}, ['binary_%d_helix_helix1'%ii for ii in range(1, helix_depth+1)])
		colours = '158AC'
		cc = 0
		for points_x, points_y in all_points:
			plt.fill(points_x, points_y, facecolor='#'+colours[4-cc]*4+'FF', edgecolor='b', linewidth=0.5, alpha=0.5)
			cc += 1
		for points_x, points_y in contour_points:
			plt.fill(points_x, points_y, facecolor='#FFBBBB', edgecolor='r', linewidth=0.5, alpha=0.7)
		points_x, points_y = all_points[-1]
		plt.fill(points_x, points_y, facecolor='#BBFFBB', edgecolor='g', linewidth=0.5, alpha=0.8)
		plt.plot([math.cos(math.pi*ii/100.) for ii in range(201)], [math.sin(math.pi*ii/100.) for ii in range(201)], 'k--', linewidth=1)
		plt.gca().set_aspect('equal', adjustable='box')
		# plt.savefig('test_circle_contour_%d'%helix_depth)
		plt.show()
		plt.close()

def test_expo_branching():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}}#, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_x'
	target2 = 'my_y'
	helix = {}#{'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	expo1 = [{'my_x': 1}, {'my_y': 1, 'my_x':2}]
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1}, exponentials = {'expo1': expo1})
	plt.style.use("seaborn-v0_8")
	for helix_depth in [3]:
		all_points = generate_branching_pattern(quad, target1, target2, 400, {"helix_depth": helix_depth}, ['binary_%d_expo_expo1'%ii for ii in range(1, helix_depth+1)])
		contour_points = generate_contour(quad, target1, target2, 40, {"helix_depth": helix_depth}, ['binary_%d_expo_expo1'%ii for ii in range(1, helix_depth+1)])
		colours = '158AC'
		cc = 0
		for points_x, points_y in all_points:
			plt.fill(points_x, points_y, facecolor='#'+colours[4-cc]*4+'FF', edgecolor='b', linewidth=0.5, alpha=0.5)
			cc += 1
		for points_x, points_y in contour_points:
			plt.fill(points_x, points_y, facecolor='#FFBBBB', edgecolor='r', linewidth=0.5, alpha=0.7)
		points_x, points_y = all_points[-1]
		plt.fill(points_x, points_y, facecolor='#BBFFBB', edgecolor='g', linewidth=0.5, alpha=0.8)
		plt.plot([ii/200. for ii in range(201)], [math.exp(ii/200.)-ii/100. for ii in range(201)], 'k--', linewidth=1, label='$y=\\exp(x)-2x$')
		plt.legend()
		# plt.gca().set_aspect('equal', adjustable='box')
		# plt.savefig('test_circle_contour_%d'%helix_depth)
		plt.show()
		plt.close()

def test_sine_error():
	linear_ineq = {}
	variables = {'my_x': {}, 'my_y': {}, 'constant_1': {'lower': 1, 'upper': 1}, 'my_angle': {'lower': -2*math.pi, 'upper': 2*math.pi}}
	target1 = 'my_angle'
	target2 = 'my_y'
	helix = {'helix1': [{'my_angle': 1}, {'constant_1': 1}, {'my_x': 1}, {'my_y': 1}, {'z_bound': 1, 'z_lower_bound': 1}]} # 'ang_max': 1, 'ang_min': -1, 
	cones = {}
	quad = cdc_formulation.QuadProblem('Helix_test', variables, {}, {}, {}, helix, {'constant_1': 1})
	return {depth: generate_error(quad, 400, {"helix_depth": depth}) for depth in [3,4,5,6,7,8]}

def error_regression(data):
	from scipy.stats import linregress
	import matplotlib.lines as mlines
	import matplotlib.patches as mpatches
	plt.style.use("seaborn-v0_8")
	handles = []
	color = 'b'
	handles.append(mpatches.Patch(color='b', label='With cuts'))
	handles.append(mpatches.Patch(color='g', label='Without cuts'))
	for mydata in data:
		x = [it for it in mydata]
		y = [mydata[it] for it in mydata]
		gradient,intercept,r_value,p_value,std_err=linregress(x,[math.log(yy) for yy in y])
		plt.scatter(x,y)
		plt.plot(x, [math.exp(gradient*t+intercept) for t in x], color, linewidth = 1)
		lab_1, lab_2, lab_3 = "r² = %.5f"%(r_value**2), "slope = $-\\log(%.3f)$"%math.exp(-gradient), "intercept = %.3f"%intercept 
		handles.append(mlines.Line2D([], [], color = color, label = lab_1.replace('-', '−')))
		handles.append(mlines.Line2D([], [], color = color, label = lab_2.replace('-', '−')))
		# mlines.Line2D([], [], color = "black", label = lab_3.replace('-', '−'))])
		color = 'g'
	plt.xlabel("Number of integer variables")
	plt.ylabel("Maximum approximation error")
	plt.legend(handles = handles)
	plt.yscale('log')
	plt.show()
	plt.close()
	# plt.savefig(archivo, bbox_inches='tight')
	# plt.clf()



if __name__ == '__main__':
	# test_circle('cdc')
	# test_circle_contour()
	# test_helix('mixed')
	# test_bilinear('mixed')
	# test_nmdt('dnmdt')
	# test_sine()
	test_sine_branching()
	# test_circle_branching()
	# test_expo_branching()
	# res = test_sine_error()
	# print(res)
	# cut = {2: 0.5392838940182695, 3: 0.07820852634514203, 4: 0.018882889348051046, 5: 0.004785630909442062}
	# uncut = {2: 0.9999378065187428, 3: 0.738823728472477, 4: 0.37093953297946713, 5: 0.19501206676142793}
	# plt.style.use("seaborn-v0_8")
	# plt.scatter(range(2,6), [uncut[ii] for ii in range(2,6)])
	# plt.scatter(range(2,6), [cut[ii] for ii in range(2,6)])
	# plt.yscale('log')
	# plt.show()
	# plt.close()


