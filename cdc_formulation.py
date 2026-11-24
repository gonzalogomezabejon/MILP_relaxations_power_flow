#!/usr/bin/env python3

import gurobipy as gp
from gurobipy import GRB
from random import random
import time
import json
import math
from matplotlib import pyplot as plt
import threading
from functools import partial
from gurobipy import nlfunc

import graph_structure # for detecting bridges in power network graphs


my_env = gp.Env(empty=True)
my_env.setParam('OutputFlag', 0)
my_env.start()

def gp_callback(model, where, data, timeout):
	data[-1] = {'memory': 0, 'nodes': 0}
	if where != GRB.Callback.POLLING:
		runtime = model.cbGet(GRB.Callback.RUNTIME)
		data[runtime] = {}
		data[runtime]['memory'] = model.cbGet(GRB.Callback.MEMUSED)
		if where == GRB.Callback.MIPNODE:
			data[runtime]['nodes'] = model.cbGet(GRB.Callback.MIPNODE_NODCNT)
			# print ('Nodes:', data[runtime]['nodes'])
		if where == GRB.Callback.MIP:
			data[runtime]['nodes'] = model.cbGet(GRB.Callback.MIP_NODCNT)
			# print ('Nodes:', data[runtime]['nodes'])
		# print (runtime, where)
		if runtime > timeout:
			model.terminate()

def cplx_inv(a, b): # inverse of a complex number a+ib
	return a/(a*a+b*b), -b/(a*a+b*b)

def cplx_prod(c1, c2): # product of two complex numbers
	return c1[0]*c2[0]-c1[1]*c2[1], c1[0]*c2[1]+c1[1]*c2[0]

class QuadProblem:
	"""
		This class models an optimization problem with continuous and binary variables,
		and several types of constraints including linear, SOC, cone shell (x² + y² = z²),
		helix (x + iy = exp(z + iw)), bilinear (xy=z) and exponential (y = exp(x)) constraints.

		For each of the different MILP relaxations of those constraints, there is a method
		that outputs it as a gurobi model, and the method solve() solves it.
	"""
	def __init__(self, name, variables, linear, inner_cones, outer_cones, helix_cones, objective, bilinears = {}, exponentials = {}):
		self.name = name
		self.variables = variables
		self.linear_ineq = linear
		self.cone_ineq = inner_cones
		self.cone_outeq = outer_cones
		self.concave_cones = {}
		self.objective = objective
		self.helix_cones = helix_cones
		self.model = None
		self.solution_data = None
		self.bilinears = bilinears
		self.exponentials = exponentials

	def display(self):
		print ("Quadratic Formulation of '%s':\n"%self.name)
		for ineq_name in self.linear_ineq:
			ineq = self.linear_ineq[ineq_name]
			sign = [it for it in ['==', '>=', '<='] if it in ineq][0]
			print (ineq_name+':', ' + '.join(['%.6f*'%ineq[var]+var for var in ineq if var != sign]), sign, ineq[sign])
		for cone_name in self.cone_ineq:
			cone = self.cone_ineq[cone_name]
			x_var = ' + '.join(['%.6f*'%cone[1][var]+var for var in cone[1]])
			y_var = ' + '.join(['%.6f*'%cone[2][var]+var for var in cone[2]])
			z_var = ' + '.join(['%.6f*'%cone[0][var]+var for var in cone[0]])
			print (cone_name+': (', x_var, ')² + (', y_var, ')² <= (', z_var, ')²')
		for var in self.variables:
			print (var, ': ', self.variables[var])

	def check_solution(self, variables, margin = 1e-8, verbose = True, fill_variables = True):
		violations, cone_violations = [], []
		if fill_variables:
			variables['constant_1'] = 1
			for var in list(variables.keys()):
				if var.startswith('pg_'):
					variables['pg_sq'+var[2:]] = variables[var]**2
		for ineq_name in self.linear_ineq:
			ineq = self.linear_ineq[ineq_name]
			sign = [it for it in ['==', '>=', '<='] if it in ineq][0]
			total = sum([ineq[var]*variables[var] for var in ineq if var != sign])
			residue = total - ineq[sign]
			if (abs(residue) > margin and sign == '==') or (residue < -margin and sign == '>=') or (residue > margin and sign == '<='):
				if verbose:
					print ('(%.2E)'%residue, ineq_name+':', ' + '.join(['%.6f*'%ineq[var]+var for var in ineq if var != sign]), sign, ineq[sign])
				violations.append(ineq_name)

		for cone_name in self.cone_ineq:
			cone = self.cone_ineq[cone_name]
			x_var = sum([cone[1][var]*variables[var] for var in cone[1]])
			y_var = sum([cone[2][var]*variables[var] for var in cone[2]])
			z_var = sum([cone[0][var]*variables[var] for var in cone[0]])
			ratio = (x_var**2 + y_var**2 - z_var**2)/(z_var**2)
			if ratio > margin:
				if verbose:
					print ('(%.2E)'%ratio, cone_name+': (', x_var, ')² + (', y_var, ')² > (', z_var, ')²')
				cone_violations.append(cone_name)
				x_var = ' + '.join(['%.6f*'%cone[1][var]+var for var in cone[1]])
				y_var = ' + '.join(['%.6f*'%cone[2][var]+var for var in cone[2]])
				z_var = ' + '.join(['%.6f*'%cone[0][var]+var for var in cone[0]])
				if verbose:
					print (cone_name+': (', x_var, ')² + (', y_var, ')² <= (', z_var, ')²')

		for cone_name in self.cone_outeq:
			cone = self.cone_outeq[cone_name]
			x_var = sum([cone[1][var]*variables[var] for var in cone[1]])
			y_var = sum([cone[2][var]*variables[var] for var in cone[2]])
			z_var = sum([cone[0][var]*variables[var] for var in cone[0]])
			ratio = (x_var**2 + y_var**2 - z_var**2)/(z_var**2)
			if abs(ratio) > margin:
				if verbose:
					print ('(%.2E)'%ratio, cone_name+': (', x_var, ')² + (', y_var, ')² != (', z_var, ')²')
				cone_violations.append(cone_name)
				x_var = ' + '.join(['%.6f*'%cone[1][var]+var for var in cone[1]])
				y_var = ' + '.join(['%.6f*'%cone[2][var]+var for var in cone[2]])
				z_var = ' + '.join(['%.6f*'%cone[0][var]+var for var in cone[0]])
				if verbose:
					print (cone_name+': (', x_var, ')² + (', y_var, ')² == (', z_var, ')²')

		for helix_name in self.helix_cones:
			cone = self.helix_cones[helix_name]
			x_var = sum([cone[2][var]*variables[var] for var in cone[2]])
			y_var = sum([cone[3][var]*variables[var] for var in cone[3]])
			z_var = sum([cone[1][var]*variables[var] for var in cone[1]])
			ang_var = sum([cone[0][var]*variables[var] for var in cone[0]])
			ratio1 = (z_var*math.cos(ang_var) - x_var)/(z_var)
			ratio2 = (z_var*math.sin(ang_var) - y_var)/(z_var)
			if abs(ratio1) > margin or abs(ratio2) > margin:
				if verbose:
					print ('(%.2E, %.2E)'%(ratio1, ratio2), helix_name+': (', x_var, ')² + (', y_var, ')² != (', z_var, ')² ~~> Angle', ang_var)
				cone_violations.append(helix_name)
				x_var = ' + '.join(['%.6f*'%cone[2][var]+var for var in cone[2]])
				y_var = ' + '.join(['%.6f*'%cone[3][var]+var for var in cone[3]])
				z_var = ' + '.join(['%.6f*'%cone[1][var]+var for var in cone[1]])
				ang_var = ' + '.join(['%.6f*'%cone[0][var]+var for var in cone[0]])
				if verbose:
					print (helix_name+': (', x_var, ')² + (', y_var, ')² == (', z_var, ')² ~~> Angle', ang_var)
		if verbose:
			print ('%d Linear Constraints fail: '%len(violations), violations)
			print ('%d Cone Constraints fail: '%len(cone_violations), cone_violations)
		return len(violations), len(cone_violations)

	def solve(self, mode='cdc', time_limit=200, **kwargs): # Formulates the model, storing it in self.model, and solves it
		if mode is not None:
			formulator = {'milp':self.formulate_milp, 'cdc':self.formulate_cdc, 'mixed':self.formulate_mixed, 'exact':self.formulate_exact,
				'nmdt':self.formulate_nmdt, 'dnmdt':self.formulate_dnmdt, 'gnb':self.formulate_gurobi_native_bilinear, 'gnp':self.formulate_gurobi_native_polar}[mode]
			formulator(**kwargs)
		if time_limit is None:
			time_limit = 200
		time0 = time.time()
		cbdata = {}
		callback_func = partial(gp_callback, data=cbdata, timeout=time_limit)
		try:
			self.model.Params.timelimit = time_limit
			self.model.Params.FeasibilityTol = 1e-9
			self.model.optimize(callback_func)
			# self.model.optimize()
		except Exception as ex:
			print ("Solver produces error:", ex)
		final_time = time.time()-time0
		self.solution_data = cbdata
		# else:
		# 	thread1 = threading.Thread(target = solver, kwargs = kwargs)
		# 	end_time = time.time() + time_limit
		# 	thread1.start()
		# 	running = True
		# 	while thread1.is_alive():
		# 		if time.time() > end_time and running:
		# 			self.model.terminate()
		# 			running = False
		# 		time.sleep(0.1)
		# 	final_time = round(time.time() - end_time + time_limit, 1)

		stats_per_iteration = []
		if self.model.status == GRB.OPTIMAL:
			print ('{*} Status is "optimal"')
			solution = {vv: self.modelvars[vv].X for vv in self.modelvars}
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 1, 'time': final_time, 'obj': self.model.ObjVal, 'sol': solution, 'mem': memory, 'nodes': nodes})

		elif self.model.status == GRB.SUBOPTIMAL:
			print ('{*} Status is "suboptimal": proceeding as optimal')
			solution = {vv: self.modelvars[vv].X for vv in self.modelvars}
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 1, 'time': final_time, 'obj': self.model.ObjVal, 'sol': solution, 'mem': memory, 'nodes': nodes})

		elif self.model.status == 3:
			self.model.computeIIS()
			print ('{*} INFEASIBLE MODEL: ', self.name)
			self.model.write('iss_case.ilp')
			return stats_per_iteration

		elif self.model.status == GRB.INTERRUPTED or self.model.status == 9:
			print ('{!!} Stopped after %.2f seconds'%final_time)
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 0, 'time': final_time, 'obj': self.model.objVal, 'bound': self.model.objBound, 'mem': memory, 'nodes': nodes})
		else:
			print ('{*} Other Gurobi Status:', self.model.status)
		return stats_per_iteration

	def add_linear_constraint(self, name, dic):
		if name in self.linear_ineq:
			print ('Cannot add linear constraint, the constraint name "%s" is already used'%name)
		else:
			self.linear_ineq[name] = dic
			if '==' in dic:
				self.model.addConstr(gp.quicksum(dic[vv]*self.modelvars[vv] for vv in dic if vv != '==') == dic['=='])
			elif '>=' in dic:
				self.model.addConstr(gp.quicksum(dic[vv]*self.modelvars[vv] for vv in dic if vv != '>=') == dic['>='])
			elif '<=' in dic:
				self.model.addConstr(gp.quicksum(dic[vv]*self.modelvars[vv] for vv in dic if vv != '<=') == dic['<='])
			else:
				raise(Exception, 'Invalid linear inequality:', dic)

	def resolve(self, time_limit=200): # Solve the model in self.model
		if time_limit is None:
			time_limit = 200
		time0 = time.time()
		cbdata = {}
		callback_func = partial(gp_callback, data=cbdata, timeout=time_limit)
		try:
			self.model.Params.timelimit = time_limit
			self.model.Params.FeasibilityTol = 1e-9
			self.model.optimize(callback_func)
		except Exception as ex:
			print ("Solver produces error:", ex)
		final_time = time.time()-time0
		self.solution_data = cbdata

		stats_per_iteration = []
		if self.model.status == GRB.OPTIMAL:
			print ('{*} Status is "optimal"')
			solution = {vv: self.modelvars[vv].X for vv in self.modelvars}
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 1, 'time': final_time, 'obj': self.model.ObjVal, 'sol': solution, 'mem': memory, 'nodes': nodes})

		elif self.model.status == GRB.SUBOPTIMAL:
			print ('{*} Status is "suboptimal": proceeding as optimal')
			solution = {vv: self.modelvars[vv].X for vv in self.modelvars}
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 1, 'time': final_time, 'obj': self.model.ObjVal, 'sol': solution, 'mem': memory, 'nodes': nodes})

		elif self.model.status == 3:
			self.model.computeIIS()
			print ('{*} INFEASIBLE MODEL: ', self.name, self.model.status)
			self.model.write('iss_case.ilp')
			return stats_per_iteration

		elif self.model.status == GRB.INTERRUPTED or self.model.status == 9:
			print ('{!!} Stopped after %.2f seconds'%final_time)
			memory = max(self.solution_data[time]['memory'] for time in self.solution_data if 'memory' in self.solution_data[time])
			nodes = max(self.solution_data[time]['nodes'] for time in self.solution_data if 'nodes' in self.solution_data[time])
			stats_per_iteration.append( {'solved': 0, 'time': final_time, 'obj': self.model.objVal, 'bound': self.model.objBound, 'mem': memory, 'nodes': nodes})
		else:
			print ('{*} Other Gurobi Status:', self.model.status)
		return stats_per_iteration

	def formulate_milp(self, cone_depth=8, helix_depth=4, interior_depth=4, concave_iterations=5, angle_precision=None):
		m = gp.Model("Quad_cdc_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add the SOC constraints by iterated lifting
		for cone in self.cone_ineq:
			z_var, x_var, y_var, params = self.cone_ineq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			etas = [m.addVar(name='eta_0_' + cone, lb=-GRB.INFINITY)]
			zetas = [m.addVar(name='zeta_0_' + cone, lb=-GRB.INFINITY)]
			ss = math.sin(ang_min)
			cc = math.cos(ang_min)
			if 'ang_min' not in params and 'ang_max' not in params:
				m.addConstr(zetas[0] >= gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
				m.addConstr(zetas[0] >= -gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
				m.addConstr(etas[0] >= gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
				m.addConstr(etas[0] >= -gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
				angle_ii = math.pi/2
			else:
				m.addConstr(zetas[0] == cc*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)+ss*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
				m.addConstr(etas[0] == cc*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)-ss*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
				angle_ii = ang_max - ang_min
			ii = 0
			while (ii < cone_depth if cone_depth is not None else angle_ii > angle_precision):
				ii += 1
				etas.append(m.addVar(name='eta_%s_'%ii + cone))
				zetas.append(m.addVar(name='zeta_%s_'%ii + cone, lb=(-GRB.INFINITY if ii < 2 else 0)))
				angle_ii = angle_ii/2.
				cc = math.cos(angle_ii)
				ss = math.sin(angle_ii)
				m.addConstr(zetas[ii] == cc*zetas[ii-1] + ss*etas[ii-1])
				m.addConstr(etas[ii] >= -ss*zetas[ii-1] + cc*etas[ii-1])
				m.addConstr(etas[ii] >= ss*zetas[ii-1] - cc*etas[ii-1])

			m.addConstr(math.cos(angle_ii/2)*zetas[-1] + math.sin(angle_ii/2)*etas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			# m.addConstr(zetas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(etas[-1] <= math.tan(angle_ii)*zetas[-1])
			if 'include_chord_cut' in params:
				m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)*math.cos((ang_max+ang_min)/2) +
					gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)*math.sin((ang_max+ang_min)/2) >= 
					gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)*math.cos((ang_max-ang_min)/2))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		# Add the helix constraints by iterated lifting
		for helix in self.helix_cones:
			ang_var, z_var, x_var, y_var, params = self.helix_cones[helix]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			etas = [m.addVar(name='eta_0_' + helix, lb=-z_bound)]
			zetas = [m.addVar(name='zeta_0_' + helix, lb=-z_bound)]
			mus = [m.addVar(name='mu_0_' + helix, lb = 0, ub = ang_max - ang_min)]
			binaries = []
			ss = math.sin(ang_min)
			cc = math.cos(ang_min)
			m.addConstr(zetas[0] == cc*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)+ss*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(etas[0] == cc*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)-ss*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			m.addConstr(mus[0] + ang_min == gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var))
			angle_ii = ang_max - ang_min
			ii = 0 # depth
			while (ii < helix_depth if angle_precision is None else angle_ii > angle_precision):
				ii += 1
				etas.append(m.addVar(name='eta_%s_'%ii + helix))
				zetas.append(m.addVar(name='zeta_%s_'%ii + helix, lb=-z_bound if ii < 2 else 0))
				mus.append(m.addVar(name='mu_%s_'%ii + helix, lb=0))
				binaries.append(m.addVar(name='binary_%s_helix_'%ii + helix, vtype = GRB.BINARY))
				self.modelvars['binary_%s_helix_'%ii + helix] = binaries[-1]
				angle_ii = angle_ii/2
				cc = math.cos(angle_ii)
				ss = math.sin(angle_ii)
				m_constant = 2*z_bound*math.sin(min(angle_ii, math.pi/2))
				m.addConstr(zetas[ii] == cc*zetas[ii-1] + ss*etas[ii-1])
				m.addConstr(etas[ii] >= -ss*zetas[ii-1] + cc*etas[ii-1])
				m.addConstr(etas[ii] >= ss*zetas[ii-1] - cc*etas[ii-1])
				m.addConstr(etas[ii] <= -ss*zetas[ii-1] + cc*etas[ii-1] + m_constant*binaries[ii-1])
				m.addConstr(etas[ii] <= ss*zetas[ii-1] - cc*etas[ii-1] + m_constant*(1-binaries[ii-1]))
				m.addConstr(mus[ii] >= mus[ii-1] - angle_ii)
				m.addConstr(mus[ii] >= angle_ii - mus[ii-1])
				m.addConstr(mus[ii] <= mus[ii-1] - angle_ii + 2*angle_ii*binaries[ii-1])
				m.addConstr(mus[ii] <= angle_ii - mus[ii-1] + 2*angle_ii*(1-binaries[ii-1]))
			m.addConstr(math.cos(angle_ii/2)*zetas[-1] + math.sin(angle_ii/2)*etas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zetas[-1] + math.tan(angle_ii/2)*etas[-1] >= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(cc*etas[-1] <= ss*zetas[-1])
			m.addConstr(mus[-1] <= angle_ii)
			# m.addConstr(etas[-1] <= mus[-1]*z_bound)
			# if 'z_lower_bound' in params:
			# 	m.addConstr(etas[-1] >= math.sin(angle_ii)/angle_ii*params['z_lower_bound']*mus[-1])

			# self.modelvars[helix] = [zetas, etas, mus, binaries]

		# Add conical surface constraints by folding
		for cone in self.cone_outeq:
			z_var, x_var, y_var, params = self.cone_outeq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			etas = [m.addVar(name='eta_0_' + cone, lb=-z_bound)]
			zetas = [m.addVar(name='zeta_0_' + cone, lb=-z_bound)]
			binaries = []
			ss = math.sin(ang_min)
			cc = math.cos(ang_min)
			m.addConstr(zetas[0] == cc*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)+ss*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(etas[0] == cc*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)-ss*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			angle_ii = ang_max - ang_min
			ii = 0 # depth
			while (ii < helix_depth if angle_precision is None else angle_ii > angle_precision):
				ii += 1
				etas.append(m.addVar(name='eta_%s_'%ii + cone))
				zetas.append(m.addVar(name='zeta_%s_'%ii + cone, lb=-z_bound if ii < 2 else 0))
				binaries.append(m.addVar(name='binary_%s_cone_'%ii + cone, vtype = GRB.BINARY))
				angle_ii = angle_ii/2
				cc = math.cos(angle_ii)
				ss = math.sin(angle_ii)
				m_constant = 2*z_bound*math.sin(min(angle_ii, math.pi/2))
				# print (ss, cc)
				m.addConstr(zetas[ii] == cc*zetas[ii-1] + ss*etas[ii-1])
				m.addConstr(etas[ii] >= -ss*zetas[ii-1] + cc*etas[ii-1])
				m.addConstr(etas[ii] >= ss*zetas[ii-1] - cc*etas[ii-1])
				m.addConstr(etas[ii] <= -ss*zetas[ii-1] + cc*etas[ii-1] + m_constant*binaries[ii-1])
				m.addConstr(etas[ii] <= ss*zetas[ii-1] - cc*etas[ii-1] + m_constant*(1-binaries[ii-1]))
			m.addConstr(math.cos(angle_ii/2)*zetas[-1] + math.sin(angle_ii/2)*etas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zetas[-1] + math.tan(angle_ii/2)*etas[-1] >= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(cc*etas[-1] <= ss*zetas[-1])
			# self.modelvars[cone] = [zetas, etas, binaries]

		check_inside_cone = True
		self.model = m
		self.model.setParam('BarHomogeneous', 1)

	def formulate_exact(self):
		m = gp.Model("Quad_exact_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			# ss = math.sin(ang_min)
			# cc = math.cos(ang_min)

			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)
			# Add angle constraints here
			if ang_max - ang_min < 2*math.pi:
				m.addConstr(math.sin(ang_max)*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) >=
					math.cos(ang_max)*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
				m.addConstr(math.sin(ang_min)*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) <=
					math.cos(ang_min)*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
				if abs(ang_min+ang_max) < 1e-3:
					m.addConstr(gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) >= math.cos(ang_max))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		# # Add helix constraints
		
		m.setParam('BarHomogeneous', 1)
		
		# m.optimize()
		self.model = m

	def formulate_mixed(self, cone_depth=8, helix_depth=4, interior_depth=4, concave_iterations=5, angle_precision=None):
		m = gp.Model("Quad_cdc_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0

			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)

			# Add angle constraints here
			if 'include_chord_cut' in params:
				m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)*math.cos((ang_max+ang_min)/2) +
					gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)*math.sin((ang_max+ang_min)/2) >= 
					gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)*math.cos((ang_max-ang_min)/2))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		# Add helix constraints by folding
		for helix in self.helix_cones:
			ang_var, z_var, x_var, y_var, params = self.helix_cones[helix]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			etas = [m.addVar(name='eta_0_' + helix, lb=-z_bound)]
			zetas = [m.addVar(name='zeta_0_' + helix, lb=-z_bound)]
			mus = [m.addVar(name='mu_0_' + helix, lb = 0, ub = ang_max - ang_min)]
			binaries = []
			ss = math.sin(ang_min)
			cc = math.cos(ang_min)
			m.addConstr(zetas[0] == cc*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)+ss*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(etas[0] == cc*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)-ss*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			m.addConstr(mus[0] + ang_min == gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var))
			angle_ii = ang_max - ang_min
			ii = 0 # depth
			while (ii < helix_depth if angle_precision is None else angle_ii > angle_precision):
				ii += 1
				etas.append(m.addVar(name='eta_%s_'%ii + helix))
				zetas.append(m.addVar(name='zeta_%s_'%ii + helix, lb=-z_bound if ii < 2 else 0))
				mus.append(m.addVar(name='mu_%s_'%ii + helix, lb=0))
				binaries.append(m.addVar(name='binary_%s_helix_'%ii + helix, vtype = GRB.BINARY))
				binaries[-1].setAttr(GRB.Attr.BranchPriority, -ii)
				self.modelvars['binary_%s_helix_'%ii + helix] = binaries[-1]
				angle_ii = angle_ii/2
				cc = math.cos(angle_ii)
				ss = math.sin(angle_ii)
				m_constant = 2*z_bound*math.sin(min(angle_ii, math.pi/2))
				m.addConstr(zetas[ii] == cc*zetas[ii-1] + ss*etas[ii-1])
				m.addConstr(etas[ii] >= -ss*zetas[ii-1] + cc*etas[ii-1])
				m.addConstr(etas[ii] >= ss*zetas[ii-1] - cc*etas[ii-1])
				m.addConstr(etas[ii] <= -ss*zetas[ii-1] + cc*etas[ii-1] + m_constant*binaries[ii-1])
				m.addConstr(etas[ii] <= ss*zetas[ii-1] - cc*etas[ii-1] + m_constant*(1-binaries[ii-1]))
				m.addConstr(mus[ii] >= mus[ii-1] - angle_ii)
				m.addConstr(mus[ii] >= angle_ii - mus[ii-1])
				m.addConstr(mus[ii] <= mus[ii-1] - angle_ii + 2*angle_ii*binaries[ii-1])
				m.addConstr(mus[ii] <= angle_ii - mus[ii-1] + 2*angle_ii*(1-binaries[ii-1]))
			m.addConstr(math.cos(angle_ii/2)*zetas[-1] + math.sin(angle_ii/2)*etas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zetas[-1] + math.tan(angle_ii/2)*etas[-1] >= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(cc*etas[-1] <= ss*zetas[-1])
			m.addConstr(mus[-1] <= angle_ii)
			#### Advanced cut constraints ####
			# m.addConstr(etas[-1] <= mus[-1]*z_bound)
			# if 'z_lower_bound' in params:
			# 	m.addConstr(etas[-1] >= math.sin(angle_ii)/angle_ii*params['z_lower_bound']*mus[-1])
			##################################

			# self.modelvars[helix] = [zetas, etas, mus, binaries]

		# Add conical surface constraints by folding
		for cone in self.cone_outeq:
			z_var, x_var, y_var, params = self.cone_outeq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			etas = [m.addVar(name='eta_0_' + cone, lb=-z_bound)]
			zetas = [m.addVar(name='zeta_0_' + cone, lb=-z_bound)]
			binaries = []
			ss = math.sin(ang_min)
			cc = math.cos(ang_min)
			m.addConstr(zetas[0] == cc*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)+ss*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(etas[0] == cc*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)-ss*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			angle_ii = ang_max - ang_min
			ii = 0 # depth
			while (ii < helix_depth if angle_precision is None else angle_ii > angle_precision):
				ii += 1
				etas.append(m.addVar(name='eta_%s_'%ii + cone))
				zetas.append(m.addVar(name='zeta_%s_'%ii + cone, lb=-z_bound if ii < 2 else 0))
				binaries.append(m.addVar(name='binary_%s_helix_'%ii + cone, vtype = GRB.BINARY))
				binaries[-1].setAttr(GRB.Attr.BranchPriority, -ii)
				angle_ii = angle_ii/2
				cc = math.cos(angle_ii)
				ss = math.sin(angle_ii)
				m_constant = 2*z_bound*math.sin(min(angle_ii, math.pi/2))
				# print (ss, cc)
				m.addConstr(zetas[ii] == cc*zetas[ii-1] + ss*etas[ii-1])
				m.addConstr(etas[ii] >= -ss*zetas[ii-1] + cc*etas[ii-1])
				m.addConstr(etas[ii] >= ss*zetas[ii-1] - cc*etas[ii-1])
				m.addConstr(etas[ii] <= -ss*zetas[ii-1] + cc*etas[ii-1] + m_constant*binaries[ii-1])
				m.addConstr(etas[ii] <= ss*zetas[ii-1] - cc*etas[ii-1] + m_constant*(1-binaries[ii-1]))
			m.addConstr(math.cos(angle_ii/2)*zetas[-1] + math.sin(angle_ii/2)*etas[-1] <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zetas[-1] + math.tan(angle_ii/2)*etas[-1] >= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(cc*etas[-1] <= ss*zetas[-1])
			# self.modelvars[cone] = [zetas, etas, binaries]

		# for expo in self.exponentials:
		# 	x_var, y_var = self.exponentials[expo]
		# 	zetas = [m.addVar(name='zeta_0_' + expo, lb=0, ub=1)]
		# 	etas = [m.addVar(name='eta_0_' + expo, lb=1, ub=3)]
		# 	binaries = []
		# 	m.addConstr(zetas[0] == gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
		# 	m.addConstr(etas[0] == gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
		# 	for ii in range(1, helix_depth+1):
		# 		m_constant = 10
		# 		zetas.append(m.addVar(name='zeta_%s_'%ii + expo, ub=0.5**ii))
		# 		etas.append(m.addVar(name='eta_%s_'%ii + expo))
		# 		binaries.append(m.addVar(name='binary_%s_expo_'%ii + expo, vtype = GRB.BINARY))
		# 		binaries[-1].setAttr(GRB.Attr.BranchPriority, -ii)
		# 		self.modelvars['binary_%s_expo_'%ii + expo] = binaries[-1]
		# 		m.addConstr(-m_constant*(1 - binaries[-1]) <= zetas[ii] - zetas[ii-1] + 0.5**ii)
		# 		m.addConstr( zetas[ii] - zetas[ii-1] + 0.5**ii <= m_constant*(1 - binaries[-1]))
		# 		m.addConstr(-m_constant*(1 - binaries[-1]) <= etas[ii] - math.exp(-0.5**ii)*etas[ii-1])
		# 		m.addConstr( etas[ii] - math.exp(-0.5**ii)*etas[ii-1] <= m_constant*(1 - binaries[-1]))
		# 		m.addConstr(-m_constant*binaries[-1] <= zetas[ii] - zetas[ii-1])
		# 		m.addConstr( zetas[ii] - zetas[ii-1] <= m_constant*binaries[-1])
		# 		m.addConstr(-m_constant*binaries[-1] <= etas[ii] - etas[ii-1])
		# 		m.addConstr(etas[ii] - etas[ii-1] <= m_constant*binaries[-1])
		# 	e_2_v = math.exp(0.5**(helix_depth))
		# 	m.addConstr(etas[-1] - 1 <= (2**helix_depth)*(e_2_v - 1)*zetas[-1])
		# 	m.addConstr(etas[-1] >= zetas[-1] + 1)
		# 	m.addConstr(etas[-1] - e_2_v >= e_2_v*(zetas[-1] - 0.5**(helix_depth)))

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()

	def formulate_cdc(self, cone_depth=8, helix_depth=4, interior_depth=4, concave_iterations=5, angle_precision=None):
		m = gp.Model("Quad_cdc_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			# ss = math.sin(ang_min)
			# cc = math.cos(ang_min)

			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)
			# Add angle constraints here

			if 'include_chord_cut' in params:
				m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)*math.cos((ang_max+ang_min)/2) +
					gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)*math.sin((ang_max+ang_min)/2) >= 
					gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)*math.cos((ang_max-ang_min)/2))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		# Add helix constraints by cdc construction
		for helix in self.helix_cones:
			ang_var, z_var, x_var, y_var, params = self.helix_cones[helix]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			eta = m.addVar(name='eta_0_' + helix, lb=-z_bound)
			zeta = m.addVar(name='zeta_0_' + helix, lb=-z_bound)
			mu = m.addVar(name='mu_0_' + helix, lb = ang_min, ub = ang_max)
			z_variable = m.addVar(name='zvar_0_' + helix, ub = z_bound)
			weights_interior = [[m.addVar(name="weight_interior_%s_%s_"%(ii,jj) + helix, ub=1) for jj in range(4)] for ii in range(2**helix_depth)]
			weights = [[m.addVar(name="weight_shared_%s_%s_"%(ii,jj) + helix, ub=1) for jj in range(3)] for ii in range(2**helix_depth+1)]
			binaries = []
			angle_total = ang_max - ang_min
			cc = math.cos(angle_total/(2**(helix_depth+1)))
			interior_points = [ # The points in only one polyhedron
			[(math.cos(ang_min+angle_total*ii/(2**helix_depth))*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))*z_bound, z_bound, ang_min+angle_total*(ii+1)/(2**helix_depth)),
			(math.cos(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound, z_bound, ang_min+angle_total*(ii+1)/(2**helix_depth)),
			(math.cos(ang_min+angle_total*(ii+1)/(2**helix_depth))*z_bound, math.sin(ang_min+angle_total*(ii+1)/(2**helix_depth))*z_bound, z_bound, ang_min+angle_total*ii/(2**helix_depth)),
			(math.cos(ang_min+angle_total*(ii+1)/(2**helix_depth))/cc*z_bound, math.sin(ang_min+angle_total*(ii+1)/(2**helix_depth))/cc*z_bound, z_bound, ang_min+angle_total*ii/(2**helix_depth))]
			for ii in range(2**helix_depth)]
			shared_points = [ # The points shared by two polyhedra
			[(math.cos(ang_min+angle_total*ii/(2**helix_depth))*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))*z_bound, z_bound, ang_min+angle_total*ii/(2**helix_depth)),
			(math.cos(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound, z_bound, ang_min+angle_total*ii/(2**helix_depth)),
			(0, 0, 0, ang_min+angle_total*ii/(2**helix_depth))]
			for ii in range(2**helix_depth + 1)]
			for ii in range(helix_depth):
				left_points = []
				right_points = []
				left_interior_points = []
				right_interior_points = []
				for jj in range(2**ii):
					left_points += list(range(1 + jj*(2**(helix_depth-ii+1)) - 2**(helix_depth-ii-1), jj*(2**(helix_depth-ii+1)) + 2**(helix_depth-ii-1)))
					left_interior_points += list(range(1 + jj*(2**(helix_depth-ii+1))-2**(helix_depth-ii-1)-1, jj*(2**(helix_depth-ii+1))+2**(helix_depth-ii-1)))
					right_points += list(range(jj*(2**(helix_depth-ii+1)) + 2**(helix_depth-ii-1) + 1, (jj+1)*(2**(helix_depth-ii+1)) - 2**(helix_depth-ii-1) ))
					right_interior_points += list(range(jj*(2**(helix_depth-ii+1))+2**(helix_depth-ii-1)+1 - 1, (jj+1)*(2**(helix_depth-ii+1))-2**(helix_depth-ii-1)))
				left_points = [it for it in left_points if 0 <= it <= 2**helix_depth]
				left_interior_points = [it for it in left_interior_points if 0 <= it < 2**helix_depth]
				right_points = [it for it in right_points if 0 <= it <= 2**helix_depth]
				right_interior_points = [it for it in right_interior_points if 0 <= it < 2**helix_depth]
				# print(left_points, right_points)
				binaries.append(m.addVar(name='binary_%s_helix_'%ii + helix, vtype = GRB.BINARY))
				m.addConstr(gp.quicksum(gp.quicksum(weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in left_interior_points) + 
						gp.quicksum(gp.quicksum(weights[jj][ii] for ii in range(len(weights[jj]))) for jj in left_points) <= binaries[-1])
				m.addConstr(gp.quicksum(gp.quicksum(weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in right_interior_points) + 
						gp.quicksum(gp.quicksum(weights[jj][ii] for ii in range(len(weights[jj]))) for jj in right_points) <= 1-binaries[-1])

			m.addConstr(zeta == gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			m.addConstr(eta == gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(mu == gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var))
			m.addConstr(z_variable == gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zeta == gp.quicksum(gp.quicksum(shared_points[jj][ii][0]*weights[jj][ii] for ii in range(len(weights[jj]))) for jj in range(2**helix_depth+1)) + 
			gp.quicksum(gp.quicksum(interior_points[jj][ii][0]*weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in range(2**helix_depth)))
			m.addConstr(eta == gp.quicksum(gp.quicksum(shared_points[jj][ii][1]*weights[jj][ii] for ii in range(len(weights[jj]))) for jj in range(2**helix_depth+1)) + 
			gp.quicksum(gp.quicksum(interior_points[jj][ii][1]*weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in range(2**helix_depth)))
			m.addConstr(z_variable == gp.quicksum(gp.quicksum(shared_points[jj][ii][2]*weights[jj][ii] for ii in range(len(weights[jj]))) for jj in range(2**helix_depth+1)) + 
			gp.quicksum(gp.quicksum(interior_points[jj][ii][2]*weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in range(2**helix_depth)))
			m.addConstr(mu == gp.quicksum(gp.quicksum(shared_points[jj][ii][3]*weights[jj][ii] for ii in range(len(weights[jj]))) for jj in range(2**helix_depth+1)) + 
			gp.quicksum(gp.quicksum(interior_points[jj][ii][3]*weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in range(2**helix_depth)))
			#m.addConstr(eta == gp.quicksum(points_a[jj][1]*weights_a[jj] + points_b[jj][1]*weights_b[jj] for jj in range(2**helix_depth+1)))
			#m.addConstr(z_variable == gp.quicksum(points_a[jj][2]*weights_a[jj] + points_b[jj][2]*weights_b[jj] for jj in range(2**helix_depth+1)))
			#m.addConstr(mu == gp.quicksum(points_a[jj][3]*weights_a[jj] + points_b[jj][3]*weights_b[jj] for jj in range(2**helix_depth+1)))
			m.addConstr(gp.quicksum(gp.quicksum(weights_interior[jj][ii] for ii in range(len(weights_interior[jj]))) for jj in range(2**helix_depth)) +
					gp.quicksum(gp.quicksum(weights[jj][ii] for ii in range(len(weights[jj]))) for jj in range(2**helix_depth+1)) == 1)


		# Add conical surface constraints by cdc construction
		for cone in self.cone_outeq:
			z_var, x_var, y_var, params = self.cone_outeq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			angle_total = ang_max - ang_min
			z_bound = params['z_bound']
			eta = m.addVar(name='eta_' + cone, lb=-z_bound)
			zeta = m.addVar(name='zeta_' + cone, lb=-z_bound)
			weight_origin = m.addVar(name='origin_' + cone, ub=1)
			# depth = ?(helix_depth, angle_precision)
			points_a = [(math.cos(ang_min+angle_total*ii/(2**helix_depth))*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))*z_bound) for ii in range(2**helix_depth+1)]
			cc = math.cos(angle_total/(2**(helix_depth+1)))
			points_b = [(math.cos(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound, math.sin(ang_min+angle_total*ii/(2**helix_depth))/cc*z_bound) for ii in range(2**helix_depth+1)]
			weights_a = [m.addVar(name="weight_a_%s_"%ii + cone, ub=1) for ii in range(2**helix_depth+1)]
			weights_b = [m.addVar(name="weight_b_%s_"%ii + cone, ub=1) for ii in range(2**helix_depth+1)]
			binaries = []

			for ii in range(helix_depth):
				left_points = []
				right_points = []
				for jj in range(2**ii):
					left_points += list(range(1 + jj*(2**(helix_depth-ii+1)) - 2**(helix_depth-ii-1), jj*(2**(helix_depth-ii+1)) + 2**(helix_depth-ii-1)))
					right_points += list(range(jj*(2**(helix_depth-ii+1)) + 2**(helix_depth-ii-1) + 1, (jj+1)*(2**(helix_depth-ii+1)) - 2**(helix_depth-ii-1) ))
				left_points = [it for it in left_points if 0 <= it <= 2**helix_depth]
				right_points = [it for it in right_points if 0 <= it <= 2**helix_depth]
				# print(left_points, right_points)
				binaries.append(m.addVar(name='binary_%s_helix_'%ii + cone, vtype = GRB.BINARY))
				m.addConstr(gp.quicksum(weights_a[jj]+weights_b[jj] for jj in left_points) <= binaries[-1])
				m.addConstr(gp.quicksum(weights_a[jj]+weights_b[jj] for jj in right_points) <= 1-binaries[-1])

			# polytopes = []
			m.addConstr(zeta == gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var))
			m.addConstr(eta == gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			m.addConstr(z_bound*(1-weight_origin) == gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(zeta == gp.quicksum(points_a[jj][0]*weights_a[jj] + points_b[jj][0]*weights_b[jj] for jj in range(2**helix_depth+1)))
			m.addConstr(eta == gp.quicksum(points_a[jj][1]*weights_a[jj] + points_b[jj][1]*weights_b[jj] for jj in range(2**helix_depth+1)))
			m.addConstr(gp.quicksum(weights_a[jj] + weights_b[jj] for jj in range(2**helix_depth+1)) + weight_origin == 1)
			# self.modelvars[cone] = [zetas, etas, binaries]

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()

	def formulate_nmdt(self, bilinear_depth=9):
		m = gp.Model("Quad_nmdt_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			# ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			# ang_min = params['ang_min'] if 'ang_min' in params else 0
			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)
			# Add angle constraints here

		decomposed_variables = {} # Stores binary decomposition of variables that appear in several constraints
		for bilinear in self.bilinears:
			z_var, y_var, x_var, params = self.bilinears[bilinear]
			new_variable = 1
			xname = list(x_var.keys())[0]
			if len([it for it in x_var]) == 1: # one variable
				if x_var[xname] == 1: # weight 1
					new_variable = 0
			if new_variable:
				betas = [m.addVar(name='beta_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
				for ii in range(bilinear_depth):
					betas[ii].setAttr(GRB.Attr.BranchPriority, -ii)
			else:
				params['xlb'] = self.variables[xname]['lower']
				params['xub'] = self.variables[xname]['upper']
				params['ylb'] = self.variables[list(y_var.keys())[0]]['lower']
				params['yub'] = self.variables[list(y_var.keys())[0]]['upper']
				if xname in decomposed_variables:
					betas = decomposed_variables[xname]
				else:
					betas = [m.addVar(name='beta_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
					for ii in range(bilinear_depth):
						betas[ii].setAttr(GRB.Attr.BranchPriority, -ii)
					decomposed_variables[xname] = betas
			mus = [m.addVar(name='mu_%d_'%ii+bilinear, lb=-GRB.INFINITY) for ii in range(bilinear_depth)]
			deltax = m.addVar(name='deltax_'+bilinear, lb=0, ub=2**(-bilinear_depth))
			deltaz = m.addVar(name='deltaz_'+bilinear, lb=-GRB.INFINITY)
			scale = 1/(params['xub']-params['xlb'])
			m.addConstr(scale*(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) - params['xlb']) ==
				gp.quicksum(2**(-1-ii)*betas[ii] for ii in range(bilinear_depth)) + deltax)
			m.addConstr(scale*(gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)) == gp.quicksum(2**(-1-ii)*mus[ii] for ii in range(bilinear_depth))
				 + deltaz + scale*params['xlb']*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
			for ii in range(bilinear_depth):
				# McCormick envelope for y*beta_j=mu_j
				m.addConstr(mus[ii] <= sum(y_var[vv]*self.modelvars[vv] for vv in y_var) + params['ylb']*betas[ii] - params['ylb'])
				m.addConstr(mus[ii] >= params['ylb']*betas[ii])
				m.addConstr(mus[ii] <= params['yub']*betas[ii])
				m.addConstr(mus[ii] >= sum(y_var[vv]*self.modelvars[vv] for vv in y_var) + params['yub']*betas[ii] - params['yub'])
			# McCormick for deltax*y = deltaz
			m.addConstr(deltaz <= 2**(-bilinear_depth)*sum(y_var[vv]*self.modelvars[vv] for vv in y_var) + params['ylb']*deltax - 2**(-bilinear_depth)*params['ylb'])
			m.addConstr(deltaz >= params['ylb']*deltax)
			m.addConstr(deltaz <= params['yub']*deltax)
			m.addConstr(deltaz >= 2**(-bilinear_depth)*sum(y_var[vv]*self.modelvars[vv] for vv in y_var) + params['yub']*deltax - 2**(-bilinear_depth)*params['yub'])

			# m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			# m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
			# 	gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)
			# Add angle constraints here

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()

	def formulate_dnmdt(self, bilinear_depth=5):
		m = gp.Model("Quad_dnmdt_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			# ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			# ang_min = params['ang_min'] if 'ang_min' in params else 0
			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)
			# Add angle constraints here

		decomposed_variables = {} # Stores binary decomposition of variables that appear in several constraints
		dnmdt_lambda = 0.5
		for bilinear in self.bilinears:
			z_var, y_var, x_var, params = self.bilinears[bilinear]
			new_variable = 1
			xname = list(x_var.keys())[0]
			yname = list(y_var.keys())[0]
			if len([it for it in x_var]) == 1 and len([it for it in y_var]): # one variable
				if x_var[xname] == 1 and y_var[yname] == 1: # weight 1
					new_variable = 0
			if new_variable:
				betas_x = [m.addVar(name='beta_x_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
				betas_y = [m.addVar(name='beta_y_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
				for ii in range(bilinear_depth):
					betas_x[ii].setAttr(GRB.Attr.BranchPriority, -ii)
					betas_y[ii].setAttr(GRB.Attr.BranchPriority, -ii)
			else:
				params['xlb'] = self.variables[xname]['lower']
				params['xub'] = self.variables[xname]['upper']
				params['ylb'] = self.variables[yname]['lower']
				params['yub'] = self.variables[yname]['upper']
				if xname in decomposed_variables:
					betas_x = decomposed_variables[xname]
				else:
					betas_x = [m.addVar(name='beta_x_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
					for ii in range(bilinear_depth):
						betas_x[ii].setAttr(GRB.Attr.BranchPriority, -ii)
					decomposed_variables[xname] = betas_x
				if yname in decomposed_variables:
					betas_y = decomposed_variables[yname]
				else:
					betas_y = [m.addVar(name='beta_y_%d_'%ii+bilinear, vtype=GRB.BINARY) for ii in range(bilinear_depth)]
					for ii in range(bilinear_depth):
						betas_y[ii].setAttr(GRB.Attr.BranchPriority, -ii)
					decomposed_variables[yname] = betas_y
			mus = [m.addVar(name='mu_%d_'%ii+bilinear) for ii in range(bilinear_depth)]
			nus = [m.addVar(name='nu_%d_'%ii+bilinear) for ii in range(bilinear_depth)]
			deltax = m.addVar(name='deltax_'+bilinear, lb=0, ub=2**(-bilinear_depth))
			deltay = m.addVar(name='deltay_'+bilinear, lb=0, ub=2**(-bilinear_depth))
			deltaz = m.addVar(name='deltaz_'+bilinear)
			scale_x = 1/(params['xub']-params['xlb'])
			scale_y = 1/(params['yub']-params['ylb'])
			m.addConstr(scale_x*(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) - params['xlb']) ==
				gp.quicksum(2**(-1-ii)*betas_x[ii] for ii in range(bilinear_depth)) + deltax)
			m.addConstr(scale_y*(gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) - params['ylb']) ==
				gp.quicksum(2**(-1-ii)*betas_y[ii] for ii in range(bilinear_depth)) + deltay)
			m.addConstr(gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var) == params['xlb']*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)
				 +(params['xub'] - params['xlb'])*(params['yub'] - params['ylb'])*(deltaz + gp.quicksum(2**(-1-ii)*(mus[ii] + nus[ii]) for ii in range(bilinear_depth)))
				 + params['ylb']*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) - params['xlb']*params['ylb'])
			ub_y_side = (1-dnmdt_lambda) + dnmdt_lambda*2**(-bilinear_depth)
			ub_x_side = dnmdt_lambda + (1-dnmdt_lambda)*2**(-bilinear_depth)
			for ii in range(bilinear_depth):
				# McCormick envelope for ((1-lambda)*y + lambda*delta_y)*beta_x_j=mu_j
				m.addConstr(mus[ii] <= dnmdt_lambda*deltay + (1-dnmdt_lambda)*scale_y*(sum(y_var[vv]*self.modelvars[vv] for vv in y_var)-params['ylb']))
				# m.addConstr(mus[ii] >= 0)
				m.addConstr(mus[ii] <= ub_y_side*betas_x[ii])
				m.addConstr(mus[ii] >= dnmdt_lambda*deltay + (1-dnmdt_lambda)*scale_y*(sum(y_var[vv]*self.modelvars[vv] for vv in y_var)-params['ylb']) + ub_y_side*betas_x[ii] - ub_y_side)
				# McCormick envelope for (lambda*x+(1-lambda)delta_x)beta_y_j=nu_j
				m.addConstr(nus[ii] <= dnmdt_lambda*deltax + (1-dnmdt_lambda)*scale_x*(sum(x_var[vv]*self.modelvars[vv] for vv in x_var)-params['xlb']))
				# m.addConstr(nus[ii] >= 0)
				m.addConstr(nus[ii] <= ub_x_side*betas_y[ii])
				m.addConstr(nus[ii] >= dnmdt_lambda*deltax + (1-dnmdt_lambda)*scale_x*(sum(x_var[vv]*self.modelvars[vv] for vv in x_var)-params['xlb']) + ub_x_side*betas_y[ii] - ub_x_side)
			# McCormick for deltax*deltay = deltaz
			m.addConstr(deltaz <= 2**(-bilinear_depth)*deltay)
			# m.addConstr(deltaz >= 0)
			m.addConstr(deltaz <= 2**(-bilinear_depth)*deltax)
			m.addConstr(deltaz >= 2**(-bilinear_depth)*deltay + 2**(-bilinear_depth)*deltax - 4**(-bilinear_depth))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()

	def formulate_gurobi_native_bilinear(self):
		m = gp.Model("Quad_gnb_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)

		# Add bilinear constraints explicitly
		for bilinear in self.bilinears:
			z_var, y_var, x_var, params = self.bilinears[bilinear]
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) ==
				gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()

	def formulate_gurobi_native_polar(self):
		m = gp.Model("Quad_gnp_"+self.name, env=my_env)
		self.modelvars = {}
		for variable in self.variables:
			bounds = self.variables[variable]
			lb = float('-inf') if 'lower' not in bounds else bounds['lower']
			ub = float('inf') if 'upper' not in bounds else bounds['upper']
			self.modelvars[variable] = m.addVar(name=variable, lb=lb, ub=ub)

		for lineq_name in self.linear_ineq:
			lineq = self.linear_ineq[lineq_name]
			if '==' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '==']) == lineq['=='])
			elif '>=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '>=']) >= lineq['>='])
			elif '<=' in lineq:
				m.addConstr(gp.quicksum([lineq[vv]*self.modelvars[vv] for vv in lineq if vv != '<=']) <= lineq['<='])

		# Add SOC constraints explicitly
		for cone in self.cone_ineq:
			z_var, y_var, x_var, params = self.cone_ineq[cone]
			m.addConstr(0 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 + 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 <= gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2)

		for cone in self.cone_outeq:
			z_var, x_var, y_var, params = self.cone_outeq[cone]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			# Make sure y = c_i-c_j and z = c_i+c_j
			# print(z_var, x_var, y_var, params)

			# auxv = m.addVar(name = 'aux_outcone_'+cone) # Auxiliary variable helps the solver?
			# m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2 == auxv)
			# m.addConstr((gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)+gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)) * 
			# 	(gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)-gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)) == auxv)

			# m.addConstr((gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)+gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)) * 
			# 	(gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)-gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)) == 
			# 	gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2)
			m.addConstr(gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)**2 == 
				gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var)**2 + gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var)**2)

			if ang_max - ang_min < 2*math.pi:
				m.addConstr(math.sin(ang_max)*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) >=
					math.cos(ang_max)*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
				m.addConstr(math.sin(ang_min)*gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) <=
					math.cos(ang_min)*gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var))
				# if abs(ang_min+ang_max) < 1e-3:
				# 	m.addConstr(gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) >=
				# 		math.cos(ang_max)*gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var))

		for helix in self.helix_cones:
			ang_var, z_var, x_var, y_var, params = self.helix_cones[helix]
			ang_max = params['ang_max'] if 'ang_max' in params else 2*math.pi
			ang_min = params['ang_min'] if 'ang_min' in params else 0
			z_bound = params['z_bound']
			m.addConstr(gp.quicksum(y_var[vv]*self.modelvars[vv] for vv in y_var) ==
				gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)*nlfunc.sin(gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var)))
			m.addConstr(gp.quicksum(x_var[vv]*self.modelvars[vv] for vv in x_var) ==
				gp.quicksum(z_var[vv]*self.modelvars[vv] for vv in z_var)*nlfunc.cos(gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var)))
			m.addConstr(gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var) <= ang_max)
			m.addConstr(gp.quicksum(ang_var[vv]*self.modelvars[vv] for vv in ang_var) >= ang_min)

		m.setObjective(gp.quicksum(self.objective[var_name]*self.modelvars[var_name] for var_name in self.objective), GRB.MINIMIZE)

		m.setParam('BarHomogeneous', 1)
		
		self.model = m
		# m.optimize()


class Power_Flow_Problem:
	"""
		This class is for instances of the Optimal Power Flow problem, defined by the
		network parameters (bus, branch and generator data).

		For each of the formulations of OPF and its relaxations, there is a method that
		returns an instance of QuadProblem
	"""
	def __init__(self, name, base_MVA, bus_data, gen_data, branch_data):
		self.name = name
		self.nodes = [it for it in bus_data]
		# print(self.nodes)
		self.shunt_g = {node: 0 for node in self.nodes}
		self.shunt_b = {node: 0 for node in self.nodes}
		self.edges, self.line_data = [], {}
		for branch in branch_data:
			aa, bb, ll = '%d'%int(branch_data[branch]['fbus']), '%d'%int(branch_data[branch]['tbus']), '%d'%int(branch_data[branch]['branch_id'])
			if (aa,bb) not in self.edges:
				self.edges.append((aa,bb))
			rr, xx, charge = branch_data[branch]['r'], branch_data[branch]['x'], branch_data[branch]['b']
			ratio, angle = branch_data[branch]['ratio'], branch_data[branch]['angle']
			if ratio == 0.0:
				ratio = 1
			self.line_data[(ll, aa, bb)] = {'charge': charge, 'ratio': ratio, 'angle': angle, 'adm_g': rr/(xx*xx+rr*rr), 'adm_b': -xx/(xx*xx+rr*rr),
				'ang_max': branch_data[branch]['angmax'], 'ang_min': branch_data[branch]['angmin'], 'rateA': branch_data[branch]['rateA']/base_MVA}
		for node in self.nodes:
			# print (bus_data[node])
			self.shunt_g[node] = bus_data[node]['Gs']/base_MVA
			self.shunt_b[node] = bus_data[node]['Bs']/base_MVA

		self.weight_constant = {it: float(gen_data[it]['c_0'].split(';')[0]) for it in gen_data}
		self.weight_linear = {it: float(gen_data[it]['c_1'])*base_MVA for it in gen_data}
		self.weight_quad = {it: float(gen_data[it]['c_2'])*base_MVA*base_MVA for it in gen_data}
		# print ('baseKV:', [bus_data[it]['baseKV'] for it in self.nodes])
		self.voltage_u = {it: bus_data[it]['Vmax'] for it in self.nodes}
		self.voltage_l = {it: bus_data[it]['Vmin'] for it in self.nodes}
		self.p_demand = {it: bus_data[it]['Pd']/base_MVA for it in self.nodes}
		self.q_demand = {it: bus_data[it]['Qd']/base_MVA for it in self.nodes}
		self.generators = {it: '%.0f'%gen_data[it]['bus'] for it in gen_data}

		self.p_max = {it: gen_data[it]['Pmax']/base_MVA for it in gen_data}
		self.p_min = {it: gen_data[it]['Pmin']/base_MVA for it in gen_data}
		self.q_max = {it: gen_data[it]['Qmax']/base_MVA for it in gen_data}
		self.q_min = {it: gen_data[it]['Qmin']/base_MVA for it in gen_data}

	def helix_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
			variables['theta_%s'%node] = {} # Set bounds to zero for reference node
		for ii, jj in self.edges:
			variables['c_%s_%s'%(ii, jj)] = {}
			variables['s_%s_%s'%(ii, jj)] = {}
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		cones = {}
		helixes = {}
		outer_cones = {}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
			constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
				f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
				f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]

		for ii, jj in self.edges:
			if (jj, ii) in self.edges:
				angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
			cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1,
				'z_bound': self.voltage_u[ii]*self.voltage_u[jj], 'z_lower_bound': self.voltage_l[ii]*self.voltage_l[jj]}
			cone2_params = {'ang_max': math.pi/2, 'ang_min': -math.pi/2, 'z_bound': self.voltage_u[ii]**2 + self.voltage_u[jj]**2}
			helixes['jabr_cone1_%s_%s'%(ii,jj)] = [{'theta_%s'%jj:1, 'theta_%s'%ii:-1}, {'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
			outer_cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
		
		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, outer_cones, helixes, objective)

	def soc_real_polar_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
			variables['s_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
		for ii, jj in self.edges:
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for jj in range(len(self.nodes)):
			for ii in range(len(self.nodes)):
				variables['ef_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
				if ii<jj:
					variables['ee_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
					variables['ff_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		cones = {}
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			if int(ii)<int(jj):
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}
			else:
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[0], f's_{jj}_{ii}': coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[1], f's_{jj}_{ii}': -coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{jj}_{ii}': -coef_ji[0], f's_{jj}_{ii}': -coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{jj}_{ii}': -coef_ji[1], f's_{jj}_{ii}': coef_ji[0], '==': 0.0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]
		for ii, jj in self.edges:
			if (jj, ii) in self.edges:
				angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
			cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1}
			cone2_params = {} #{'ang_max': math.pi/2, 'ang_min': -math.pi/2}
			if int(ii) < int(jj):
				cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
				cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
			else:
				cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(jj,ii):1}, {'s_%s_%s'%(jj,ii):1}, cone1_params]
				cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]

		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, {}, {}, objective)

	def soc_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
		for ii, jj in self.edges:
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for jj in range(len(self.nodes)):
			for ii in range(jj):
				variables['c_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
				variables['s_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		cones = {}
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			if int(ii)<int(jj):
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}
				# Angle constraints, should be implicit in cone constraints
				# constraints['angle_max_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': math.sin(self.line_data[line]['ang_max']*math.pi/180), 
				# 	f's_{ii}_{jj}': -math.cos(self.line_data[line]['ang_max']*math.pi/180), '>=':0}
				# constraints['angle_min_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': -math.sin(self.line_data[line]['ang_min']*math.pi/180), 
				# 	f's_{ii}_{jj}': math.cos(self.line_data[line]['ang_min']*math.pi/180), '>=':0}
				# constraints['chord_constraint_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': 1, 
				# 	'aux_line_%s_%s'%(ii,jj): -math.cos(math.pi/6), '>=':0}
			else:
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[0], f's_{jj}_{ii}': coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[1], f's_{jj}_{ii}': -coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{jj}_{ii}': -coef_ji[0], f's_{jj}_{ii}': -coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{jj}_{ii}': -coef_ji[1], f's_{jj}_{ii}': coef_ji[0], '==': 0.0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]
		for ii, jj in self.edges:
			if (jj, ii) in self.edges:
				angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
			cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1}
			cone2_params = {} #{'ang_max': math.pi/2, 'ang_min': -math.pi/2}
			if int(ii) < int(jj):
				cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
				cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
			else:
				cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(jj,ii):1}, {'s_%s_%s'%(jj,ii):1}, cone1_params]
				cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]

		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, {}, {}, objective)

	def bounds_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
		for ii, jj in self.edges:
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for jj in range(len(self.nodes)):
			for ii in range(jj):
				variables['c_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
				variables['s_%s_%s'%(self.nodes[ii], self.nodes[jj])] = {}
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		cones = {}
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			if int(ii)<int(jj):
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}
				# Angle constraints, should be implicit in cone constraints
				constraints['angle_max_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': math.sin(angle_bounds[(ii,jj)][0]*math.pi/180), 
					f's_{ii}_{jj}': -math.cos(self.line_data[line]['ang_max']*math.pi/180), '>=':0}
				constraints['angle_min_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': -math.sin(angle_bounds[(ii,jj)][1]*math.pi/180), 
					f's_{ii}_{jj}': math.cos(self.line_data[line]['ang_min']*math.pi/180), '>=':0}
				constraints['chord_constraint_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': 1, 
					'aux_line_%s_%s'%(ii,jj): -math.cos(math.pi/6), '>=':0}
			else:
				constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[0], f's_{jj}_{ii}': coef_ij[1], f'==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
					f'c_{jj}_{ii}': -coef_ij[1], f's_{jj}_{ii}': -coef_ij[0], '==': 0.0}
				constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
					f'c_{jj}_{ii}': -coef_ji[0], f's_{jj}_{ii}': -coef_ji[1], '==': 0.0}
				constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
					f'c_{jj}_{ii}': -coef_ji[1], f's_{jj}_{ii}': coef_ji[0], '==': 0.0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]

		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, {}, {}, objective)

	def tree_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
			variables['theta_%s'%node] = {} # Set bounds to zero for reference node
		for ii, jj in self.edges:
			variables['c_%s_%s'%(ii, jj)] = {}
			variables['s_%s_%s'%(ii, jj)] = {}
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		cones = {}
		helixes = {}
		outer_cones = {}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
			constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
				f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
				f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]
		for ii, jj in self.edges:
			if (jj, ii) in self.edges:
				angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
			cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1,
				'z_bound': self.voltage_u[ii]*self.voltage_u[jj], 'z_lower_bound': self.voltage_l[ii]*self.voltage_l[jj]}
			cone2_params = {'ang_max': math.pi/2, 'ang_min': -math.pi/2, 'z_bound': self.voltage_u[ii]**2 + self.voltage_u[jj]**2}
			outer_cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
			outer_cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
		
		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, outer_cones, helixes, objective)

	def jabr_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
		for ii, jj in self.edges:
			variables['c_%s_%s'%(ii, jj)] = {}
			variables['s_%s_%s'%(ii, jj)] = {}
			variables['aux_line_%s_%s'%(ii, jj)] = {'lower': 0} # Artificial variable for the cone decomposition
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic
		cones = {}
		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
			constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
				f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
				f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}
			# Angle constraints, should be implicit in cone constraints
			# constraints['angle_max_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': math.sin(self.line_data[line]['ang_max']*math.pi/180), 
			# 	f's_{ii}_{jj}': -math.cos(self.line_data[line]['ang_max']*math.pi/180), '>=':0}
			# constraints['angle_min_line_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': -math.sin(self.line_data[line]['ang_min']*math.pi/180), 
			# 	f's_{ii}_{jj}': math.cos(self.line_data[line]['ang_min']*math.pi/180), '>=':0}
			# constraints['chord_constraint_%s_%s_%s'%(ll,jj,ii)] = {f'c_{ii}_{jj}': 1, 
			# 	'aux_line_%s_%s'%(ii,jj): -math.cos(math.pi/6), '>=':0}

			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]
		for ii, jj in self.edges:
			if (jj, ii) in self.edges:
				angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
			cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1}
			cone2_params = {} #{'ang_max': math.pi/2, 'ang_min': -math.pi/2}
			cones['jabr_cone1_%s_%s'%(ii,jj)] = [{'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
			cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
		outer_cones = {}
		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, outer_cones, objective)

	def rectangular_formulation(self):
		variables = {'constant_1': {'lower': 1, 'upper': 1}}
		objective = {'constant_1': sum(self.weight_constant[it] for it in self.weight_constant)}

		for gen in self.generators:
			variables['pg_%s'%gen] = {'lower':self.p_min[gen], 'upper':self.p_max[gen]}
			objective['pg_%s'%gen] = self.weight_linear[gen]
			variables['qg_%s'%gen] = {'lower':self.q_min[gen], 'upper':self.q_max[gen]}
		for node in self.nodes:
			variables['c_%s'%node] = {'lower':self.voltage_l[node]**2, 'upper':self.voltage_u[node]**2}
			# variables['theta_%s'%node] = {} # Set bounds to zero for reference node
			variables['e_%s'%node] = {'lower': -self.voltage_l[node], 'upper':self.voltage_l[node]}
			variables['f_%s'%node] = {'lower': -self.voltage_l[node], 'upper':self.voltage_l[node]}
		for ii, jj in self.edges:
			variables['c_%s_%s'%(ii, jj)] = {}
			variables['s_%s_%s'%(ii, jj)] = {}
			variables['ee_%s_%s'%(ii, jj)] = {}
			variables['ef_%s_%s'%(ii, jj)] = {}
			variables['ef_%s_%s'%(jj, ii)] = {}
			variables['ff_%s_%s'%(ii, jj)] = {}
		for ll, ii, jj in self.line_data:
			variables['flow_p_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_q_%s_%s_%s'%(ll, ii, jj)] = {}
			variables['flow_p_%s_%s_%s'%(ll, jj, ii)] = {}
			variables['flow_q_%s_%s_%s'%(ll, jj, ii)] = {}
		constraints = {}
		for ii in self.nodes:
			p_dic = {'c_%s'%ii: -self.shunt_g[ii]}
			q_dic = {'c_%s'%ii: self.shunt_b[ii]}
			for gen in self.generators:
				if self.generators[gen] == ii:
					p_dic['pg_%s'%gen] = 1
					q_dic['qg_%s'%gen] = 1
			for ll, aa, bb in self.line_data:
				if aa == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, aa, bb)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, aa, bb)] = -1
				if bb == ii:
					p_dic['flow_p_%s_%s_%s'%(ll, bb, aa)] = -1
					q_dic['flow_q_%s_%s_%s'%(ll, bb, aa)] = -1

			p_dic['=='] = self.p_demand[ii]
			q_dic['=='] = self.q_demand[ii]
			constraints['p_demand_%s'%ii] = p_dic
			constraints['q_demand_%s'%ii] = q_dic

		angle_bounds = {(ii,jj): (180, -180) for ii, jj in self.edges}
		cones = {}
		helixes = {}
		outer_cones = {}
		bilinears = {}
		for line in self.line_data:
			ll, ii, jj = line
			bb, ratio, angle = self.line_data[line]['charge'], self.line_data[line]['ratio'], self.line_data[line]['angle']
			t_real, t_im = ratio*math.cos(angle/180*math.pi), ratio*math.sin(angle/180*math.pi)
			adm_g, adm_b = self.line_data[line]['adm_g'], self.line_data[line]['adm_b']
			coef_ij = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, t_im))
			coef_ji = cplx_prod((adm_g, -adm_b), cplx_inv(t_real, -t_im))
			angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], self.line_data[line]['ang_max']), max(angle_bounds[(ii,jj)][1], self.line_data[line]['ang_min']))

			constraints['p_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_p_{ll}_{ii}_{jj}': -1, f'c_{ii}': adm_g/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[0], f's_{ii}_{jj}': -coef_ij[1], f'==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,ii,jj)] = {f'flow_q_{ll}_{ii}_{jj}': -1, f'c_{ii}': (-adm_b - bb/2)/(ratio*ratio),
				f'c_{ii}_{jj}': -coef_ij[1], f's_{ii}_{jj}': coef_ij[0], '==': 0.0}
			constraints['p_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_p_{ll}_{jj}_{ii}': -1, f'c_{jj}': adm_g,
				f'c_{ii}_{jj}': -coef_ji[0], f's_{ii}_{jj}': coef_ji[1], '==': 0.0}
			constraints['q_flow_line_%s_%s_%s'%(ll,jj,ii)] = {f'flow_q_{ll}_{jj}_{ii}': -1, f'c_{jj}': -adm_b - bb/2,
				f'c_{ii}_{jj}': -coef_ji[1], f's_{ii}_{jj}': -coef_ji[0], '==': 0.0}

			constraints[f'polar_c_{ll}_{ii}_{jj}'] = {f'c_{ii}_{jj}': -1, 'ff_%s_%s'%(ii, jj): 1, 'ee_%s_%s'%(ii, jj): 1, '==': 0.0}
			constraints[f'polar_s_{ll}_{ii}_{jj}'] = {f's_{ii}_{jj}': -1, 'ef_%s_%s'%(ii, jj): -1, 'ef_%s_%s'%(jj, ii): 1, '==': 0.0}
			constraints[f'min_ang_{ll}_{ii}_{jj}'] = {f'c_{ii}_{jj}': -math.sin(angle_bounds[(ii, jj)][1]*math.pi/180),
				f's_{ii}_{jj}': math.cos(angle_bounds[(ii, jj)][1]*math.pi/180), '>=': 0.0}
			constraints[f'max_ang_{ll}_{ii}_{jj}'] = {f'c_{ii}_{jj}': -math.sin(angle_bounds[(ii, jj)][0]*math.pi/180),
				f's_{ii}_{jj}': math.cos(angle_bounds[(ii, jj)][0]*math.pi/180), '<=': 0.0}
			cones['line_capacity_%s_%s_%s'%(ll,jj,ii)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{jj}_{ii}': 1}, {f'flow_q_{ll}_{jj}_{ii}': 1}, {}]
			cones['line_capacity_%s_%s_%s'%(ll,ii,jj)] = [{'constant_1': self.line_data[line]['rateA']}, {f'flow_p_{ll}_{ii}_{jj}': 1}, {f'flow_q_{ll}_{ii}_{jj}': 1}, {}]
		for ii, jj in self.edges:
			bilinears['bil1_%s_%s'%(ii, jj)] = [{'ee_%s_%s'%(ii, jj): 1}, {'e_%s'%ii: 1}, {'e_%s'%jj: 1}, {}]
			bilinears['bil2_%s_%s'%(ii, jj)] = [{'ef_%s_%s'%(ii, jj): 1}, {'e_%s'%ii: 1}, {'f_%s'%jj: 1}, {}]
			bilinears['bil3_%s_%s'%(ii, jj)] = [{'ef_%s_%s'%(jj, ii): 1}, {'e_%s'%jj: 1}, {'f_%s'%ii: 1}, {}]
			bilinears['bil4_%s_%s'%(ii, jj)] = [{'ff_%s_%s'%(ii, jj): 1}, {'f_%s'%ii: 1}, {'f_%s'%jj: 1}, {}]
		for ii in self.nodes: # use univariate parabola instead
			variables['ee_%s'%ii] = {}
			variables['ff_%s'%ii] = {}
			bilinears['sq1_%s'%ii] = [{'ee_%s'%ii: 1}, {'e_%s'%ii: 1}, {'e_%s'%ii: 1}, {}]
			bilinears['sq2_%s'%ii] = [{'ff_%s'%ii: 1}, {'f_%s'%ii: 1}, {'f_%s'%ii: 1}, {}]
			constraints[f'polar_c_{ii}'] = {f'c_{ii}': -1, f'ee_{ii}': 1, f'ff_{ii}': 1, '==': 0}
		# 	if (jj, ii) in self.edges:
		# 		angle_bounds[(ii,jj)] = (min(angle_bounds[(ii,jj)][0], -angle_bounds[(jj,ii)][1]), max(angle_bounds[(ii,jj)][1], -angle_bounds[(jj,ii)][0]))
		# 	cone1_params = {'ang_max': angle_bounds[(ii, jj)][0]*math.pi/180, 'ang_min': angle_bounds[(ii, jj)][1]*math.pi/180, 'include_chord_cut': 1,
		# 		'z_bound': self.voltage_u[ii]*self.voltage_u[jj], 'z_lower_bound': self.voltage_l[ii]*self.voltage_l[jj]}
		# 	cone2_params = {'ang_max': math.pi/2, 'ang_min': -math.pi/2, 'z_bound': self.voltage_u[ii]**2 + self.voltage_u[jj]**2}
		# 	helixes['jabr_cone1_%s_%s'%(ii,jj)] = [{'theta_%s'%jj:1, 'theta_%s'%ii:-1}, {'aux_line_%s_%s'%(ii,jj):1}, {'c_%s_%s'%(ii,jj):1}, {'s_%s_%s'%(ii,jj):1}, cone1_params]
		# 	outer_cones['jabr_cone2_%s_%s'%(ii,jj)] = [{'c_%s'%ii:1, 'c_%s'%jj:1}, {'aux_line_%s_%s'%(ii,jj):2}, {'c_%s'%ii:1, 'c_%s'%jj:-1}, cone2_params]
		
		par_sq_cone = {'ang_max': math.pi/2, 'ang_min': -math.pi/2} # Parameters for the cone used for y=x^2, here x >= 0 so we can ignore half the plane
		for gen in self.weight_quad:
			if self.weight_quad[gen] > 0.0:
				variables['pg_sq_%s'%gen] = {'lower': 0}
				cones['square_objective_%s'%gen] = [{'pg_sq_%s'%gen:1, 'constant_1':1}, {'pg_%s'%gen:2}, {'pg_sq_%s'%gen:1, 'constant_1':-1}, par_sq_cone]
				objective['pg_sq_%s'%gen] = self.weight_quad[gen]
			if self.weight_quad[gen] < 0.0:
				print ('Objective function is not convex! c_2[%s] = '%str(gen), self.weight_quad[gen])
		return QuadProblem(self.name, variables, constraints, cones, outer_cones, helixes, objective, bilinears = bilinears)

def read_pglib(filename, include_0_status = False):
	# include_0_status: remove "offline" generators (with status=0) if set to False
	with open(filename, 'r') as file:
		lines = file.read().split('\n')

	bus_data = {}
	ii = lines.index('mpc.bus = [')
	base_MVA = float([line for line in lines[:ii] if line[:14] == 'mpc.baseMVA = '][0][14:-2])
	bus_lines = lines[ii+1:ii+lines[ii:].index('];')]
	bus_entries = 'bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin'.split(' ')

	for line in bus_lines:
		dic = {aa:float(bb) for aa,bb in zip(bus_entries, line[1:line.find(';')].replace(' ', '').split('\t'))}
		bus_data['%.0f'%dic['bus_i']] = dic

	gen_data = {}
	ii = lines.index('mpc.gen = [')
	gen_entries = 'bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin'.split(' ')
	gen_lines = lines[ii+1:ii+lines[ii:].index('];')]

	for ii in range(len(gen_lines)):
		line = gen_lines[ii]
		generator_data = {aa:float(bb) for aa,bb in zip(gen_entries, line[1:line.find(';')].replace(' ', '').split('\t'))}
		if include_0_status or generator_data['status'] == float('1'):
			gen_data[ii] = generator_data
	ii = lines.index('mpc.gencost = [')
	gen_cost_lines = lines[ii+1:ii+lines[ii:].index('];')]
	gen_cost_entries = 'constant_2 startup shutdown n c_2 c_1 c_0'.split(' ')

	for ii in gen_data:
		line = gen_cost_lines[ii]
		for aa,bb in zip(gen_cost_entries, line[1:line.find(';')].replace(' ', '').split('\t')):
			gen_data[ii][aa] = bb
	gen_data = {'%.0f'%ii :gen_data[ii] for ii in gen_data}

	branch_data = {}
	ii = lines.index('mpc.branch = [')
	branch_lines = lines[ii+1:ii+lines[ii:].index('];')]
	branch_entries = 'fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax'.split(' ')

	for line, idx in zip(branch_lines, range(len(branch_lines))):
		dic = {aa:float(bb) for aa,bb in zip(branch_entries, line[1:line.find(';')].replace(' ', '').split('\t'))}
		if include_0_status or dic['status'] == float('1'):
			branch_data[idx] = {'branch_id': idx+1, **dic}
			
	return base_MVA, bus_data, gen_data, branch_data

