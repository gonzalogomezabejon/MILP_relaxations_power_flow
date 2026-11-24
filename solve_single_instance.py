
import cdc_formulation as form
import time, sys
import gurobipy as gp
from gurobipy import GRB

# This script solves a single instance of the OPF problem using the implemented relaxations in gurobi

def callback(model, where):
	if where != GRB.Callback.POLLING:
		runtime = model.cbGet(GRB.Callback.RUNTIME)
		if runtime > 60:
			model.terminate()

if __name__ == '__main__':
	# model = gp.read('model.mps')
	# model.optimize(callback)
	# sys.exit()
	filenm = 'pglib-opf/sad/pglib_opf_case24_ieee_rts__sad.m'
	filenm = 'pglib-opf/sad/pglib_opf_case30_ieee__sad.m'
	filenm = 'pglib-opf/sad/pglib_opf_case60_c__sad.m'
	filenm = 'nesta-0.7.0/opf/sad/nesta_case29_edin__sad.m'
	filenm = 'pglib-opf/sad/pglib_opf_case3_lmbd__sad.m'
	filenm = 'pglib-opf/pglib_opf_case5_pjm.m'
	base_MVA, bus_data, gen_data, branch_data = form.read_pglib(filenm)
	obj = form.Power_Flow_Problem("Instance", base_MVA, bus_data, gen_data, branch_data)
	# quad_problem = obj.helix_formulation()
	quad_problem = obj.rectangular_formulation()
	time0 = time.time()
	# lp_results = quad_problem.solve(mode='mixed', time_limit=3600, cone_depth=7, interior_depth=1, concave_iterations=1, helix_depth=12)
	lp_results = quad_problem.solve(mode='gnb', time_limit=3600)
	quad_problem.model.write('model.mps')
	# print(lp_results)
	print ('Time: ', time.time()-time0)
	for it in lp_results[0]:
		if it != 'sol':
			print (it, '->', lp_results[0][it])
	case_type = 'sad' if '__sad' in filenm else ('api' if '__api' in filenm else 'typ')
	baseline = 'results/results_julia_'+case_type+'.txt'
	with open(baseline, 'r') as f:
		lines = f.read().split('\n')
	for line in lines:
		if line.split('\t')[0] in filenm:
			gap = 100. - 100*lp_results[0]['obj']/float(line.split('\t')[1]) if lp_results[0]['solved'] else 100. - 100*lp_results[0]['bound']/float(line.split('\t')[1])
			print ('Gap:', '%.2f%%'%gap)
	sys.exit()