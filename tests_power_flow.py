
import cdc_formulation as form
from tqdm import tqdm
import iss_checker

# This script tests basic functionality like reading OPF instances from pglib-opf
# and that their respective IPOPT solutions are feasible in the polar formulation

if __name__ == '__main__':

	lib_dir = 'pglib-opf/'

	# printing the model
	# base_MVA, gen_data, bus_data, branch_data = read_pglib(lib_dir+'pglib_opf_case60_c.m')
	# obj = Power_Flow_Problem('Case3', base_MVA, gen_data, bus_data, branch_data)
	# quad_problem = obj.jabr_formulation()
	# quad_problem.display()

	import os
	test_cases = []
	for a,b,c in os.walk(lib_dir):
		test_cases = c
		break
	results = {}
	files_read = []
	files_to_use = []
	test_cases = [it for it in test_cases if it[-2:] == '.m']
	for filenm in tqdm(test_cases): # ['pglib_opf_case179_goc.m']:#
		try:
			base_MVA, bus_data, gen_data, branch_data = form.read_pglib(lib_dir+filenm)
			files_read.append(filenm)
			if len(branch_data) < 5000:
				files_to_use.append(filenm)
		except Exception as ex:
			print ('{*} Failed to read file:', filenm)
			continue
	errors = len(test_cases)-len(files_read)
	print ('\033[91m{-} Test 1 failed' if errors else '\033[92m{+} Test 1 Passed', '(%d/%d) errors\033[0m'%(errors, len(test_cases)))

	import subprocess

	relaxation_correct = []
	for filenm in tqdm(files_to_use):
		try:
			subprocess.run(['julia', 'benchmark_solution.jl', lib_dir+filenm, 'julia_sol.txt'], stdout = subprocess.DEVNULL)
			base_MVA, bus_data, gen_data, branch_data = form.read_pglib(lib_dir+filenm, include_0_status = False)
			obj = form.Power_Flow_Problem('Case_n', base_MVA, bus_data, gen_data, branch_data)
			variables = iss_checker.parse_file('julia_sol.txt', bus_data, gen_data, branch_data)
			jabr_formulation = obj.jabr_formulation()
			violations, cone_vio = jabr_formulation.check_solution(variables, margin = 1e-5, verbose = False)
			helix_formulation = obj.helix_formulation()
			h_violations, h_cone_vio = helix_formulation.check_solution(variables, margin = 1e-5, verbose = False)
			if violations == 0 and cone_vio == 0 and h_violations == 0 and h_cone_vio == 0:
				relaxation_correct.append(filenm)
			else:
				print('{*} IPOpt solution unfeasible in compact model:', filenm)
		except Exception as ex:
			print ('{*} Error running IPOPT feasibility:', filenm)
			print ('Exception:', ex)

	errors = len(files_to_use)-len(relaxation_correct)
	print ('\033[91m{-} Test 2 failed' if errors else '\033[92m{+} Test 2 Passed', '(%d/%d) errors\033[0m'%(errors, len(files_to_use)))
