
import cdc_formulation as form
import time
import datetime

if __name__ == '__main__':
	# lib_dir = 'pglib-opf/api/'
	# lib_dir = 'pglib-opf/sad/'
	lib_dir = 'pglib-opf/'
	filename = 'results/cdc_compact_delete.txt'
	# filename = 'results/bs.txt'

	import os
	test_cases = []

	for lib_dir in ['pglib-opf/api/', 'pglib-opf/sad/', 'pglib-opf/']:
	# for lib_dir in ['nesta-0.7.0/opf/sad/', 'nesta-0.7.0/opf/', 'nesta-0.7.0/opf/api/']:
		for a,b,c in os.walk(lib_dir):
			test_cases = test_cases + c
			break
	cases = {}
	results = {}
	target_depth = [7,9]#[6,8]# [6,7,8,9]
	flag=1
	lib_dir = 'nesta-0.7.0/opf/'
	lib_dir = 'pglib-opf/'

	for filenm in test_cases: # ['pglib_opf_case179_goc.m']:#
		folder = 'api/' if '__api' in filenm else ('sad/' if '__sad' in filenm else '')
		if filenm[-2:] != '.m':
			continue
		try:
			base_MVA, bus_data, gen_data, branch_data = form.read_pglib(lib_dir+folder+filenm)
		except Exception as ex:
			print ('{*} Failed to read file:', filenm)
			continue
		case_name = filenm[10:-2]
		if len(bus_data) < 99 or len(bus_data) > 300:# or len(branch_data)<5:
			continue
		# if flag:
		# 	if case_name == 'case200_activ':
		# 		flag = 0
		# 	# else:
		# 	continue
		print ('\nCase ', case_name)
		print (len(bus_data), 'nodes   ', len(gen_data), 'generators &  ', len(branch_data), 'branches')
		cases[case_name] = {'nodes': len(bus_data), 'generators': len(gen_data), 'branches': len(branch_data)}

		case_type = 'sad' if '__sad' in filenm else ('api' if '__api' in filenm else 'typ')
		baseline = 'results/results_julia_'+case_type+'.txt'
		with open(baseline, 'r') as f:
			lines = f.read().split('\n')
		for line in lines:
			if line.split('\t')[0] in filenm:
				ipopt_objective = float(line.split('\t')[1])

		results[case_name] = {'objective': {}, 'time': {}, 'memory': {}, 'nodes': {}, 'objective': {}, 'gap': {}}
		obj = form.Power_Flow_Problem(filenm[:-2], base_MVA, bus_data, gen_data, branch_data)
		helix_problem = obj.helix_formulation()
		quad_problem = obj.rectangular_formulation()
		for depth in target_depth:
			# print ('\n############## Case %s, MIP relaxation (compact), depth=%d #############'%(case_name, depth))
			# time0 = time.time()
			# # print(helix_problem.cone_ineq)
			# lp_results = helix_problem.solve(mode='mixed', time_limit=3600, cone_depth=depth, interior_depth=1, concave_iterations=1, helix_depth=depth)
			# results[case_name]['time']['Compact_%d'%depth] = time.time()-time0
			# print ('Finished in %f'%results[case_name]['time']['Compact_%d'%depth], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			# print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
			# obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
			# results[case_name]['objective']['Compact_%d'%depth] = obj
			# results[case_name]['memory']['Compact_%d'%depth] = lp_results[0]['mem']
			# results[case_name]['nodes']['Compact_%d'%depth] = lp_results[0]['nodes']
			# gap = (100.-100*obj/ipopt_objective)
			# results[case_name]['gap']['Compact_%d'%depth] = gap
			# print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
			# # helix_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)

			# print ('\n############## Case %s, MIP relaxation (CDC), depth=%d #############'%(case_name, depth))
			# time0 = time.time()
			# # print(helix_problem.cone_ineq)
			# lp_results = helix_problem.solve(mode='cdc', time_limit=3600, cone_depth=depth, interior_depth=1, concave_iterations=1, helix_depth=depth)
			# results[case_name]['time']['CDC_%d'%depth] = time.time()-time0
			# if lp_results:
			# 	print ('Finished in %f'%results[case_name]['time']['CDC_%d'%depth], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			# 	print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
			# 	obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
			# 	results[case_name]['objective']['CDC_%d'%depth] = obj
			# 	results[case_name]['memory']['CDC_%d'%depth] = lp_results[0]['mem']
			# 	results[case_name]['nodes']['CDC_%d'%depth] = lp_results[0]['nodes']
			# 	gap = (100.-100*obj/ipopt_objective)
			# 	results[case_name]['gap']['CDC_%d'%depth] = gap
			# 	print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
			# 	# helix_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)
			# else:
			# 	results[case_name]['memory']['CDC_%d'%depth] = results[case_name]['nodes']['CDC_%d'%depth] = 0
			# 	results[case_name]['objective']['CDC_%d'%depth] = results[case_name]['gap']['CDC_%d'%depth] = 0

			bilinear_depth = depth // 2

			# print ('\n############## Case %s, MIP relaxation (DNMDT), depth=%d #############'%(case_name, bilinear_depth))
			# time0 = time.time()
			# # print(quad_problem.cone_ineq)
			# lp_results = quad_problem.solve(mode='dnmdt', time_limit=3600, bilinear_depth = bilinear_depth)
			# if lp_results:
			# 	results[case_name]['time']['DNMDT_%d'%bilinear_depth] = time.time()-time0
			# 	print ('Finished in %f'%results[case_name]['time']['DNMDT_%d'%bilinear_depth], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			# 	print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
			# 	obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
			# 	results[case_name]['objective']['DNMDT_%d'%bilinear_depth] = obj
			# 	results[case_name]['memory']['DNMDT_%d'%bilinear_depth] = lp_results[0]['mem']
			# 	results[case_name]['nodes']['DNMDT_%d'%bilinear_depth] = lp_results[0]['nodes']
			# 	gap = (100.-100*obj/ipopt_objective)
			# 	results[case_name]['gap']['DNMDT_%d'%bilinear_depth] = gap
			# 	print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
			# else:
			# 	results[case_name]['memory']['DNMDT_%d'%bilinear_depth] = results[case_name]['nodes']['DNMDT_%d'%bilinear_depth] = 0
			# 	results[case_name]['objective']['DNMDT_%d'%bilinear_depth] = results[case_name]['gap']['DNMDT_%d'%bilinear_depth] = 0
			# # quad_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)

			# print ('\n############## Case %s, MIP relaxation (NMDT), depth=%d #############'%(case_name, bilinear_depth))
			# time0 = time.time()
			# # print(quad_problem.cone_ineq)
			# lp_results = quad_problem.solve(mode='nmdt', time_limit=3600, bilinear_depth=bilinear_depth)
			# if lp_results:
			# 	results[case_name]['time']['NMDT_%d'%bilinear_depth] = time.time()-time0
			# 	print ('Finished in %f'%results[case_name]['time']['NMDT_%d'%bilinear_depth], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			# 	print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
			# 	obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
			# 	results[case_name]['objective']['NMDT_%d'%bilinear_depth] = obj
			# 	results[case_name]['memory']['NMDT_%d'%bilinear_depth] = lp_results[0]['mem']
			# 	results[case_name]['nodes']['NMDT_%d'%bilinear_depth] = lp_results[0]['nodes']
			# 	gap = (100.-100*obj/ipopt_objective)
			# 	results[case_name]['gap']['NMDT_%d'%bilinear_depth] = gap
			# 	print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
			# else:
			# 	results[case_name]['memory']['NMDT_%d'%bilinear_depth] = results[case_name]['nodes']['NMDT_%d'%bilinear_depth] = 0
			# 	results[case_name]['objective']['NMDT_%d'%bilinear_depth] = results[case_name]['gap']['NMDT_%d'%bilinear_depth] = 0
			# # quad_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)

		# print ('\n############## Case %s, Full Rectangular (gurobi), depth=infty #############'%(case_name))
		# time0 = time.time()
		# # print(quad_problem.cone_ineq)
		# lp_results = quad_problem.solve(mode='gnb', time_limit=3600)
		# if lp_results:
		# 	results[case_name]['time']['Gurobi_'] = time.time()-time0
		# 	print ('Finished in %f'%results[case_name]['time']['Gurobi_'], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
		# 	print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
		# 	obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
		# 	results[case_name]['objective']['Gurobi_'] = obj
		# 	results[case_name]['memory']['Gurobi_'] = lp_results[0]['mem']
		# 	results[case_name]['nodes']['Gurobi_'] = lp_results[0]['nodes']
		# 	gap = (100.-100*obj/ipopt_objective)
		# 	results[case_name]['nodes']['Gurobi_'] = gap
		# 	print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
		# else:
		# 	results[case_name]['memory']['Gurobi_'] = results[case_name]['nodes']['Gurobi_'] = 0
		# 	results[case_name]['objective']['Gurobi_'] = results[case_name]['gap']['Gurobi_'] = 0
		# # quad_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)


		print ('\n############## Case %s, Full Polar (gurobi), depth=infty #############'%(case_name))
		time0 = time.time()
		# print(quad_problem.cone_ineq)
		lp_results = helix_problem.solve(mode='gnp', time_limit=3600)
		if lp_results:
			results[case_name]['time']['GurobiP_'] = time.time()-time0
			print ('Finished in %f'%results[case_name]['time']['GurobiP_'], 'seconds on', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
			print ("Memory: ", lp_results[0]['mem'], "       Nodes: ", lp_results[0]['nodes'])
			obj = lp_results[0]['obj'] if lp_results[0]['solved'] else lp_results[0]['bound']
			results[case_name]['objective']['GurobiP_'] = obj
			results[case_name]['memory']['GurobiP_'] = lp_results[0]['mem']
			results[case_name]['nodes']['GurobiP_'] = lp_results[0]['nodes']
			gap = (100.-100*obj/ipopt_objective)
			results[case_name]['nodes']['GurobiP_'] = gap
			print ('Objective: ', obj, '  Gap: %.02f%%'%gap)
		else:
			results[case_name]['memory']['GurobiP_'] = results[case_name]['nodes']['GurobiP_'] = 0
			results[case_name]['objective']['GurobiP_'] = results[case_name]['gap']['GurobiP_'] = 0
		# quad_problem.check_solution(lp_results[0]['sol'], margin = 1e-5, fill_variables = False)
		
	for it in results:
		print ('Stats:', it, results[it])
	models = list(results[list(results.keys())[0]]['objective'].keys())
	with open(filename, 'w') as f:
		f.write('case\t'+'\t'.join(['%s objective\t%s time\t%s memory\t%s nodes\t%s objective\t%s gap'%(model, model, model, model, model, model) for model in models]))
		for case_name in results:
			f.write('\n'+case_name+'\t'+'\t'.join(['%.2f\t%.2f\t%.3f\t%d\t%.3f\t%.3f'%(results[case_name]['objective'][model],
				results[case_name]['time'][model], results[case_name]['memory'][model], results[case_name]['nodes'][model],
				results[case_name]['objective'][model], results[case_name]['gap'][model]) for model in models]))
			# print ('Cone percentage:', '%.2f'%results[it][-1][0], it)

	# print ('\n|    Case    | '+ ' | '.join(['{:^9d}'.format(ii) for ii in range(1,7)]), '\n|-'+'|'.join(['-'*11]*7) )
	# for it in results:
	# 	print ('| %10s | '%it[4:13] + ' | '.join(['{:9.2f}'.format(results[it]['lp_violations'][ii][0]['obj']) for ii in range(1,7)]))
	with open('ipopt_solutions.txt', 'r') as f:
		lines = f.read().split('\n')

	obj_ipopt = lines[0].split('\t').index('objective_total')
	obj_soc = lines[0].split('\t').index('objective_soc')
	obj_sdp = lines[0].split('\t').index('objective_sdp')
	baseline = {line.split('\t')[0]: (float(line.split('\t')[obj_ipopt]), float(line.split('\t')[obj_soc]), float(line.split('\t')[obj_sdp])) for line in lines[1:-1]}
	# baseline = [["case3_lmbd", "case5_pjm", "case14_ieee", "case24_ieee_rts", "case30_as", "case30_ieee", "case39_epri", "case57_ieee", "case60_c", "case73_ieee_rts", "case89_pegase", "case118_ieee", "case162_ieee_dtc", "case179_goc", "case197_snem", "case200_activ", "case240_pserc", "case300_ieee", "case500_goc", "case588_sdet", "case793_goc", "case1354_pegase", "case1803_snem", "case1888_rte", "case1951_rte", "case2000_goc", "case2312_goc", "case2383wp_k", "case2736sp_k", "case2737sop_k", "case2742_goc", "case2746wp_k", "case2746wop_k", "case2848_rte", "case2853_sdet", "case2868_rte", "case2869_pegase", "case3012wp_k", "case3022_goc", "case3120sp_k", "case3375wp_k", "case3970_goc", "case4020_goc", "case4601_goc", "case4619_goc", "case4661_sdet", "case4837_goc", "case4917_goc", "case5658_epigrids", "case6468_rte", "case6470_rte", "case6495_rte", "case6515_rte", "case7336_epigrids", "case8387_pegase", "case9241_pegase", "case9591_goc", "case10000_goc", "case10192_epigrids", "case10480_goc", "case13659_pegase", "case19402_goc", "case20758_epigrids", "case24464_goc", "case30000_goc", "case78484_epigrids"], 
	# 	[0.019678115844726562, 0.4022819995880127, 0.41129589080810547, 0.09222197532653809, 0.03761696815490723, 0.08340620994567871, 0.07701301574707031, 0.0755760669708252, 0.25785398483276367, 0.1507129669189453, 0.3699021339416504, 0.27025699615478516, 0.3873891830444336, 0.5916929244995117, 0.32007288932800293, 0.25044798851013184, 3.418570041656494, 0.6651160717010498, 1.2731740474700928, 1.0408480167388916, 1.4705429077148438, 4.060840129852295, 9.101161003112793, 9.819712162017822, 7.176493167877197, 8.129146099090576, 6.157668828964233, 7.599519968032837, 6.180613994598389, 6.639733076095581, 24.092157125473022, 6.060914993286133, 6.694109916687012, 14.379122018814087, 8.48564100265503, 14.424613952636719, 10.67537784576416, 9.336115837097168, 12.29654598236084, 9.117607831954956, 10.056735038757324, 19.30965495109558, 41.89067220687866, 22.090655088424683, 20.42939591407776, 17.942494869232178, 22.693082094192505, 18.6300208568573, 24.253079891204834, 53.65568709373474, 34.947628021240234, 54.93614387512207, 42.633119106292725, 28.32610011100769, 55.689955949783325, 53.415639877319336, 57.43434000015259, 43.28433012962341, 51.6811900138855, 56.159862995147705, 70.01711201667786, 119.25563716888428, 103.88664388656616, 1240.0700190067291, 390.25721502304077, 1347.7599740028381],
	# 	[5812.64293503618, 17551.890838594758, 2178.080410819641, 63352.201086064175, 803.1273103853115, 8208.515427944541, 138415.5625287722, 37589.33820419316, 92693.6695979166, 189764.0771673735, 107285.67306998349, 97213.60693895924, 108075.64369843062, 754266.4141684102, 1.5016545569315944, 27557.570465648296, 3.329670055972561e6, 565219.9718650873, 454945.978311204, 313139.77996788063, 260197.8478786532, 1.258843985100568e6, 98334.55691408037, 1.4025308246576267e6, 2.0855818360669948e6, 973432.4715402615, 441330.3289176641, 1.8681915835349108e6, 1.3080149931456242e6, 777727.6837669068, 275705.45262737107, 1.631707930248126e6, 1.2082585006475565e6, 1.286608191272638e6, 2.052386718513573e6, 2.0096053046873603e6, 2.4627904324824684e6, 2.6008427211516215e6, 601383.8430137298, 2.14796907483637e6, 7.438169412911984e6, 960985.2533687336, 822247.2841290073, 826241.5226459657, 476703.7249870737, 2.251344065334859e6, 872255.3156134635, 1.3877910158801272e6, 1.2073137799883753e6, 2.069730145296015e6, 2.237570580215981e6, 3.067825301481965e6, 2.82549964118996e6, 1.882390221232304e6, 2.771392261880956e6, 6.243090383619432e6, 1.0616835533567932e6, 1.3540313143350312e6, 1.686922671847054e6, 2.3146479995697928e6, 8.948048917265542e6, 1.9778153972977179e6, 2.618635680815579e6, 2.6295312892210935e6, 1.1423315687882868e6, 1.5315885662626693e7]]

	# julia_time_result = {case: (time, obj) for case, time, obj in zip(baseline[0], baseline[1], baseline[2])}
	# soc_gap = {'case3_lmbd': 1.32, 'case5_pjm': 14.55, 'case14_ieee': 0.11, 'case24_ieee_rts': 0.02, 'case30_as': 0.06, 'case30_ieee': 18.84, 'case39_epri': 0.56, 'case57_ieee': 0.16, 'case60_c': 0.07, 'case73_ieee_rts': 0.04, 'case89_pegase': 0.75, 'case118_ieee': 0.91, 'case162_ieee_dtc': 5.95, 'case179_goc': 0.16, 'case197_snem': 0.05, 'case200_activ': 0.01, 'case240_pserc': 2.78, 'case300_ieee': 2.63, 'case500_goc': 0.25, 'case588_sdet': 2.14, 'case793_goc': 1.33, 'case1354_pegase': 1.57, 'case1803_snem': 8.03, 'case1888_rte': 2.05, 'case1951_rte': 0.14, 'case2000_goc': 0.31, 'case2312_goc': 1.90, 'case2383wp_k': 1.04, 'case2736sp_k': 0.31, 'case2737sop_k': 0.27, 'case2742_goc': 1.35, 'case2746wop_k': 0.37, 'case2746wp_k': 0.33, 'case2848_rte': 0.13, 'case2853_sdet': 0.91, 'case2868_rte': 0.10, 'case2869_pegase': 1.01, 'case3012wp_k': 1.03, 'case3022_goc': 2.77, 'case3120sp_k': 0.56, 'case3375wp_k': 0.55, 'case3970_goc': 0.344}

	print('\\hline Case & Nodes & Edges & SOC gap & SDP gap & '+ ' & '.join([model+' Gap & '+model+' Time' for model in models]) + ' \\\\\\hline')
	for case in sorted(results, key=lambda x: cases[x]['nodes']):
		# gap = 100*(julia_time_result[case][1] - results[case]['lp_violations'][target_depth][0]['obj'])/julia_time_result[case][1]
		soc_gap = 100*(baseline[case][0]-baseline[case][1])/baseline[case][0]
		sdp_gap = 100*(baseline[case][0]-baseline[case][2])/baseline[case][0]
		gap_and_time = ' & '.join(['%.2f\\%% & %.2f s'%(100*(baseline[case][0]-results[case]['objective'][model])/baseline[case][0], results[case]['time'][model]) for model in models])
		print (case[4:13].replace('_', '\\_'), '&', cases[case]['nodes'], '&', cases[case]['branches'], '& %.2f\\%% & %.2f\\%% &'%(soc_gap, sdp_gap), gap_and_time, '\\\\')
