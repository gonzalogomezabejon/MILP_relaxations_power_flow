
using PowerModels
using Ipopt
using SDPLR
using SCS
# using SDPT3
using Mosek
using MosekTools

# test_cases = ["case3_lmbd", "case3_lmbd", "case5_pjm", "case14_ieee", "case24_ieee_rts", "case30_as", "case30_ieee", "case39_epri", "case57_ieee", "case60_c", "case73_ieee_rts", "case89_pegase", "case118_ieee"]
# test_cases = ["case162_ieee_dtc", "case162_ieee_dtc", "case179_goc", "case197_snem", "case200_activ", "case240_pserc", "case300_ieee"]#, "case500_goc", "case588_sdet", "case793_goc", "case1354_pegase", "case1803_snem", "case1888_rte", "case1951_rte", "case2000_goc", "case2312_goc", "case2383wp_k", "case2736sp_k", "case2737sop_k", "case2742_goc", "case2746wp_k", "case2746wop_k", "case2848_rte", "case2853_sdet", "case2868_rte", "case2869_pegase", "case3012wp_k", "case3022_goc", "case3120sp_k", "case3375wp_k", "case3970_goc", "case4020_goc", "case4601_goc", "case4619_goc", "case4661_sdet", "case4837_goc", "case4917_goc", "case5658_epigrids", "case6468_rte", "case6470_rte", "case6495_rte", "case6515_rte", "case7336_epigrids", "case8387_pegase", "case9241_pegase", "case9591_goc", "case10000_goc", "case10192_epigrids", "case10480_goc", "case13659_pegase", "case19402_goc", "case20758_epigrids", "case24464_goc", "case30000_goc", "case78484_epigrids"]
test_cases = ["case6_c", "case6_c", "case9_wscc", "case6_ww"]
exec_times = zeros(0)
objectives = zeros(0)
suffix = ""
prefix = "pglib_opf_"
prefix = "nesta_"
folder = ""
filename = "nesta_typ"
library_folder = "pglib-opf/"
library_folder = "nesta-0.7.0/opf/"

# test_cases = ["case9mod", "case9mod", "case22loop", "case39mod1", "case39mod2"] #, "case118mod", "case300mod", "LMBM3", "WB2", "WB3", "WB5"]
# suffix = ""
# prefix = ""
# folder = "local_optima/"
# filename = "local_optima"


open("results/results_julia_"*filename*".txt","w") do io
	println(io, "Case\tobjective_search_ipopt\ttime_search_ipopt\tobjective_qc_ipopt\ttime_qc_ipopt\tobjective_soc_ipopt\ttime_soc_ipopt\tobjective_sdp_mosek\ttime_sdp_mosek")
	# \tobjective_sdp_sdplr\ttime_sdp_sdplr\tobjective_sdp_scs\ttime_sdp_scs
	for test_case in test_cases

		println("Running solver for ", test_case)
		# data = PowerModels.parse_file("/home/gonzalo/.jqsoft/pglib-opf/pglib_opf_"*test_case*".m")
		data = PowerModels.parse_file(library_folder*folder*prefix*test_case*suffix*".m")
		time0 = time()
		result = solve_model(data, ACPPowerModel, Ipopt.Optimizer, build_opf)
		# result = solve_ac_opf("/home/gonzalo/.jqsoft/pglib-opf/pglib_opf_"*test_case*".m", Ipopt.Optimizer)
		exec_time_0 = time() - time0
		println("Execution time", exec_time_0)
		append!( exec_times, exec_time_0)
		#solve_opf(file, ACPPowerModel, optimizer; kwargs...)
		println("Objective fields:", result["objective"])
		append!( objectives, result["objective"] )

		time0 = time()
		result_qc_ipopt = solve_model(data, QCRMPowerModel, Ipopt.Optimizer, build_opf)
		# result = solve_ac_opf("pglib-opf/pglib_opf_"*test_case*".m", Ipopt.Optimizer)
		exec_time_1 = time() - time0
		# time0 = time()
		# result_sdp_sdplr = solve_model(data, SDPWRMPowerModel, SDPLR.Optimizer, build_opf)
		# # result_sdp = solve_model(data, SDPWRMPowerModel, SDPT3.Optimizer, build_opf)
		# exec_time_2 = time() - time0

		# time0 = time()
		# result_sdp_scs = solve_model(data, SDPWRMPowerModel, optimizer_with_attributes(SCS.Optimizer, "verbose"=>0), build_opf)
		# exec_time_3 = time() - time0

		time0 = time()
		result_soc_ipopt = solve_model(data, SOCWRPowerModel, Ipopt.Optimizer, build_opf)
		exec_time_4 = time() - time0

		time0 = time()
		result_sdp_mosek = solve_model(data, SDPWRMPowerModel, Mosek.Optimizer, build_opf)
		exec_time_5 = time() - time0

		println(io, test_case, "\t", result["objective"], "\t", exec_time_0, "\t", result_qc_ipopt["objective"], "\t", exec_time_1, 
			# "\t", result_sdp_sdplr["objective"], "\t", exec_time_2, "\t", result_sdp_scs["objective"], "\t", exec_time_3,
			"\t", result_soc_ipopt["objective"], "\t", exec_time_4, "\t", result_sdp_mosek["objective"], "\t", exec_time_5)
		# println(io, test_case, "\t", result["objective"], "\t", exec_time_0, "\t", result_sdp_mosek["objective"], "\t", exec_time_5)
		println("\n\n\n Solved ", test_case, " IPOPT: ", result["objective"], " MOSEK: ", result_sdp_mosek["objective"], "\n\n\n")


	# nodes = keys(result["solution"]["bus"])
	# println("Nodes:", nodes)
	# e_var = nodes
	# f_var = nodes
	# for node in nodes
	# 	println(io, node, " voltage: ", result["solution"]["bus"][string(node)])
	# end

	#println("Solution: ", result["solution"])
	# println("Solution fields:", keys(result["solution"]))
	# println("Result fields:", keys(result))
	# sys.exit()

	# gen_nodes = keys(result["solution"]["gen"])
	# branches = keys(result["solution"]["branch"])

	# open("julia_sol.txt","w") do io

	# 	for generator in gen_nodes
	# 		println(io, "aux_pg_", generator, " = ", result["solution"]["gen"][string(generator)]["pg"])
	# 		println(io, "aux_qg_", generator, " = ", result["solution"]["gen"][string(generator)]["qg"])
	# 	end

	# 	for kk in nodes
	# 		va = result["solution"]["bus"][string(kk)]["va"]
	# 		vm = result["solution"]["bus"][string(kk)]["vm"]
	# 		println(io, "va_", kk, " = ", va)
	# 		println(io, "vm_", kk, " = ", vm)
	# 	end

	# 	for branch in branches
	# 		println(io, "p_flow_from_", branch, " = ", result["solution"]["branch"][string(branch)]["pf"])
	# 		println(io, "q_flow_from_", branch, " = ", result["solution"]["branch"][string(branch)]["qf"])
	# 		println(io, "p_flow_to_", branch, " = ", result["solution"]["branch"][string(branch)]["pt"])
	# 		println(io, "q_flow_to_", branch, " = ", result["solution"]["branch"][string(branch)]["qt"])
	# 	end
	# end
	end
end

println(test_cases)
println(exec_times)
println(objectives)

# ["case3_lmbd", "case5_pjm", "case14_ieee", "case24_ieee_rts", "case30_as", "case30_ieee", "case39_epri", "case57_ieee", "case60_c", "case73_ieee_rts", "case89_pegase", "case118_ieee", "case162_ieee_dtc", "case179_goc", "case197_snem", "case200_activ", "case240_pserc", "case300_ieee", "case500_goc", "case588_sdet", "case793_goc", "case1354_pegase", "case1803_snem", "case1888_rte", "case1951_rte", "case2000_goc", "case2312_goc", "case2383wp_k", "case2736sp_k", "case2737sop_k", "case2742_goc", "case2746wp_k", "case2746wop_k", "case2848_rte", "case2853_sdet", "case2868_rte", "case2869_pegase", "case3012wp_k", "case3022_goc", "case3120sp_k", "case3375wp_k", "case3970_goc", "case4020_goc", "case4601_goc", "case4619_goc", "case4661_sdet", "case4837_goc", "case4917_goc", "case5658_epigrids", "case6468_rte", "case6470_rte", "case6495_rte", "case6515_rte", "case7336_epigrids", "case8387_pegase", "case9241_pegase", "case9591_goc", "case10000_goc", "case10192_epigrids", "case10480_goc", "case13659_pegase", "case19402_goc", "case20758_epigrids", "case24464_goc", "case30000_goc", "case78484_epigrids"]
# [0.019678115844726562, 0.4022819995880127, 0.41129589080810547, 0.09222197532653809, 0.03761696815490723, 0.08340620994567871, 0.07701301574707031, 0.0755760669708252, 0.25785398483276367, 0.1507129669189453, 0.3699021339416504, 0.27025699615478516, 0.3873891830444336, 0.5916929244995117, 0.32007288932800293, 0.25044798851013184, 3.418570041656494, 0.6651160717010498, 1.2731740474700928, 1.0408480167388916, 1.4705429077148438, 4.060840129852295, 9.101161003112793, 9.819712162017822, 7.176493167877197, 8.129146099090576, 6.157668828964233, 7.599519968032837, 6.180613994598389, 6.639733076095581, 24.092157125473022, 6.060914993286133, 6.694109916687012, 14.379122018814087, 8.48564100265503, 14.424613952636719, 10.67537784576416, 9.336115837097168, 12.29654598236084, 9.117607831954956, 10.056735038757324, 19.30965495109558, 41.89067220687866, 22.090655088424683, 20.42939591407776, 17.942494869232178, 22.693082094192505, 18.6300208568573, 24.253079891204834, 53.65568709373474, 34.947628021240234, 54.93614387512207, 42.633119106292725, 28.32610011100769, 55.689955949783325, 53.415639877319336, 57.43434000015259, 43.28433012962341, 51.6811900138855, 56.159862995147705, 70.01711201667786, 119.25563716888428, 103.88664388656616, 1240.0700190067291, 390.25721502304077, 1347.7599740028381]
# [5812.64293503618, 17551.890838594758, 2178.080410819641, 63352.201086064175, 803.1273103853115, 8208.515427944541, 138415.5625287722, 37589.33820419316, 92693.6695979166, 189764.0771673735, 107285.67306998349, 97213.60693895924, 108075.64369843062, 754266.4141684102, 1.5016545569315944, 27557.570465648296, 3.329670055972561e6, 565219.9718650873, 454945.978311204, 313139.77996788063, 260197.8478786532, 1.258843985100568e6, 98334.55691408037, 1.4025308246576267e6, 2.0855818360669948e6, 973432.4715402615, 441330.3289176641, 1.8681915835349108e6, 1.3080149931456242e6, 777727.6837669068, 275705.45262737107, 1.631707930248126e6, 1.2082585006475565e6, 1.286608191272638e6, 2.052386718513573e6, 2.0096053046873603e6, 2.4627904324824684e6, 2.6008427211516215e6, 601383.8430137298, 2.14796907483637e6, 7.438169412911984e6, 960985.2533687336, 822247.2841290073, 826241.5226459657, 476703.7249870737, 2.251344065334859e6, 872255.3156134635, 1.3877910158801272e6, 1.2073137799883753e6, 2.069730145296015e6, 2.237570580215981e6, 3.067825301481965e6, 2.82549964118996e6, 1.882390221232304e6, 2.771392261880956e6, 6.243090383619432e6, 1.0616835533567932e6, 1.3540313143350312e6, 1.686922671847054e6, 2.3146479995697928e6, 8.948048917265542e6, 1.9778153972977179e6, 2.618635680815579e6, 2.6295312892210935e6, 1.1423315687882868e6, 1.5315885662626693e7]
