
using PowerModels
using Ipopt

# A simple script to solve OPF for one test case using IPOPT

function solution2file(problemfile, solutionfile)
	# println("Running solver for ", test_case)
	# result = solve_ac_opf("pglib-opf/"*test_case, Ipopt.Optimizer)
	result = solve_ac_opf(problemfile, Ipopt.Optimizer)

	#solve_opf(file, ACPPowerModel, optimizer; kwargs...)

	nodes = keys(result["solution"]["bus"])
	# println("Nodes:", nodes)
	e_var = nodes
	f_var = nodes
	# for node in nodes
	# 	println(io, node, " voltage: ", result["solution"]["bus"][string(node)])
	# end

	# println("Solution: ", result["solution"])
	# println("Solution fields:", keys(result["solution"]))
	# println("Generators", keys(result["solution"]["gen"]), result["solution"]["gen"])
	# sys.exit()

	gen_nodes = keys(result["solution"]["gen"])
	branches = keys(result["solution"]["branch"])

	# println(gen_nodes)

	open(solutionfile,"w") do io

		for generator in gen_nodes
			println(io, "aux_pg_", generator, " = ", result["solution"]["gen"][string(generator)]["pg"])
			println(io, "aux_qg_", generator, " = ", result["solution"]["gen"][string(generator)]["qg"])
		end

		for kk in nodes
			va = result["solution"]["bus"][string(kk)]["va"]
			vm = result["solution"]["bus"][string(kk)]["vm"]
			println(io, "va_", kk, " = ", va)
			println(io, "vm_", kk, " = ", vm)
		end
		# for ii in nodes
		# 	va = result["solution"]["bus"][string(ii)]["va"]
		# 	vm = result["solution"]["bus"][string(ii)]["vm"]
		# 	e_ii = vm*cos(va)
		# 	f_ii = vm*sin(va)
		# 	for jj in nodes
		# 		va = result["solution"]["bus"][string(jj)]["va"]
		# 		vm = result["solution"]["bus"][string(jj)]["vm"]
		# 		e_jj = vm*cos(va)
		# 		f_jj = vm*sin(va)
		# 		if ii != jj # parse(Int64, ii) < parse(Int64, jj)
		# 			println(io, "s_",ii,"_",jj," = ", e_ii*f_jj - e_jj*f_ii)
		# 			println(io, "c_",ii,"_",jj," = ", e_ii*e_jj + f_ii*f_jj)
		# 		end
		# 		if ii == jj
		# 			println(io, "c_",ii," = ", e_ii*e_ii + f_ii*f_ii)
		# 		end
		# 	end
		# end

		for branch in branches
			println(io, "p_flow_from_", branch, " = ", result["solution"]["branch"][string(branch)]["pf"])
			println(io, "q_flow_from_", branch, " = ", result["solution"]["branch"][string(branch)]["qf"])
			println(io, "p_flow_to_", branch, " = ", result["solution"]["branch"][string(branch)]["pt"])
			println(io, "q_flow_to_", branch, " = ", result["solution"]["branch"][string(branch)]["qt"])
		end
	end
end

if length(ARGS) == 0

	test_case = "pglib_opf_case179_goc.m"
	test_case = "pglib_opf_case60_c.m"
	test_case = "pglib_opf_case2383wp_k.m"
	test_case = "pglib_opf_case200_activ.m"
	test_case = "pglib_opf_case1803_snem.m"
	test_case = "pglib_opf_case24_ieee_rts.m"
	test_case = "pglib_opf_case3_lmbd.m"

	solution2file("pglib-opf/"*test_case, "julia_sol.txt")

end

if length(ARGS) == 2
	solution2file(ARGS[1], ARGS[2])
end




