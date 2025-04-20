import gurobipy as gp
from gurobipy import GRB
from reader import parse_cor, parse_tim, parse_sto_dep, parse_sto_indep, split_stages, build_extensive_form

def solve_extensive_form(A_ext, b_ext, c_ext, row_sense):
    m, n = A_ext.shape

    model = gp.Model("AirliftExtensiveForm")
    model.setParam("OutputFlag", 0)

    # Add variables (all continuous and nonnegative by default)
    x = model.addVars(n, lb=0.0, vtype=GRB.CONTINUOUS, name="x")

    # Add constraints
    for i in range(m):
        expr = gp.LinExpr()
        for j in range(n):
            if A_ext[i, j] != 0.0:
                expr.add(A_ext[i, j] * x[j])

        sense = row_sense[i]
        rhs = b_ext[i]

        if sense == "L":
            model.addConstr(expr <= rhs, name=f"c{i}")
        elif sense == "E":
            model.addConstr(expr == rhs, name=f"c{i}")
        elif sense == "G":
            model.addConstr(expr >= rhs, name=f"c{i}")
        else:
            raise ValueError(f"Unknown constraint sense '{sense}' at row {i}")

    # Set objective
    obj_expr = gp.LinExpr()
    for j in range(n):
        if c_ext[j] != 0.0:
            obj_expr.add(c_ext[j] * x[j])
    model.setObjective(obj_expr, GRB.MINIMIZE)

    # Solve
    model.optimize()

    # Extract results
    if model.status == GRB.OPTIMAL:
        print("Model solved to optimality.")
        solution = model.getAttr("x", x)
        objective = model.objVal
        return solution, objective
    else:
        print("Model not solved to optimality.")
        return None, None


cor_res = parse_cor("airlift/AIRL.cor")
tim_res = parse_tim("airlift/AIRL.tim")

split = split_stages(
    A=cor_res["A"],
    b=cor_res["b"],
    c=cor_res["c"],
    var_names=cor_res["var_names"],
    variable_stage=tim_res
)

sto_first = parse_sto_dep("airlift/AIRL.sto.first", cor_res["row_names"], cor_res["b"])
sto_second = parse_sto_indep("airlift/AIRL.sto.second", cor_res["row_names"], cor_res["b"])

ext_form_first = build_extensive_form(
    A1=split["A1"],
    A2=split["A2"],
    b=cor_res["b"],
    c1=split["c1"],
    c2=split["c2"],
    sto_scenarios=sto_first,
    row_sense=cor_res["row_sense"]
)

ext_form_second = build_extensive_form(
    A1=split["A1"],
    A2=split["A2"],
    b=cor_res["b"],
    c1=split["c1"],
    c2=split["c2"],
    sto_scenarios=sto_second,
    row_sense=cor_res["row_sense"]
)


solution_first, obj_first = solve_extensive_form(
    A_ext=ext_form_first["A_ext"],
    b_ext=ext_form_first["b_ext"],
    c_ext=ext_form_first["c_ext"],
    row_sense=ext_form_first["row_sense"]
)

solution_second, obj_second = solve_extensive_form(
    A_ext=ext_form_second["A_ext"],
    b_ext=ext_form_second["b_ext"],
    c_ext=ext_form_second["c_ext"],
    row_sense=ext_form_second["row_sense"]
)

print("---sto.first LP---\n")
print("Optimal objective value:", obj_first)

# First-stage decision variables
n1_first = ext_form_first["n1"]
print("First-stage decision:")
for i in range(n1_first):
    print(f"x{i} = {solution_first[i]}")
    
print("\n---sto.second LP---\n")
print("Optimal objective value:", obj_second)

# Second-stage decision variables
n2_first = ext_form_second["n1"]
print("Second-stage decision:")
for i in range(n2_first):
    print(f"x{i} = {solution_second[i]}")
    