import numpy as np
import gurobipy as gp
from gurobipy import GRB
from reader import parse_cor, parse_tim, parse_sto_dep, parse_sto_indep, split_stages, build_extensive_form
from itertools import product

def solve_subproblem(x_val, A1, A2, b_s, c2, row_sense):
    m, n2 = A2.shape
    y_model = gp.Model("BendersSubproblem")
    y_model.setParam("OutputFlag", 0)

    # Variables
    y = y_model.addVars(n2, lb=0.0, vtype=GRB.CONTINUOUS, name="y")

    # RHS: b_s - A1 x
    rhs_resid = b_s - A1 @ x_val

    # Constraints: A2 y <= rhs_resid
    constraints = []
    for i in range(m):
        expr = gp.LinExpr()
        for j in range(n2):
            if A2[i, j] != 0:
                expr.add(A2[i, j] * y[j])
        sense = row_sense[i]
        if sense == "L":
            constraints.append(y_model.addConstr(expr <= rhs_resid[i]))
        elif sense == "E":
            constraints.append(y_model.addConstr(expr == rhs_resid[i]))
        elif sense == "G":
            constraints.append(y_model.addConstr(expr >= rhs_resid[i]))
        else:
            raise ValueError(f"Unknown constraint sense: {sense}")

    # Objective
    obj = gp.LinExpr()
    for j in range(n2):
        if c2[j] != 0:
            obj.add(c2[j] * y[j])
    y_model.setObjective(obj, GRB.MINIMIZE)

    # Solve
    y_model.optimize()

    # Results
    if y_model.status == GRB.OPTIMAL:
        obj_val = y_model.ObjVal
        duals = [constr.Pi for constr in constraints]
        return {
            'status': 'optimal',
            'value': obj_val,
            'duals': np.array(duals)
        }
    elif y_model.status == GRB.INFEASIBLE:
        return {
            'status': 'infeasible'
        }
    else:
        raise RuntimeError(f"Unexpected subproblem status: {y_model.status}")

def benders_slim(A1, A2, b_scenarios, c1, c2, row_sense, tol=1e-6, max_iters=1000):
    m, n1 = A1.shape
    _, n2 = A2.shape
    K = len(b_scenarios)

    # Step 1: Initialize master model
    master = gp.Model("SlimBendersMaster")
    master.setParam("OutputFlag", 0)

    # Variables: x (first-stage), theta (estimate of recourse cost)
    x = master.addVars(n1, lb=0.0, vtype=GRB.CONTINUOUS, name="x")
    theta = master.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="theta")

    # Objective: c1^T x + theta
    master.setObjective(
        gp.quicksum(c1[i] * x[i] for i in range(n1)) + theta, GRB.MINIMIZE)

    # Add first-stage constraints
    for i in range(m):
        if row_sense[i] == "L":
            master.addConstr(
                gp.quicksum(A1[i, j] * x[j] for j in range(n1)) <= b_scenarios[0][1][i],
                name=f"master_constr_L_{i}"
            )
        elif row_sense[i] == "E":
            master.addConstr(
                gp.quicksum(A1[i, j] * x[j] for j in range(n1)) == b_scenarios[0][1][i],
                name=f"master_constr_E_{i}"
            )
        elif row_sense[i] == "G":
            master.addConstr(
                gp.quicksum(A1[i, j] * x[j] for j in range(n1)) >= b_scenarios[0][1][i],
                name=f"master_constr_G_{i}"
            )
        else:
            raise ValueError(f"Invalid row sense {row_sense[i]}")

    iteration = 0
    UB = 1e10     # large but finite
    LB = -1e10    # small but finite
    while iteration < max_iters and abs((UB - LB) / (abs(UB) + 1e6)) > tol:
        iteration += 1
        if iteration % 50 == 0:
            print(f"Iteration {iteration}:")
            print(f"Upper Bound: {UB:.4f}")
            print(f"Lower Bound: {LB:.4f}")

        master.optimize()

        if master.status != GRB.OPTIMAL:
            raise RuntimeError("Master problem not solved to optimality.")

        x_val = np.array([x[i].X for i in range(n1)])
        theta_val = theta.X

        expected_rec_cost = 0.0
        cuts_added = 0

        for prob, b_s in b_scenarios:
            sub_res = solve_subproblem(x_val, A1, A2, b_s, c2, row_sense)

            if sub_res["status"] == "optimal":
                expected_rec_cost += prob * sub_res["value"]
                pi = sub_res["duals"]

                # Optimality cut: θ ≥ πᵀ(b_s - A1 x)
                rhs = np.dot(pi, b_s)
                lhs_expr = gp.quicksum(-pi[i] * gp.quicksum(A1[i, j] * x[j] for j in range(n1)) for i in range(m))
                master.addConstr(theta >= rhs + lhs_expr)
                cuts_added += 1

            elif sub_res["status"] == "infeasible":
                raise NotImplementedError("Feasibility cuts not yet implemented.")
            else:
                raise RuntimeError("Unexpected subproblem status.")

        UB = master.objVal
        LB = sum(c1[i] * x_val[i] for i in range(n1)) + expected_rec_cost

    return {
        "x": x_val,
        "theta": theta_val,
        "objective": master.objVal,
        "iterations": iteration
    }
    
def benders_fat(A1, A2, b_scenarios, c1, c2, row_sense, row_names, tol=1e-4, max_iters=50):
    m, n1 = A1.shape
    _, n2 = A2.shape
    K = len(b_scenarios)

    master = gp.Model("FatBenders")
    master.setParam("OutputFlag", 0)

    # === VARIABLES ===
    x = master.addVars(n1, lb=0.0, vtype=GRB.CONTINUOUS, name="x")
    y = {
        k: master.addVars(n2, lb=0.0, vtype=GRB.CONTINUOUS, name=f"y_{k}")
        for k in range(K)
    }

    # === OBJECTIVE ===
    obj = gp.quicksum(c1[i] * x[i] for i in range(n1))
    for k, (prob, _) in enumerate(b_scenarios):
        obj += prob * gp.quicksum(c2[j] * y[k][j] for j in range(n2))
    master.setObjective(obj, GRB.MINIMIZE)

    # === FIRST-STAGE CONSTRAINTS ===
    # Identify first-stage constraints (those that do NOT depend on demand)
    first_stage_idx = [i for i, r in enumerate(row_sense) if "HOURS" in row_names[i]]

    for i in first_stage_idx:
        sense = row_sense[i]
        lhs = gp.quicksum(A1[i, j] * x[j] for j in range(n1))
        rhs_val = b_scenarios[0][1][i]

        if sense == "L":
            master.addConstr(lhs <= rhs_val, name=f"first_L_{i}")
        elif sense == "E":
            master.addConstr(lhs == rhs_val, name=f"first_E_{i}")
        elif sense == "G":
            master.addConstr(lhs >= rhs_val, name=f"first_G_{i}")

    for k in range(K):
        dummy_constr = gp.quicksum(y[k][j] for j in range(n2)) <= 1e6
        master.addConstr(dummy_constr, name=f"dummy_bound_{k}")

    # === FAT BENDERS LOOP ===
    scenario_constrs = {k: set() for k in range(K)}
    converged = False
    iteration = 0

    while not converged and iteration < max_iters:
        iteration += 1
        print(f"\nIteration {iteration}")
        master.optimize()

        if master.status != GRB.OPTIMAL:
            print("Master problem not solved to optimality.")
            print("Model Status:", master.status)
            break

        x_val = np.array([x[i].X for i in range(n1)])
        converged = True

        for k, (prob, b_k) in enumerate(b_scenarios):
            y_val = np.array([y[k][j].X for j in range(n2)])
            lhs_val = A1 @ x_val + A2 @ y_val
            violated = False

            for i in range(m):
                if i in scenario_constrs[k]:
                    continue  # constraint already added

                lhs = lhs_val[i]
                rhs = b_k[i]
                sense = row_sense[i]

                if (
                    (sense == "L" and lhs > rhs + tol)
                    or (sense == "G" and lhs < rhs - tol)
                    or (sense == "E" and abs(lhs - rhs) > tol)
                ):
                    violated = True
                    scenario_constrs[k].add(i)

                    # Add this scenario constraint
                    lhs_expr = (
                        gp.quicksum(A1[i, j] * x[j] for j in range(n1)) +
                        gp.quicksum(A2[i, j] * y[k][j] for j in range(n2))
                    )

                    if sense == "L":
                        master.addConstr(lhs_expr <= rhs, name=f"fat_scen_{k}_L_{i}")
                    elif sense == "E":
                        master.addConstr(lhs_expr == rhs, name=f"fat_scen_{k}_E_{i}")
                    elif sense == "G":
                        master.addConstr(lhs_expr >= rhs, name=f"fat_scen_{k}_G_{i}")

            if violated:
                converged = False

    return {
        "x": np.array([x[i].X for i in range(n1)]),
        "y": {k: np.array([y[k][j].X for j in range(n2)]) for k in range(K)},
        "objective": master.objVal,
        "iterations": iteration
    }


cor_data = parse_cor("airlift/AIRL.cor")
tim_data = parse_tim("airlift/AIRL.tim")
sto_first = parse_sto_dep("airlift/AIRL.sto.first", cor_data["row_names"], cor_data["b"])
sto_second = parse_sto_indep("airlift/AIRL.sto.second", cor_data["row_names"], cor_data["b"])

split = split_stages(
    A=cor_data["A"],
    b=cor_data["b"],
    c=cor_data["c"],
    var_names=cor_data["var_names"],
    variable_stage=tim_data
)

# slim_benders_result_first = benders_slim(
#     A1=split["A1"],
#     A2=split["A2"],
#     b_scenarios=sto_first,  # or sto_second
#     c1=split["c1"],
#     c2=split["c2"],
#     row_sense=cor_data["row_sense"]
# )

# slim_benders_result_second = benders_slim(
#     A1=split["A1"],
#     A2=split["A2"],
#     b_scenarios=sto_second,
#     c1=split["c1"],
#     c2=split["c2"],
#     row_sense=cor_data["row_sense"]
# )


fat_benders_result_first = benders_fat(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_first,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"]
)

fat_benders_result_second = benders_fat(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_second,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"]
)


# print("\n---Slim Benders decomposition results for .sto.first---")
# print("\nOptimal objective:", slim_benders_result_first["objective"])
# print("First-stage solution:", slim_benders_result_first["x"])
# print('Number of iterations:', slim_benders_result_first["iterations"])

# print("\n---Slim Benders decomposition results for .sto.second---")
# print("\nOptimal objective:", slim_benders_result_second["objective"])
# print("First-stage solution:", slim_benders_result_second["x"])
# print('Number of iterations:', slim_benders_result_second["iterations"])

print("\n---Fat Benders decomposition results for .sto.first---")
print("\nOptimal objective:", fat_benders_result_first["objective"])
print("First-stage solution:", fat_benders_result_first["x"])
print('Number of iterations:', fat_benders_result_first["iterations"])

print("\n---Fat Benders decomposition results for .sto.second---")
print("\nOptimal objective:", fat_benders_result_second["objective"])
print("First-stage solution:", fat_benders_result_second["x"])
print('Number of iterations:', fat_benders_result_second["iterations"])