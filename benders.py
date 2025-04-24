import numpy as np
import gurobipy as gp
from gurobipy import GRB

def solve_subproblem(x_val, T2, A2, b2, c2, sense2):
    """
    Solve the second-stage subproblem:
        minimize c2^T y
        s.t. A2 y {<=,=,>=} b2 - T2 x_val
    """
    m2, n2 = A2.shape
    y_model = gp.Model("BendersSubproblem")
    y_model.setParam("OutputFlag", 0)

    # Variables
    y = y_model.addVars(n2, lb=0.0, vtype=GRB.CONTINUOUS, name="y")

    # Residual RHS: b2 - T2 @ x_val
    rhs_resid = b2 - T2.dot(x_val)

    # Recourse constraints
    constraints = []
    for i in range(m2):
        expr = gp.LinExpr()
        for j in range(n2):
            if A2[i, j] != 0:
                expr.add(A2[i, j]*y[j])
        sense = sense2[i]
        if sense == "L":
            constraints.append(y_model.addConstr(expr <= rhs_resid[i]))
        elif sense == "E":
            constraints.append(y_model.addConstr(expr == rhs_resid[i]))
        elif sense == "G":
            constraints.append(y_model.addConstr(expr >= rhs_resid[i]))
        else:
            raise ValueError(f"Unknown sense: {sense}")

    # Objective: minimize c2^T y
    obj = gp.quicksum(c2[j] * y[j] for j in range(n2) if c2[j] != 0)
    y_model.setObjective(obj, GRB.MINIMIZE)
    y_model.optimize()

    if y_model.status == GRB.OPTIMAL:
        return {
            'status': 'optimal',
            'value': y_model.ObjVal,
            'duals': np.array([con.Pi for con in constraints])
        }
    elif y_model.status == GRB.INFEASIBLE:
        return {'status': 'infeasible'}
    else:
        raise RuntimeError(f"Subproblem status: {y_model.status}")

def benders_slim(A1, A2, b_scenarios, c1, c2, row_sense, row_names=None, tol=1e-6, max_iters=50):
    """
    Solve the Benders decomposition problem using the slim Benders approach.
    
    Arguments: 
        A1: Coefficient matrix for first-stage variables (dense)
        A2: Coefficient matrix for second-stage variables (dense),
        b_scenarios: List of tuples (probability, b^k) for each scenario. b1 is the first m rows of b^k.
        c1: Coefficients for first-stage objective function
        c2: Coefficients for second-stage objective function
        row_sense: List of constraint senses ('L', 'E', 'G')
        row_names: List of constraint names
        tol: Tolerance for convergence
        max_iters: Maximum number of iterations
    Returns:
        Mapping with optimal values of x and theta, objective value, number of iterations,
        and the number of constraints added to the master problem.
    """

    # Dimensions
    m, n1 = A1.shape
    _, n2 = A2.shape

    # identify first-stage rows
    first_rows = [i for i,name in enumerate(row_names) if 'HOURS' in name.upper()]
    second_rows = [i for i in range(m) if i not in first_rows]

    m1, m2 = len(first_rows), len(second_rows)

    # slice matrices
    A1_master = A1[first_rows, :]
    T2_sub   = A1[second_rows, :]
    A2_sub   = A2[second_rows, :]
    c2_sub   = c2

    sense1 = [row_sense[i] for i in first_rows]
    sense2 = [row_sense[i] for i in second_rows]

    # first-stage RHS (b1)
    prob0, full_b0 = b_scenarios[0]
    b1 = np.array([full_b0[i] for i in first_rows])
    
    print(f"Building Benders decomposition model with {m1} first-stage and {m2} second-stage constraints, {n1} first-stage variables and {n2} second-stage variables.")
    
    # Step 1: Initialize master model
    master = gp.Model("SlimBendersMaster")
    master.setParam("OutputFlag", 0)

    # Variables: x (first-stage), theta (estimate of recourse cost)
    x = master.addVars(n1, lb=0.0, vtype=GRB.CONTINUOUS, name="x")
    theta = master.addVar(lb=0,vtype=GRB.CONTINUOUS, name="theta")

    # Objective: c1^T x + theta
    master.setObjective(gp.quicksum(c1[i] * x[i] for i in range(n1)) + theta, GRB.MINIMIZE)

    # First-stage constraints: A1 x <= b1
    for idx, i in enumerate(first_rows):
        expr = gp.quicksum(A1_master[idx,j] * x[j] for j in range(n1))
        s = sense1[idx]
        if s=='L': 
            master.addConstr(expr <= b1[idx])
        elif s=='E': 
            master.addConstr(expr == b1[idx])
        elif s=='G': 
            master.addConstr(expr >= b1[idx])
        else:
            raise ValueError(f"Invalid row sense {row_sense[i]}")
    
    iteration = 0
    cuts_added = 0
    UB = 1e10     # large but finite
    LB = -1e10    # small but finite
    
    print("Starting Benders decomposition loop")
    while iteration < max_iters and abs((UB - LB) / (abs(UB) + 1e6)) > tol:
        iteration += 1
        if iteration % 25 == 0:
            print(f"Iteration {iteration}:")
            print(f"Upper Bound: {UB:.4f}")
            print(f"Lower Bound: {LB:.4f}")

        master.optimize()

        if master.status != GRB.OPTIMAL:
            print("Master problem not solved to optimality.")
            print("Model Status:", master.status)
            if int(master.status) == 3:
                print("Infeasibility detected. Check model constraints.")
                master.computeIIS()
                master.write("infeasible.ilp")
                master.write("benders_slim_model.lp")
                print("Infesability report written to infeasible.ilp, model written to benders_slim_model.lp")
            return None
        
        x_val = np.array([x[i].X for i in range(n1)])
        theta_val = theta.X
        expected_rec_cost = 0.0
        pi_bar = np.zeros(m2)
        rhs_bar = 0

        for prob, full_b in b_scenarios:
            b2 = np.array([full_b[i] for i in second_rows])
            sub_res = solve_subproblem(x_val, T2_sub, A2_sub, b2, c2, sense2)
            if sub_res["status"] == "optimal":
                expected_rec_cost += prob * sub_res["value"]
                pi_bar += sub_res["duals"] * prob  # dual variables for the subproblem constraints
                # compute the constant term    
                rhs_bar += prob * (sub_res["duals"] @ b2)
            elif sub_res["status"] == "infeasible":
                print(f"Subproblem infeasible for scenario with probability {prob}.")
                return None
            else:
                print('I don\'t know what the problem is')
        lhs_bar = gp.quicksum(-pi_bar[i] * gp.quicksum(T2_sub[i,j] * x[j] for j in range(n1)) for i in range(m2))
        master.addConstr(theta >= rhs_bar + lhs_bar, name=f"benders_cut_{iteration}")
        cuts_added += 1

        LB = master.objVal
        UB = sum(c1[i] * x_val[i] for i in range(n1)) + expected_rec_cost

    return {
        "x": x_val,
        "theta": theta_val,
        "objective": master.objVal,
        "iterations": iteration,
        "constraints_added": cuts_added
    }
    
def benders_fat(A1, A2, b_scenarios, c1, c2, row_sense, row_names, tol=1e-4, max_iters=50):
    """
    Solve the Benders decomposition problem using the fat Benders approach.
    
    Arguments:
        A1: Coefficient matrix for first-stage variables (dense)
        A2: Coefficient matrix for second-stage variables (dense)
        b_scenarios: List of tuples (probability, b^k) for each scenario
        c1: Coefficients for first-stage objective function
        c2: Coefficients for second-stage objective function
        row_sense: List of constraint senses ('L', 'E', 'G')
        row_names: List of constraint names
        tol: Tolerance for convergence
        max_iters: Maximum number of iterations
    Returns:
        Mapping with optimal values of x and y, objective value, number of iterations,
        and the number of constraints added to the master problem.
    """
    
    # Reset Gurobi model
    # gp.Model.reset()
    
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
        if iteration % 10 == 0:
            print(f"Iteration {iteration}:")
            print(f"Objective: {master.objVal:.4f}")
        master.optimize()

        if master.status != GRB.OPTIMAL:
            print("Master problem not solved to optimality.")
            print("Model Status:", master.status)
            if int(master.status) == 3:
                print("Infeasibility detected. Check model constraints.")
                master.computeIIS()
                master.write("infeasible.ilp")
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
        "iterations": iteration,
        "constraints_added": sum(len(v) for v in scenario_constrs.values())
    }