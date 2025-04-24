import numpy as np
import gurobipy as gp
from gurobipy import GRB

def solve_subproblem_slim(x_val, T2, A2, b2, c2, sense2):
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

def solve_subproblem_fat(x_val, A1, A2, b_s, c2, row_sense, lb2, ub2):
    """
    Solve the second stage subproblem:
        min  c2^T y
        s.t. A1 x + A2 y {<=,=,>=}  b_s
    """
    m, n2 = A2.shape
    mdl = gp.Model("BendersSub")
    mdl.Params.OutputFlag = 0
    mdl.Params.Threads     = 1

    if lb2 is None: lb2 = np.zeros(n2)
    if ub2 is None: ub2 = np.full (n2,  GRB.INFINITY)
    y = mdl.addVars(n2, lb=lb2, ub=ub2, name="y")
    
    resid = b_s - A1 @ x_val
    cons  = []
    for i in range(m):
        expr = gp.quicksum(A2[i, j] * y[j] for j in range(n2))
        s    = row_sense[i]
        if   s == "L": cons.append(mdl.addConstr(expr <= resid[i]))
        elif s == "E": cons.append(mdl.addConstr(expr == resid[i]))
        elif s == "G": cons.append(mdl.addConstr(expr >= resid[i]))
        else: raise ValueError(f"unknown sense {s}")

    mdl.setObjective(gp.quicksum(c2[j] * y[j] for j in range(n2)), GRB.MINIMIZE)
    mdl.optimize()

    if mdl.Status == GRB.OPTIMAL:
        return {"status": "optimal",
                "value" : mdl.ObjVal,
                "duals" : np.array([c.Pi for c in cons])}
    elif mdl.Status == GRB.INFEASIBLE:
        return {"status": "infeasible"}
    raise RuntimeError(f"sub-problem status = {mdl.Status}")

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
            sub_res = solve_subproblem_slim(x_val, T2_sub, A2_sub, b2, c2, sense2)
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
        "cuts_added": cuts_added
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
    
    m, n1 = A1.shape
    K = len(b_scenarios)

    master = gp.Model("FatMaster")
    master.setParam("OutputFlag", 0)

    # Variables
    x = master.addVars(n1, lb=0.0, name="x")
    theta = {k: master.addVar(lb=0.0, name=f"theta_{k}") for k in range(K)}

    # Objective: c1^T x + sum_k p_k * theta_k
    master.setObjective(
        gp.quicksum(c1[j]*x[j] for j in range(n1)) +
        gp.quicksum(b_scenarios[k][0] * theta[k] for k in range(K)),
        GRB.MINIMIZE
    )

    # First-stage constraints (from scenario theta)
    for i in range(m):
        rhs0 = b_scenarios[0][1][i]
        s = row_sense[i]
        lhs = gp.quicksum(A1[i, j]*x[j] for j in range(n1))
        if s == "L": 
            master.addConstr(lhs <= rhs0)
        elif s == "E": 
            master.addConstr(lhs == rhs0)
        elif s == "G": 
            master.addConstr(lhs >= rhs0)

    # Main Benders loop
    UB, LB, iters = 1e12, -1e12, 0
    while iters < max_iters and UB - LB > tol * (1 + abs(UB)):
        iters += 1
        master.optimize()
        if master.status != GRB.OPTIMAL:
            raise RuntimeError("Fat master not optimal")

        x_val = np.array([x[j].X for j in range(n1)])
        expected = 0.0

        for k, (p_k, b_k) in enumerate(b_scenarios):
            res = solve_subproblem_fat(x_val, A1, A2, b_k, c2, row_sense, None, None)
            if res["status"] == "infeasible":
                raise RuntimeError("Feasibility cut generation not implemented.")

            pi = res["duals"]
            expected += p_k * res["value"]

            const_term = float(pi @ b_k)
            coeffs_j = pi @ A1  # ∇(dual) · A1 — one coefficient per x_j
            cut_expr = const_term - gp.quicksum(coeffs_j[j] * x[j] for j in range(n1))
            master.addConstr(theta[k] >= cut_expr, name=f"cut_k{k}_it{iters}")

        UB = master.ObjVal
        LB = sum(c1[j] * x_val[j] for j in range(n1)) + expected

    return {
        "x" : np.array([x[j].X for j in range(n1)]),
        "theta" : {k: theta[k].X for k in range(K)},
        "objective" : master.ObjVal,
        "iterations": iters,
        "cuts_added" : iters * K
    }
