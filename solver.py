import gurobipy as gp
from gurobipy import GRB

def solve_extensive_form(A_ext, b_ext, c_ext, row_sense):
    """
    Solves the full problem to the airlift problem using Gurobi.
    
    Parameters:
        A_ext (numpy.ndarray): The extended constraint matrix.
        b_ext (numpy.ndarray): The right-hand side vector for the constraints.
        c_ext (numpy.ndarray): The coefficients for the objective function.
        row_sense (list): The list of constraint senses ('L', 'E', 'G').
    Returns: 
        solution (list): The optimal solution vector.
        objective (float): The optimal objective value.
    """
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