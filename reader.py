from collections import defaultdict
import numpy as np
from itertools import product

def parse_cor(cor_file_path):
    with open(cor_file_path, 'r') as f:
        lines = f.readlines()

    # Initialize structures
    rows = []
    row_types = {}      # map from row name → type (N, L, E)
    columns = defaultdict(list)  # var → list of (row, coef)
    obj_row = None
    rhs = {}
    var_names = []
    reading_section = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith("*"):
            continue

        if line.startswith("ROWS"):
            reading_section = "ROWS"
            continue
        elif line.startswith("COLUMNS"):
            reading_section = "COLUMNS"
            continue
        elif line.startswith("RHS"):
            reading_section = "RHS"
            continue
        elif line.startswith("ENDATA"):
            break
        elif line.startswith("NAME"):
            continue

        tokens = line.split()
        if reading_section == "ROWS":
            row_type, row_name = tokens
            row_types[row_name] = row_type
            rows.append(row_name)
            if row_type == "N":
                obj_row = row_name

        elif reading_section == "COLUMNS":
            var = tokens[0]
            var_names.append(var) if var not in var_names else None
            for i in range(1, len(tokens), 2):
                row = tokens[i]
                val = float(tokens[i+1])
                columns[var].append((row, val))

        elif reading_section == "RHS":
            rhs_name = tokens[0]  # usually 'RIGHT'
            for i in range(1, len(tokens), 2):
                row = tokens[i]
                val = float(tokens[i+1])
                rhs[row] = val

    # Final list of variables and constraints
    num_vars = len(var_names)
    num_rows = len([r for r in rows if r != obj_row])

    # Map row/var to indices
    var_to_idx = {var: i for i, var in enumerate(var_names)}
    row_names = [r for r in rows if r != obj_row]
    row_to_idx = {row: i for i, row in enumerate(row_names)}

    # Initialize matrix A, vectors c, b
    A = np.zeros((num_rows, num_vars))
    c = np.zeros(num_vars)
    b = np.zeros(num_rows)
    row_sense = []

    for var, coeffs in columns.items():
        j = var_to_idx[var]
        for row, val in coeffs:
            if row == obj_row:
                c[j] = val
            else:
                i = row_to_idx[row]
                A[i, j] = val

    for row in row_names:
        idx = row_to_idx[row]
        b[idx] = rhs.get(row, 0.0)
        sense = row_types[row]
        row_sense.append(sense)

    return {
        'A': A,
        'b': b,
        'c': c,
        'row_names': row_names,
        'var_names': var_names,
        'row_to_idx': row_to_idx,
        'var_to_idx': var_to_idx,
        'row_sense': row_sense
    }
    
def parse_tim(tim_file_path):
    with open(tim_file_path, 'r') as f:
        lines = f.readlines()

    variable_stage = {}
    in_periods_section = False

    for line in lines:
        line = line.strip()
        if line.startswith("PERIODS"):
            in_periods_section = True
            continue
        elif line.startswith("ENDATA"):
            break
        elif not in_periods_section or not line:
            continue

        tokens = line.split()
        if len(tokens) != 3:
            continue  # skip malformed lines

        var, _, period = tokens
        if period.upper() == "PERIOD1":
            variable_stage[var] = 1
        elif period.upper() == "PERIOD2":
            variable_stage[var] = 2
        else:
            raise ValueError(f"Unknown period: {period}")

    return variable_stage

def split_stages(A, b, c, var_names, variable_stage):
    first_stage_vars = [v for v in var_names if variable_stage.get(v, 2) == 1]
    second_stage_vars = [v for v in var_names if variable_stage.get(v, 2) == 2]

    first_stage_indices = [var_names.index(v) for v in first_stage_vars]
    second_stage_indices = [var_names.index(v) for v in second_stage_vars]

    A1 = A[:, first_stage_indices]   # First-stage columns of A
    A2 = A[:, second_stage_indices] # Second-stage columns of A
    c1 = c[first_stage_indices]
    c2 = c[second_stage_indices]

    return {
        'A1': A1,
        'A2': A2,
        'c1': c1,
        'c2': c2,
        'b': b,
        'first_stage_vars': first_stage_vars,
        'second_stage_vars': second_stage_vars,
        'first_stage_indices': first_stage_indices,
        'second_stage_indices': second_stage_indices
    }

def parse_sto_dep(sto_file_path, row_names, default_rhs):
    scenarios = []
    current_rhs = default_rhs.copy()
    prob = None
    in_block = False

    with open(sto_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("*"):
                continue
            if line.startswith("ENDATA"):
                break

            tokens = line.split()

            if tokens[0] == "STOCH" or tokens[0] == "BLOCKS":
                continue

            elif tokens[0] == "BL":
                # Save previous block
                if prob is not None:
                    scenarios.append((prob, current_rhs.copy()))
                # Start new block
                prob = float(tokens[3])
                current_rhs = default_rhs.copy()

            elif tokens[0] == "RIGHT":
                for i in range(1, len(tokens), 2):
                    row = tokens[i]
                    val = float(tokens[i + 1])
                    if row in row_names:
                        idx = row_names.index(row)
                        current_rhs[idx] = val
                    else:
                        raise ValueError(f"Unknown row name {row} in .sto file")

    # Final scenario
    if prob is not None:
        scenarios.append((prob, current_rhs.copy()))

    return scenarios  # list of tuples: (probability, rhs vector)

def parse_sto_indep(sto_file_path, row_names, default_rhs):
    rhs_values = defaultdict(list)  # row → list of (value, prob)

    with open(sto_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("STOCH", "INDEP", "*")):
                continue
            if line.startswith("ENDATA"):
                break

            tokens = line.split()
            if tokens[0] == "RIGHT":
                row = tokens[1]
                val = float(tokens[2])
                prob = float(tokens[4])
                rhs_values[row].append((val, prob))

    # Cartesian product of independent RHS values

    demand_rows = list(rhs_values.keys())
    scenarios = []

    for value_combo in product(*[rhs_values[r] for r in demand_rows]):
        # value_combo = ((val1, prob1), (val2, prob2), ...)
        b_s = default_rhs.copy()
        total_prob = 1.0

        for (val, prob), row in zip(value_combo, demand_rows):
            idx = row_names.index(row)
            b_s[idx] = val
            total_prob *= prob

        scenarios.append((total_prob, b_s.copy()))

    return scenarios

def build_extensive_form(A1, A2, b, c1, c2, sto_scenarios, row_sense):
    m, n1 = A1.shape
    _, n2 = A2.shape
    K = len(sto_scenarios)

    total_vars = n1 + K * n2
    total_rows = K * m

    # Big A matrix: each block is [A1 | 0 ... A2 ... 0]
    A_ext = np.zeros((total_rows, total_vars))
    b_ext = np.zeros(total_rows)
    c_ext = np.zeros(total_vars)
    row_sense_ext = []

    # First-stage cost
    c_ext[:n1] = c1

    for k, (prob, b_k) in enumerate(sto_scenarios):
        row_start = k * m
        row_end = (k + 1) * m

        # Place A1 block
        A_ext[row_start:row_end, :n1] = A1
        # Place A2 block in scenario's y^k slot
        col_start = n1 + k * n2
        col_end = col_start + n2
        A_ext[row_start:row_end, col_start:col_end] = A2

        # RHS
        b_ext[row_start:row_end] = b_k

        # Scenario-weighted second-stage cost
        c_ext[col_start:col_end] = prob * c2

        # Extend row senses
        row_sense_ext.extend(row_sense)

    return {
        'A_ext': A_ext,
        'b_ext': b_ext,
        'c_ext': c_ext,
        'row_sense': row_sense_ext,
        'n1': n1,
        'n2': n2,
        'K': K
    }
