import numpy as np
import pandas as pd
from itertools import product
import scipy.sparse as sp
import os

def lognormal_params(mean, std):
    """
    Converts mean and standard deviation to log-normal parameters.
    The log-normal distribution is defined by its mean and standard deviation
    in the log space.
    
    Parameters:
        mean (float): Desired mean of the log-normal distribution.
        std (float): Desired standard deviation of the log-normal distribution.
    Returns:
        mu_n (float): Mean of the underlying normal distribution.
        sigma_n (float): Standard deviation of the underlying normal distribution.
    """
    sigma_n_sq = np.log(1 + (std**2) / (mean**2))
    mu_n = np.log(mean) - sigma_n_sq / 2
    sigma_n = np.sqrt(sigma_n_sq)
    return mu_n, sigma_n


def generate_joint_log_normal_samples(mu_vec, sigma_vec, S):
    """
    Generates S joint log-normal samples for multiple independent demands.

    Parameters:
        mu_vec (array-like): Vector of means for each demand.
        sigma_vec (array-like): Vector of standard deviations for each demand.
        S (int): Number of scenarios to generate.
    Returns:
        samples: list of mappings with keys 'DEMAND1', 'DEMAND2', ..., 'prob'
    """
    
    J = len(mu_vec)
    mu_n, sigma_n = [], []
    samples = []
    for mu, sigma in zip(mu_vec, sigma_vec):
        mu_n_i, sigma_n_i = lognormal_params(mu, sigma)
        mu_n.append(mu_n_i)
        sigma_n.append(sigma_n_i)
        
    # For each scenario, sample all demands to happen at the same time
    for _ in range(S):
        demand_vals = np.random.normal(mu_n, sigma_n)
        demand = np.exp(demand_vals)
        prob = 1.0 / S
        sample = {"prob": prob}
        for j in range(J):
            sample[f"DEMAND{j+1}"] = demand[j]
        samples.append(sample)
    return samples

def generate_indep_log_normal_samples(mu_vec, sigma_vec, N_vec):
    """
    Generates Ni independent samples for DEMANDi for i=1, 2, ..., n
    using log-normal distributions with specified means and stds.

    Parameters:
        mu_vec (array-like): Means for each demand.
        sigma_vec (array-like): Standard deviations for each demand.
        N_vec (array-like): Number of samples for each demand.
    Returns:
        scenarios: list of dicts with keys 'DEMAND1', 'DEMAND2', ..., 'prob'
    """
    J = len(mu_vec)
    samples = []
    for i in range(J):
        mu_n_i, sigma_n_i = lognormal_params(mu_vec[i], sigma_vec[i])
        demand_vals = np.random.normal(mu_n_i, sigma_n_i, N_vec[i])
        demand = np.exp(demand_vals)
        samples.append(demand)
    
    combos = list(product(*samples))
    prob = 1.0 / len(combos)
    
    scenarios = []
    for combo in combos:
        scenario = {"prob": prob}
        for j in range(J):
            scenario[f"DEMAND{j+1}"] = combo[j]
        scenarios.append(scenario)
    
    return scenarios

def write_smps(I, J, S, indep=False, folder=".", base_name="airlift_large", seed=42):
    np.random.seed(seed)
    n1 = I * J
    nrec = I * J * (J - 1)
    nyplusminus = 2 * J
    n2 = nrec + nyplusminus
    m1 = I
    m2 = I * J + J

    # Demand generation
    mu_vec = np.random.uniform(1000, 1500, size=J)
    sigma_vec = np.random.uniform(50, 300, size=J)
    if indep:
        b_scenarios = generate_indep_log_normal_samples(mu_vec, sigma_vec, [S] * J)
    else:
        b_scenarios = generate_joint_log_normal_samples(mu_vec, sigma_vec, S)

    c1 = np.random.uniform(2000, 10000, size=n1)
    switching_costs = [np.random.uniform(2000, 15000) for _ in range(nrec)]
    c_plus = np.random.uniform(300, 600, size=J)
    c_minus = np.zeros(J)
    c2 = np.concatenate([switching_costs, c_plus, c_minus])
    
    # --- Structured matrix and RHS construction ---
    
    a = np.random.uniform(20, 60, size=(I, J)) 
    b = np.random.uniform(20, 75, size=(I, J))
    A1 = np.zeros((I + I * J + J, I * J))
    for i in range(I):
        for j in range(J):
            A1[i, i * J + j] = a[i, j]  # a[i][j] is flight hours needed for aircraft i on route j
            
     # 2. Negative a_ij coefficients (used in recourse availability constraints)
    row = I
    for i in range(I):
        for j in range(J):
            A1[row, i * J + j] = -a[i, j]
            row += 1

    # 3. Positive b_ij coefficients (used in demand constraints)
    for j in range(J):
        for i in range(I):
            A1[I + I * J + j, i * J + j] = b[i, j]

    nrec = I * J * (J - 1)  # number of x_ijk variables
    nyplusminus = 2 * J     # y^+ and y^-
    n2 = nrec + nyplusminus
    A2 = np.zeros((I * J + J, n2))  # I*J availability + J demand constraints

    aijk = {
        (i, j, k): np.random.uniform(15, 60)
        for i in range(I)
        for j in range(J)
        for k in range(J) if k != j
    }
    
    # --- Availability Constraints (I * J) ---
    row = 0
    col = 0
    for i in range(I):
        for j in range(J):
            for k in range(J):
                if k != j:
                    A2[row, col] = aijk[(i, j, k)]  # aijk[i][j][k] is flight hours needed for aircraft i on route j switched to k
                    col += 1
            row += 1
        
    # --- Demand Constraints (J) ---
    # Continuing col from above (should now be nrec)
    # Add negative b_ij * aijk/aij for flights switched away
    # Add positive b_ij for flights switched to
    # Add y+ as +1, y- as -1

    for j in range(J):  # for each demand constraint
        demand_row = I * J + j
        col = 0
        for i in range(I):
            for j2 in range(J):
                for k in range(J):
                    if k == j2:
                        continue
                    # flights switched AWAY from route j
                    if j2 == j:
                        aij = a[i][j2]
                        if aij == 0:
                            continue  # prevent division by zero
                        A2[demand_row, col] -= b[i][j] * aijk[(i, j2, k)] / aij
                    # flights switched TO route j
                    if k == j:
                        A2[demand_row, col] += b[i][j]
                    col += 1
        # Add y+ and y- coefficients
        A2[demand_row, nrec + j] = 1     # y_j^+
        A2[demand_row, nrec + J + j] = -1  # y_j^-

    c_full = np.concatenate([c1, c2])

    # Structured RHS
    flight_hours = np.random.uniform(7000, 8000, size=I)
    rhs_vec = np.concatenate([
        flight_hours,
        np.zeros(I * J),  # AVAIL rows
        np.ones(J) # DEMAND rows
    ])

    # Structured row names and senses
    
    # First-stage variable names and rows
    x1_names = [f"X{i+1}{j+1}" for i in range(I) for j in range(J)]
    row_names_1 = (
    [f"HOURS{i+1}" for i in range(I)] +                     # I rows
    [f"AVAIL{i+1}{j+1}" for i in range(I) for j in range(J)] +  # I*J rows
    [f"DEMAND{j+1}" for j in range(J)]                      # J rows)
    )

    # Second-stage variable names and rows
    xijk_names = [f"X{i+1}{j+1}{k+1}" for i in range(I) for j in range(J) for k in range(J) if j != k]
    yplus_names = [f"YPLUS{j+1}" for j in range(J)]
    yminus_names = [f"YMINUS{j+1}" for j in range(J)]

    x2_names = xijk_names + yplus_names + yminus_names
    row_names_2 = [f"AVAIL{i+1}{j+1}" for i in range(I) for j in range(J)] + [f"DEMAND{j+1}" for j in range(J)]

    # Full variable and row lists
    var_names = x1_names + x2_names
    row_names = row_names_1
    
    row_types = ['L'] * m1 + ['L'] * (I * J) + ['E'] * J
    row_index = {name: i for i, name in enumerate(row_names)}
    c_full = np.concatenate([c1, c2])  # Cost vector must match var_names order

    with open(os.path.join(folder, f"{base_name}.cor"), 'w') as f:
        f.write(f"NAME          {base_name.upper()}\n")
        f.write("ROWS\n")
        f.write(" N  OBJ\n")
        for name, sense in zip(row_names, row_types):
            f.write(f" {sense}  {name}\n")

        f.write("COLUMNS\n")

        # A1 columns
        for j, var in enumerate(x1_names):
            entries = []
            for i, row in enumerate(row_names_1):
                val = A1[i, j]
                if abs(val) > 1e-10:
                    entries.append((row, val))
            if abs(c1[j]) > 1e-10:
                entries.insert(0, ("OBJ", c1[j]))
            for k in range(0, len(entries), 2):
                line = f"    {var:8}"
                for row, val in entries[k:k+2]:
                    line += f" {row:8} {val: .6f}"
                f.write(line + "\n")

        # A2 columns
        for j, var in enumerate(x2_names):
            entries = []
            for i, row in enumerate(row_names_2):
                val = A2[i, j]
                if abs(val) > 1e-10:
                    entries.append((row, val))
            if abs(c2[j]) > 1e-10:
                entries.insert(0, ("OBJ", c2[j]))
            for k in range(0, len(entries), 2):
                line = f"    {var:8}"
                for row, val in entries[k:k+2]:
                    line += f" {row:8} {val: .6f}"
                f.write(line + "\n")

        f.write(" RHS           RIGHT\n")
        for row, val in zip(row_names, rhs_vec):
            if abs(val) > 1e-10:
                f.write(f"    RIGHT     {row:8} {val:.6f}\n")

        f.write("ENDATA\n")

    # Write .tim
    with open(os.path.join(folder, f"{base_name}.tim"), 'w') as f:
        f.write(f"TIME          {base_name.upper()}\n")
        f.write("PERIODS\n")

        # First-stage Xij → HOURS#
        for i in range(I):
            for j in range(J):
                var = f"X{i+1}{j+1}"
                activity = f"HOURS{i+1}"
                f.write(f"    {var:8} {activity:10} PERIOD1\n")

        # Second-stage Xijk → AVAILij
        for i in range(I):
            for j in range(J):
                for k in range(J):
                    if k != j:
                        var = f"X{i+1}{j+1}{k+1}"
                        activity = f"AVAIL{i+1}{j+1}"
                        f.write(f"    {var:8} {activity:10} PERIOD2\n")

        # YPLUSj and YMINUSj → DEMANDj
        for j in range(J):
            f.write(f"    YPLUS{j+1:<4}  DEMAND{j+1:<5} PERIOD2\n")
        for j in range(J):
            f.write(f"    YMINUS{j+1:<3}  DEMAND{j+1:<5} PERIOD2\n")

        f.write("ENDATA\n")


    # Write .sto.first
    with open(os.path.join(folder, f"{base_name}.sto.first"), 'w') as f:
        f.write("STOCH        AIRL\nBLOCKS       DISCRETE\n")
        for s, scenario in enumerate(b_scenarios):
            prob = scenario["prob"]
            f.write(f"BL  B{s+1}     PERIOD2       {prob:.6f}\n")
            tempstr = "    RIGHT"
            for j in range(J):
                row_idx = m1 + I * J + j
                tempstr += f"    {row_names[row_idx]:8} {scenario[f'DEMAND{j+1}']:.6f}"
            tempstr += "\n"
            f.write(tempstr)
        f.write("ENDATA\n")

    # Write .sto.second
    with open(os.path.join(folder, f"{base_name}.sto.second"), 'w') as f:
        f.write("STOCH     AIRL\nINDEP     DISCRETE\n")
        for j in range(J):
            row_idx = m1 + I * J + j
            seen = set()
            for scenario in b_scenarios:
                val = scenario[f"DEMAND{j+1}"]
                prob = scenario["prob"]
                key = (round(val, 6), round(prob, 6))
                if key not in seen:
                    seen.add(key)
                    f.write(f"    RIGHT    {row_names[row_idx]:8} {val:.6f}       PERIOD2 {prob:.6f}\n")
        f.write("ENDATA\n")