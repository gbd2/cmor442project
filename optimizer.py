import numpy as np
from gurobipy import Model, GRB

class AirlineOptimizer:
    def __init__(self, cost_matrix, aircraft_to_hours, capacity_matrix, max_flight_hours, demand_prob_tuple):
        self.cost_matrix = cost_matrix # c_{ij}
        self.aircraft_to_hours = aircraft_to_hours # a_{ij}
        self.capacity_matrix = capacity_matrix # b_{ij}
        self.max_flight_hours = max_flight_hours # F_i
        self.demand_prob_tuple = demand_prob_tuple # (d_j, p(d_{j} = X)) (random variable representing demand for route j)
        self.num_aircraft, self.num_routes = aircraft_to_hours.shape # (I, J)


    # ========= Original Integer Program ============= #
    def Solve_ZI(self):
        I = self.num_aircraft
        J = self.num_routes
        model = Model("Stochastic-Airlies-ZI")


        # ======== First Stage ============ #
        # Decision vector
        x1 = model.addVars(I, J, vtype=GRB.BINARY, name="x1")


        # c1 vector: stacked c_ij and zeros
        c1 = np.concatenate([
            self.cost_matrix.reshape(I * J),  # c_ij coefficients
            np.zeros(I)             # for s_i slack variables
        ])
        
        # b1 vector: All max flight times per aircraft
        b1 = self.max_flight_hours

        # A1 matrix: [A | I]
        A_Left = np.empty()
        for i in range(I*J):
            A_left = np.concatenate((np.zeros(i), self.aircraft_to_hours[i]), np.zeros(I*J - J) axis=None)
            A_right = np.identity(I*I)



        # ======== Second Stage ============ #
        # Binary variables x_{ijk}, where k ∈ {1, 2, ..., J} \ {J}
        x2 = model.addVars(I, J, J-1, vtype=GRB.BINARY, name="x2")

        # Slack variables s_{ij} (assumed continuous and non-negative)
        s = model.addVars(I, J, vtype=GRB.CONTINUOUS, lb=0.0, name="s")

        # Positive and negative deviation variables y^+_j and y^-_j
        y_plus = model.addVars(J, vtype=GRB.CONTINUOUS, lb=0.0, name="y_plus")
        y_minus = model.addVars(J, vtype=GRB.CONTINUOUS, lb=0.0, name="y_minus")

        # Now concatenate into a single list in the order from the image
        # x_ijk with k=1 is skipped (assumed from context: starts from x_{i1,2} onward)
        x2_vars = []

        # Add x_{ijk} for k=2 only, across i and j (x_{ij2})
        for i in range(I):
            for j in range(J):
                for k in range(J):
                    if k != j:
                        x2_vars.append(x[i, j, 1]) # k=1 means index 2 (0-based)

        # Add slack variables s_{ij}
        for i in range(I):
            for j in range(J):
                x2_vars.append(s[i, j])

        # Add y^+_j
        for j in range(J):
            x2_vars.append(y_plus[j])

        # Add y^-_j
        for j in range(J):
            x2_vars.append(y_minus[j])

        # Optional: print variable names in x2 vector
        print("x2 vector variables:")
        for var in x2_vars:
            print(var.VarName)

        
        

    

    # b1 vector: F_i values
    b1 = np.random.rand(I)

    

    # ======== Second Stage ============ #


    
    def Benders(self, model):
        