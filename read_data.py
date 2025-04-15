import numpy as np

def parse_cor_file(file_path):
    """
    Parses a .cor file in SMPS format and returns the model components.
    
    Args:
        file_path (str): Path to the .cor file.
        
    Returns:
        dict: A dictionary containing the model name, rows, columns, and right-hand side (RHS) values.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    name = ""
    rows = {}
    columns = {}
    rhs = {}
    current_section = None

    # Parse the file line by line
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue  # Skip empty lines and comments

        if line.startswith('NAME'):
            name = line.split()[1]
        elif line.startswith('ROWS'):
            current_section = 'ROWS'
            continue
        elif line.startswith('COLUMNS'):
            current_section = 'COLUMNS'
            continue
        elif line.startswith('RHS'):
            current_section = 'RHS'
            continue
        elif line.startswith('ENDATA'):
            break  # End of data

        if current_section == 'ROWS':
            parts = line.split()
            row_type = parts[0]
            row_name = parts[1]
            rows[row_name] = row_type
        elif current_section == 'COLUMNS':
            parts = line.split()
            if parts[0] not in columns:
                columns[parts[0]] = []
            if len(parts) < 4:
                col_name = parts[0]
                col_type = parts[1]
                col_value = float(parts[2])
                columns[col_name].append((col_type, col_value))
            else:
                col_name = parts[0]
                col_type1 = parts[1]
                col_value1 = float(parts[2])
                col_type2 = parts[3]
                col_value2 = float(parts[4])
                columns[col_name].append((col_type1, col_value1))
                columns[col_name].append((col_type2, col_value2))
        elif current_section == 'RHS':
            parts = line.split()
            rhs_name = parts[1]
            rhs_value = float(parts[2])
            rhs[rhs_name] = rhs_value

    return {
        'name': name,
        'rows': rows,
        'columns': columns,
        'rhs': rhs
    }
    
def parse_tim_file(file_path):
    """
    Parses a .tim file in SMPS format and returns the time periods.
    
    Args:
        file_path (str): Path to the .tim file.
        
    Returns:
        dict: A dictionary containing the model name and periods.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    name = ""
    periods = {}
    current_section = None

    # Parse the file line by line
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue  # Skip empty lines and comments

        if line.startswith('TIME'):
            name = line.split()[1]
        elif line.startswith('PERIODS'):
            current_section = 'PERIODS'
            continue
        elif line.startswith('ENDATA'):
            break  # End of data

        if current_section == 'PERIODS':
            parts = line.split()
            period = parts[-1]
            if period not in periods:
                periods[period] = {}
            col_name = parts[0]
            if col_name not in periods[period]:
                periods[period][col_name] = []
            period_name = parts[1]
            periods[period][col_name].append(period_name)

    return {
        'name': name,
        'periods': periods
    }
    
def parse_sto_first_file(file_path):
    """
    Parses a .sto.first file in SMPS format and returns the stochastic data.
    
    Args:
        file_path (str): Path to the .sto file.
        
    Returns:
        tuple: A tuple containing the stochastic information and blocks.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    blocks = {}
    info = {}
    current_block = None

    # Parse the file line by line
    
    for line in lines:
        line = line.strip().split()
        if line[0] == 'STOCH':
            info['STOCH'] = line[1]
        elif line[0] == 'BLOCKS':
            info['BLOCKS'] = line[1]
        elif line[0] == 'BL':
            current_block = line[1]
            if current_block not in blocks:
                blocks[current_block] = {}
            period = line[2]
            if period not in blocks[current_block]:
                blocks[current_block][period] = []
            prob = float(line[3])
        elif line[0] == 'RIGHT':
            demand1name = line[1]
            demand1value = float(line[2])
            demand2name = line[3]
            demand2value = float(line[4])
            blocks[current_block][period].append({demand1name: demand1value, demand2name: demand2value, 'prob': prob})
        elif line[0] == 'ENDATA':
            break
    return info, blocks

def parse_sto_second_file(file_path):
    """
    Parses a .sto.second file in SMPS format and returns the stochastic data.
    
    Args:
        file_path (str): Path to the .sto file.
        
    Returns:
        tuple: A tuple containing the stochastic information and blocks.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    blocks = {}
    info = {}
    
    for line in lines:
        line = line.strip().split()
        if line[0] == 'STOCH':
            info['STOCH'] = line[1]
        elif line[0] == 'INDEP':
            info['BLOCKS'] = line[1]
        elif line[0] == 'ENDATA':
            break
        else:
            demandname = line[1]
            demandvalue = float(line[2])
            period = line[3]
            prob = float(line[4])
            if period not in blocks:
                blocks[period] = []
            blocks[period].append({demandname: demandvalue, 'prob': prob})
    return info, blocks 

def clean_data(cor_file, tim_file, sto_first_file, sto_second_file):
    """
    Cleans the data and puts it in the format that is ready
    for modeling.
    
    Inputs:
        cor_file: path to the .cor file
        tim_file: path to the .tim file
        sto_first_file: path to the .sto.first file
        sto_second_file: path to the .sto.second file
        
    Outputs: cost_matrix: cost matrix of the aircraft-to-route represented as a numpy array
             aircraft_to_hours: aircraft-to-route hours represented as a numpy array
             capacity_matrix: capacity matrix of the aircraft-to-route represented as a numpy array
             max_flight_hours: maximum flight hours of each aircraft represented as a numpy array
             demand_prob_tuple: tuple of demand and probability for each route
    """
    # Parse the files
    cor_data = parse_cor_file(cor_file)
    sto_first_data = parse_sto_first_file(sto_first_file)
    sto_second_data = parse_sto_second_file(sto_second_file)
    cost_matrix = np.array([])
    aircraft_to_hours = np.array([])
    capacity_matrix = np.array([])
    max_flight_hours = np.array([])
    demand_prob_tuple = []
    
    # Extract data from parsed files
    for colname, i in enumerate(cor_data['columns'].keys()):
        if len(colname == 3):
            # We know that the column is X_{ij} (i.e. aircraft i to route j)
            cost_matrix = np.append(cost_matrix, cor_data['columns'][colname][0][1])
            capacity_matrix = np.append(capacity_matrix, cor_data['columns'][colname][-1][1])
            aircraft_to_hours = np.append(aircraft_to_hours, cor_data['columns'][colname][1][1])
    for colname, i in enumerate(cor_data['rhs'].keys()):
        if len(colname == 6):
            # We know that the column is F_i (i.e. maximum flight hours of aircraft i)
            max_flight_hours = np.append(max_flight_hours, cor_data['rhs'][colname])
            
            
    demand_prob_tuple = []
    for block, periods in sto_first_data[1].items():
        for period, demands in periods.items():
            for demand in demands:
                demand_prob_tuple.append((demand, demand['prob']))

    return cost_matrix, aircraft_to_hours, capacity_matrix, max_flight_hours, demand_prob_tuple

# print(clean_data('airlift/AIRL.cor', 'airlift/AIRL.tim', 'airlift/AIRL.sto.first', 'airlift/AIRL.sto.second'))

cor_results = parse_cor_file('airlift/AIRL.cor')

print(cor_results['columns'].keys())