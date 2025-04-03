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
