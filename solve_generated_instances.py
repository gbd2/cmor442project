from benders import benders_fat, benders_slim
from reader import parse_cor, parse_tim, parse_sto_dep, parse_sto_indep, split_stages
from generator import write_smps
import time
import os

os.makedirs('test', exist_ok=True)

SMALL_NAME="small_test_instance"
MEDIUM_NAME="medium_test_instance"
LARGE_NAME = "large_test_instance"

print("Generating small test instance...")
write_smps(I=10, J=10, S=25, folder ='test', indep=False, base_name=SMALL_NAME)
print("Generating medium test instance...")
write_smps(I=20, J=20, S=50, folder ='test', indep=False, base_name=MEDIUM_NAME)
print("Generating large test instance...")
write_smps(I=30, J=30, S=100, folder ='test', indep=False, base_name=LARGE_NAME)

small_cor_data = parse_cor('test/' + SMALL_NAME + ".cor")
small_tim_data = parse_tim('test/' + SMALL_NAME + ".tim")
small_sto_first = parse_sto_dep('test/' + SMALL_NAME + ".sto.first", small_cor_data["row_names"], small_cor_data["b"])

small_split = split_stages(
    A=small_cor_data["A"],
    b=small_cor_data["b"],
    c=small_cor_data["c"],
    var_names=small_cor_data["var_names"],
    variable_stage=small_tim_data
)

# Write a1, a2, b_scenarios, c1, c2 to files
# with open("ex.data", "w") as f:
#     f.write("A1\n")
#     for row in small_split["A1"]:
#         f.write(" ".join(map(str, row)) + "\n")
#     f.write("\nA2\n")
#     for row in small_split["A2"]:
#         f.write(" ".join(map(str, row)) + "\n")
#     f.write("\nb_scenarios\n")
#     for prob, b in small_sto_first:
#         f.write(f"{prob} " + " ".join(map(str, b)) + "\n")
#     f.write("\nc1\n")
#     f.write(" ".join(map(str, small_split["c1"])) + "\n")
#     f.write("\nc2\n")
#     f.write(" ".join(map(str, small_split["c2"])) + "\n")
        
medium_cor_data = parse_cor('test/' + MEDIUM_NAME + ".cor")
medium_tim_data = parse_tim('test/' + MEDIUM_NAME + ".tim")
medium_sto_first = parse_sto_dep('test/' + MEDIUM_NAME + ".sto.first", medium_cor_data["row_names"], medium_cor_data["b"])

medium_split = split_stages(
    A=medium_cor_data["A"],
    b=medium_cor_data["b"],
    c=medium_cor_data["c"],
    var_names=medium_cor_data["var_names"],
    variable_stage=medium_tim_data
)

large_cor_data = parse_cor('test/' + LARGE_NAME + ".cor")
large_tim_data = parse_tim('test/' + LARGE_NAME + ".tim")
large_sto_first = parse_sto_dep('test/' + LARGE_NAME + ".sto.first", large_cor_data["row_names"], large_cor_data["b"])


large_split = split_stages(
    A=large_cor_data["A"],
    b=large_cor_data["b"],
    c=large_cor_data["c"],
    var_names=large_cor_data["var_names"],
    variable_stage=large_tim_data
)

time_slim_small = time.time()
slim_benders_result_small = benders_slim(
    A1=small_split["A1"],
    A2=small_split["A2"],
    b_scenarios=small_sto_first,
    c1=small_split["c1"],
    c2=small_split["c2"],
    row_sense=small_cor_data["row_sense"],
    row_names=small_cor_data["row_names"]
)
time_slim_small = time.time() - time_slim_small
time_fat_small = time.time()
fat_benders_result_small = benders_fat(
    A1=small_split["A1"],
    A2=small_split["A2"],
    b_scenarios=small_sto_first,
    c1=small_split["c1"],
    c2=small_split["c2"],
    row_sense=small_cor_data["row_sense"],
    row_names=small_cor_data["row_names"]
)
time_fat_small = time.time() - time_fat_small

print("\n---Slim Benders decomposition results for small test instance")
print("\nOptimal objective:", slim_benders_result_small["objective"])
print("First-stage solution:", slim_benders_result_small["x"])
print('Number of iterations:', slim_benders_result_small["iterations"])
print('Time taken:', time_slim_small, 'seconds')

print("\n---Fat Benders decomposition results for the small test instance")
print("\nOptimal objective:", fat_benders_result_small["objective"])
print("First-stage solution:", fat_benders_result_small["x"])
print('Number of iterations:', fat_benders_result_small["iterations"])
print('Time taken:', time_fat_small, 'seconds')

time_slim_medium = time.time()
slim_benders_result_medium = benders_slim(
    A1=medium_split["A1"],
    A2=medium_split["A2"],
    b_scenarios=medium_sto_first,
    c1=medium_split["c1"],
    c2=medium_split["c2"],
    row_sense=medium_cor_data["row_sense"],
    row_names=medium_cor_data["row_names"],
)
time_slim_medium = time.time() - time_slim_medium
time_fat_medium = time.time()
fat_benders_result_medium = benders_fat(
    A1=medium_split["A1"],
    A2=medium_split["A2"],
    b_scenarios=medium_sto_first,
    c1=medium_split["c1"],
    c2=medium_split["c2"],
    row_sense=medium_cor_data["row_sense"],
    row_names=medium_cor_data["row_names"],
)
time_fat_medium = time.time() - time_fat_medium


print("\n---Slim Benders decomposition results for medium test instance")
print("\nOptimal objective:", slim_benders_result_medium["objective"])
print("First-stage solution:", slim_benders_result_medium["x"])
print('Number of iterations:', slim_benders_result_medium["iterations"])
print('Time taken:', time_slim_medium, 'seconds')

print("\n---Fat Benders decomposition results for the medium test instance")
print("\nOptimal objective:", fat_benders_result_medium["objective"])
print("First-stage solution:", fat_benders_result_medium["x"])
print('Number of iterations:', fat_benders_result_medium["iterations"])
print('Time taken:', time_fat_medium, 'seconds')

time_slim_large = time.time()
slim_benders_result_large = benders_slim(
    A1=large_split["A1"],
    A2=large_split["A2"],
    b_scenarios=large_sto_first,
    c1=large_split["c1"],
    c2=large_split["c2"],
    row_sense=large_cor_data["row_sense"],
    row_names=large_cor_data["row_names"],
)
time_slim_large = time.time() - time_slim_large
time_fat_large = time.time()
fat_benders_result_large = benders_fat(
    A1=large_split["A1"],
    A2=large_split["A2"],
    b_scenarios=large_sto_first,
    c1=large_split["c1"],
    c2=large_split["c2"],
    row_sense=large_cor_data["row_sense"],
    row_names=large_cor_data["row_names"]
)
time_fat_large = time.time() - time_fat_large


print("\n---Slim Benders decomposition results for large test instance")
print("\nOptimal objective:", slim_benders_result_large["objective"])
print("First-stage solution:", slim_benders_result_large["x"])
print('Number of iterations:', slim_benders_result_large["iterations"])
print('Time taken:', time_slim_large, 'seconds')

print("\n---Fat Benders decomposition results for the large test instance")
print("\nOptimal objective:", fat_benders_result_large["objective"])
print("First-stage solution:", fat_benders_result_large["x"])
print('Number of iterations:', fat_benders_result_large["iterations"])
print('Time taken:', time_fat_large, 'seconds')

# Write all the data to a file

with open("generated_instance_results.txt", "w") as f:
    f.write("Small test instance results:\n")
    f.write(f"Optimal objective (slim): {slim_benders_result_small['objective']}\n")
    f.write(f"First-stage solution (slim): {slim_benders_result_small['x']}\n")
    f.write(f"Number of iterations (slim): {slim_benders_result_small['iterations']}\n")
    f.write(f"Time taken (slim): {time_slim_small} seconds\n\n")
    f.write(f"Cuts added (slim): {slim_benders_result_small['cuts_added']}\n")
    f.write(f"Optimal objective (fat): {fat_benders_result_small['objective']}\n")
    f.write(f"First-stage solution (fat): {fat_benders_result_small['x']}\n")
    f.write(f"Number of iterations (fat): {fat_benders_result_small['iterations']}\n")
    f.write(f"Time taken (fat): {time_fat_small} seconds\n\n")
    f.write(f"Cuts added (fat): {fat_benders_result_small['cuts_added']}\n\n")

    f.write("Medium test instance results:\n")
    f.write(f"Optimal objective (slim): {slim_benders_result_medium['objective']}\n")
    f.write(f"First-stage solution (slim): {slim_benders_result_medium['x']}\n")
    f.write(f"Number of iterations (slim): {slim_benders_result_medium['iterations']}\n")
    f.write(f"Time taken (slim): {time_slim_medium} seconds\n\n")
    f.write(f"Cuts added (slim): {slim_benders_result_medium['cuts_added']}\n")
    f.write(f"Optimal objective (fat): {fat_benders_result_medium['objective']}\n")
    f.write(f"First-stage solution (fat): {fat_benders_result_medium['x']}\n")
    f.write(f"Number of iterations (fat): {fat_benders_result_medium['iterations']}\n")
    f.write(f"Time taken (fat): {time_fat_medium} seconds\n\n")
    f.write(f"Cuts added (fat): {fat_benders_result_medium['cuts_added']}\n\n")
    
    f.write("Large test instance results:\n")
    f.write(f"Optimal objective (slim): {slim_benders_result_large['objective']}\n")
    f.write(f"First-stage solution (slim): {slim_benders_result_large['x']}\n")
    f.write(f"Number of iterations (slim): {slim_benders_result_large['iterations']}\n")
    f.write(f"Time taken (slim): {time_slim_large} seconds\n\n")
    f.write(f"Cuts added (slim): {slim_benders_result_large['cuts_added']}\n")
    f.write(f"Optimal objective (fat): {fat_benders_result_large['objective']}\n")
    f.write(f"First-stage solution (fat): {fat_benders_result_large['x']}\n")
    f.write(f"Number of iterations (fat): {fat_benders_result_large['iterations']}\n")
    f.write(f"Time taken (fat): {time_fat_large} seconds\n\n")
    f.write(f"Cuts added (fat): {fat_benders_result_large['cuts_added']}\n\n")
    
