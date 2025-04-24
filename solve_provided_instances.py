from benders import benders_fat, benders_slim
from solver import solve_extensive_form
from reader import parse_cor, parse_tim, parse_sto_dep, parse_sto_indep, split_stages, build_extensive_form
import time

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

x1_vars = [v for v in cor_data["var_names"] if len(v) == 3 and v.startswith("X")]


# Write a1, a2, b_scenarios, c1, c2 to files

# with open("airlift/AIRL.data", "w") as f:
#     f.write("A1\n")
#     for row in split["A1"]:
#         f.write(" ".join(map(str, row)) + "\n")
#     f.write("\nA2\n")
#     for row in split["A2"]:
#         f.write(" ".join(map(str, row)) + "\n")
#     f.write("\nb_scenarios\n")
#     for prob, b in sto_first:
#         f.write(f"{prob} " + " ".join(map(str, b)) + "\n")
#     f.write("\nc1\n")
#     f.write(" ".join(map(str, split["c1"])) + "\n")
#     f.write("\nc2\n")
#     f.write(" ".join(map(str, split["c2"])) + "\n")
        
full_model_form_first = build_extensive_form(
    A1=split["A1"],
    A2=split["A2"],
    b=cor_data["b"],
    c1=split["c1"],
    c2=split["c2"],
    sto_scenarios=sto_first,
    row_sense=cor_data["row_sense"])

full_model_form_second = build_extensive_form(
    A1=split["A1"],
    A2=split["A2"],
    b=cor_data["b"],
    c1=split["c1"],
    c2=split["c2"],
    sto_scenarios=sto_second,
    row_sense=cor_data["row_sense"])

full_model_first_time = time.time()
full_model_first_solution, full_model_first_objective = solve_extensive_form(
    A_ext=full_model_form_first["A_ext"],
    b_ext=full_model_form_first["b_ext"],
    c_ext=full_model_form_first["c_ext"],
    row_sense=full_model_form_first["row_sense"]
)
full_model_first_time = time.time() - full_model_first_time

full_model_second_time = time.time()
full_model_second_solution, full_model_second_objective = solve_extensive_form(
    A_ext=full_model_form_second["A_ext"],
    b_ext=full_model_form_second["b_ext"],
    c_ext=full_model_form_second["c_ext"],
    row_sense=full_model_form_second["row_sense"]
)
full_model_second_time = time.time() - full_model_second_time

slim_benders_first_time = time.time()
slim_benders_result_first = benders_slim(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_first,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"],
    max_iters=1000
)
slim_benders_first_time = time.time() - slim_benders_first_time

slim_benders_second_time = time.time()
slim_benders_result_second = benders_slim(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_second,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"],
    max_iters=1000
)
slim_benders_second_time = time.time() - slim_benders_second_time

fat_benders_result_first_time = time.time()
fat_benders_result_first = benders_fat(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_first,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"],
    max_iters=1000
)
fat_benders_result_first_time = time.time() - fat_benders_result_first_time

fat_benders_second_time = time.time()
fat_benders_result_second = benders_fat(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_second,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"],
    max_iters=1000
)
fat_benders_second_time = time.time() - fat_benders_second_time

print("\n---Full model results for .sto.first---")
print("\nOptimal objective:", full_model_first_objective)
print("First-stage solution:", full_model_first_solution)
print('Time taken:', full_model_first_time, 'seconds')

print("\n---Full model results for .sto.second---")
print("\nOptimal objective:", full_model_second_objective)
print("First-stage solution:", full_model_second_solution)
print('Time taken:', full_model_second_time, 'seconds')

print("\n---Slim Benders decomposition results for .sto.first---")
print("\nOptimal objective:", slim_benders_result_first["objective"])
print("First-stage solution:", slim_benders_result_first["x"])
print('Number of iterations:', slim_benders_result_first["iterations"])
print('Time taken:', slim_benders_first_time, 'seconds')
print('Cuts added:', slim_benders_result_first["cuts_added"])

print("\n---Slim Benders decomposition results for .sto.second---")
print("\nOptimal objective:", slim_benders_result_second["objective"])
print("First-stage solution:", slim_benders_result_second["x"])
print('Number of iterations:', slim_benders_result_second["iterations"])
print('Time taken:', slim_benders_second_time, 'seconds')
print('Cuts added:', slim_benders_result_second["cuts_added"])

print("\n---Fat Benders decomposition results for .sto.first---")
print("\nOptimal objective:", fat_benders_result_first["objective"])
print("First-stage solution:", fat_benders_result_first["x"])
print('Number of iterations:', fat_benders_result_first["iterations"])
print('Time taken:', fat_benders_result_first_time, 'seconds')
print('Cuts added:', fat_benders_result_first["cuts_added"])

print("\n---Fat Benders decomposition results for .sto.second---")
print("\nOptimal objective:", fat_benders_result_second["objective"])
print("First-stage solution:", fat_benders_result_second["x"])
print('Number of iterations:', fat_benders_result_second["iterations"])
print('Time taken:', fat_benders_second_time, 'seconds')
print('Cuts added:', fat_benders_result_second["cuts_added"])