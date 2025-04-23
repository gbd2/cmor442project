from benders import benders_fat, benders_slim
from reader import parse_cor, parse_tim, parse_sto_dep, parse_sto_indep, split_stages
from generator import write_smps

base_name="test_instance"

print("Generating test instance...")
write_smps(I=10, J=10, S=25, indep=False, base_name=base_name)
cor_data = parse_cor("test_instance.cor")
tim_data = parse_tim("test_instance.tim")
print("Test instance generated.")
sto_first = parse_sto_dep("test_instance.sto.first", cor_data["row_names"], cor_data["b"])

split = split_stages(
    A=cor_data["A"],
    b=cor_data["b"],
    c=cor_data["c"],
    var_names=cor_data["var_names"],
    variable_stage=tim_data
)

# Write a1, a2, b_scenarios, c1, c2 to files
with open("ex.data", "w") as f:
    f.write("A1\n")
    for row in split["A1"]:
        f.write(" ".join(map(str, row)) + "\n")
    f.write("\nA2\n")
    for row in split["A2"]:
        f.write(" ".join(map(str, row)) + "\n")
    f.write("\nb_scenarios\n")
    for prob, b in sto_first:
        f.write(f"{prob} " + " ".join(map(str, b)) + "\n")
    f.write("\nc1\n")
    f.write(" ".join(map(str, split["c1"])) + "\n")
    f.write("\nc2\n")
    f.write(" ".join(map(str, split["c2"])) + "\n")
        

slim_benders_result_first = benders_slim(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_first,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"]
)
fat_benders_result_first = benders_fat(
    A1=split["A1"],
    A2=split["A2"],
    b_scenarios=sto_first,
    c1=split["c1"],
    c2=split["c2"],
    row_sense=cor_data["row_sense"],
    row_names=cor_data["row_names"]
)


print("\n---Slim Benders decomposition results for .sto.first---")
print("\nOptimal objective:", slim_benders_result_first["objective"])
print("First-stage solution:", slim_benders_result_first["x"])
print('Number of iterations:', slim_benders_result_first["iterations"])

print("\n---Fat Benders decomposition results for .sto.first---")
print("\nOptimal objective:", fat_benders_result_first["objective"])
print("First-stage solution:", fat_benders_result_first["x"])
print('Number of iterations:', fat_benders_result_first["iterations"])
