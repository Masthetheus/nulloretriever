"Allows testing of trivial extensions functions."

from nulloretriever.analysis.trivial_extensions import *

v1 = 0
v2s = [0,1,2,3,4]
case_choice = int(input("Which case to test?"))

if case_choice == 0:
    k1 = 4
    k2 = 5
    half_k = k1//2
    print("DEBUG: Obter lista v1")
    v1_extensions = v1_possible_extensions(v1, half_k)
    print(v1_extensions)
    v1_bits = obtain_sequence_list_from_v1(v1_extensions[1],3)
    print(v1_extensions[1],v1_bits)
    v2_extensions = []
    print("DEBUG: Obter lista v2's")
    for v2 in v2s:
        obtain_v2_extensions(v2 ,v2_extensions, half_k)
    print(v2_extensions)
    v2_bits = obtain_sequence_list_from_v1(v2_extensions[6],half_k)
    print(v2_extensions[6],v2_bits) 
    print("DEBUG: assign v2 to v1")
    v2s_extensions = assign_v2_to_v1(v2s,half_k)
    print(v2s_extensions)
elif case_choice == 1:
    k1 = 5
    k2 = 6
    half_k = k1//2
    print("DEBUG: Obter lista v1")
    v1_extensions = v1_possible_extensions(v1, half_k)
    print(v1_extensions)
    v1_bits = obtain_sequence_list_from_v1(v1_extensions[1],3)
    print(v1_extensions[1],v1_bits)

