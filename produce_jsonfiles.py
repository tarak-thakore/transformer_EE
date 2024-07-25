import os
import json


# Define the directory containing CSV files
csv_directory = "/exp/dune/app/users/rrichi/MLProject/"

# Define the path to save JSON configuration files
output_directory = "/home/rrichi/MLProject/transformer_EE/transformer_ee/config/"

maxprong_dict={}
with open('maxprong_dict.txt','r') as file:
    for line in file:
        line=line.strip()
        if line:
            parts=line.rsplit(' ',1)
            if len(parts)==2:
                key,value_str=parts
                value=int(value_str)
                maxprong_dict[key]=value



# Define different types for vector, scalar, and target variables
vector_type = [
        "Final_State_Particles_PDG",
        "Final_State_Particles_Mass",
        "Final_State_Particles_Energy",
        "Final_State_Particles_Momentum_X",
        "Final_State_Particles_Momentum_Y",
        "Final_State_Particles_Momentum_Z",
        "Final_State_Particles_CosTheta",
        "Final_State_Particles_Theta"
]
scalar_types = [
    [
        "tot_fKE",
        "p_tot",
        "P_miss"
    ],
    [  
        "tot_fKE",
        "p_tot",
        "P_miss"
    ],
    [
        
        "tot_fKE",
        "p_tot",
        "P_miss"
    ],
    [
        
        "tot_fKE",
        "p_tot"
    ],
    [
        
        "tot_fKE",
        "p_tot",
        "P_miss"

    ],
    [
        
        "tot_fKE",
        "p_tot",
        "P_miss"
    ],
    [
        
        "tot_fKE",
        "p_tot",
        "P_miss"

    ],
    [
        "tot_fKE",
        "p_tot",
        "P_miss"
    ],
    [
        "tot_fKE",
        "p_tot",
        "P_miss"

    ],
    [
        "tot_fKE",
        "p_tot",
        "P_miss"

    ]

]

target_types = [
    [
        "Nu_Theta",
        "Nu_CosTheta"
    ],
    [
        "Nu_Energy"
    ],
    [
        "Nu_Energy",
        "Nu_Mom_X",
        "Nu_Mom_Y",
        "Nu_Mom_Z"     
    ],

    [
        "Nu_Energy",
        "Nu_Mom_X",
        "Nu_Mom_Y",
        "Nu_Mom_Z",
        "P_miss"
    ],
    [
        "Nu_Energy",
        "Nu_Mom_X",
        "Nu_Mom_Y",
        "Nu_Mom_Z",
        "Nu_Baseline"      
    ],
    [
        "Nu_Energy",
        "Nu_Baseline"    
    ],
    [
        "Nu_Energy",
        "Nu_Theta"

    ],
    [
        "Nu_Theta"
    ],
    [
        "Nu_CosTheta"
    ],
    [
        "Nu_Mom_X",
        "Nu_Mom_Y",
        "Nu_Mom_Z"     
    ]


]

loss_functions=["mean squared error","mean absolute error","mean absolute percentage error"]
acronyms=["MSE","MAE","MAPE"]

# Function to generate a configuration for a given CSV file and type combination
def generate_configuration(csv_file, scalar_type, target_type,loss_function,acronym,num_prong):

    # Join all target types into a single string separated by underscores
    target_type_str = '_'.join(target_type)
    # Create the JSON configuration dictionary
    config = {
        "data_path": csv_file,
        "vector": vector_type,
        "scalar": scalar_type,
        "target": target_type,
        "max_num_prongs": num_prong,
        "batch_size_train": 1024,
        "batch_size_valid": 256,
        "batch_size_test": 3000,
        "test_size": 0.6,
        "valid_size": 0.05,
        "seed": 0,
        "loss": {
            "kwargs": {
                "coefficients": [0.5]*len(target_type),
                "base_loss_names": [loss_function] * len(target_type)
            }
        },
        "optimizer": {
            "name": "Adam",
            "kwargs": {
                "lr": 0.001
            }
        },
        "model": {
            "name": "Transformer_EE_MV",
            "kwargs": {}
        },
        "save_path": f"/home/rrichi/MLProject/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}",
        "model_phys_name": f"GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"
    }

    return config

# List all CSV files in the directory
csv_files = [f for f in os.listdir(csv_directory) if f.endswith('.csv')]

# Generate and save JSON configuration files for each CSV file, scalar type, target type, and loss function combination
for csv_file in csv_files:
    for scalar_type, target_type in zip(scalar_types, target_types):
        for loss_function, acronym in zip(loss_functions, acronyms):
            prong_val=maxprong_dict[csv_file]
            config = generate_configuration(csv_file, scalar_type, target_type, loss_function, acronym,prong_val)
            target_type_str = '_'.join(target_type)
            json_filetojoin = f"input_GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{target_type_str}_{acronym}.json"
            json_filename = os.path.join(output_directory, json_filetojoin)
            with open(json_filename, "w") as json_file:
                json.dump(config, json_file, indent=4)
