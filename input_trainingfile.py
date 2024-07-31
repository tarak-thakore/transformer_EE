import json
import os
import argparse
from transformer_ee.train import MVtrainer

csv_directory = "/exp/dune/app/users/rrichi/MLProject/"
json_directory="/home/rrichi/MLProject/transformer_EE/transformer_ee/config/"
acronyms = ["MSE", "MAE", "MAPE"]

def training(csv_file, input_dict, acronym, target_type_str):
    input_dict["data_path"] = os.path.join(csv_directory, csv_file)
    input_dict["model"]["kwargs"]["nhead"] = 2
    input_dict["model"]["epochs"] = 100
    input_dict["model"]["kwargs"]["num_layers"] = 5
    input_dict["optimizer"]["name"] = "Adam"
    input_dict["optimizer"]["kwargs"]["lr"] = 0.001
    input_dict["save_path"] = f"/home/rrichi/MLProject/save/model/LossVars_{target_type_str}/GENIEv3-0-6-Honda-Truth-hA-LFG_{os.path.splitext(csv_file)[0]}_{acronym}"
    
    return input_dict

def main(json_file, csv_file):
    base_jsonfile = os.path.splitext(json_file)[0]
    json_filename=os.path.join(json_directory,json_file)
    for acronym in acronyms:
        if acronym in base_jsonfile:
            with open(json_filename, encoding="UTF-8", mode="r") as f:
                input_d = json.load(f)
                target_string = base_jsonfile.split('_' + acronym)[0].split('pi_', 1)[-1]
                train_dict = training(csv_file, input_d, acronym, target_string)
                my_trainer = MVtrainer(train_dict)
                my_trainer.train()
                my_trainer.eval()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train a model using specified JSON and CSV files.")
    parser.add_argument("json_file", type=str, help="The path to the JSON file.")
    parser.add_argument("csv_file", type=str, help="The path to the CSV file.")
    args = parser.parse_args()

    main(args.json_file, args.csv_file)
