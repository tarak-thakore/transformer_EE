import json

from transformer_ee.train import MVtrainer

#with open("transformer_ee/config/input_DUNE_atmo.json", encoding="UTF-8", mode="r") as f:
#with open("transformer_ee/config/input_DUNE_atmo-4m.json", encoding="UTF-8", mode="r") as f:
with open("/home/jbarrow/MLProject2/transformer_EE/transformer_ee/config/NEW_FullyVectorizedFinalStates/input_GENIEv3-0-6-Honda-Truth-hA-LFG_AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi_MAE.json", encoding="UTF-8", mode="r") as f:
    input_d = json.load(f)


input_d["data_path"]="/exp/dune/app/users/rrichi/MLProject/AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi.csv"
# input_d["model"]["name"] = "Transformer_EE_v4"
input_d["model"]["kwargs"]["nhead"] = 2
input_d["model"]["epochs"] = 100
input_d["model"]["kwargs"]["num_layers"] = 5
#input_d["optimizer"]["name"] = "sgd"
input_d["optimizer"]["name"] = "Adam"
input_d["optimizer"]["kwargs"]["lr"] = 0.001
#input_d["optimizer"]["kwargs"]["momentum"] = 0.9
input_d["save_path"] = "/home/jbarrow/MLProject2/save/model/NEW_FullyVectorizedFinalStates/GENIEv3-0-6-Honda-Truth-hA-LFG_AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi_MAE/"
# input_d["weight"] = {"name": "FlatSpectraWeights", "kwargs": {"maxweight": 5, "minweight": 0.2}}


my_trainer = MVtrainer(input_d)
my_trainer.train()
my_trainer.eval()
