#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <variant>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include "TCanvas.h"
#include "TStyle.h"
#include <tuple>
#include <stdexcept>
#include <filesystem>
#include <algorithm>
// Function to convert a string to the appropriate type
// Function to convert a string to the appropriate type
const double PI = 3.141592653589793;
const double R=6371000;
const double h=15000;

TH1D *Init_Nu_Mom;
TH1D *Init_Nu_Theta; 
TH1D *Init_Nu_CosTheta; 
TH1D*Init_Nu_Phi;
TH2D *Oscillogram;
TH1D *Fin_CC_Lept_Mom;
TH1D *Fin_CC_Lept_Theta;
TH1D *Fin_CC_Lept_CosTheta;
TH1D *Fin_Prot_Mom;
TH1D *Fin_PiPlus_Mom;
TH1D *Fin_PiMinus_Mom;
TH1D *Fin_PiZero_Mom;
TH1D *Fin_Prot_Mult;
TH1D *Fin_PiPlus_Mult;
TH1D *Fin_PiMinus_Mult;
TH1D *Fin_PiZero_Mult;

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// using DictionaryType = std::map<std::string, std::variant<int, float, std::string>>;
// std::variant<int, float, std::string> getValue(const DictionaryType& dictionary, const std::string& key) {
//     const auto& value = dictionary.at(key);

//     if (std::holds_alternative<int>(value)) {
//         return std::get<int>(value);
//     } else if (std::holds_alternative<float>(value)) {
//         return std::get<float>(value);
//     } else if (std::holds_alternative<std::string>(value)) {
//         return std::get<std::string>(value);
//     } else {
//         throw std::runtime_error("Unexpected type encountered for key '" + key + "'.");
//     }
// }


std::variant<int, float, std::string> convert_value(const std::string& value) {
    // Try to convert to integer
    std::istringstream intStream(value);
    int intValue;
    if (intStream >> intValue && intStream.eof()) {
        // //std::cout << "conversion, returning int: " << value << std::endl;
        return intValue;
    }

    // If conversion to integer fails, try to convert to float
    std::istringstream floatStream(value);
    float floatValue;
    if (floatStream >> floatValue && floatStream.eof()) {
        // //std::cout << "conversion, returning float: " << value << std::endl;
        return floatValue;
    }

    // If conversion to both integer and float fails, return the original string
    // //std::cout << "No conversion, returning string: " << value << std::endl;
    return value;
}

// Function to process the input file and return a dictionary
std::map<std::string, std::variant<int, float, std::string>> process_file(const std::string& input_file) {
    std::ifstream infile(input_file);
    std::map<std::string, std::variant<int, float, std::string>> dictionary;

    std::string line;
    while (std::getline(infile, line)) {
        size_t delimiter_pos = line.find(':');
        if (delimiter_pos != std::string::npos) {
            std::string key = line.substr(0, delimiter_pos);
            std::string value = line.substr(delimiter_pos + 1);

            // Trim leading and trailing spaces from key and value
            size_t start = key.find_first_not_of(" \t\r\n");
            size_t end = key.find_last_not_of(" \t\r\n");
            if (start != std::string::npos && end != std::string::npos) {
                key = key.substr(start, end - start + 1);
            }
            start = value.find_first_not_of(" \t\r\n");
            end = value.find_last_not_of(" \t\r\n");
            if (start != std::string::npos && end != std::string::npos) {
                value = value.substr(start, end - start + 1);
            }

            // Convert value and store in dictionary
            dictionary[key] = convert_value(value);
        }
    }

    infile.close();
    return dictionary;
}

std::string variantToString(const std::variant<int, float, std::string>& value) {
    if (std::holds_alternative<int>(value)) {
        return std::to_string(std::get<int>(value));
    } else if (std::holds_alternative<float>(value)) {
        return std::to_string(std::get<float>(value));
    } else if (std::holds_alternative<std::string>(value)) {
        return std::get<std::string>(value);
    }
    throw std::runtime_error("Unsupported type for conversion to string");
}


std::tuple<double, double, double> kinematics(const double tStdHepP4[], int j) {
    double fAbsoluteParticleMomentum = sqrt(pow(tStdHepP4[4 * j], 2) + pow(tStdHepP4[4 * j + 1], 2) + pow(tStdHepP4[4 * j + 2], 2));
    double fInvMass = sqrt(pow(tStdHepP4[4 * j + 3], 2) - pow(fAbsoluteParticleMomentum, 2));
    double fKE = tStdHepP4[4 * j + 3] - fInvMass;

    return std::make_tuple(fAbsoluteParticleMomentum, fInvMass, fKE);
}

//function to print the particle level info to outfile
void writeVectorToFile(std::ofstream& outfile, const std::vector<double>& values) {
    for (int j = 0; j < values.size(); j++) {
        // Be sure to end the vector with a quote
        if (j < values.size() - 1) {
            outfile << std::setprecision(3) << values[j] << ",";
        }
        if (j == values.size() - 1) {
            outfile << std::setprecision(3) << values[j] << "\"";
        }
    }
    outfile << ",\"";
}
void getValue(const std::map<std::string, std::variant<int, float, std::string>>& dictionary,
              const std::string& key,
              std::variant<int, float, std::string>& name) {

    const auto& value = dictionary.at(key);

    // Assigning value to name based on its type
    if (std::holds_alternative<int>(value)) {
        name = std::get<int>(value);
    } else if (std::holds_alternative<float>(value)) {
        name = std::get<float>(value);
    } else if (std::holds_alternative<std::string>(value)) {
        name = std::get<std::string>(value);
    }
}

double getDoubleValue(const std::variant<int, float, std::string>& value) {
    if (std::holds_alternative<int>(value)) {
        return static_cast<double>(std::get<int>(value));
    } else if (std::holds_alternative<float>(value)) {
        return static_cast<double>(std::get<float>(value));
    } 
    throw std::runtime_error("Unsupported type for conversion to double");
}
std::tuple<double, double> kinematics_massless(const double tStdHepP4[], int j) {
    double fAbsoluteParticleMomentum = sqrt(pow(tStdHepP4[4 * j], 2) + pow(tStdHepP4[4 * j + 1], 2) + pow(tStdHepP4[4 * j + 2], 2));
    double fKE = fAbsoluteParticleMomentum;

    return std::make_tuple(fAbsoluteParticleMomentum, fKE);
}

double calc_baseline(const double tStdHepP4[],double fAbsoluteParticleMomentum,int j){
    double theta_z=acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
    double eta=PI-theta_z;
    double baseline = (sqrt(pow((R + h), 2) - pow((R * sin(eta)), 2)) + (R * cos(eta))) / 1000.0;

    return baseline;
}

double phi_nu(const double tStdHepP4[],double fAbsoluteParticleMomentum, int j){
    double num= sqrt(pow(tStdHepP4[4*j],2)+pow(tStdHepP4[4*j+2],2));
    double phi= (180./PI)*acos(num/fAbsoluteParticleMomentum);

    return phi;

}


void lepton_kinematics(const double tStdHepP4[],int j,const int tStdHepPdg[],std::vector<int>& pdgs, std::vector<double>& masses,
                    std::vector<double>& energies, std::vector<double>& pxs,
                    std::vector<double>& pys, std::vector<double>& pzs,std::vector<double>& costheta_arr,std::vector<double>& theta_arr,double& tot_fKE, double& tot_fpx, double& tot_fpy, double& tot_fpz,
                    const std::map<std::string, std::variant<int,float, std::string>>& dictionary){
    
    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
    

    double muon_ke=getDoubleValue(dictionary.at("Muon_KE"));
    double electron_ke=getDoubleValue(dictionary.at("Electron_KE"));

    //std::cout<<"Lepton_mom is "<<fAbsoluteParticleMomentum<<" and fKE_lept is  "<<fKE<< endl;
   

    if ((tStdHepPdg[j] == 13 || tStdHepPdg[j] == -13) && fKE >= muon_ke || 
        (tStdHepPdg[j] == 11 || tStdHepPdg[j] == -11) && fKE >= electron_ke) {
        // //std::cout << "lepton loop for j: " << j << " and status: " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;


        pdgs.push_back(tStdHepPdg[j]);
        masses.push_back(fInvMass);
        energies.push_back(fKE);
        pxs.push_back(tStdHepP4[4 * j]);
        pys.push_back(tStdHepP4[4 * j + 1]);
        pzs.push_back(tStdHepP4[4 * j + 2]);
        costheta_arr.push_back(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
        theta_arr.push_back((180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));

        tot_fKE += fKE;
        tot_fpx += tStdHepP4[4 * j];
        tot_fpy += tStdHepP4[4 * j + 1];
        tot_fpz += tStdHepP4[4 * j + 2];

        Fin_CC_Lept_Mom->Fill(1000. * fAbsoluteParticleMomentum);
        Fin_CC_Lept_CosTheta->Fill(tStdHepP4[4 * j + 1] / fAbsoluteParticleMomentum);
        Fin_CC_Lept_Theta->Fill(acos(tStdHepP4[4 * j + 1] / fAbsoluteParticleMomentum));
    }

    else {


        pdgs.push_back(tStdHepPdg[j]);
        masses.push_back(fInvMass);
        energies.push_back(fKE);
        pxs.push_back(tStdHepP4[4 * j]);
        pys.push_back(tStdHepP4[4 * j + 1]);
        pzs.push_back(tStdHepP4[4 * j + 2]);
        costheta_arr.push_back(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
        theta_arr.push_back((180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));

        tot_fKE += fKE;
        tot_fpx += tStdHepP4[4 * j];
        tot_fpy += tStdHepP4[4 * j + 1];
        tot_fpz += tStdHepP4[4 * j + 2];

        Fin_CC_Lept_Mom->Fill(1000. * fAbsoluteParticleMomentum);
        Fin_CC_Lept_CosTheta->Fill(tStdHepP4[4 * j + 1] / fAbsoluteParticleMomentum);
        Fin_CC_Lept_Theta->Fill(acos(tStdHepP4[4 * j + 1] / fAbsoluteParticleMomentum));
    
    }
}

void finalparticles_info(const double tStdHepP4[],int j,const int tStdHepPdg[],std::vector<int>& pdgs, std::vector<double>& masses,
                    std::vector<double>& energies, std::vector<double>& pxs,
                    std::vector<double>& pys, std::vector<double>& pzs,std::vector<double>& costheta_arr,std::vector<double>& theta_arr,double& tot_fKE, double& tot_fpx, double& tot_fpy, double& tot_fpz,
                    const std::map<std::string, std::variant<int,float, std::string>>& dictionary){

    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);

    //define costheta and theta array


    double prot_ke= getDoubleValue(dictionary.at("Proton_KE"));
    double kaon_ke=getDoubleValue(dictionary.at("K+-_KE"));
    double pion_ke=getDoubleValue(dictionary.at("Pi+-_KE"));
    
        
    if ((tStdHepPdg[j] == 2212 && fKE>=prot_ke) ||
        (tStdHepPdg[j]==211 && fKE>=pion_ke) || 
        (tStdHepPdg[j]==-211 && fKE>=pion_ke) ||
        (tStdHepPdg[j]==111) || 
        ((tStdHepPdg[j]== 321 || tStdHepPdg[j]== -321 || tStdHepPdg[j]== 130 || tStdHepPdg[j]== 310) && fKE>=kaon_ke)) {        
        
        
        
        pdgs.push_back(tStdHepPdg[j]);
        masses.push_back(fInvMass);
        energies.push_back(fKE);
        pxs.push_back(tStdHepP4[4 * j]);
        pys.push_back(tStdHepP4[4 * j + 1]);
        pzs.push_back(tStdHepP4[4 * j + 2]);        
        costheta_arr.push_back(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
        theta_arr.push_back((180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));

        tot_fKE += fKE;
        tot_fpx += tStdHepP4[4 * j];
        tot_fpy += tStdHepP4[4 * j + 1];
        tot_fpz += tStdHepP4[4 * j + 2];
            
        if (tStdHepPdg[j]==2212){
            Fin_Prot_Mom->Fill(1000.*fAbsoluteParticleMomentum);
        }
    
        if (tStdHepPdg[j]==211){
            Fin_PiPlus_Mom->Fill(1000.*fAbsoluteParticleMomentum);
        }

        if (tStdHepPdg[j]==-211){
            Fin_PiMinus_Mom->Fill(1000.*fAbsoluteParticleMomentum);
        }

        if (tStdHepPdg[j]==111){
            Fin_PiZero_Mom->Fill(1000.*fAbsoluteParticleMomentum);
        }
    }
    
}

int max_num=0;

void do_resize(std::vector<int>& pdgs, std::vector<double>& masses,
    std::vector<double>& energies, std::vector<double>& pxs,
    std::vector<double>& pys, std::vector<double>& pzs,std::vector<double>& costheta_arr,std::vector<double>& theta_arr){
    if (pdgs.size() < max_num) 
    {
    pdgs.resize(max_num, 0);
    }

    if (energies.size() < max_num) 
    {
    energies.resize(max_num, 0);
    }
    if (masses.size() < max_num) 
    {
    masses.resize(max_num, 0);
    }
    if (pxs.size() < max_num) 
    {
    pxs.resize(max_num, 0);
    }
    if (pys.size() < max_num) 
    {
    pys.resize(max_num, 0);
    }
    if (pzs.size() < max_num) 
    {
    pzs.resize(max_num, 0);
    }
    if (costheta_arr.size() < max_num) 
    {
    costheta_arr.resize(max_num, 0);
    }
    if (theta_arr.size() < max_num) 
    {
    theta_arr.resize(max_num, 0);
    }


}

   
void testcopy(const std::string& input_file) {
      //Set a maximum event particle list length
    const int NMaxParticlesPerEvent = 1000;


    std::map<std::string, std::variant<int,float, std::string>> dictionary = process_file(input_file);

    // Print the contents of the dictionary in the order they were read
// Print the contents of the dictionary in the order they were read
    // for (const auto& pair : dictionary) {
    //     //std::cout << pair.first << ": ";
    //     if (std::holds_alternative<int>(pair.second)) {
    //         //std::cout << std::get<int>(pair.second);
    //     } else if (std::holds_alternative<float>(pair.second)) {
    //         //std::cout << std::get<float>(pair.second);
    //     } else {
    //         //std::cout << std::get<std::string>(pair.second);
    //     }
    //     //std::cout << std::endl;
    // }


    // Load ROOT file (example usage)
    std::string Input_Root_file = std::get<std::string>(dictionary["INPUT_ROOT_FILE"]);
    TFile *SignalFile = TFile::Open(Input_Root_file.c_str(), "READ");
    TTree *SignalTree = (TTree*)SignalFile->Get("gRooTracker");

    int NumberEntries = SignalTree->GetEntries();

    // Handle Num_events
    auto Num_events_var = dictionary["EventNum"];
    int Num_events;

    if (std::holds_alternative<std::string>(Num_events_var) && std::get<std::string>(Num_events_var) == "All") {
        Num_events = NumberEntries;
    } else if (std::holds_alternative<int>(Num_events_var)) {
        Num_events = std::get<int>(Num_events_var);
        if (Num_events > NumberEntries) {
            //std::cout << "Max reached" << std::endl;
            Num_events = NumberEntries;
        }
    }

    std::vector<std::variant<int, std::string>> Init_PDG_arr;
    const auto& Init_PDG_Values = dictionary.at("Initial_PDG");
    
    if (std::holds_alternative<int>(Init_PDG_Values)) {
        Init_PDG_arr.push_back(std::get<int>(Init_PDG_Values));
    } else if (std::holds_alternative<std::string>(Init_PDG_Values)) {
        const std::string& str_value = std::get<std::string>(Init_PDG_Values);
        if (str_value == "Any") {
            Init_PDG_arr.push_back("Any");
        } else {
            // Split by commas and add each token as an integer
            std::vector<std::string> tokens = split(str_value, ',');
            for (const std::string& token : tokens) {
                Init_PDG_arr.push_back(std::stoi(token));
            }
        }
    }

    // // Output the result
    // //std::cout << "Init_PDG_arr: ";
    // for (const auto& value : Init_PDG_arr) {
    //     if (std::holds_alternative<int>(value)) {
    //         //std::cout << std::get<int>(value) << " ";
    //     }
    //     // } else if (std::holds_alternative<std::string>(value)) {
    //     //     //std::cout << std::get<std::string>(value) << " ";
    //     // }
    // }
    // //std::cout << std::endl;

    // double Initial_Nu_Energy_min;
    // double Initial_Nu_Energy_max;

    // std::variant<int, float, std::string> Initial_Nu_Energy_min;
    // std::variant<int, float, std::string> Initial_Nu_Energy_max;
    std::variant<int, float, std::string> num_proton_value;
    std::variant<int, float, std::string> num_pion_value;
    std::variant<int, float, std::string> Final_lepton_PDG;

    // getValue(dictionary,"Initial_Nu_Energy_min", Initial_Nu_Energy_min);
    // getValue(dictionary,"Initial_Nu_Energy_max",Initial_Nu_Energy_max);

    double Min_Nu_energy= getDoubleValue(dictionary.at("Initial_Nu_Energy_min"));
    double Max_Nu_energy=getDoubleValue(dictionary.at("Initial_Nu_Energy_max"));

    // //std::cout << "Min_Nu: " << Min_Nu_energy << std::endl;
    // //std::cout << "Max_Nu:" << Max_Nu_energy<< std::endl;

    // if (std::holds_alternative<int>(dictionary.at("Initial_Nu_Energy_min"))) {
    //     Initial_Nu_Energy_min = static_cast<double>(std::get<int>(dictionary.at("Initial_Nu_Energy_min")));
    // } else if (std::holds_alternative<float>(dictionary.at("Initial_Nu_Energy_min"))) {
    //     Initial_Nu_Energy_min = static_cast<double>(std::get<float>(dictionary.at("Initial_Nu_Energy_min")));
    // }

    // if (std::holds_alternative<int>(dictionary.at("Initial_Nu_Energy_max"))) {
    //     Initial_Nu_Energy_max = static_cast<double>(std::get<int>(dictionary.at("Initial_Nu_Energy_max")));
    // } else if (std::holds_alternative<float>(dictionary.at("Initial_Nu_Energy_max"))) {
    //     Initial_Nu_Energy_max = static_cast<double>(std::get<float>(dictionary.at("Initial_Nu_Energy_max")));
    // }


    // if (std::holds_alternative<int>(Initial_Nu_Energy_min)) {
    //    //std::cout << "Min_nu(int): " << std::get<int>(Initial_Nu_Energy_min) << std::endl;
    // } else if (std::holds_alternative<float>(Initial_Nu_Energy_min)) {
    //     //std::cout << "Min_nu(float) " << std::get<float>(Initial_Nu_Energy_min) << std::endl;
    // }

    // if (std::holds_alternative<int>(Initial_Nu_Energy_max)) {
    //    //std::cout << "Max_nu(int): " << std::get<int>(Initial_Nu_Energy_max) << std::endl;
    // } else if (std::holds_alternative<float>(Initial_Nu_Energy_max)) {
    //     //std::cout << "Max_nu(float) " << std::get<float>(Initial_Nu_Energy_max) << std::endl;
    // }


    getValue(dictionary,"num_proton",num_proton_value);
    getValue(dictionary,"num_pion",num_pion_value);
    getValue(dictionary,"Final_state_lepton_PDG_code_or_process",Final_lepton_PDG);


    int int_num_prot = -1;
    std::string str_num_prot;
    bool is_num_prot_int = false;

    if (std::holds_alternative<int>(num_proton_value)) {
        int_num_prot = std::get<int>(num_proton_value);
        is_num_prot_int = true;
    } else if (std::holds_alternative<std::string>(num_proton_value)) {
        str_num_prot = std::get<std::string>(num_proton_value);
    }

    int int_num_pion = -1;
    std::string str_num_pion;
    bool is_num_pion_int = false;

    if (std::holds_alternative<int>(num_pion_value)) {
        int_num_pion = std::get<int>(num_pion_value);
        is_num_pion_int = true;
    } else if (std::holds_alternative<std::string>(num_pion_value)) {
        str_num_pion = std::get<std::string>(num_pion_value);
    }

    int int_Final_lepton_PDG = -1;
    std::string str_Final_lepton_PDG;
    bool is_Final_lepton_PDG_int = false;

    if (std::holds_alternative<int>(Final_lepton_PDG)) {
        int_Final_lepton_PDG = std::get<int>(Final_lepton_PDG);
        is_Final_lepton_PDG_int = true;
    } else if (std::holds_alternative<std::string>(Final_lepton_PDG)) {
        str_Final_lepton_PDG = std::get<std::string>(Final_lepton_PDG);
    }

    std::string num_events_str = variantToString(dictionary["EventNum"]);
    std::string proton_num_str = variantToString(dictionary["num_proton"]);
    std::string pion_num_str = variantToString(dictionary["num_pion"]);
    
    std::string add_file = std::get<std::string>(dictionary["OUTPUT_ROOT_FILE"]);
    std::string Output_Root_file=add_file+num_events_str+"_"+proton_num_str+"p"+pion_num_str+"pi.root";

    TFile *treefile = new TFile(Output_Root_file.c_str(),"RECREATE");

    std::string output_name = std::get<std::string>(dictionary["OUTPUT_NAME"]);
    std::string outfile_name= output_name+num_events_str+"_"+ proton_num_str+"p"+pion_num_str+"pi.csv";

    ofstream outfile(outfile_name);

    //define the header for the CSV file
    outfile << "\"Event_Index\",\"Nu_PDG\",\"Nu_Energy\",\"Nu_Mom_X\",\"Nu_Mom_Y\",\"Nu_Mom_Z\",\"Nu_CosTheta\",\"Nu_Theta\",\"Nu_Phi\",\"Nu_Baseline\",\"Final_State_Particles_PDG\",\"Final_State_Particles_Mass\",\"Final_State_Particles_Energy\",\"Final_State_Particles_Momentum_X\",\"Final_State_Particles_Momentum_Y\",\"Final_State_Particles_Momentum_Z\",\"Final_State_Particles_CosTheta\",\"Final_State_Particles_Theta\",\"tot_fKE\",\"p_tot\",\"P_miss\"\n";
    //outfile << "\"Event_Index\",\"Initial_State_Neutrino_PDG\",\"Initial_State_Neutrino_Energy\",\"Initial_State_Neutrino_Momentum_X\",\"Initial_State_Neutrino_Momentum_Y\",\"Initial_State_Neutrino_Momentum_Z\",\"Initial_Neutrino_CosTheta\",\"Initial_Neutrino_Theta\",\"Final_State_Particles_PDG\",\"Final_State_Particles_Mass\",\"Final_State_Particles_Energy\",\"Final_State_Particles_Momentum_X\",\"Final_State_Particles_Momentum_Y\",\"Final_State_Particles_Momentum_Z\",\"Final_State_Particles_CosTheta\",\"Final_State_Particles_Theta\",\"tot_fKE\",\"p_tot\",\"P_miss\"\n";


    // //print the output file
    // //std::cout << "Output_Root_file: " << Output_Root_file << std::endl;
    // //std::cout << "Outfile_name: " << outfile_name << std::endl;

    //directory here

    std::string output_directory=std::get<std::string>(dictionary["OUTPUT_DIR"]);
    std::string directory=output_directory+num_events_str+"_"+ proton_num_str+"p"+pion_num_str+"pi";

    try {
        if (std::filesystem::create_directories(directory)) {
            std::cout << "Directory created successfully: " << directory << std::endl;
        } else {
            std::cout << "Directory already exists or could not be created: " << directory << std::endl;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "An error occurred while creating the directory: " << e.what() << std::endl;
    }

    std::string last_name=output_name+num_events_str+proton_num_str+"p"+pion_num_str+"pi_";

    gStyle->SetStatY(0.9); gStyle->SetStatX(0.9);

    TCanvas *c1 = new TCanvas("canvas","canvas",1500,1000);

    Init_Nu_Mom = new TH1D((last_name + "Init_Nu_Mom").c_str(), (last_name + "Initial Neutrino Momentum, hA 2018").c_str(), 1000, 0., 16000.);
    Init_Nu_Mom->GetXaxis()->SetTitle("Initial Neutrino Momentum [MeV/c]");
    Init_Nu_Mom->GetYaxis()->SetTitle("Counts");
    Init_Nu_Mom->SetLineColor(kGreen);
    Init_Nu_Mom->SetLineWidth(3);

    //Initial neutrino angle with respect to zenith (y-axis in DUNE geometry/coordinate system)
    Init_Nu_Theta = new TH1D((last_name+"Init_Nu_Theta").c_str(),(last_name+"Initial Neutrino #theta, hA 2018").c_str(),1000, 0.,0.);
    Init_Nu_Theta->GetXaxis()->SetTitle("Initial Neutrino #theta [rad]");
    Init_Nu_Theta->GetYaxis()->SetTitle("Counts");
    Init_Nu_Theta->SetLineColor(kGreen);
    Init_Nu_Theta->SetLineWidth(3);

    Init_Nu_CosTheta = new TH1D((last_name+"Init_Nu_CosTheta").c_str(),(last_name+"Initial Neutrino cos(#theta), hA 2018").c_str(),1000,0.,0.);//,16.0);
    Init_Nu_CosTheta->GetXaxis()->SetTitle("Initial Neutrino Cos(#theta)");
    Init_Nu_CosTheta->GetYaxis()->SetTitle("Counts");
    Init_Nu_CosTheta->SetLineColor(kGreen);
    Init_Nu_CosTheta->SetLineWidth(3);

    //for phi_nu
    Init_Nu_Phi = new TH1D((last_name+"Init_Nu_Phi").c_str(),(last_name+"Initial Neutrino #phi, hA 2018").c_str(),1000, 0.,0.);
    Init_Nu_Phi->GetXaxis()->SetTitle("Initial Neutrino #phi[deg])");
    Init_Nu_Phi->GetYaxis()->SetTitle("Counts");
    Init_Nu_Phi->SetLineColor(kGreen);
    Init_Nu_Phi->SetLineWidth(3);
    

    //Make oscillogram plot
    Oscillogram = new TH2D((last_name+"Oscillogram_hA_2018").c_str(),(last_name+"Oscillogram, hA 2018").c_str(),1000,0.,0.,1000,0.,0.);
    Oscillogram->GetXaxis()->SetTitle("cos(#theta_{#nu}) from Homestake Zenith");
    Oscillogram->GetYaxis()->SetTitle("Initial Neutrino Energy [GeV]");

    Fin_CC_Lept_Mom = new TH1D((last_name + "Fin_CC_Lept_Mom").c_str(), (last_name + "Final State CC-Produced Lepton Momentum, hA 2018").c_str(), 1000, 0., 16000.);
    Fin_CC_Lept_Mom->GetXaxis()->SetTitle("Final State CC-Produced Lepton Momentum [MeV/c]");
    Fin_CC_Lept_Mom->GetYaxis()->SetTitle("Counts");
    Fin_CC_Lept_Mom->SetLineColor(kViolet);
    Fin_CC_Lept_Mom->SetLineWidth(3);

    Fin_CC_Lept_Theta = new TH1D((last_name + "Fin_CC_Lept_Theta").c_str(), (last_name + "Final CC Lepton #theta, hA 2018").c_str(), 1000, 0.,0.);
    Fin_CC_Lept_Theta->GetXaxis()->SetTitle("Final CC Lepton #theta [rad]");
    Fin_CC_Lept_Theta->GetYaxis()->SetTitle("Counts");
    Fin_CC_Lept_Theta->SetLineColor(kGreen);
    Fin_CC_Lept_Theta->SetLineWidth(3);

    Fin_CC_Lept_CosTheta = new TH1D((last_name + "Fin_CC_Lept_CosTheta").c_str(), (last_name + "Final CC Lepton cos(#theta), hA 2018").c_str(), 1000, 0., 0.);
    Fin_CC_Lept_CosTheta->GetXaxis()->SetTitle("Final CC Lepton Cos(#theta)");
    Fin_CC_Lept_CosTheta->GetYaxis()->SetTitle("Counts");
    Fin_CC_Lept_CosTheta->SetLineColor(kGreen);
    Fin_CC_Lept_CosTheta->SetLineWidth(3);

    Fin_Prot_Mom = new TH1D((last_name + "Fin_Prot_Mom").c_str(), (last_name + "Final State Proton Momentum, hA 2018").c_str(), 1000, 0., 4500.);

    // Fin_Prot_Mom = new TH1D((last_name + "Fin_Prot_Mom_low").c_str(), (last_name + "Final State Proton Momentum_low, hA 2018").c_str(), 1000, -100, 500);
    Fin_Prot_Mom->GetXaxis()->SetTitle("Final State Proton Momentum [MeV/c]");
    Fin_Prot_Mom->GetYaxis()->SetTitle("Counts");
    Fin_Prot_Mom->SetLineColor(kRed);
    Fin_Prot_Mom->SetLineWidth(3);

    Fin_PiPlus_Mom = new TH1D((last_name + "Fin_PiPlus_Mom").c_str(), (last_name + "Final State #pi^{+} Momentum, hA 2018").c_str(), 1000, 0., 4500.);
    // Fin_PiPlus_Mom = new TH1D((last_name + "Fin_PiPlus_Mom_low").c_str(), (last_name + "Final State #pi^{+} Momentum_low, hA 2018").c_str(), 1000, -100, 500);
    Fin_PiPlus_Mom->GetXaxis()->SetTitle("Final State #pi^{+} Momentum [MeV/c]");
    Fin_PiPlus_Mom->GetYaxis()->SetTitle("Counts");
    Fin_PiPlus_Mom->SetLineColor(kRed+1);
    Fin_PiPlus_Mom->SetLineWidth(3);

    Fin_PiMinus_Mom = new TH1D((last_name + "Fin_PiMinus_Mom").c_str(), (last_name + "Final State #pi^{-} Momentum, hA 2018").c_str(), 1000, 0., 4500.);
    // Fin_PiMinus_Mom = new TH1D((last_name + "Fin_PiMinus_Mom_low").c_str(), (last_name + "Final State #pi^{-} Momentum_Low, hA 2018").c_str(), 1000, -100, 500);
    Fin_PiMinus_Mom->GetXaxis()->SetTitle("Final State #pi^{-} Momentum [MeV/c]");
    Fin_PiMinus_Mom->GetYaxis()->SetTitle("Counts");
    Fin_PiMinus_Mom->SetLineColor(kBlue+1);
    Fin_PiMinus_Mom->SetLineWidth(3);

    Fin_PiZero_Mom = new TH1D((last_name + "Fin_PiZero_Mom").c_str(), (last_name + "Final State #pi^{0} Momentum, hA 2018").c_str(), 1000, 0., 4500.);
    // Fin_PiZero_Mom = new TH1D((last_name + "Fin_PiZero_Mom_low").c_str(), (last_name + "Final State #pi^{0} Momentum_Low, hA 2018").c_str(), 1000, -100,500);
    Fin_PiZero_Mom->GetXaxis()->SetTitle("Final State #pi^{0} Momentum [MeV/c]");
    Fin_PiZero_Mom->GetYaxis()->SetTitle("Counts");
    Fin_PiZero_Mom->SetLineColor(kGreen+1);
    Fin_PiZero_Mom->SetLineWidth(3);

    Fin_Prot_Mult = new TH1D((last_name + "Fin_Prot_Mult").c_str(), (last_name + "Final State Proton Multiplicity, hA 2018").c_str(), 20, -0.5, 19.5);
    Fin_Prot_Mult->GetXaxis()->SetTitle("Final State Proton Multiplicity");
    Fin_Prot_Mult->GetYaxis()->SetTitle("Events");
    Fin_Prot_Mult->SetLineColor(kRed);
    Fin_Prot_Mult->SetLineWidth(3);

    Fin_PiPlus_Mult = new TH1D((last_name + "Fin_PiPlus_Mult").c_str(), (last_name + "Final State #pi^{+} Multiplicity, hA 2018").c_str(), 7, -0.5, 6.5);
    Fin_PiPlus_Mult->GetXaxis()->SetTitle("Final State #pi^{+} Multiplicity");
    Fin_PiPlus_Mult->GetYaxis()->SetTitle("Events");
    Fin_PiPlus_Mult->SetLineColor(kRed+1);
    Fin_PiPlus_Mult->SetLineWidth(3);

    Fin_PiMinus_Mult = new TH1D((last_name + "Fin_PiMinus_Mult").c_str(), (last_name + "Final State #pi^{-} Multiplicity, hA 2018").c_str(), 7, -0.5, 6.5);
    Fin_PiMinus_Mult->GetXaxis()->SetTitle("Final State #pi^{-} Multiplicity");
    Fin_PiMinus_Mult->GetYaxis()->SetTitle("Events");
    Fin_PiMinus_Mult->SetLineColor(kBlue+1);
    Fin_PiMinus_Mult->SetLineWidth(3);

    Fin_PiZero_Mult = new TH1D((last_name + "Fin_PiZero_Mult").c_str(), (last_name + "Final State #pi^{0} Multiplicity, hA 2018").c_str(), 7, -0.5, 6.5);
    Fin_PiZero_Mult->GetXaxis()->SetTitle("Final State #pi^{0} Multiplicity");
    Fin_PiZero_Mult->GetYaxis()->SetTitle("Events");
    Fin_PiZero_Mult->SetLineColor(kGreen+1);
    Fin_PiZero_Mult->SetLineWidth(3);

    int tStdHepN = 0; //the num of particles in an event
    int tStdHepStatus[NMaxParticlesPerEvent] = {0}; // an array with the all the number of elements in the tStdHepStatus set to 0
    int tStdHepPdg[NMaxParticlesPerEvent] = {0};
    double tStdHepP4[4*NMaxParticlesPerEvent] = {0.};
    double tStdHepX4[4*NMaxParticlesPerEvent] = {0.};

    SignalTree->SetBranchAddress("StdHepN",&tStdHepN);
    SignalTree->SetBranchAddress("StdHepStatus",&tStdHepStatus);
    SignalTree->SetBranchAddress("StdHepPdg",&tStdHepPdg);
    SignalTree->SetBranchAddress("StdHepP4",&tStdHepP4);
    SignalTree->SetBranchAddress("StdHepX4",&tStdHepX4);

    std::vector<double> energies, masses, pxs, pys, pzs,costheta_arr,theta_arr;
    std::vector<int> pdgs;
    

    // Reading different kinetic energies from the dictionary
    double muon_ke=getDoubleValue(dictionary.at("Muon_KE"));
    double electron_ke=getDoubleValue(dictionary.at("Electron_KE"));
    double prot_ke= getDoubleValue(dictionary.at("Proton_KE"));
    double kaon_ke=getDoubleValue(dictionary.at("K+-_KE"));
    double pion_ke=getDoubleValue(dictionary.at("Pi+-_KE"));

    //For printing different KEs:
    // //std::cout << "Proton KE: " << prot_ke << std::endl;
    // //std::cout << "Pi+- KE: " << pion_ke<< std::endl;
    // //std::cout << "K+- KE: " << kaon_ke << std::endl;

    // //std::cout << "Muon KE: " << muon_ke << std::endl;
    // //std::cout << "Electron KE: " << electron_ke << std::endl;
    // int max_num=0;
    for(int i=0; i<Num_events; i++)
    {   //Initializing multiplicity variables
        SignalTree->GetEntry(i);
        
        if ( abs( tStdHepPdg[0] - tStdHepPdg[4] ) > 1 ) continue;

        int n_prot=0, n_piplus=0, n_piminus=0, n_pizero=0, n_pi=0, visible_count=0;

        energies.clear(); masses.clear(); pxs.clear(); pys.clear(); pzs.clear();costheta_arr.clear();theta_arr.clear();
        pdgs.clear();

        double tot_fpx=0,tot_fpy=0,tot_fpz=0;
        double tot_fKE=0;

        bool skip_event = false;
        
        //Print out the current event number to screen 
        // //std::cout << "The current event number is " << i << " of " << NumberEntries << endl;
        //Get the event entry information
        
        

        // if(i%100000==0){std::cout << "The current event number is " << i << " of " << NumberEntries << endl;}  
        // //std::cout << "The current event number is " << i << " of " << NumberEntries << endl;            
        // //std::cout << "secondly, Number of elements in the array: " << tStdHepN/*leng*/ << std::endl; 
        int tStdHepPdg_nu, lepton_pdg_instore;
        double fAbsoluteParticleMomentum_nu,fKE_nu,fAbsoluteParticleMomentum_lept,fKE_lept;


        for (int j=0;j<tStdHepN;j++){

            if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                tStdHepPdg_nu=tStdHepPdg[j];
                auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                fAbsoluteParticleMomentum_nu=fAbsoluteParticleMomentum;
                fKE_nu=fKE;

            }

            if (tStdHepStatus[j] == 1){
                if (tStdHepPdg[j] == 11 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -15 || tStdHepPdg[j] == 16 || tStdHepPdg[j] == -16){                
                    lepton_pdg_instore=tStdHepPdg[j];
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    fAbsoluteParticleMomentum_lept=fAbsoluteParticleMomentum;
                    fKE_lept=fKE;
                    if ((tStdHepPdg[j]==13 && fKE_lept>=muon_ke) || (tStdHepPdg[j]==-13 && fKE_lept>=muon_ke) || (tStdHepPdg[j]==11 && fKE_lept>=electron_ke) || (tStdHepPdg[j]==-11 && fKE_lept>=electron_ke)){
                        visible_count++;
                    }

                }

                if (tStdHepPdg[j] == 2212){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    if (fKE>=prot_ke){
                        visible_count++;
                        n_prot++;
                    }

                }
                if (tStdHepPdg[j]==211 || tStdHepPdg[j]==-211){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    // fAbsoluteParticleMomentum_pion=fAbsoluteParticleMomentum;
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    if (fKE>=pion_ke){
                        visible_count++;
                        n_pi++;

                        if (tStdHepPdg[j]==211){
                            n_piplus++;
                        }
                        else if (tStdHepPdg[j]==-211){
                            n_piminus++;
                        }
                    }
                }

                if (tStdHepPdg[j]==111){
                    visible_count++, n_pi++,n_pizero++;
                }

                if (tStdHepPdg[j]==321|| tStdHepPdg[j]==-321 || tStdHepPdg[j]==130 || tStdHepPdg[j]==310){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    // fAbsoluteParticleMomentum_kaon=fAbsoluteParticleMomentum;
                    if (fKE>=kaon_ke){
                        visible_count++;

                    }
                }
            }

        }
        
        if (visible_count>max_num){
            max_num=visible_count;   
            //std::cout << "for i:" << i <<"  max num:"<<max_num<< endl;        
        }

    }

    //////////////START OF THE 2ND EVENT LOOP//////
    for(int i=0; i<Num_events; i++)
    {   //Initializing multiplicity variables
        SignalTree->GetEntry(i);
        
        if ( abs( tStdHepPdg[0] - tStdHepPdg[4] ) > 1 ) continue;

        int n_prot=0, n_piplus=0, n_piminus=0, n_pizero=0, n_pi=0, visible_count=0;

        energies.clear(); masses.clear(); pxs.clear(); pys.clear(); pzs.clear();costheta_arr.clear();theta_arr.clear();
        pdgs.clear();

        double tot_fpx=0,tot_fpy=0,tot_fpz=0;
        double tot_fKE=0;

        bool skip_event = false;
        
        //Print out the current event number to screen 
        // //std::cout << "The current event number is " << i << " of " << NumberEntries << endl;
        //Get the event entry information
        
        

        if(i%100000==0){std::cout << "The current event number is " << i << " of " << NumberEntries << endl;}  
        // //std::cout << "The current event number is " << i << " of " << NumberEntries << endl;            
        // //std::cout << "secondly, Number of elements in the array: " << tStdHepN/*leng*/ << std::endl; 
        int tStdHepPdg_nu, lepton_pdg_instore;
        double fAbsoluteParticleMomentum_nu,fKE_nu,fAbsoluteParticleMomentum_lept,fKE_lept;


        for (int j=0;j<tStdHepN;j++){

            if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                tStdHepPdg_nu=tStdHepPdg[j];
                auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                fAbsoluteParticleMomentum_nu=fAbsoluteParticleMomentum;
                fKE_nu=fKE;

            }

            if (tStdHepStatus[j] == 1){
                if (tStdHepPdg[j] == 11 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -15 || tStdHepPdg[j] == 16 || tStdHepPdg[j] == -16){                
                    lepton_pdg_instore=tStdHepPdg[j];
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    fAbsoluteParticleMomentum_lept=fAbsoluteParticleMomentum;
                    fKE_lept=fKE;
                    if ((tStdHepPdg[j]==13 && fKE_lept>=muon_ke) || (tStdHepPdg[j]==-13 && fKE_lept>=muon_ke) || (tStdHepPdg[j]==11 && fKE_lept>=electron_ke) || (tStdHepPdg[j]==-11 && fKE_lept>=electron_ke)){
                        visible_count++;
                    }

                }

                if (tStdHepPdg[j] == 2212){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    if (fKE>=prot_ke){
                        visible_count++;
                        n_prot++;
                    }

                }
                if (tStdHepPdg[j]==211 || tStdHepPdg[j]==-211){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    // fAbsoluteParticleMomentum_pion=fAbsoluteParticleMomentum;
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    if (fKE>=pion_ke){
                        visible_count++;
                        n_pi++;

                        if (tStdHepPdg[j]==211){
                            n_piplus++;
                        }
                        else if (tStdHepPdg[j]==-211){
                            n_piminus++;
                        }
                    }
                }

                if (tStdHepPdg[j]==111){
                    visible_count++, n_pi++,n_pizero++;
                }

                if (tStdHepPdg[j]==321|| tStdHepPdg[j]==-321 || tStdHepPdg[j]==130 || tStdHepPdg[j]==310){
                    auto [fAbsoluteParticleMomentum, fInvMass,fKE] = kinematics(tStdHepP4, j);
                    if (fAbsoluteParticleMomentum<=0.){
                        skip_event = true;
                        break;
                    }
                    // fAbsoluteParticleMomentum_kaon=fAbsoluteParticleMomentum;
                    if (fKE>=kaon_ke){
                        visible_count++;

                    }
                }
            }

        }
        

            
        
        if (visible_count==0 || skip_event){
            // //std::cout << "I am skipping for i: " << i << endl;
            continue;
        }



        bool found = false;
        bool string_found=false;
        for (const auto& value : Init_PDG_arr) {
            if (std::holds_alternative<int>(value) && std::get<int>(value) == tStdHepPdg_nu) {
                found = true;
                break;
            }

            else if (std::holds_alternative<std::string>(value)){
                string_found=true;
                break;
            }
        }

        // if (found) {
        //     //std::cout << "tStdHepPdg_nu is in Init_PDG_arr" <<"for i  "<<i<< std::endl;
        // } else {
        //     //std::cout << "tStdHepPdg_nu is not in Init_PDG_arr" << "for i  "<<i<<std::endl;
        // }

        if (found && (is_Final_lepton_PDG_int && int_Final_lepton_PDG == lepton_pdg_instore) && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {

            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);

                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);
                    
                    
                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";

                    // Output the results
                    // //std::cout << "i AM HEre for i:"<<i<<"and j:"<< j << std::endl;

                    // //std::cout << "Absolute Particle Momentum: " << fAbsoluteParticleMomentum << std::endl;
                    // // //std::cout << "Invariant Mass: " << fInvMass << std::endl;
                    // //std::cout << "Kinetic Energy: " << fKE << std::endl;
                    
                }

                if(tStdHepStatus[j] == 1){
                    //Charged current cases
                    if(tStdHepPdg[j] == 11 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == -15){
                        lepton_kinematics(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }

                    else if (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310){                       
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }
                }

            }   
            
            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //NumuInclusiveNpNpi/NueInclusiveNpNpi

        if (found && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "Inclusive") && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {
            
            // //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";
                    // Output the results
                    // //std::cout << "i AM HEre for i:"<<i<<"and j:"<< j << std::endl;

                    // //std::cout << "Absolute Particle Momentum: " << fAbsoluteParticleMomentum << std::endl;
                    // // //std::cout << "Invariant Mass: " << fInvMass << std::endl;
                    // //std::cout << "Kinetic Energy: " << fKE << std::endl;
                    
                }

                if(tStdHepStatus[j] == 1){
                    //Charged current cases
                    if(tStdHepPdg[j] == 11 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == -15){
                        lepton_kinematics(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }

                    else if (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310){                    
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //NumuCCNpNpi/NueCCNpNpi
        if (found && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "CC") && (lepton_pdg_instore==11 || lepton_pdg_instore==-11 || lepton_pdg_instore==13 || lepton_pdg_instore==-13 || lepton_pdg_instore==15 || lepton_pdg_instore==-15) && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {

            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                    
                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";
                    // Output the results
                    // //std::cout << "i AM HEre for i:"<<i<<"and j:"<< j << std::endl;

                    // //std::cout << "Absolute Particle Momentum: " << fAbsoluteParticleMomentum << std::endl;
                    // // //std::cout << "Invariant Mass: " << fInvMass << std::endl;
                    // //std::cout << "Kinetic Energy: " << fKE << std::endl;
                    
                }

                if(tStdHepStatus[j] == 1){
                    //Charged current cases
                    if(tStdHepPdg[j] == 11 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == -15){
                        lepton_kinematics(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }

                    else if (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310){              
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //NumuNCNpNpi/NueNCNpNpi
        if (found && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "NC") && (lepton_pdg_instore==12 || lepton_pdg_instore==-12 || lepton_pdg_instore==14 || lepton_pdg_instore==-14 || lepton_pdg_instore==16 || lepton_pdg_instore==-16) && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {

            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";                  
                    
                }

                if(tStdHepStatus[j] == 1 && (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310)){
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //AnyNuCCNpNpi
        if (string_found && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "CC") && (lepton_pdg_instore==11 || lepton_pdg_instore==-11 || lepton_pdg_instore==13 || lepton_pdg_instore==-13 || lepton_pdg_instore==15 || lepton_pdg_instore==-15) && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {

            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){

                
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);
                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";                  

                    // Output the results
                    // //std::cout << "i AM HEre for i:"<<i<<"and j:"<< j << std::endl;

                    // //std::cout << "Absolute Particle Momentum: " << fAbsoluteParticleMomentum << std::endl;
                    // // //std::cout << "Invariant Mass: " << fInvMass << std::endl;
                    // //std::cout << "Kinetic Energy: " << fKE << std::endl;
                    
                }

                if(tStdHepStatus[j] == 1){
                    //Charged current cases
                    if(tStdHepPdg[j] == 11 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == -15){
                        lepton_kinematics(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }

                    else if (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310){
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //AnyNuNCNpNpi

        if (string_found  && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "NC") && (lepton_pdg_instore==12 || lepton_pdg_instore==-12 || lepton_pdg_instore==14 || lepton_pdg_instore==-14 || lepton_pdg_instore==16 || lepton_pdg_instore==-16) && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {

            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                        
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);

                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);



                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";
                    
                }

                if(tStdHepStatus[j] == 1 && (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310)){
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        }

        //AnyNuInclusiveNpNpi(Basically All)
        if (string_found && (!is_Final_lepton_PDG_int && str_Final_lepton_PDG == "Inclusive") && Min_Nu_energy < fKE_nu && fKE_nu < Max_Nu_energy &&
        (
            (is_num_prot_int && is_num_pion_int && n_prot == int_num_prot && n_pi == int_num_pion) ||
            (!is_num_prot_int && str_num_prot == "N" && is_num_pion_int && n_pi == int_num_pion) ||
            (is_num_prot_int && n_prot == int_num_prot && !is_num_pion_int && str_num_pion == "N") ||
            (!is_num_prot_int && str_num_prot == "N" && !is_num_pion_int && str_num_pion == "N")
        )) {
            
            //std::cout << "I am here yay for i  "<<i<< std::endl;
            
            for(int j=0; j<tStdHepN; j++){
            // //std::cout <<" i issss "<<i<< " and j isssss " << j << " and status is " << tStdHepStatus[j] <<" and pdg " << tStdHepPdg[j]<< std::endl;

            
                if(tStdHepStatus[j] == 0 && (tStdHepPdg[j] == -16 || tStdHepPdg[j] == -14 || tStdHepPdg[j] == -12 || tStdHepPdg[j] == 12 || tStdHepPdg[j] == 14 || tStdHepPdg[j] == 16)){
                    auto [fAbsoluteParticleMomentum, fKE] = kinematics_massless(tStdHepP4, j);

                    double baseline= calc_baseline(tStdHepP4, fAbsoluteParticleMomentum, j);
                    double phi=phi_nu(tStdHepP4,fAbsoluteParticleMomentum, j);

                    Init_Nu_Mom->Fill(1000.*fAbsoluteParticleMomentum);
                    Init_Nu_CosTheta->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum);
                    Init_Nu_Theta->Fill(acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum));
                    Init_Nu_Phi->Fill(phi);
                    Oscillogram->Fill(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum,fAbsoluteParticleMomentum);

                    //outfile << "\"" << i << "\",";                  
                    outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" << phi << "\",\""<< baseline << "\",\"";

                    //outfile << "\"" << i << "\"," << std::setprecision(3)  << "\"" << tStdHepPdg[j] << "\",\"" << fAbsoluteParticleMomentum << "\",\"" << tStdHepP4[4*j] << "\",\"" << tStdHepP4[4*j+1] << "\",\"" << tStdHepP4[4*j+2] << "\",\"" << tStdHepP4[4*j+1]/fAbsoluteParticleMomentum << "\",\"" << (180./PI)*acos(tStdHepP4[4*j+1]/fAbsoluteParticleMomentum) << "\",\"" ;
                }

                if(tStdHepStatus[j] == 1){
                    //Charged current cases
                    if(tStdHepPdg[j] == 11 || tStdHepPdg[j] == 13 || tStdHepPdg[j] == 15 || tStdHepPdg[j] == -11 || tStdHepPdg[j] == -13 || tStdHepPdg[j] == -15){
                        lepton_kinematics(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }

                    else if (tStdHepPdg[j] == 2212 || tStdHepPdg[j] == 211 || tStdHepPdg[j] == -211 || tStdHepPdg[j] == 111 || tStdHepPdg[j] == 321 || tStdHepPdg[j] == -321 || tStdHepPdg[j] == 130 || tStdHepPdg[j] == 310){
                        // if (i == 74 && tStdHepPdg[j]==2212) {
                        //     //std::cout<<"pz is"<<tStdHepP4[4*j+3]<<std::endl;
                        //     auto [fAbsoluteParticleMomentum, fInvMass, fKE] = kinematics(tStdHepP4, j);
                        //     //std::cout << "mass is " << fInvMass << std::endl;
                        //     //std::cout << "KE is " << fKE << std::endl;
                        // }
                        //std::cout<<"I am here for i:"<<i<<endl;
                        finalparticles_info(tStdHepP4,j,tStdHepPdg,pdgs,masses,
                        energies, pxs, pys, pzs,costheta_arr,theta_arr, tot_fKE, tot_fpx, tot_fpy, tot_fpz,dictionary);
                    }
                }

            }   

            do_resize(pdgs,masses,energies,pxs,pys,pzs,costheta_arr,theta_arr);

            Fin_Prot_Mult->Fill(n_prot); n_prot=0;
            Fin_PiPlus_Mult->Fill(n_piplus); n_piplus=0;
            Fin_PiMinus_Mult->Fill(n_piminus); n_piminus=0;
            Fin_PiZero_Mult->Fill(n_pizero); n_pizero=0;

            for(int j=0; j<pdgs.size(); j++)
            {
            //Be sure to end the vector with a quote
            if ( j <  pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << ",";}
            if ( j == pdgs.size()-1 ) { outfile << setprecision(3) << pdgs[j] << "\"";}
            }

            outfile << ",\"";

            writeVectorToFile(outfile, masses);
            writeVectorToFile(outfile, energies);
            writeVectorToFile(outfile, pxs);
            writeVectorToFile(outfile, pys);
            writeVectorToFile(outfile, pzs);
            writeVectorToFile(outfile,costheta_arr);
            writeVectorToFile(outfile,theta_arr);

            double p_tot=sqrt(pow(tot_fpx,2)+pow(tot_fpy,2)+pow(tot_fpz,2));
            double p_miss=tot_fKE-p_tot;
            outfile << tot_fKE << "\",\"" << p_tot << "\",\""<< p_miss << "\"\n";

        } 

        // //std::cout << "Contents of Vectors:" << std::endl;
        // //std::cout << "PDGs: ";
        // for (const auto& pdg : pdgs) {
        //     //std::cout << pdg << " ";
        // }
        // //std::cout << std::endl;

        // //std::cout << "Masses: ";
        // for (const auto& mass : masses) {
        //     //std::cout << mass << " ";
        // }
        // //std::cout << std::endl;

        // //std::cout << "Energies: ";
        // for (const auto& energy : energies) {
        //     //std::cout << energy << " ";
        // }
        // //std::cout << std::endl;

        // //std::cout << "Px: ";
        // for (const auto& px : pxs) {
        //     //std::cout << px << " ";
        // }
        // //std::cout << std::endl;

        // //std::cout << "Py: ";
        // for (const auto& py : pys) {
        //     //std::cout << py << " ";
        // }
        // //std::cout << std::endl;

        // //std::cout << "Pz: ";
        // for (const auto& pz : pzs) {
        //     //std::cout << pz << " ";
        // }
        // //std::cout << std::endl;

    } // End of event loop
    
    std::cout <<"final max num: "<<max_num<< endl;

    // Create a map to store strings as keys and integers as values
    std::map<std::string, int> maxprong_dict;

    // Insert key-value pairs into the map
    maxprong_dict[outfile_name] = max_num;

    // Open a file in write mode
    std::ofstream file("maxprong_dict.txt", std::ios_base::app);


    // Write key-value pairs to the file
    for (const auto& pair :maxprong_dict) {
        file << pair.first << " " << pair.second << std::endl;
    }




    //Write out histograms to file
    Init_Nu_Mom->Write();
    Init_Nu_CosTheta->Write();
    Init_Nu_Theta->Write();
    Init_Nu_Phi->Write();
    Oscillogram->Write();
    Fin_CC_Lept_Mom->Write();
    Fin_CC_Lept_Theta->Write();
    Fin_CC_Lept_CosTheta->Write();
    Fin_Prot_Mom->Write();
    Fin_PiPlus_Mom->Write();
    Fin_PiMinus_Mom->Write();
    Fin_PiZero_Mom->Write();
    

    //std::cout << "Initial Neutrino Momentum, hA 2018" << endl;
    //Init_Nu_Mom->Print("all");
    Init_Nu_Mom->Draw("hist");
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->SetLogy(1);
    c1->Print((directory + "/" + last_name + "Init_Nu_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Init_Nu_Mom.root").c_str()); c1->Clear();

    //Now for the Nu angle histograms
    Init_Nu_CosTheta->Draw("hist");
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Init_Nu_CosTheta.png").c_str());
    c1->Print((directory + "/" + last_name + "Init_Nu_CosTheta.root").c_str());
    c1->Clear();


    //And now for the initial neutrino oscillogram
    Init_Nu_Theta->Draw("hist");
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Init_Nu_Theta.png").c_str());
    c1->Print((directory + "/" + last_name + "Init_Nu_Theta.root").c_str());
    c1->Clear();

    Init_Nu_Phi->Draw("hist");
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Init_Nu_Phi.png").c_str());
    c1->Print((directory + "/" + last_name + "Init_Nu_Phi.root").c_str());
    c1->Clear();

    Oscillogram->Draw("colz");
    c1->SetLogy(1);
    c1->SetLogz(1);
    c1->Print((directory + "/" + last_name + "Oscillogram.png").c_str());
    c1->Print((directory + "/" + last_name + "Oscillogram.root").c_str());
    c1->Clear();

    Fin_CC_Lept_Mom->Draw("hist");
    c1->SetLogy(1);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_Mom.root").c_str());
    c1->Clear();

    Fin_CC_Lept_Theta->Draw("hist");
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_Theta.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_Theta.root").c_str());
    c1->Clear();

    Fin_CC_Lept_CosTheta->Draw("hist");
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_CosTheta.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_CC_Lept_CosTheta.root").c_str());
    c1->Clear();

    Fin_Prot_Mom->Draw("hist");
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->SetLogy(1);
    c1->Print((directory + "/" + last_name + "Fin_Prot_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_Prot_Mom.root").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_Prot_Mom_low.png").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_Prot_Mom_low.root").c_str());
    c1->Clear();

    Fin_Prot_Mult->Draw("hist");
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->SetLogy(0);
    c1->Print((directory + "/" + last_name + "Fin_Prot_Mult.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_Prot_Mult.root").c_str());
    c1->Clear();

    Fin_PiPlus_Mom->Draw("hist");
    c1->SetLogy(1);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mom.root").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mom_low.png").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mom_low.root").c_str());
    c1->Clear();

    Fin_PiPlus_Mult->Draw("hist");
    c1->SetLogy(0);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mult.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiPlus_Mult.root").c_str());
    c1->Clear();

    Fin_PiMinus_Mom->Draw("hist");
    c1->SetLogy(1);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mom.root").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mom_low.png").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mom_low.root").c_str());
    c1->Clear();

    Fin_PiMinus_Mult->Draw("hist");
    c1->SetLogy(0);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mult.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiMinus_Mult.root").c_str());
    c1->Clear();
    
    //pi0 histograms
    Fin_PiZero_Mom->Draw("hist");
    c1->SetLogy(1);
    c1->BuildLegend(0.5,0.3,0.9,0.7);  
    c1->Print((directory + "/" + last_name + "Fin_PiZero_Mom.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiZero_Mom.root").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiZero_Mom_low.png").c_str());
    // c1->Print((directory + "/" + last_name + "Fin_PiZero_Mom_Low.root").c_str());
    c1->Clear();

    
    Fin_PiZero_Mult->Draw("hist");
    c1->SetLogy(0);
    c1->BuildLegend(0.5,0.3,0.9,0.7);
    c1->Print((directory + "/" + last_name + "Fin_PiZero_Mult.png").c_str());
    c1->Print((directory + "/" + last_name + "Fin_PiZero_Mult.root").c_str());
    c1->Clear();

    treefile->Write();
    c1->Close();
}





int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1; // Exit with error code 1 if incorrect usage
    }

    std::string input_file = argv[1];
    testcopy(input_file);

    return 0; // Exit successfully
}

void run_test(const char* input_file) {
    std::string file(input_file);
    testcopy(file);
}




    // double fAbsoluteParticleMomentum, fInvMass, fKE;
    // kinematics(tStdHepP4, j, fAbsoluteParticleMomentum, fInvMass, fKE);

  