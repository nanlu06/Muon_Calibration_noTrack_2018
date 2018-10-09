#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>

#include "NtuplerReader.h"

int main(int argc, const char* argv[]) {

       TChain *mc = new TChain("hcalHBHEMuon/TREE");

       TString listname = argv[1]; //"preselection_Data_Run2018A_PromptReco_v1_json1.list";
       TString prefix = argv[2];
       TString nvtx_sel_str = argv[3];
       TString path = argv[4];
       int nvtx_sel = nvtx_sel_str.Atoi();
       ifstream listfile;
       listfile.open(listname);
       if (! listfile.good()) {
          std::cerr << "Error: " << listname << " cannot be opened." << std::endl;
       return -1;
       }
       else std::cerr << "Opened: " << listname << std::endl;

       std::string filename;
       while (getline(listfile, filename)) {
       mc->Add(filename.c_str());
       std::cerr << "Added file: " << filename << std::endl;
       }
       NtuplerReader *readMC = new NtuplerReader(mc);
       readMC->Loop(path, prefix, nvtx_sel);
       listfile.close();

        //char file[100];
        //        //sprintf(file, "%s", argv[1]);
        //                //cout << "read file: " << file << " ... " << endl;
        //
       return 1;
}
