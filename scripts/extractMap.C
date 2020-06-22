#include <iostream>
#include <vector>

#pragma link C++ class vector < vector < short>> + ;

#include "Math/PositionVector3D.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TTree.h"

#include "DepositedCharge.hpp"
#include "PixelCharge.hpp"
#include "PropagatedCharge.hpp"

using namespace std;

int test() {
  gSystem->Load("/home/mpetric/allpix/lib/libAllpixObjects.so");

  TFile* file0                  = TFile::Open("conf/output/data.root");
  TTree* propagated_charge_tree = static_cast<TTree*>(file0->Get("PropagatedCharge"));
  TTree* deposited_charge_tree  = static_cast<TTree*>(file0->Get("DepositedCharge"));
  TTree* pixel_charge_tree      = static_cast<TTree*>(file0->Get("PixelCharge"));

  TBranch* propagated_charge_branch = propagated_charge_tree->FindBranch("clicpix2");
  TBranch* deposited_charge_branch  = deposited_charge_tree->FindBranch("clicpix2");
  TBranch* pixel_charge_branch      = pixel_charge_tree->FindBranch("clicpix2");

  std::vector<allpix::PixelCharge*> pixelCharges;
  pixel_charge_branch->SetObject(&pixelCharges);

  std::vector<allpix::DepositedCharge*> depositedCharges;
  deposited_charge_branch->SetObject(&depositedCharges);

  std::vector<allpix::PropagatedCharge*> propagatedCharges;
  propagated_charge_branch->SetObject(&propagatedCharges);

  TH2F* pixel      = new TH2F("pixel", "pixel", 20, 1.55, 1.63, 20, 1.55, 1.63);
  TH2F* deposited  = new TH2F("deposited", "deposited", 200, 1.55, 1.63, 200, 1.55, 1.63);
  TH2F* propagated = new TH2F("propagated", "propagated", 200, 1.55, 1.63, 200, 1.55, 1.63);

  vector<vector<Short_t>> store(1000, vector<Short_t>(0));
  cout << store.size();
  for (int i = 0; i < pixel_charge_branch->GetEntries(); i++) {
    pixel_charge_branch->GetEntry(i);
    deposited_charge_branch->GetEntry(i);
    propagated_charge_branch->GetEntry(i);
    cout << depositedCharges[0]->getLocalPosition().X() << "\t" << depositedCharges[0]->getLocalPosition().Y() << "\t"
         << depositedCharges[0]->getLocalPosition().Z() << "\t" << endl;

    if (pixelCharges.size() > 0) {
      for (auto charge : pixelCharges) {
        Short_t collectedCharge = ((double)charge->getCharge()) / ((double)depositedCharges[0]->getCharge()) * 100;
        Short_t xpos            = int(charge->getIndex().X()) - 64;
        Short_t ypos            = int(charge->getIndex().Y()) - 64;
        store[i].emplace_back(collectedCharge);
        store[i].emplace_back(xpos);
        store[i].emplace_back(ypos);
      }
    } else {
      //if there is no propagation we leave the vector empty
    }
  }

  TFile* fout = TFile::Open("map.root", "RECREATE");
  fout->WriteObject(&store, "store");
  fout->Close();

  TFile*                   fin = TFile::Open("map.root", "READ");
  vector<vector<Short_t>>* reread;
  fin->GetObject("store", reread);  // I try to retrieve the vector
  for (auto entry : *reread) {
    for (auto vec : entry) {
      cout << vec << "\t";
    }
    cout << endl;
  }

  return 0;
}
