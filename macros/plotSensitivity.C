#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"

struct Rate{
  double mass;
  double eps;
  double rate1;
  double rate2;
  double vtxEff;
};

static bool customSort(const Rate & a,
		       const Rate & b)
{
  if (fabs(a.mass-b.mass)<0.01) return a.eps <= b.eps;
  else return a.mass <= b.mass;
};

static bool sortPair(const std::pair<double,double> & a,
		     const std::pair<double,double> & b)
  
{
  return a.first <= b.first;
};

void readInput(const std::string & prodfile,
	       std::vector<Rate> & ratevec){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Rate tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    tmp.rate1 = 0;
    tmp.rate2 = 0;
    tmp.vtxEff = 1;
    lInput>>tmp.mass>>tmp.eps>>tmp.rate1>>tmp.rate2>>tmp.vtxEff;

    //protect against nan
    if (tmp.rate1!=tmp.rate1) tmp.rate1=0;
    if (tmp.rate2!=tmp.rate2) tmp.rate2=0;

    //protect against 0 values
    if (tmp.mass>0 && tmp.rate2>0){
      ratevec.push_back(tmp);
    } else {
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps << " " << tmp.rate1 << " " << tmp.rate2 << std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};

int plotSensitivity(){
  
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);

  const double clsval = 2.3;
  
  const std::string prod="/afs/cern.ch/user/t/takmete/ShipDPAnalysis-master/data/180608/";
  
  TLatex lat;
  char lbuf[500];
  
  const unsigned nP = 4;
  std::string proc[nP] = {"meson","pbrem","qcd","dp_qcd"};

  std::ofstream foutCLS[nP];
  
  std::vector<Rate>lRate[nP];
  TH2F *hSens[nP];

  TH2F *hvtxEff[nP];
  
  TCanvas *mycC[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "myc" << iP ;
    mycC[iP] = new TCanvas(lname.str().c_str(),"",1);
  }

  for (unsigned iP(0); iP<nP; ++iP){

    std::cout << " - Processing " << proc[iP] << std::endl;
    
    std::ostringstream lname;
    lname << prod << proc[iP] << ".dat";
    readInput(lname.str(),lRate[iP]);

    std::sort(lRate[iP].begin(), lRate[iP].end(), customSort);

    lname.str("");
    lname << "ForNico" << proc[iP] << ".txt";
    foutCLS[iP].open(lname.str().c_str());

    
    //    for (unsigned iN(0); iN<lRate[iP].size(); ++iN){
    //std::cout << lRate[iP][iN].mass << " "
    //		<< lRate[iP][iN].eps << " "
    //		<< lRate[iP][iN].rate2 << std::endl;
    //}

    lname.str("");
    lname << "hSensitivity_" << proc[iP];
    hSens[iP] = new TH2F(lname.str().c_str(),";mass (GeV);#gamma' coupling to SM log_{10}(#varepsilon);rate for 2.10^{10} p.o.t.",120,0,9.2,100,-9,-5);
    if (!hSens[iP]) return 1;

    lname.str("");
    lname << "hvtxEff_" << proc[iP];
    hvtxEff[iP] = new TH2F(lname.str().c_str(),";mass (GeV);#gamma' coupling to SM log_{10}(#varepsilon);vessel efficiency (%)",120,0,9.2,100,-9,-5);

    double mass = 0;
    unsigned nM = 0;
    unsigned nE = 0;
    TGraph *gr = 0;
    TGraph *grTest = 0;
    mycC[iP]->cd();
    gPad->SetRightMargin(0.05);
    TLine *line = 0;
    TLine *line1 = 0;
    TLine *line2 = 0;
    
    
    
    mycC[iP]->Print(("CheckRate_"+proc[iP]+".pdf[").c_str());
    
    for (unsigned iN(0); iN<lRate[iP].size(); ++iN){
      
      hSens[iP]->Fill(lRate[iP][iN].mass,log10(lRate[iP][iN].eps),lRate[iP][iN].rate2);
      hvtxEff[iP]->Fill(lRate[iP][iN].mass,log10(lRate[iP][iN].eps),lRate[iP][iN].vtxEff*100);
      //std::cout << lRate[iP][iN].mass << " "
      //	<< lRate[iP][iN].eps << " "
      //	<< lRate[iP][iN].rate2/(lRate[iP][iN].eps*lRate[iP][iN].eps) << " "
      //	<< (lRate[iP][iN].rate1-lRate[iP][iN].rate2)/lRate[iP][iN].rate1
      //	<< std::endl;

      double tmpRate = lRate[iP][iN].rate2;
      //if (tmpRate>50) tmpRate=50;
      //if (tmpRate<0.05) tmpRate=0.05;
      if (fabs(lRate[iP][iN].mass-mass)>0.01 || iN==lRate[iP].size()-1) {

	if (iN==lRate[iP].size()-1){
	  gr->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate));
	  nE += 1;
	}
	if (gr) {
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  gPad->SetGridx(1);
	  gr->SetMarkerStyle(20);
	  gr->SetTitle(";#gamma' coupling to SM log(#varepsilon);log(rate) for 2.10^{10} p.o.t.");
	  gr->Draw("AP");
	  gr->GetYaxis()->SetRangeUser(-4,4);
	  sprintf(lbuf,"mDP = %3.3f GeV",mass);
	  lat.DrawLatexNDC(0.2,0.8,lbuf);
	  line = new TLine(-9,log10(clsval),-5,log10(clsval));
	  line->SetLineColor(6);
	  line->Draw("same");

	  
	  std::cout << " doing interpolation" << std::endl;
	  //interpolate and get 2.3 event points.
	  std::vector<std::pair<double,double> > limitVec;
	  lname.str("");
	  lname << "hTest_" << proc[iP] << "_M" << mass;
	  grTest=new TGraph();
	  grTest->SetName(lname.str().c_str());
	  for (unsigned ieps(0); ieps<1000000;++ieps){
	    double tmpeps=-8.+ieps*3./1000000;
	    double limit = gr->Eval(tmpeps);
	    if (ieps%1000==0) grTest->SetPoint(ieps/1000,tmpeps,limit);
	    
	    double diff = fabs(limit-log10(clsval));
	    if (diff<0.01) limitVec.push_back(std::pair<double,double>(diff,tmpeps));
	    //if (ieps>700&&ieps<800&&iP==2&&mass==1.0) std::cout << "check " << mass << " " << tmpeps << " " << limit << std::endl;
	  }
	  grTest->Draw("Psame");
	  if (limitVec.size()>0){
	    std::sort(limitVec.begin(),limitVec.end(),sortPair);
	    //for (unsigned iele(0); iele<limitVec.size();++iele){
	      //if (iP==1) std::cout << mass << " " << iele << " " << limitVec[iele].second << " " << limitVec[iele].first << std::endl;
	    //}
	    std::cout << mass << " " << limitVec[0].second << " " << limitVec[0].first << std::endl;
	    unsigned inext = 1;
	    while (fabs(limitVec[inext].second-limitVec[0].second)<0.01) {
	      inext++;
	      if (inext>limitVec.size()-1) break;
	    }
	    std::cout << mass << " " << inext << " " << pow(10,limitVec[inext].second) << " " << pow(10,limitVec[inext].first) << std::endl;
	    foutCLS[iP] << mass << ", " << pow(10,limitVec[0].second) << std::endl;
	    foutCLS[iP] << mass << ", " << pow(10,limitVec[inext].second) << std::endl;
	    line1 = new TLine(limitVec[0].second,-4,limitVec[0].second,4);
	    line1->SetLineColor(6);
	    line1->Draw("same");
	    line2 = new TLine(limitVec[inext].second,-4,limitVec[inext].second,4);
	    line2->SetLineColor(6);
	    line2->Draw("same");
	  } else {
	    std::cout << " -- did not find any point next to " << clsval << std::endl;
	  }
	  
	  mycC[iP]->Update();
	  mycC[iP]->Print(("CheckRate_"+proc[iP]+".pdf").c_str());
	}
	nE=0;
	mass = lRate[iP][iN].mass;
	lname.str("");
	lname << "hCheckRate_" << proc[iP] << "_M" << mass;
	std::cout << lname.str() << std::endl;
	gr = new TGraph();
	//gr->SetBit(TGraph::kIsSortedX);
	gr->SetName(lname.str().c_str());
	gr->SetPoint(0,log10(lRate[iP][iN].eps),log10(tmpRate));
	
	nE += 1;
	nM += 1;
      }//if different mass
      else {
	gr->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate));
	nE += 1;
      }
      
    }//loop on rates
    std::cout << " --- Found " << nM << " mass points." << std::endl;
    mycC[iP]->Print(("CheckRate_"+proc[iP]+".pdf]").c_str());

    foutCLS[iP].close();
    
  }//loop on proc mode
  
  
  TCanvas *myc[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "myc" << iP ;
    myc[iP] = new TCanvas(lname.str().c_str(),"",1);
  }
  TCanvas *mycV[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "mycV" << iP ;
    mycV[iP] = new TCanvas(lname.str().c_str(),"",1);
  }
  
  for (unsigned iP(0); iP<nP; ++iP){
    myc[iP]->cd();
    gPad->SetLogx(0);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    if (iP==0) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==1) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==2) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    hSens[iP]->GetZaxis()->SetRangeUser(0.01,100);
    hSens[iP]->Draw("colz");
    lat.DrawLatexNDC(0.1,1.0,proc[iP].c_str());
    myc[iP]->Update();
    myc[iP]->Print(("Sensitivity_"+proc[iP]+".pdf").c_str());

    mycV[iP]->cd();
    if (iP==0) hvtxEff[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==1) hvtxEff[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==2) hvtxEff[iP]->GetXaxis()->SetRangeUser(0,9.1);
    hvtxEff[iP]->GetZaxis()->SetRangeUser(1,100);
    hvtxEff[iP]->Draw("colz");
    lat.DrawLatexNDC(0.1,1.0,proc[iP].c_str());
    mycV[iP]->Update();
    mycV[iP]->Print(("VtxEfficiency_"+proc[iP]+".pdf").c_str());

    gPad->SetRightMargin(0.05);
    hvtxEff[iP]->Draw("text");
    lat.DrawLatexNDC(0.2,0.8,proc[iP].c_str());
    mycV[iP]->Update();
    mycV[iP]->Print(("VtxEfficiency_"+proc[iP]+"_text.pdf").c_str());

    
  }
  
    
  return 0;
}
