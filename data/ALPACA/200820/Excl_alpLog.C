void Excl_alpLog()
{
//=========Macro generated from canvas: c1/
//=========  (Tue Aug 25 01:20:53 2020) by ROOT version 6.20/04
   TCanvas *c1 = new TCanvas("c1", "",0,45,1920,1035);
   c1->SetHighLightColor(2);
   c1->Range(0.6470792,-8.75,3.207187,-1.25);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetLogy();
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetRightMargin(0.05);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1[30] = {
   8,
   20,
   40,
   80,
   100,
   200,
   300,
   400,
   500,
   600,
   700,
   800,
   900,
   1000,
   1100,
   1100,
   1000,
   900,
   800,
   700,
   600,
   500,
   400,
   300,
   200,
   100,
   80,
   40,
   20,
   8};
   Double_t Graph0_fy1[30] = {
   0.006075835,
   0.002090682,
   0.0007604734,
   0.0001869812,
   0.0001158601,
   2.822494e-05,
   1.223061e-05,
   6.548231e-06,
   4.025576e-06,
   2.667995e-06,
   1.90723e-06,
   1.416224e-06,
   1.090009e-06,
   8.616443e-07,
   6.89325e-07,
   3.359023e-08,
   3.563149e-08,
   3.805844e-08,
   4.124393e-08,
   4.558228e-08,
   5.137833e-08,
   5.947932e-08,
   7.210079e-08,
   9.405199e-08,
   1.410376e-07,
   3.065173e-07,
   4.01176e-07,
   9.832434e-07,
   2.469039e-06,
   8.975589e-06};
   TGraph *graph = new TGraph(30,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle(";m_{#alpha}(MeV/c^{2});g_{#alpha#gamma}[GeV^{-1}]");
   graph->SetFillStyle(1000);
   graph->SetLineColor(4);
   graph->SetMarkerColor(4);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","",100,8,1200);
   Graph_Graph01->SetMinimum(1e-08);
   Graph_Graph01->SetMaximum(0.01);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetTitle("m_{#alpha}(MeV/c^{2})");
   Graph_Graph01->GetXaxis()->SetLabelFont(42);
   Graph_Graph01->GetXaxis()->SetLabelSize(0.04);
   Graph_Graph01->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph01->GetXaxis()->SetTickLength(0.02);
   Graph_Graph01->GetXaxis()->SetTitleOffset(1.15);
   Graph_Graph01->GetXaxis()->SetTitleFont(42);
   Graph_Graph01->GetYaxis()->SetTitle("g_{#alpha#gamma}[GeV^{-1}]");
   Graph_Graph01->GetYaxis()->SetLabelFont(42);
   Graph_Graph01->GetYaxis()->SetLabelSize(0.04);
   Graph_Graph01->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph01->GetYaxis()->SetTickLength(0.02);
   Graph_Graph01->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph01->GetYaxis()->SetTitleFont(42);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.04);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.04);
   Graph_Graph01->GetZaxis()->SetTickLength(0.02);
   Graph_Graph01->GetZaxis()->SetTitleOffset(1);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph01);
   
   graph->Draw("al");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
