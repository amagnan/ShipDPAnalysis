void Excl_alpLog()
{
//=========Macro generated from canvas: c1/
//=========  (Fri Jul 31 12:42:06 2020) by ROOT version 6.20/04
   TCanvas *c1 = new TCanvas("c1", "",0,45,1920,1035);
   c1->SetHighLightColor(2);
   c1->Range(1.766901,-8.625,3.128792,-2.375);
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
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle(";m_{a}(MeV/c^{2});g_{a#gamma}[GeV^{-1}]");
   
   Double_t Graph_fx1[22] = {
   2,
   4,
   8,
   10,
   20,
   80,
   300,
   500,
   600,
   700,
   800,
   900,
   1000,
   1100,
   900,
   600,
   500,
   400,
   300,
   200,
   100,
   80};
   Double_t Graph_fy1[22] = {
   0.0001970552,
   2.998859e-05,
   9.744331e-06,
   7.026839e-06,
   2.69081e-06,
   0.0004466281,
   2.865095e-05,
   3.916706e-06,
   2.497014e-06,
   1.845517e-06,
   1.371172e-06,
   9.740853e-07,
   8.336543e-07,
   6.647388e-07,
   3.330481e-08,
   2.305876e-08,
   4.533201e-08,
   7.024542e-08,
   1.025543e-07,
   1.536169e-07,
   3.335285e-07,
   4.360166e-07};
   TGraph *graph = new TGraph(22,Graph_fx1,Graph_fy1);
   graph->SetName("Graph");
   graph->SetTitle("ALP");
   graph->SetFillStyle(1000);
   graph->SetLineColor(4);
   graph->SetLineWidth(3);
   graph->SetMarkerColor(4);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","ALP",100,1.8,1209.8);
   Graph_Graph1->SetMinimum(2.075289e-08);
   Graph_Graph1->SetMaximum(0.0004912886);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.04);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph1->GetXaxis()->SetTickLength(0.02);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.65);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.04);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph1->GetYaxis()->SetTickLength(0.02);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.04);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.04);
   Graph_Graph1->GetZaxis()->SetTickLength(0.02);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   multigraph->Add(graph,"");
   multigraph->Draw("AL");
   multigraph->GetXaxis()->SetLimits(80, 1150);
   multigraph->GetXaxis()->SetTitle("m_{a}(MeV/c^{2})");
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelOffset(0.007);
   multigraph->GetXaxis()->SetLabelSize(0.04);
   multigraph->GetXaxis()->SetTitleSize(0.04);
   multigraph->GetXaxis()->SetTickLength(0.02);
   multigraph->GetXaxis()->SetTitleOffset(0.65);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("g_{a#gamma}[GeV^{-1}]");
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelOffset(0.007);
   multigraph->GetYaxis()->SetLabelSize(0.04);
   multigraph->GetYaxis()->SetTitleSize(0.04);
   multigraph->GetYaxis()->SetTickLength(0.02);
   multigraph->GetYaxis()->SetTitleOffset(0.95);
   multigraph->GetYaxis()->SetTitleFont(42);
   multigraph->SetMinimum(1e-08);
   multigraph->SetMaximum(0.001);
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
