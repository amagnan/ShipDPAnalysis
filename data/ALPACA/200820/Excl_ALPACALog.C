void Excl_ALPACALog()
{
//=========Macro generated from canvas: c1/
//=========  (Tue Aug 25 01:02:15 2020) by ROOT version 6.20/04
   TCanvas *c1 = new TCanvas("c1", "",0,45,1920,1035);
   c1->SetHighLightColor(2);
   c1->Range(0,0,1,1);
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
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle(";m_{#gamma^{D}}(MeV/c^{2});#varepsilon^{2}");
   multigraph->Draw("AL");
   multigraph->GetXaxis()->SetLimits(1, 3500);
   multigraph->GetXaxis()->SetTitle("m_{#gamma^{D}}(MeV/c^{2})");
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetLabelOffset(0.007);
   multigraph->GetXaxis()->SetLabelSize(0.04);
   multigraph->GetXaxis()->SetTitleSize(0.04);
   multigraph->GetXaxis()->SetTickLength(0.02);
   multigraph->GetXaxis()->SetTitleOffset(0.65);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("#varepsilon^{2}");
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetLabelOffset(0.007);
   multigraph->GetYaxis()->SetLabelSize(0.04);
   multigraph->GetYaxis()->SetTitleSize(0.04);
   multigraph->GetYaxis()->SetTickLength(0.02);
   multigraph->GetYaxis()->SetTitleOffset(0.95);
   multigraph->GetYaxis()->SetTitleFont(42);
   multigraph->SetMinimum(1e-18);
   multigraph->SetMaximum(1e-05);
   
   TLegend *leg = new TLegend(0.75,0.7,0.95,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->Draw();
   
   leg = new TLegend(0.75,0.7,0.95,0.9,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
