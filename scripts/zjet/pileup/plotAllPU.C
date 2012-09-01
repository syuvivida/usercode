{
  
  h_input_nint_data_central->SetLineColor(1);
  h_true_nint_mc_before->SetLineColor(2);
  h_true_nint_mc_after_standalone_central->SetLineColor(4);
  h_true_nint_mc_after_standalone_central->SetLineStyle(4);

  h_input_nint_mc->SetLineColor(2);
  h_input_nint_mc->SetLineStyle(2);
 
  h_input_nint_data_central->SetXTitle("Number of true interactions");
  h_input_nint_data_central->DrawNormalized();
  h_true_nint_mc_before->DrawNormalized("hesame");
  h_input_nint_mc->DrawNormalized("hesame");
  h_true_nint_mc_after_standalone_central->DrawNormalized("hesame");
  gPad->SetLogy(1);
  float x1NDC = 0.647;
  float y1NDC = 0.572;
  float x2NDC = 0.849;
  float y2NDC = 0.909;

  TLegend* leg = new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  leg->SetHeader("2011 central: 68 mb");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->AddEntry(h_input_nint_data_central,"Data input");
  leg->AddEntry(h_input_nint_mc, "MC input");
  leg->AddEntry(h_true_nint_mc_before, "MC output");
  leg->AddEntry(h_true_nint_mc_after_standalone_central, "Reweighed MC");
  leg->Draw("same");
  
  c1->Print("AllPUDist_central.eps");
  c1->Print("AllPUDist_central.gif");
  c1->Print("AllPUDist_central.pdf");
}
