void plotFEWZ(std::string inputFile, std::string xtitle="0.5(y_{Z}-y_{jet1})",
	      std::string outputPrefix="test")
{
  

  float ey[50];
  float y[50];
  float cvalue[50];
  float evalue[50];
  float pdfperror[50];
  float pdfmerror[50];

  float maxheight=-9999;

  ifstream fin;
  fin.open(inputFile.data());
  int nBIN=1;
  fin >> y[0] >> cvalue[0] >> evalue[0] >> pdfperror[0] >> pdfmerror[0];
  while(!fin.eof())
    {
      fin >> y[nBIN] >> cvalue[nBIN] >> evalue[nBIN] >> 
	pdfperror[nBIN] >> pdfmerror[nBIN];
      ey[nBIN-1] = (y[nBIN]-y[nBIN-1])*0.5;

      float thismax = cvalue[nBIN]+evalue[nBIN];
      if(thismax > maxheight)
	maxheight = thismax;

      nBIN++;

    }

  cout << "There are " << nBIN << " bins" << endl;
  cout << "maxheight = " << maxheight << endl;
  ey[nBIN-1] = ey[nBIN-2];

  for(int i=0;i<nBIN;i++)
    {
      cout << ey[i] << "\t" << y[i] << "\t" << cvalue[i] << "\t"
	   << evalue[i] << "\t" << pdfperror[i] << "\t" 
	   << pdfmerror[i] << endl;      
    }

  TGraphErrors *gr_central= new TGraphErrors(nBIN, y, cvalue, ey, evalue); 
  gr_central->GetXaxis()->SetTitle(xtitle.data());
  gr_central->GetYaxis()->SetTitle("Weight");
  gr_central->GetYaxis()->SetTitleOffset(1.5);
  gr_central->SetTitle("");
  gr_central->SetMaximum(1.6*maxheight);
  gr_central->SetLineWidth(2);

  TGraphAsymmErrors *gr_pdf = new TGraphAsymmErrors(nBIN, y, cvalue, ey, ey, 
			      pdfmerror, pdfperror);
  
  gr_pdf->GetXaxis()->SetTitle(xtitle.data());
  gr_pdf->GetYaxis()->SetTitle("Weight");
  gr_pdf->GetYaxis()->SetTitleOffset(1.2);
  gr_pdf->SetTitle("");

  gr_pdf->SetFillColor(5);
  gr_pdf->SetLineColor(5);
  
  TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);
  
  gr_central->Draw("ap");
  gr_pdf->Draw("pe2,same");
  gr_central->Draw("p,same");
  
  TLegend *legr = new TLegend(0.481322,0.775298,0.787356,0.906250);
  legr->SetBorderSize(0);
  legr->SetFillColor(0);
  legr->SetFillStyle(0);
  legr->SetTextSize(0.03);
  legr->AddEntry(gr_central, "FEWZ+MSTW2008NNLO", "pl");
  legr->AddEntry(gr_pdf, "PDF uncertainty", "f");
  
  legr->Draw("same");

  c3->Print(Form("%s.eps",outputPrefix.data()));
  c3->Print(Form("%s.gif",outputPrefix.data()));
  c3->Print(Form("%s.pdf",outputPrefix.data()));


}
