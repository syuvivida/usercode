#include <TCut.h>
#include "tightSelection.h"

void overlayData(TChain* t1,
		 std::string var1,
		 int decCode, 
		 std::string psname,
		 int nbin, float xmin, float xmax, 
		 std::string xtitle="", std::string ytitle="",
		 TCut cut1="",bool logy=false)
{

  // for cuts the same for barrel and endcap
  TCut allCut = cut1 + myCut;

  // only for data to choose 900 GeV or 2360GeV runs
  TCut dataCut;
  if(psname.find("2360") != std::string::npos)
    dataCut= Cut_2360GeV;
  else
    dataCut= Cut_900GeV;

  allCut += dataCut;
  // if we are scaling relative to all data, uncommented out the following line
  //  TCut scaleCut = allCut; // doesnt separate endcap and barrel

  // now separate barrel and endcap
  if(decCode==1)
    allCut += barrelCut;
  else if(decCode==2)
    allCut += endcapCut;
  // if we are scaling barrel and endcap separately, uncommented the following line
  TCut scaleCut = allCut;  


  float binwidth = (xmax-xmin)/(float)nbin;

  char tempName[300];
  sprintf(tempName, "Candidates per %1.3lf", binwidth);
  if(binwidth < 0.001)
  sprintf(tempName, "Candidates per %1.4lf", binwidth);
  if(binwidth < 0.0001)
  sprintf(tempName, "Candidates per %1.5lf", binwidth);
  if(ytitle == "")ytitle = tempName;
  
  // first read in the runs that need to be checked
  int TOTALRUNS = 0;
  int runNo[100]={0};
  ifstream fin;
  fin.open("runs.txt");
  int temp; 
  fin >> temp;
  if(!fin.eof())runNo[TOTALRUNS]=temp;
  while(!fin.eof() && TOTALRUNS < 100)
    {
      TOTALRUNS ++;
      fin >> temp;
      if(!fin.eof())runNo[TOTALRUNS]=temp;
    }

  for(int k=0; k< TOTALRUNS; k++)
    cout << runNo[k] << endl;

  if(TOTALRUNS < 1)return;
 
  const int numRuns = TOTALRUNS;

  TH1F* h[numRuns];
  TH1F* hs[numRuns];
  char CutName[300];
  char HisName[300];
  char ScaleName[300];
  char title[300];

  int MIN = 1e9;  
  int MINRUN = -1;
  int firstNonZero = -1;

  // fill histograms for each run
  int NHISTOS_ENTRIES = 0;
  for(int i=0; i < numRuns; i++)
    {

      sprintf(title,"Run %d", runNo[i]);
      sprintf(CutName,"run == %d", runNo[i]);
      TCut thisRunCut(CutName);

      // fill histogram for showing
      sprintf(HisName,"h%02i",i);
      h[i] = new TH1F(HisName,title,nbin,xmin,xmax);
      h[i]->SetXTitle(xtitle.data());  
      h[i]->SetYTitle(ytitle.data());  
      std::string histo = HisName;
      std::string tempDraw = var1 + ">>" + histo;
      t1->Draw(tempDraw.data(), allCut && thisRunCut);

      // fill histograms for scaling
      sprintf(ScaleName,"hs%02i",i);
      sprintf(title,"Run %d", runNo[i]);
      hs[i] = new TH1F(ScaleName,title,nbin,xmin,xmax);
      histo = ScaleName;
      std::string tempDraw2 = var1 + ">>" + histo;
      t1->Draw(tempDraw2.data(), scaleCut && thisRunCut);

      int nScaleEve = hs[i]->GetEntries();
      int nEve = h[i]->GetEntries();
      
      // find out the first non-zero histogram
//       if(nEve >0 && firstNonZero < 0)firstNonZero = i;
      if(runNo[i] == 124023)firstNonZero = i;
      if(nEve >0)NHISTOS_ENTRIES++;

      if(nScaleEve < MIN && nScaleEve>0)
	{
	  MIN= nScaleEve;
	  MINRUN = runNo[i];
	}

    }
 
  cout << "Run " << MINRUN << " has minimum entries = " << MIN << endl;
  MIN = hs[firstNonZero]->GetEntries();
  cout << "Run 124023 has minimum entries = " << MIN << endl;

  // find out the maximum and scale every histogram to have the same entries;
  float MAX = -999.;
  for(int i=0; i < numRuns; i++)
    {

      int nEve = h[i]->GetEntries();
      int nScaleEve = hs[i]->GetEntries();
      if( nEve < 1 || nScaleEve < 1)continue;
      float scale = float(MIN)/float(nScaleEve);
      cout << "nScaleEve = " << nScaleEve << ", nEve = " << nEve 
	   << ", scale " << i << " = " << scale << endl;
      h[i]->Sumw2();
      h[i]->Scale(scale);

      float thisMax = h[i]->GetMaximum()  
// 	+h[i]->GetBinError(h[i]->GetMaximumBin())
	;
      if(thisMax > MAX)MAX = thisMax;

      cout << "h" << i << " integral = " <<  h[i]->Integral() << endl;
    }


  MAX *= 1.2;
  cout << "MAX = " << MAX << " and firstNonZero = " << firstNonZero << endl;

  for(int i= 0; i < numRuns; i++)
    h[i]->SetMaximum(MAX);
      

  // making plots now
  std::string decName;
  if(decCode==0)decName = "all";
  else if(decCode==1)decName = "barrel";
  else if(decCode==2)decName = "endcap";
  
  TLegend* leg4 = new TLegend(0.652529,0.364407,0.859425,0.898305);
  leg4->SetFillColor(0);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.04);
  leg4->SetBorderSize(0);
 
  std::string dirname = 
    "$CMSSW_BASE/src/CRAB/preprod_figures/";
  psname = psname + "_" + decName;
  
  TCanvas* c1 = new TCanvas("c1","",500,500);
  const int HistosPerCanvas = 3;
  h[firstNonZero]->SetTitle("");
  h[firstNonZero]->SetLineColor(1);
  h[firstNonZero]->SetLineStyle(1);
  h[firstNonZero]->SetLineWidth(2);
  h[firstNonZero]->SetMarkerSize(0);
  h[firstNonZero]->GetXaxis()->SetNdivisions(5);
  char name[300];
  cout << "Run " << runNo[firstNonZero] << " included" << endl;
 
  int nHistos = 0;
  int nCanvas = 0;
  int colorCODE[]={kRed, kBlue, kViolet};
  int markerCODE[]={24,21,3};

  for(int i=0; i < numRuns; i++)
    {
      if(i == firstNonZero)continue;
      int nEve = h[i]->GetEntries();
      sprintf(title,"Run %d", runNo[i]);
      if(nEve<1)continue;
      nHistos++;

      // setting color and line style
      int tempNumber =(nHistos-1)%HistosPerCanvas;
      int COLOR = colorCODE[tempNumber] ; //((int)(tempNumber/2)+1)*2;
      int MARKER = markerCODE[tempNumber];
      int STYLE = nHistos%2 == 1? 1: 3;
      cout << COLOR << endl;
      h[i]->SetLineColor(COLOR);
      h[i]->SetMarkerStyle(MARKER);
      h[i]->SetLineStyle(STYLE);
      h[i]->SetMarkerColor(COLOR);
      h[i]->SetMarkerSize(1);
      h[i]->SetLineWidth(2);

      if(nHistos%HistosPerCanvas==1)
	{
	  if(logy)
	    c1->SetLogy(1);
	  else
	    c1->SetLogy(0);
	  h[firstNonZero]->SetMarkerSize(0);
	  h[firstNonZero]->Draw("histe");
	  leg4->Clear();
	  leg4->SetHeader(decName.data());
	  sprintf(name,"Run %d", runNo[firstNonZero]);
	  leg4->AddEntry(h[firstNonZero],name);
	}
      h[i]->Draw("same");
      leg4->AddEntry(h[i],title);

      cout << "Run " << runNo[i] << " included" << endl;
      if(nHistos%HistosPerCanvas== 0 || nHistos== (NHISTOS_ENTRIES-1))
	{
	  nCanvas ++;
	  leg4->Draw("same");
	  sprintf(name,"%d",nCanvas);
	  std::string countName(name);
	  std::string epsfilename = dirname + "dataOnly_" + psname + countName + ".eps";
	  std::string giffilename = dirname + "dataOnly_" + psname + countName + ".gif";

	  c1->Print(epsfilename.data());
// 	  c1->Print(giffilename.data());
	}
      
    }



//   TFile* outFile = new TFile("test.root","recreate");

//   for(int i=0; i < numRuns; i++)
//     h[i]->Write();
//   outFile->Close();


  for(int i=0; i < numRuns; i++)
    { delete h[i]; delete hs[i];}
  delete c1;
  delete leg4;
}
		     
