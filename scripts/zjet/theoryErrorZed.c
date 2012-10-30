TGraph * getParams(TString PLOT ,TString TYPE, float offset){
  //
  // takes the scale and pdf errors parameterized by 
  // two lines and created a tgraph
  //
  // DIST=0(Ysum),1(Ydif),2(Yzed),3(Yjet)
  // VAR=0(SCALE) ,1(PDF), 3(correlated sytematic)
  //
  float p0=999;
  float p1=999;
  float q0=999;
  float q1=999;
  float xmax=2.0;
  cout << "PLOTTING       " << PLOT << endl;

  if( TYPE == "scale" ){
    if(PLOT=="Yzed"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      14.1376
      //Ndf                       =            8
      p0                        =      1.02187   ;//+/-   0.00252102  
      p1                        =   0.00222095   ;//+/-   0.0028169   
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      14.4821
      //Ndf                       =            8
      q0                        =     0.952714   ;//+/-   0.00264769  
      q1                        =  -0.00625964   ;//+/-   0.00281455  

      xmax = 2.2;
            
    }
    else if(PLOT=="Yjet"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      18.8519
      //Ndf                       =            8
      p0                        =      1.02254   ;//+/-   0.00411505  
      p1                        =  -0.00100786   ;//+/-   0.00433872  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      5.44619
      //Ndf                       =            8
      q0                        =     0.949235   ;//+/-   0.00351431  
      q1                        = -0.000704664   ;//+/-   0.0040205   

      xmax = 2.4;
            
    }
    else if(PLOT=="Ysum"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      8.17304
      //Ndf                       =            8
      p0                        =      1.03688   ;//+/-   0.0027989   
      p1                        =    -0.018442   ;//+/-   0.00356851  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      5.84382
      //Ndf                       =            8
      q0                        =     0.948807   ;//+/-   0.00280022  
      q1                        =  0.000426846   ;//+/-   0.00358723  

      xmax = 2.2;
            
    }
    else if(PLOT=="Ydif"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      31.9518
      //Ndf                       =            8
      p0                        =      0.99725   ;//+/-   0.00273221  
      p1                        =    0.0550617   ;//+/-   0.00456508  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      9.15572
      //Ndf                       =            8
      q0                        =     0.953604   ;//+/-   0.00232232  
      q1                        =   -0.0124717   ;//+/-   0.00405408  
      xmax = 2.0;
    }        
  }else{//PDF
    if(PLOT=="Yjet"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =       12.616
      //Ndf                       =            8
      p0                        =      1.04991   ;//+/-   0.00398275  
      p1                        =  -0.00283503   ;//+/-   0.00388914  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      9.59712
      //Ndf                       =            8
      q0                        =     0.981345   ;//+/-   0.00412055  
      q1                        =   0.00712834   ;//+/-   0.00397635  

      xmax = 2.4;
            
    }
    else if(PLOT=="Yzed"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      13.9214
      //Ndf                       =            8
      p0                        =      1.04328   ;//+/-   0.00234591  
      p1                        =   0.00553852   ;//+/-   0.00261301  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =       9.5016
      //Ndf                       =            8
      q0                        =     0.987812   ;//+/-   0.00411731  
      q1                        =   0.00230216   ;//+/-   0.00360298  
      xmax = 2.2;
            
    }
    else if(PLOT=="Ysum"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =        11.37
      //Ndf                       =            8
      p0                        =      1.05204   ;//+/-   0.00266447  
      p1                        =  -0.00490591   ;//+/-   0.00363166  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =       4.6955
      //Ndf                       =            8
      q0                        =     0.990225   ;//+/-   0.00387649  
      q1                        =  -0.00216067   ;//+/-   0.0045875   
      xmax = 2.2;
            
    }
    else if(PLOT=="Ydif"){
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      43.5413
      //Ndf                       =            8
      p0                        =      1.04224   ;//+/-   0.002366    
      p1                        =    0.0202095   ;//+/-   0.00402052  
            
      //****************************************
      //Minimizer is Linear
      //Chi2                      =      11.8057
      //Ndf                       =            8
      q0                        =     0.990368   ;//+/-   0.0030726   
      q1                        =   0.00135132   ;//+/-   0.00475005  
      xmax = 2.0;
    }
  }
    
  p0=1.0;
  q0=1.0;
    
  cout <<" p0="<<p0<<" p1="<<p1<<" q0="<<q0<<" q1="<<q1<<endl;
    
  TGraph *TGa  =makeGraph( p0, p1, q0, q1, offset,xmax);
  //TGra->Print("all");
  return TGa;
}


void colorIt(TH1F *hMidA, int kCyan , int kStyle){
  hMidA->SetMarkerColor(kCyan);
  hMidA->SetLineColor(kCyan);
  hMidA->SetMarkerStyle(0);//24
  hMidA->SetFillColor(kCyan);
  hMidA->SetFillStyle(kStyle);
}

void colorJt(TGraph *hMidA, int kCyan , int kStyle){
  hMidA->SetMarkerColor(kCyan);
  hMidA->SetLineColor(kCyan);
  hMidA->SetMarkerStyle(24);//24
  hMidA->SetFillColor(kCyan);
  hMidA->SetFillStyle(kStyle);
}


TGraph * makeGraph(float a1, float b1, float a2, float b2, float offset, 
		   float xmax){
  //
  // takes two line params as input outputs a graph with
  // asym errors
  //
  float x[300]; float y[300]; float dy[300];float dx[300];
  float y1,y2;
  int npoints=300;
    
  float dx1=xmax/npoints;
  for ( int i=0; i<npoints;i++){
    x[i]=(i+0.5)*dx1;
    y1=a1+b1*x[i];
    y2=a2+b2*x[i];
    y[i]=((y1+y2)/2.)+ offset;
    dy[i]=(y1-y2)/2;
    dx[i]=0.0;
  }
  TGraph *gr1 = new TGraphErrors(npoints,x,y,dx,dy);
  //gr1->Print("all");
  return gr1;
    
}


void theoryErrorZed(TString PLOT="Ydif", float OFFSET=0 ) 
{
  /*
   *  To be used to draw systematics on an existing plot
   *
   *  Created by linn on 6/2/12.
   *  Copyright 2012 __FIU/CMS__. All rights reserved.
   *
   */
    
  //this is for testing
    
  gROOT->Reset();     
  gStyle->SetOptStat(0);
    
  TH1F *scale1 = new TH1F("","", 150, 0.0,3.0);
  scale1->FillRandom("gaus", 1);
/*   TCanvas *canvas = new TCanvas("canvas","",600,600); */
/*   canvas->cd();  */
  scale1->SetMaximum(1.2);
  scale1->SetMinimum(0.8);
  colorIt(scale1, kBlack, 3000);
/*   scale1->Draw("pe1");// works */
    
  // end of test
    
  // put errors for PLOT="YSum"
  // OFFSET default is Y=0 on a linear plot
    
  TGraph * TGa = getParams( PLOT ,"scale", OFFSET);
/*   colorJt(TGa,kGreen-2, 3006); */
  colorJt(TGa,kGreen-2, 3004);
  TGa->Draw("4 same");// works
    
  TGraph * TGb = getParams( PLOT ,"pdfs", OFFSET);
  colorJt(TGb,kPink-8, 3005);
  TGb->Draw("4 same");
    
        
  TLegend *legend1 = new TLegend(0.713,0.813,0.871,0.954);
  legend1->SetFillColor(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.045);
  legend1->SetBorderSize(0);
  legend1->AddEntry(TGa,"Scale","f");
  legend1->AddEntry(TGb,"PDF","f");
  //legend1->AddEntry(TGc,"Correlated","f");
  legend1->Draw("same");
    
  //TGa->Print("all");
    
    
}





