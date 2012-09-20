TGraph * getParams(TString PLOT ,TString TYPE, float offset){
	//
	// takes the scale and pdf errors parameterized by 
	// two lines and created a tgraph
	//
	// DIST=0(Ysum),1(Ydif),2(Yzed),3(Yjet)
	// VAR=0(SCALE) ,1(PDF), 3(correlated sytematic)
	//
	float p0=999;
	float p1=000;
	float q0=999;
	float q1=999;
	cout << "PLOTTING       " << PLOT << endl;
    if( TYPE == "scale" ){
        
        if( PLOT == "Yzed" ){ 
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      14.1378
            //Ndf                       =            8
            p0                        =     0.998474   ;//+/-   0.00246329  
            p1                        =      0.00217   ;//+/-   0.0027524   
            
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      14.4823
            //Ndf                       =            8
            q0                        =      1.00467   ;//+/-   0.00279207  
            q1                        =    -0.006601   ;//+/-   0.00296803  
            
        }else if (PLOT == "Yjet"){
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      18.8518
            //Ndf                       =            8
            p0                        =     0.999123   ;//+/-   0.00402082  
            p1                        = -0.000984803   ;//+/-   0.00423938  
            
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      5.44617
            //Ndf                       =            8
            q0                        =        1.001   ;//+/-   0.00370595  
            q1                        = -0.000743084   ;//+/-   0.00423974  
            
        }else if(PLOT== "Ysum"){ 
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      8.17309
            //Ndf                       =            8
            p0                        =      1.01314   ;//+/-   0.00273481  
            p1                        =   -0.0180198   ;//+/-   0.0034868   
            
            //****************************************
            //MInimizer is Linear
            //Chi2                      =      5.84382
            //Ndf                       =            8
            q0                        =      1.00055   ;//+/-   0.00295293  
            q1                        =  0.000450039   ;//+/-   0.00378285  
            
        }else if( PLOT=="Ydif"){ 
            //****************************************
            //MInimizer is Linear
            //Chi2                      =       31.952
            //Ndf                       =            8
            p0                        =     0.974416   ;//+/-   0.00266965  
            p1                        =    0.0538008   ;//+/-   0.00446055  
            
            //****************************************
            //MInimizer is Linear
            //Chi2                      =       9.1557
            //Ndf                       =            8
            q0                        =      1.00561   ;//+/-   0.00244896  
            q1                        =   -0.0131519   ;//+/-   0.00427516  
        }
    }else{//PDF
        if (PLOT="Yjet"){ 
            //****************************************
            //Minimizer is Linear
            //Chi2                      =      177.093
            //NDf                       =            8
            p0                        =     1.00208   ;//+/-   0.00480356  
            p1                        =   0.00340742   ;//+/-   0.00508689  
            
            //****************************************
            //Minimizer is Linear
            //Chi2                      =      67.8396
            //NDf                       =            8
            q0                        =     0.991653   ;//+/-   0.0037726   
            q1                        =   -0.0260927   ;//+/-   0.0041261   
            
        }else if(PLOT=="YZed"){ 
            //****************************************
            //Minimizer is Linear
            //Chi2                      =      126.619
            //NDf                       =            8
            p0                        =      1.00457   ;//+/-   0.00318095  
            p1                        =  -0.00197225   ;//+/-   0.00337169  
            
            //****************************************
            //Minimizer is Linear
            //Chi2                      =      71.1142
            //NDf                       =            8
            q0                        =     0.990407   ;//+/-   0.0032956   
            q1                        =  -0.00354171   ;//+/-   0.00373747  
            
        }else if(PLOT="Ysum"){ 
            //****************************************
            //Minimizer is Linear
            //Chi2                      =       8.4352
            //NDf                       =            8
            p0                        =      1.01113   ;//+/-   0.00339403  
            p1                        =  -0.00490665   ;//+/-   0.00426406  
            
            //****************************************
            //Minimizer is Linear
            //Chi2                      =       3.6631
            //NDf                       =            8
            q0                        =      1.00028   ;//+/-   0.00309168  
            q1                        =   -0.0082579   ;//+/-   0.00385568  
            
        } else if(PLOT=="Ydif"){ 
            //****************************************
            //Minimizer is Linear
            //Chi2                      =        28.29
            //NDf                       =            8
            p0                        =     0.996068   ;//+/-   0.00296783  
            p1                        =    0.0178861   ;//+/-   0.00434411  
            
            //****************************************
            //Minimizer is Linear
            //Chi2                      =      21.7771
            //NDf                       =            8
            q0                        =     0.995272   ;//+/-   0.00233179  
            q1                        =   0.00307972   ;//+/-   0.00416143  
        }
    }
    
    cout <<" p0="<<p0<<" p1="<<p1<<" q0="<<q0<<" q1="<<q1<<endl;
    
    TGraph *TGa  =makeGraph( p0, p1, q0, q1, offset);
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


TGraph * makeGraph(float a1, float b1, float a2, float b2, float offset){
    //
    // takes two line params as input outputs a graph with
    // asym errors
    //
    float x[200]; float y[200]; float dy[200];float dx[200];
    float y1,y2,dy1;	
    for ( int i=0; i<200;i++){
        x[i]=(i+0.5)*0.02;
        y1=a1+b1*x[i];
        y2=a2+b2*x[i];
        y[i]=((y1+y2)/2.)+ offset;
        dy[i]=(y1-y2)/2;
        dx[i]=0.0;
    }
    TGraph *gr1 = new TGraphErrors(200,x,y,dx,dy);
    //gr1->Print("all");
    return gr1;
    
}


void theoryError(TString PLOT="Ydif", float OFFSET=0 ) 
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
/*     TCanvas *canvas = new TCanvas("canvas","",600,600); */
/*     canvas->cd();  */
    scale1->SetMaximum(1.2);
    scale1->SetMinimum(0.8);
    colorIt(scale1, kBlack, 3000);
/*     scale1->Draw("pe1same");// works */
    
    // end of test
    
    // put errors for PLOT="YSum"
    // OFFSET default is Y=0 on a linear plot
    
    TGraph * TGa = getParams( PLOT ,"scale", OFFSET);
    colorJt(TGa,kGreen-2, 3006);
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





