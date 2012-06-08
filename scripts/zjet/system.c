TGraph * getParams(int VAR, int DIST , float offset){
	//
	// takes the scale and pdf errors parameterized by 
	// two lines and created a tgraph
	//
	// DIST=0(YSUM),1(YDIF),2(ELSE)
	// VAR=0(SCALE) ,1(PDF), 3(correlated sytematic)
	//
	float p0=999;
	float p1=000;
	float q0=999;
	float q1=999;
	if(VAR==0){
		if ( DIST==0 ){
			cout << "Inside getParams: SCALE, YSUM:";
			//VAR Ysum GAM=0
			//****************************************
			//Minimizer is Linear
			//Chi2                      =      9.04276
			//NDf                       =           10
			p0                        =       1.0369;    
			p1                        =    -0.018478;     
			//****************************************
			//Minimizer is Linear
			//Chi2                      =      6.50297
			//NDf                       =           10
			q0                        =     0.948965 ;   
			q1                        =  0.000120763;    
		}elseif( DIST==1){
			cout << "Inside getParams: SCALE, YDIF:";
			//VAR Ydif GAM=0
			///****************************************
			//Minimizer is Linear
			//Chi2                      =      32.1349
			//NDf                       =           10
			p0                        =     0.997267;    
			p1                        =    0.0550207;    
			//****************************************
			//Minimizer is Linear
			//Chi2                      =      9.67188
			//NDf                       =           10
			q0                        =      0.95367 ;    
			q1                        =   -0.0126505 ;   
			//****************************************
		}else {
			cout << "Inside getParams: UNDEFINED";
			//other distributions
		}

		//****************************************
		
	}
	elseif(VAR==1){
		
		if(DIST ==1){
			cout << "Inside getParams: PDFVAR,YDIF:";
			//PDFVAR YDIF GAM=0
			//****************************************
			////Minimizer is Linear
			//Chi2                      =      33.9173
			//NDf                       =           10
			p0                        =      1.03905;     
			p1                        =    0.0215708 ;    
			//****************************************
			//Minimizer is Linear
			//Chi2                      =      19.5871
			//NDf                       =           10
			q0                        =     0.991011;     
			q1                        = -0.000168945 ;   
		}else if(DIST==0){
			cout << "Inside getParams: PDFVAR,YSUM:";
			//PDFVAR Ysum GAM=0
			//****************************************
			//Minimizer is Linear
			//Chi2                      =       10.773
			//NDf                       =           10
			p0                        =      1.05264;    
			p1                        =  -0.00800108 ;    
			//****************************************
			//Minimizer is Linear
			//Chi2                      =      3.18639
			//NDf                       =           10
			q0                        =     0.991236 ;    
			q1                        =  -0.00444716 ;  
		}else {
			cout << "Inside getParams: UNDEFINED";
			// other dostributions
		}

	}
	else {// VAR=2
		
		if ( DIST==0 ){	
			cout << "Inside getParams: UNCORR YSUM:" ;
			// this is where correlated systematics go
			// four parameters
			float p0=1.02;
			float p1=0.00;
			float q0=0.99;
			float q1=0.00;
			
		}elseif( DIST==1){
			cout << "Inside getParams: UNCORR YDIF:";
		    float p0=1.02;
		    float p1=0.00;
		    float q0=0.99;
		    float q1=0.00;
		}else {
		cout << "Inside getParams: UNDEFINED";
		// other dostributions
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
	float x[20]; float y[20]; float dy[20];float dx[20];
	float y1,y2,dy1;	
	for ( int i=0; i<20;i++){
		x[i]=(i+0.5)*0.2;
		y1=a1+b1*x[i];
		y2=a2+b2*x[i];
		y[i]=((y1+y2)/2.)+ offset;
		dy[i]=(y1-y2)/2;
        dx[i]=0.0;
	}
	TGraph *gr1 = new TGraphErrors(15,x,y,dx,dy);
	//gr1->Print("all");
	return gr1;
	
  }


void system(int DIST=0, float OFFSET=0 ) 
{
	/*
	 *  To be used to draw systematics on an existing plot
	 *
	 *  Created by linn on 6/2/12.
	 *  Copyright 2012 __FIU/CMS__. All rights reserved.
	 *
	 */
	
	//this is for testing
	/*
	TH1F *scale1 = new TH1F("scale1","scale1", 15, 0.0,3.0);
	scale1->FillRandom("gaus", 10);
	TCanvas *canvas = new TCanvas("canvas","",600,600);
	canvas->cd(); 
	scale1->SetMaximum(1.5);
	scale1->SetMinimum(0.5);
	colorIt(scale1, kYellow, 3000);
    scale1->Draw("");// works
	 */
	// end of test

	// put errors for Ydif (DIST=1), YSum(DIST=0)
	// OFFSET default is Y=1 on a linear plot
	
	TGraph * TGa = getParams( 0     ,DIST  ,OFFSET);
	colorJt(TGa,kBlue, 3006);
	TGa->Draw("4 same");// works
	
	TGraph * TGb = getParams( 1     ,DIST  ,OFFSET);
	colorJt(TGb,kRed, 3005);
    TGb->Draw("4 same");
	
	TGraph * TGc = getParams( 2     ,DIST  ,OFFSET);
	colorJt(TGc,kGreen, 3004);
    TGc->Draw("4 same");
	
	
	TLegend *legend1 = new TLegend(0.179412,0.234282,0.457653,0.483077);
	legend1->SetFillStyle(0);
	legend1->SetFillColor(0);
	legend1->SetBorderSize(0);
	legend1->SetTextSizePixels(20);
	legend1->AddEntry(TGa,"Scale","f");
	legend1->AddEntry(TGb,"PDF","f");
	legend1->AddEntry(TGc,"Correlated","f");
	legend1->Draw("same");
	
	
	//TGa->Print("all");

		
}





