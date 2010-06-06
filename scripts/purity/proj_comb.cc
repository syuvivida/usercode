#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "format.h"
#include "TMath.h"
#include "TRandom2.h"

#define W_MASS          80.403
#define Z_MASS          91.1876
#define MUON_MASS       0.105658
#define ELECTRON_MASS   0.0005109989
#define SB_TYPE         2
#define SPLITSAMPLE     0 //0 not split; 1 odd number for test; 2 even number for template
 
double dPhi(double p1,double p2)
{
        double dp = p1 - p2;
        if (TMath::Abs(dp+3.14159265358979323846*2.) < TMath::Abs(dp)) dp += 3.14159265358979323846*2.;
        else
        if (TMath::Abs(dp-3.14159265358979323846*2.) < TMath::Abs(dp)) dp -= 3.14159265358979323846*2.;

        return TMath::Abs(dp);
}

double dR(double e1, double e2, double p1, double p2)
{
        return sqrt(pow(e1-e2,2)+pow(dPhi(p1,p2),2));
}

struct sample_t {
char filename[128];
char tag[128];
double xsec; // in unit of pb
double ngen; // no. of events in generation(in unit of K)
}; 

// luminosity in pb^-1
#define TARGET_LUMINOSITY 0.1

enum{
_PhoJet15,	// 00
_PhoJet30,	// 01
_PhoJet80,	// 02
// _PhoJet170,	// 03
_QCD15,         // 04
_QCD30,         // 05
_QCD80,         // 06
_QCD170,         // 07
// _QCD300,         // 08
// _QCD470,         // 09
// _EGdata05,         // 09
_EGdata,           // 09
_sig_sum,		// 10
_bkg_sum,	        // 11
N_SAMPLES		// 12
};

enum{
_EB,			// 00
_EE,			// 01
//_ALL,			// 02
N_CATEGORIES		// 03
};

char tag_category[N_CATEGORIES][32] = {
"EB",
"EE",
};

char type_category[SB_TYPE][32] = {
"SIG",
"BKG",
};

struct sample_t sample[N_SAMPLES] = {
{"PhotonJet_Pt15.root",      "PhoJet15",        192200.  , 1143.390  },
{"PhotonJet_Pt30.root",      "PhoJet30",        20070.   ,  786.794  },
{"PhotonJet_Pt80.root",      "PhoJet80",        556.5    , 1207.711  },
// {"PhotonJet_Pt170.root",     "PhoJet170",       2.437    , 1139.400  },
// {"PhotonJet_Pt300.root",     "PhoJet300",       1.636    , 954.266  },
{"QCD_Pt15.root",               "QCD15",      876200000. , 6190.500  },
{"QCD_Pt30.root",               "QCD30",       60410000. , 5269.664  },
{"QCD_Pt80.root",               "QCD80",         923800. , 3221.800  },
{"QCD_Pt170.root",              "QCD170",         25470. , 3171.950  },
// {"QCD_Pt300.root",              "QCD300",          1256. , 3282.660  },
// {"QCD_Pt470.root",              "QCD470",         87.98  , 2000.  },
// {"EG_data_May05_2010.root",     "EGdata05",             1. , 0.0001  },
{"EGdata_132440-135735.root",     "EGdata",             1. , 0.0001  },
{"signal sum",			     "sig_sum", 	     -1.     ,-1    },
{"background sum",		     "bkg_sum", 	     -1.     ,-1    },
};

void proj_comb(bool sumw2=false)
{      

    TFile *fout = new TFile("proj_comb.root","recreate");
	
    TH1F *h_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_et_noIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_eta[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_phi[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_r9[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_r9_noIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_ecalIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_ecalIsoSB[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_hcalIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_trkIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_comb3Iso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_comb3IsoSB[N_CATEGORIES][N_SAMPLES][SB_TYPE];

    TH2F *h_comb3Iso_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_comb3IsoSB_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_ecalIso_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_ecalIsoSB_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_hcalIso_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_trkIso_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_r9_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_SIEIE_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_HoverE_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH2F *h_ESRatio_et[N_CATEGORIES][N_SAMPLES][SB_TYPE];

    TH1F *h_ecalIsoEt[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_hcalIsoEt[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_comb3IsoEt[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_trkIsoEt[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_HoverE[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_SEE[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_SIEIE[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_ESRatio[N_CATEGORIES][N_SAMPLES][SB_TYPE];    

    TH1F *h_HoverE_noIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_SIEIE_noIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];
    TH1F *h_ESRatio_noIso[N_CATEGORIES][N_SAMPLES][SB_TYPE];    

    
    double evtcount[N_CATEGORIES][N_SAMPLES][SB_TYPE][2];
    double sigcount[N_SAMPLES][SB_TYPE][2];
     
    char buffer[128];
    
    for(int file=0; file<N_SAMPLES; file++) { 
    	for(int cate=0; cate<N_CATEGORIES; cate++) {
	    for(int type=0; type<SB_TYPE; type++) {
	    
	    	sprintf(buffer,"h_%s_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_et[cate][file][type] = new TH1F(buffer,buffer,  60, 0., 300.);
		if(sumw2)h_et[cate][file][type]->Sumw2();
		  
	    	sprintf(buffer,"h_%s_et_noIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_et_noIso[cate][file][type] = new TH1F(buffer,buffer,  60, 0., 300.);	       
		if(sumw2)h_et_noIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_eta_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_eta[cate][file][type] = new TH1F(buffer,buffer,  60, -3.,3.);	       
		if(sumw2)h_eta[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_phi_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_phi[cate][file][type] = new TH1F(buffer,buffer,  60, TMath::Pi()*-1.1, TMath::Pi()*1.1);	       
		if(sumw2)h_phi[cate][file][type]->Sumw2();
		
	    	sprintf(buffer,"h_%s_r9_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_r9[cate][file][type] = new TH1F(buffer,buffer,  110, 0.,1.1);		    
		if(sumw2)h_r9[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_r9_noIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_r9_noIso[cate][file][type] = new TH1F(buffer,buffer,  110, 0.,1.1);		    
		if(sumw2)h_r9_noIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ecalIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ecalIso[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_ecalIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ecalIsoSB_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ecalIsoSB[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_ecalIsoSB[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ecalIso_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ecalIso_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		    
		if(sumw2)h_ecalIso_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ecalIsoSB_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ecalIsoSB_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		    	
		if(sumw2)h_ecalIsoSB_et[cate][file][type]->Sumw2();
		
	    	sprintf(buffer,"h_%s_comb3IsoSB_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_comb3IsoSB_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		   
		if(sumw2)h_comb3IsoSB_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_hcalIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_hcalIso[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_hcalIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_trkIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_trkIso[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_trkIso[cate][file][type]->Sumw2();

		// added by Eiko
	    	sprintf(buffer,"h_%s_hcalIso_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_hcalIso_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		    
		if(sumw2)h_hcalIso_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_trkIso_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_trkIso_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		    
		if(sumw2)h_trkIso_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_r9_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_r9_et[cate][file][type] = new TH2F(buffer,buffer,  110, 0., 1.1, 60, 0., 300.);		    
		if(sumw2)h_r9_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_SIEIE_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_SIEIE_et[cate][file][type] = new TH2F(buffer,buffer,  50, 0., 0.05, 60, 0., 300.);		    
		if(sumw2)h_SIEIE_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ESRatio_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ESRatio_et[cate][file][type] = new TH2F(buffer,buffer, 48, -0.1000000001,1.0999999999, 60, 0., 300.);		    
		if(sumw2)h_ESRatio_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_HoverE_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_HoverE_et[cate][file][type] = new TH2F(buffer,buffer, 40, -0.02,0.18, 60, 0., 300.);		    
		if(sumw2)h_HoverE_et[cate][file][type]->Sumw2();

		// 
	    	sprintf(buffer,"h_%s_comb3Iso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_comb3Iso[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_comb3Iso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_comb3IsoSB_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_comb3IsoSB[cate][file][type] = new TH1F(buffer,buffer,  120, -1., 11.);		    
		if(sumw2)h_comb3IsoSB[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_comb3Iso_et_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_comb3Iso_et[cate][file][type] = new TH2F(buffer,buffer,  120, -1., 11., 60, 0., 300.);		    
		if(sumw2)h_comb3Iso_et[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ecalIsoEt_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ecalIsoEt[cate][file][type] = new TH1F(buffer,buffer,  80, -0.05,0.35);		    
		if(sumw2)h_ecalIsoEt[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_hcalIsoEt_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_hcalIsoEt[cate][file][type] = new TH1F(buffer,buffer,  80, -0.05,0.35);		    
		if(sumw2)h_hcalIsoEt[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_trkIsoEt_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_trkIsoEt[cate][file][type] = new TH1F(buffer,buffer,  80, -0.05,0.35);		    
		if(sumw2)h_trkIsoEt[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_comb3IsoEt_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_comb3IsoEt[cate][file][type] = new TH1F(buffer,buffer,  120, -0.1,1.1);		    
		if(sumw2)h_comb3IsoEt[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_HoverE_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_HoverE[cate][file][type] = new TH1F(buffer,buffer,  40, -0.02,0.18);		    
		if(sumw2)h_HoverE[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_SEE_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_SEE[cate][file][type] = new TH1F(buffer,buffer,  50, 0.,0.1);		    
		if(sumw2)h_SEE[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_SIEIE_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_SIEIE[cate][file][type] = new TH1F(buffer,buffer,  50, 0.,0.05);		    
		if(sumw2)h_SIEIE[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ESRatio_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ESRatio[cate][file][type] = new TH1F(buffer,buffer,  48, -0.1000000001,1.0999999999);		    
		if(sumw2)h_ESRatio[cate][file][type]->Sumw2();

		//noiso
	    	sprintf(buffer,"h_%s_HoverE_noIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_HoverE_noIso[cate][file][type] = new TH1F(buffer,buffer,  40, -0.02,0.18);		    
		if(sumw2)h_HoverE_noIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_SIEIE_noIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_SIEIE_noIso[cate][file][type] = new TH1F(buffer,buffer,  50, 0.,0.05);		    
		if(sumw2)h_SIEIE_noIso[cate][file][type]->Sumw2();

	    	sprintf(buffer,"h_%s_ESRatio_noIso_%s_%s",tag_category[cate],sample[file].tag, type_category[type]);
    	    	h_ESRatio_noIso[cate][file][type] = new TH1F(buffer,buffer,  48, -0.1000000001,1.0999999999);		    
		if(sumw2)h_ESRatio_noIso[cate][file][type]->Sumw2();

	    	evtcount[cate][file][type][0] = evtcount[cate][file][type][1] = 0.; 
	    }	
	    sigcount[file][0][0] = sigcount[file][0][1] = 0.;
	    sigcount[file][1][0] = sigcount[file][1][1] = 0.;
	}
    }
    
    printf(" %d samples need to be projected \n", N_SAMPLES );
    for(int file=0; file<N_SAMPLES; file++) {
    
    	if (sample[file].xsec<0.) continue;
		
        TFile *f1 = new TFile(sample[file].filename);
        TTree *root = new TTree();
	
	double scaling_factor = TARGET_LUMINOSITY*sample[file].xsec/sample[file].ngen/1000.;
	if ( SPLITSAMPLE > 0 ) scaling_factor *= 2.; //*2. if only read half
	
	int _isSigMC=0; int _isBkgMC=0; int _isData=0;
	if ( strcmp(sample[file].filename,"PhotonJet_Pt15.root")==0 ||
	     strcmp(sample[file].filename,"PhotonJet_Pt30.root")==0 ||
	     strcmp(sample[file].filename,"PhotonJet_Pt80.root")==0 ||
	     strcmp(sample[file].filename,"PhotonJet_Pt170.root")==0 ||
	     strcmp(sample[file].filename,"PhotonJet_Pt300.root")==0 )		     
	  { _isSigMC=1; 
	    printf(" loading Signal MC \n"); 
	    root = (TTree*)f1->FindObjectAny("Analysis");
	  }
	if ( strcmp(sample[file].filename,"QCD_Pt15.root")==0 ||
	     strcmp(sample[file].filename,"QCD_Pt30.root")==0 ||
	     strcmp(sample[file].filename,"QCD_Pt80.root")==0 ||
	     strcmp(sample[file].filename,"QCD_Pt170.root")==0 ||
	     strcmp(sample[file].filename,"QCD_Pt300.root")==0 ||
	     strcmp(sample[file].filename,"QCD_Pt470.root")==0 )		     
	  { _isBkgMC=1; printf(" loading QCD MC \n"); 
	    root = (TTree*)f1->FindObjectAny("Analysis");
	  }

	if ( strcmp(sample[file].filename,"EGdata_132440-135735.root")==0)
	  { _isData=1; 
	    scaling_factor=1.;
	    root = (TTree*)f1->FindObjectAny("Analysis");
	  }

	printf("Now processing %s, scaling factor = %g\n",sample[file].filename,scaling_factor);

	double pthat_min=0; double pthat_max=0.;
	if (strcmp(sample[file].filename,"PhotonJet_Pt15.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt15.root")==0 ) {
	  pthat_min=15; pthat_max=30; printf("only take pthat %5.1f ~ %5.1f \n", pthat_min,pthat_max);}
	if (strcmp(sample[file].filename,"PhotonJet_Pt30.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt30.root")==0 ) {
	  pthat_min=30; pthat_max=80; printf("only take pthat  %5.1f ~ %5.1f \n", pthat_min,pthat_max);}
	if (strcmp(sample[file].filename,"PhotonJet_Pt80.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt80.root")==0 ) {
	  pthat_min=80; pthat_max=170; printf("only take pthat %5.1f ~ %5.1f \n", pthat_min,pthat_max);}
	if (strcmp(sample[file].filename,"PhotonJet_Pt170.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt170.root")==0 ) {
	  pthat_min=170; pthat_max=300; printf("only take pthat %5.1f ~ %5.1f \n", pthat_min,pthat_max);}
	if (strcmp(sample[file].filename,"PhotonJet_Pt300.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt300.root")==0 ) {
	  pthat_min=300; pthat_max=470; printf("only take pthat %5.1f ~ %5.1f \n", pthat_min,pthat_max);}
	if (strcmp(sample[file].filename,"PhotonJet_Pt470.root")==0 ||
	    strcmp(sample[file].filename,"QCD_Pt470.root")==0 ) {
	  pthat_min=470; pthat_max=1000; printf("only take pthat %5.1f ~ %5.1f \n", pthat_min,pthat_max);}

	EvtInfoBranches EvtInfo;
        
	EvtInfo.Register(root);
		
	FILE *ffISO;
	if (_isData==1)  ffISO = fopen("EGdata.dat","w");

	int Naccept=0; int reject=0;
	TRandom2 *trd = new TRandom2();
	printf(" %d entries to be loaded \n", root->GetEntries());
   	for(int entry=0;entry<root->GetEntries();entry++) {
//       	for(int entry=0;entry<100000;entry++) {
                root->GetEntry(entry);
 		double tr = trd->Uniform(0,1);
		if (entry%100==0)
		printf("%4.1f%% done.\r",(float)entry/(float)root->GetEntriesFast()*100.);		

		double FACTORS[20];

		FACTORS [ 0] =  5.;//2.1;
		FACTORS [ 1] =  0.0 ;
		FACTORS [ 2] =  0.0 ;
		FACTORS [ 3] =  0.0500 ;
		FACTORS [ 4] =  0.0100 ;

		FACTORS [ 5] =  5.;//1.45;
		FACTORS [ 6] =  0.00 ;
		FACTORS [ 7] =  0.00 ;
		FACTORS [ 8] =  0.0500 ;
		FACTORS [ 9] =  0.1000 ;

		double EB_ECALISO_MAX       = FACTORS[0];	
		double EB_HCALISO_MAX       = FACTORS[1];	
		double EB_TRKISO_MAX        = FACTORS[2];	
		double EB_HOVERE_MAX        = FACTORS[3];	
		double EB_SIEIE_MAX         = FACTORS[4];	

		double EE_ECALISO_MAX       = FACTORS[5];	
		double EE_HCALISO_MAX       = FACTORS[6];	
		double EE_TRKISO_MAX        = FACTORS[7];	
		double EE_HOVERE_MAX        = FACTORS[8];	
		double EE_ESRATIO_MIN       = FACTORS[9];	
		
		// cut bits, reserved for histogram booking
		enum {
		_ECALISO,
		_HCALISO,
		_TRKISO,
		_HOVERE,
		_SIEIE,
		_ESRATIO
		};
		
		// Trigger bit, single electron or single muon
// 		bool accept = false;
// 		if (EvtInfo.TrgBook[0]==1) accept = true;
// 		if (EvtInfo.TrgBook[2]==1) accept = true;
// 		if (!accept) cut_bits |= (1<<_HLT);
 		if ( EvtInfo.HLT_Photon15_L1R != true ) continue;

 		if ( _isData==1 && !(!EvtInfo.TTBit[36] && !EvtInfo.TTBit[37] && !EvtInfo.TTBit[38] && !EvtInfo.TTBit[39] && !EvtInfo.vtxIsFake && EvtInfo.vtxNdof > 4 && TMath::Abs(EvtInfo.vtxZ) <= 15) ) continue;

// 		if ( (_isSigMC==1 || _isBkgMC==1) && (EvtInfo.ptHat < pthat_min || EvtInfo.ptHat > pthat_max) ) continue;
		int type=0; // for signal or background type	
		
		for ( int ipho=0; ipho<EvtInfo.nPhotons; ipho++) {

		  int cut_bits = 0;
		  
		  if ( (_isSigMC==1 || _isBkgMC==1) && (EvtInfo.et[ipho] < pthat_min || EvtInfo.et[ipho] > pthat_max) ) continue;
// 		  if ( _isSigMC==1 && !(EvtInfo.isGenMatched[ipho]==1&&EvtInfo.genCalIsoDR04[ipho]<5.0) ) continue;
// 		  if ( _isBkgMC==1 && EvtInfo.genCalIsoDR04[ipho]<=5.0 ) continue;			  		  
  		  if ( (_isSigMC==1 || _isBkgMC==1) && (EvtInfo.isGenMatched[ipho]==1&& TMath::Abs(EvtInfo.genMomId[ipho])<=22&&EvtInfo.genCalIsoDR04[ipho]<5.0) ) type=0;
//   		  if ( (_isSigMC==1 || _isBkgMC==1) && (EvtInfo.isGenMatched[ipho]==1&& TMath::Abs(EvtInfo.genMomId[ipho])<=22) ) type=0;
 		  else if(_isSigMC==1 || _isBkgMC==1) type=1;
		  
		  int good_LS=0;
//		  if (
//		      ){
//		    good_LS = 1;
//		  }
// 		  if ( _isData==1 && good_LS!=1 ) continue;
		  
		  if ( EvtInfo.et[ipho]<15. ) continue;
//  		  if ( EvtInfo.isLoose[ipho] != 1 ) continue;
// 		  if ( EvtInfo.isEBGap[ipho]==1 || EvtInfo.isEEGap[ipho]==1 || EvtInfo.isEBEEGap[ipho]==1 ) continue;
		  
		  if ( TMath::Abs(EvtInfo.eta[ipho]) > 1.45 && TMath::Abs(EvtInfo.eta[ipho]) < 1.7 ) continue;
		  if ( TMath::Abs(EvtInfo.eta[ipho] > 2.5) ) continue;

		  //do selection after selection
		  Naccept++;
		  //if ( _isData==0 && SPLITSAMPLE==1 && tr>0.5) {
		  if ( _isData==0 && SPLITSAMPLE==1 && Naccept%2==1) {
		    reject++;
		    continue; // read even events for testing}
		  }
		  //if ( _isData==0 && SPLITSAMPLE==2 && tr<0.5) {
		  if ( _isData==0 && SPLITSAMPLE==2 && Naccept%2==0) {
		    reject++;
		    continue; // read odd events for template
		  }

		  double comb3Iso = EvtInfo.ecalRecHitSumEtConeDR04[ipho] + 
		    EvtInfo.hcalTowerSumEtConeDR04[ipho] + EvtInfo.trkSumPtHollowConeDR04[ipho] ;

		  int cate = -1;
		  if ( TMath::Abs(EvtInfo.eta[ipho])<1.45 ) {
		    cate = _EB;
		    //remove EB spike by Swiss-cross or 1st bin of sigmaietaieta
// 		    if ( _isData==1 && 1-((EvtInfo.eRight[ipho]+EvtInfo.eLeft[ipho]+EvtInfo.eTop[ipho]+EvtInfo.eBottom[ipho])/EvtInfo.eMax[ipho]) > 0.95 ) continue;
		    
		    if ( _isData==1 && EvtInfo.isEB[ipho] && 
			 !(EvtInfo.seedSeverity[ipho]!=3&&EvtInfo.seedSeverity[ipho]!=4&&EvtInfo.seedRecoFlag[ipho]!=2) ) continue;
 		    if ( EvtInfo.sigmaIetaIeta[ipho]<0.002 ) continue;

		    //only take spike
// 		    if ( _isData==1 && EvtInfo.isEB[ipho] && 
// 			 (EvtInfo.seedSeverity[ipho]!=3&&EvtInfo.seedSeverity[ipho]!=4&&EvtInfo.seedRecoFlag[ipho]!=2) ) continue;
// 		    if ( EvtInfo.sigmaIetaIeta[ipho]>0.002 ) continue;

		    
//  		    if ( EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] > EB_ECALISO_MAX) cut_bits |= (1<<_ECALISO);
//  		    if ( EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] > EB_HCALISO_MAX) cut_bits |= (1<<_HCALISO);
//    		    if ( EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho] > EB_TRKISO_MAX) cut_bits |= (1<<_TRKISO);
 
		    if ( comb3Iso > EB_ECALISO_MAX ) cut_bits |= (1<<_ECALISO);
  		    if ( EvtInfo.hadronicOverEm[ipho] > EB_HOVERE_MAX) cut_bits |= (1<<_HOVERE);
   		    if ( EvtInfo.sigmaIetaIeta[ipho] > EB_SIEIE_MAX) cut_bits |= (1<<_SIEIE);


		  }


		  if ( TMath::Abs(EvtInfo.eta[ipho])>1.7 && TMath::Abs(EvtInfo.eta[ipho])<2.5) {
		    cate = _EE;
//   		    if ( EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] > EE_ECALISO_MAX) cut_bits |= (1<<_ECALISO);
//   		    if ( EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] > EE_HCALISO_MAX) cut_bits |= (1<<_HCALISO);
//    		    if ( EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho] > EE_TRKISO_MAX) cut_bits |= (1<<_TRKISO);
 		    if ( comb3Iso > EE_ECALISO_MAX ) cut_bits |= (1<<_ECALISO);
   		    if ( EvtInfo.hadronicOverEm[ipho] > EE_HOVERE_MAX) cut_bits |= (1<<_HOVERE);
   		    if ( TMath::Abs(EvtInfo.ESRatio[ipho]) < EE_ESRATIO_MIN ) cut_bits |= (1<<_ESRATIO);
		  }

		  if ( cate==-1 ) continue; 

//  		  float ecalhcal_fisher = 0.324 +
//  		    EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] * (-1.541) +
//  		    EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] *( -0.943 );
//  		  if  ( TMath::Abs(EvtInfo.eta[ipho])>1.55 ) {
//  		    ecalhcal_fisher = 0.635 +
//  		      EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] * (-6.803) +
//  		      EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] *( -5.264 );
//  		  }
// 		  float ecalhcal_fisher = 0.607 +
// 		    EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] * (-2.490) +
// 		    EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] *(-1.181) +
// 		    EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho] *(-0.948);
// 		  if  ( TMath::Abs(EvtInfo.eta[ipho])>1.55 ) {
// 		    ecalhcal_fisher = 0.557 +
// 		      EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho] * (-5.010) +
// 		      EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho] *( -4.202 ) +
// 		      EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho] *(-0.206);
// 		  }


		  if ((cut_bits & (~(1<<_ECALISO)))==0){
		    h_ecalIso[cate][file][type]->Fill(EvtInfo.ecalRecHitSumEtConeDR04[ipho],scaling_factor);
		    h_ecalIsoEt[cate][file][type]->Fill(EvtInfo.ecalRecHitSumEtConeDR04[ipho]/EvtInfo.et[ipho],scaling_factor);
		    h_ecalIso_et[cate][file][type]->Fill(EvtInfo.ecalRecHitSumEtConeDR04[ipho], EvtInfo.et[ipho],scaling_factor);

 		    h_hcalIso[cate][file][type]->Fill(EvtInfo.hcalTowerSumEtConeDR04[ipho],scaling_factor); 
 		    h_hcalIsoEt[cate][file][type]->Fill(EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho],scaling_factor);
		    h_hcalIso_et[cate][file][type]->Fill(EvtInfo.hcalTowerSumEtConeDR04[ipho], EvtInfo.et[ipho],scaling_factor);

  		    h_trkIso[cate][file][type]->Fill(EvtInfo.trkSumPtHollowConeDR04[ipho],scaling_factor);  		
  		    h_trkIsoEt[cate][file][type]->Fill(EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho],scaling_factor);
		    h_trkIso_et[cate][file][type]->Fill(EvtInfo.trkSumPtHollowConeDR04[ipho], EvtInfo.et[ipho],scaling_factor);

		    h_comb3Iso[cate][file][type]->Fill(comb3Iso, scaling_factor); 
		    h_comb3IsoEt[cate][file][type]->Fill(comb3Iso, scaling_factor); 
		    h_et_noIso[cate][file][type]->Fill(EvtInfo.et[ipho],scaling_factor);  
		    h_comb3Iso_et[cate][file][type]->Fill(comb3Iso, EvtInfo.et[ipho], scaling_factor); 
		    h_r9_noIso[cate][file][type]->Fill(EvtInfo.r9[ipho],scaling_factor);  
		    h_r9_et[cate][file][type]->Fill(EvtInfo.r9[ipho], EvtInfo.et[ipho],scaling_factor);

		    if(_isData==1)fprintf(ffISO,"%f  %f %f \n",comb3Iso, EvtInfo.et[ipho], EvtInfo.eta[ipho]);
		  }

		  if ((cut_bits & (~(1<<_ECALISO | 1<<_HOVERE)))==0 && (cut_bits & 1<<_HOVERE)){
		    h_ecalIsoSB[cate][file][type]->Fill(EvtInfo.ecalRecHitSumEtConeDR04[ipho],scaling_factor);
		    h_ecalIsoSB_et[cate][file][type]->Fill(EvtInfo.ecalRecHitSumEtConeDR04[ipho], EvtInfo.et[ipho],scaling_factor);
		  }

// 		  if ((cut_bits & (~(1<<_HCALISO | 1<<_ECALISO | 1<<_TRKISO)))==0) {
// 		    h_comb3Iso[cate][file][type]->Fill(comb3Iso, scaling_factor); 
// 		    h_comb3IsoEt[cate][file][type]->Fill(comb3Iso, scaling_factor); 
// 		    h_et_noIso[cate][file][type]->Fill(EvtInfo.et[ipho],scaling_factor);  
// 		    h_comb3Iso_et[cate][file][type]->Fill(comb3Iso, EvtInfo.et[ipho], scaling_factor); 
// 		    h_r9_noIso[cate][file][type]->Fill(EvtInfo.r9[ipho],scaling_factor);  
// 		  }

		  if ( cate == _EB ) {
		    if (  (cut_bits & (~(1<<_ECALISO | 1<<_SIEIE)))==0 ) {
		      if ( EvtInfo.sigmaIetaIeta[ipho] > 0.011 && EvtInfo.sigmaIetaIeta[ipho] < 0.012 ) {
			h_comb3IsoSB[cate][file][type]->Fill(comb3Iso, scaling_factor); 
			h_comb3IsoSB_et[cate][file][type]->Fill(comb3Iso, EvtInfo.et[ipho], scaling_factor); 
		      }		 
		    }
		  }
		  if ( cate == _EE ) {
		    if (  (cut_bits & (~(1<<_ECALISO | 1<<_ESRATIO)))==0 ) {
		      if ( TMath::Abs(EvtInfo.ESRatio[ipho])<0.05 ) {
			h_comb3IsoSB[cate][file][type]->Fill(comb3Iso, scaling_factor); 
			h_comb3IsoSB_et[cate][file][type]->Fill(comb3Iso, EvtInfo.et[ipho], scaling_factor); 
		      }
		    }		 
		  }
		  
// 		  if ((cut_bits & (~(1<<_HCALISO)))==0){
// 		    h_hcalIso[cate][file][type]->Fill(EvtInfo.hcalTowerSumEtConeDR04[ipho],scaling_factor); 
// 		    h_hcalIsoEt[cate][file][type]->Fill(EvtInfo.hcalTowerSumEtConeDR04[ipho]/EvtInfo.et[ipho],scaling_factor);
//  		  }

//  		  if ((cut_bits & (~(1<<_TRKISO)))==0){
//  		    h_trkIso[cate][file][type]->Fill(EvtInfo.trkSumPtHollowConeDR04[ipho],scaling_factor);  		
//  		    h_trkIsoEt[cate][file][type]->Fill(EvtInfo.trkSumPtHollowConeDR04[ipho]/EvtInfo.et[ipho],scaling_factor);  		
//  		  }
		  if ((cut_bits & (~(1<<_HOVERE)))==0)
		    {
		      h_HoverE[cate][file][type]->Fill(EvtInfo.hadronicOverEm[ipho],scaling_factor);  		
		      h_HoverE_et[cate][file][type]->Fill(EvtInfo.hadronicOverEm[ipho], EvtInfo.et[ipho],scaling_factor);
		    }

		  if ((cut_bits & (~(1<<_HOVERE |1<<_HCALISO | 1<<_ECALISO | 1<<_TRKISO )))==0)
		    h_HoverE_noIso[cate][file][type]->Fill(EvtInfo.hadronicOverEm[ipho],scaling_factor);  		

		  if ((cut_bits & (~(1<<_SIEIE)))==0)
		    {
		      h_SIEIE[cate][file][type]->Fill(EvtInfo.sigmaIetaIeta[ipho],scaling_factor);  		
		      h_SIEIE_et[cate][file][type]->Fill(EvtInfo.sigmaIetaIeta[ipho], EvtInfo.et[ipho],scaling_factor);
		    }
		  if ((cut_bits & (~(1<<_SIEIE |1<<_HCALISO | 1<<_ECALISO | 1<<_TRKISO )))==0)
		    h_SIEIE_noIso[cate][file][type]->Fill(EvtInfo.sigmaIetaIeta[ipho],scaling_factor);  		

		  if ((cut_bits & (~(1<<_ESRATIO)))==0)
		    {
		      h_ESRatio[cate][file][type]->Fill(TMath::Abs(EvtInfo.ESRatio[ipho]),scaling_factor);  		
		      h_ESRatio_et[cate][file][type]->Fill(TMath::Abs(EvtInfo.ESRatio[ipho]), EvtInfo.et[ipho],scaling_factor);
		    }
		  if ((cut_bits & (~(1<<_ESRATIO |1<<_HCALISO | 1<<_ECALISO | 1<<_TRKISO )))==0)
		    h_ESRatio_noIso[cate][file][type]->Fill(TMath::Abs(EvtInfo.ESRatio[ipho]),scaling_factor);  		
		  
		  if (cut_bits!=0) continue; //do the full cut, book et, and do counting
		
		  h_et[cate][file][type]->Fill(EvtInfo.et[ipho],scaling_factor);  
		  h_eta[cate][file][type]->Fill(EvtInfo.eta[ipho],scaling_factor);  
		  h_phi[cate][file][type]->Fill(EvtInfo.phi[ipho],scaling_factor);  
		  h_r9[cate][file][type]->Fill(EvtInfo.r9[ipho],scaling_factor);  
		  h_SEE[cate][file][type]->Fill(EvtInfo.sigmaEtaEta[ipho],scaling_factor);  
		
		  evtcount[cate][file][type][0]+=1.;		
		  if (cate==_EB || cate==_EE) sigcount[file][type][0]+=1.;	
		}//end of ipho loop
        }	
	
	printf("100.0%% done.\n");
	
	for(int cate=0; cate<N_CATEGORIES; cate++) {
	  for(int type=0; type<SB_TYPE; type++) {
		evtcount[cate][file][type][1] = sqrt(evtcount[cate][file][type][0])*scaling_factor;
		if (evtcount[cate][file][type][1]<scaling_factor) evtcount[cate][file][type][1] = scaling_factor;
        	evtcount[cate][file][type][0] =      evtcount[cate][file][type][0] *scaling_factor;
	  }
	}
	
	for(int type=0; type<SB_TYPE; type++) {
	  sigcount[file][type][1] = sqrt(sigcount[file][type][0])*scaling_factor;
	  if (sigcount[file][type][1]<scaling_factor) sigcount[file][type][1] = scaling_factor;
	  sigcount[file][type][0] =      sigcount[file][type][0] *scaling_factor;
	}
	printf("%s processed %d, rejected %d", sample[file].filename, Naccept, reject);

	if(_isData==1) fclose(ffISO);
	delete f1;	
    }
	
    for(int file=_PhoJet15; file<=_QCD170; file++) {		
      for(int type=0; type<SB_TYPE; type++) {	  
        for(int cate=0; cate<N_CATEGORIES; cate++) {
	  h_et           [cate][_sig_sum][type]->Add(h_et[cate][file][type]);
	  h_et_noIso     [cate][_sig_sum][type]->Add(h_et_noIso[cate][file][type]);
	  h_eta          [cate][_sig_sum][type]->Add(h_eta[cate][file][type]);
	  h_phi          [cate][_sig_sum][type]->Add(h_phi[cate][file][type]);
	  h_r9           [cate][_sig_sum][type]->Add(h_r9[cate][file][type]);
	  h_r9_noIso     [cate][_sig_sum][type]->Add(h_r9_noIso     [cate][file][type]);
	  h_ecalIso[cate][_sig_sum][type]->Add(h_ecalIso[cate][file][type]);
	  h_ecalIsoSB[cate][_sig_sum][type]->Add(h_ecalIsoSB[cate][file][type]);
	  h_ecalIso_et[cate][_sig_sum][type]->Add(h_ecalIso_et[cate][file][type]);

	  // added by Eiko
	  h_hcalIso_et[cate][_sig_sum][type]->Add(h_hcalIso_et[cate][file][type]);
	  h_trkIso_et[cate][_sig_sum][type]->Add(h_trkIso_et[cate][file][type]);
	  h_SIEIE_et[cate][_sig_sum][type]->Add(h_SIEIE_et[cate][file][type]);
	  h_r9_et[cate][_sig_sum][type]->Add(h_r9_et[cate][file][type]);
	  h_HoverE_et[cate][_sig_sum][type]->Add(h_HoverE_et[cate][file][type]);
	  h_ESRatio_et[cate][_sig_sum][type]->Add(h_ESRatio_et[cate][file][type]);

	  h_ecalIsoSB_et[cate][_sig_sum][type]->Add(h_ecalIsoSB_et[cate][file][type]);
	  h_hcalIso[cate][_sig_sum][type]->Add(h_hcalIso[cate][file][type]);
	  h_comb3Iso[cate][_sig_sum][type]->Add(h_comb3Iso[cate][file][type]);
	  h_comb3IsoSB[cate][_sig_sum][type]->Add(h_comb3IsoSB[cate][file][type]);
	  h_comb3IsoSB_et[cate][_sig_sum][type]->Add(h_comb3IsoSB_et[cate][file][type]);
	  h_comb3Iso_et[cate][_sig_sum][type]->Add(h_comb3Iso_et[cate][file][type]);
	  h_trkIso [cate][_sig_sum][type]->Add(h_trkIso [cate][file][type]);
	  h_ecalIsoEt[cate][_sig_sum][type]->Add(h_ecalIsoEt[cate][file][type]);
	  h_hcalIsoEt[cate][_sig_sum][type]->Add(h_hcalIsoEt[cate][file][type]);
	  h_comb3IsoEt[cate][_sig_sum][type]->Add(h_comb3IsoEt[cate][file][type]);
	  h_trkIsoEt [cate][_sig_sum][type]->Add(h_trkIsoEt [cate][file][type]);
	  h_HoverE [cate][_sig_sum][type]->Add(h_HoverE [cate][file][type]);
	  h_SEE    [cate][_sig_sum][type]->Add(h_SEE    [cate][file][type]);
	  h_SIEIE  [cate][_sig_sum][type]->Add(h_SIEIE  [cate][file][type]);
	  h_ESRatio[cate][_sig_sum][type]->Add(h_ESRatio[cate][file][type]);
	  h_HoverE_noIso [cate][_sig_sum][type]->Add(h_HoverE_noIso [cate][file][type]);
	  h_SIEIE_noIso  [cate][_sig_sum][type]->Add(h_SIEIE_noIso  [cate][file][type]);
	  h_ESRatio_noIso[cate][_sig_sum][type]->Add(h_ESRatio_noIso[cate][file][type]);

	  evtcount[cate][_sig_sum][type][0] += evtcount[cate][file][type][0];
	  evtcount[cate][_sig_sum][type][1] =  sqrt(pow(evtcount[cate][_sig_sum][type][1],2) + pow(evtcount[cate][file][type][1],2));
	  }
	sigcount[_sig_sum][type][0] += sigcount[file][type][0];
	sigcount[_sig_sum][type][1] =  sqrt(pow(sigcount[_sig_sum][type][1],2) + pow(sigcount[file][type][1],2));
      }
    }
    
    for(int file=_PhoJet15; file<=_QCD170; file++) {		
	for(int type=0; type<SB_TYPE; type++) {	  
	  for(int cate=0; cate<N_CATEGORIES; cate++) {	
	  h_et     [cate][_bkg_sum][type]->Add(h_et     [cate][file][type]);
	  h_et_noIso     [cate][_bkg_sum][type]->Add(h_et_noIso     [cate][file][type]);
	  h_eta    [cate][_bkg_sum][type]->Add(h_eta    [cate][file][type]);
	  h_phi    [cate][_bkg_sum][type]->Add(h_phi    [cate][file][type]);
	  h_r9     [cate][_bkg_sum][type]->Add(h_r9     [cate][file][type]);
	  h_r9_noIso     [cate][_bkg_sum][type]->Add(h_r9_noIso     [cate][file][type]);
	  h_ecalIso[cate][_bkg_sum][type]->Add(h_ecalIso[cate][file][type]);
	  h_ecalIsoSB[cate][_bkg_sum][type]->Add(h_ecalIsoSB[cate][file][type]);
	  h_ecalIso_et[cate][_bkg_sum][type]->Add(h_ecalIso_et[cate][file][type]);

	  // added by Eiko
	  h_hcalIso_et[cate][_bkg_sum][type]->Add(h_hcalIso_et[cate][file][type]);
	  h_trkIso_et[cate][_bkg_sum][type]->Add(h_trkIso_et[cate][file][type]);
	  h_SIEIE_et[cate][_bkg_sum][type]->Add(h_SIEIE_et[cate][file][type]);
	  h_r9_et[cate][_bkg_sum][type]->Add(h_r9_et[cate][file][type]);
	  h_HoverE_et[cate][_bkg_sum][type]->Add(h_HoverE_et[cate][file][type]);
	  h_ESRatio_et[cate][_bkg_sum][type]->Add(h_ESRatio_et[cate][file][type]);

	  h_ecalIsoSB_et[cate][_bkg_sum][type]->Add(h_ecalIsoSB_et[cate][file][type]);
	  h_hcalIso[cate][_bkg_sum][type]->Add(h_hcalIso[cate][file][type]);
	  h_comb3Iso[cate][_bkg_sum][type]->Add(h_comb3Iso[cate][file][type]);
	  h_comb3IsoSB[cate][_bkg_sum][type]->Add(h_comb3IsoSB[cate][file][type]);
	  h_comb3IsoSB_et[cate][_bkg_sum][type]->Add(h_comb3IsoSB_et[cate][file][type]);
	  h_comb3Iso_et[cate][_bkg_sum][type]->Add(h_comb3Iso_et[cate][file][type]);
	  h_trkIso [cate][_bkg_sum][type]->Add(h_trkIso [cate][file][type]);
	  h_ecalIsoEt[cate][_bkg_sum][type]->Add(h_ecalIsoEt[cate][file][type]);
	  h_hcalIsoEt[cate][_bkg_sum][type]->Add(h_hcalIsoEt[cate][file][type]);
	  h_comb3IsoEt[cate][_bkg_sum][type]->Add(h_comb3IsoEt[cate][file][type]);
	  h_trkIsoEt [cate][_bkg_sum][type]->Add(h_trkIsoEt [cate][file][type]);
	  h_HoverE [cate][_bkg_sum][type]->Add(h_HoverE [cate][file][type]);
	  h_SEE    [cate][_bkg_sum][type]->Add(h_SEE    [cate][file][type]);
	  h_SIEIE  [cate][_bkg_sum][type]->Add(h_SIEIE  [cate][file][type]);
	  h_ESRatio[cate][_bkg_sum][type]->Add(h_ESRatio[cate][file][type]);
	  h_HoverE_noIso [cate][_bkg_sum][type]->Add(h_HoverE_noIso [cate][file][type]);
	  h_SIEIE_noIso  [cate][_bkg_sum][type]->Add(h_SIEIE_noIso  [cate][file][type]);
	  h_ESRatio_noIso[cate][_bkg_sum][type]->Add(h_ESRatio_noIso[cate][file][type]);

	  evtcount[cate][_bkg_sum][type][0] += evtcount[cate][file][type][0];
	  evtcount[cate][_bkg_sum][type][1] =  sqrt(pow(evtcount[cate][_bkg_sum][type][1],2) + pow(evtcount[cate][file][type][1],2));
	  }
	  sigcount[_bkg_sum][type][0] += sigcount[file][type][0];
	  sigcount[_bkg_sum][type][1] =  sqrt(pow(sigcount[_bkg_sum][type][1],2) + pow(sigcount[file][type][1],2));
	}
    }
    
    //dump the event counting:    
    printf("===========================================================================================\n");	 
    printf("Event count summary (EB|EE|Sig. Sum):\n");
    printf("===========================================================================================\n");	 
    for(int type=0; type<SB_TYPE; type++) {	  
      for(int file=0; file<N_SAMPLES; file++) {
    	printf("%10s ",sample[file].tag);
    	printf("| %10.3f +-%6.3f ",evtcount[_EB][file][type][0],evtcount[_EB][file][type][1]);
    	printf("| %10.3f +-%6.3f ",evtcount[_EE][file][type][0],evtcount[_EE][file][type][1]);
    	printf("| %10.3f +-%6.3f\n",sigcount[file][type][0],sigcount[file][type][1]);
	
      }
    }
    printf("===========================================================================================\n");	     

    for(int file=_PhoJet15; file<=_PhoJet80; file++) {
        for(int type=0; type<SB_TYPE; type++) {	  
	double S  = sigcount[file][type][0];
	double dS = sigcount[file][type][1];
	double N  = sigcount[_bkg_sum][type][0];
	double dN = N*0.2;
	double fitness = 2.*(sqrt(S+N) - sqrt(N+dN))*sqrt((N+dN)/(N+N+dN));

	printf("S: %20s: %10.3f +-%6.3f, S/N = %6.3f, S/sqrt(S+N) = %6.3f (signif: %6.3f)\n",
	sample[file].tag,S,dS,
	S/N,S/sqrt(S+N),fitness);
	}	
    }
    		
    fout->cd();
    fout->Write();

    TH1F* hcombIsoSpike = (TH1F*)h_comb3Iso[0][7][0]->Clone();
    TH1F* hecalIsoSpike = (TH1F*)h_ecalIso[0][7][0]->Clone();
    TH1F* hhcalIsoSpike = (TH1F*)h_hcalIso[0][7][0]->Clone();
    TH1F* htrkIsoSpike  = (TH1F*)h_trkIso[0][7][0]->Clone();

    TFile *fout2 = new TFile("spike.root","recreate");
    
    hcombIsoSpike ->Write() ;
    hecalIsoSpike ->Write() ;
    hhcalIsoSpike ->Write() ;
    htrkIsoSpike  ->Write() ;

    fout2->Close();
}
