const double minZPt  = 0.0;
const double minMll = 70.0;
const double maxMll =110.0;

const double minJetPt=30.0;
const double maxJetEta=2.4;
const double mindR=0.5;
const double minLepPt = 20.0;

const double minEleBarrelEta = 0.0;
const double maxEleBarrelEta = 1.442;
const double minEleEndcapEta = 1.566;
const double maxEleEndcapEta = 2.4;

const double maxMuoEta = 2.4;

const double fBinsPt01[]= {15,20,25,30,40,55,75,105,150,210,315,500};
const double fBinsPt02[]= {15,20,25,30,40,55,75,105,150,210,315,450};
const double fBinsPt03[]= {15,20,25,30,40,55,75,105,150,300};
const double fBinsPt04[]= {15,20,25,30,40,55,75,105,200};

const int nPtBins01 = sizeof(fBinsPt01)/sizeof(fBinsPt01[0])-1;
const int nPtBins02 = sizeof(fBinsPt02)/sizeof(fBinsPt02[0])-1;
const int nPtBins03 = sizeof(fBinsPt03)/sizeof(fBinsPt03[0])-1;
const int nPtBins04 = sizeof(fBinsPt04)/sizeof(fBinsPt04[0])-1;

const double fBinsY[]={0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7};
const int nYBins = sizeof(fBinsY)/sizeof(fBinsY[0])-1;

