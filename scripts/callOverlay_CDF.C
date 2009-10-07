{

  std::string pythiaF = "egamma_pthat0_CDF_masshisto.root";
  std::string wgnloF  = "wgamma_nlo_wpg_CDF_CTEQ5L.root";
  std::string wgF     = "wpg_CDF_baur.root";

 gROOT->ProcessLine(".L overLay.C");


 overLay("mw_CDF",pythiaF,"h_mW2",wgnloF,"h300",
	 wgF,"h_mW2","M(W) [GeV]",-1,-1,false,1.1);

 overLay("mwg_CDF",pythiaF,"h_mWg2",wgnloF,"h23",
	 wgF,"h_mWg2","M(W#gamma) [GeV]",50,250,false,0.35);

 overLay("ctm_CDF",pythiaF,"hmT3",wgnloF,"h24",
	 wgF,"hmT3","Cluster mass [GeV]",-1,-1,false,0.1);

 overLay("gpt_CDF",pythiaF,"h_gpta",wgnloF,"h12",
	 wgF,"h_gpta","p_{T}(#gamma) [GeV]",0,150,true);

 overLay("ept_CDF",pythiaF,"h_lpta",wgnloF,"h13",
	 wgF,"h_lpta","p_{T}(electron) [GeV]",-1,-1,true);

 overLay("npt_CDF",pythiaF,"h_npta",wgnloF,"h14",
	 wgF,"h_npta","p_{T}(#nu) [GeV]",0,200,true);

 
 overLay("geta_CDF",pythiaF,"h_getaa",wgnloF,"h17",
	 wgF,"h_getaa","#eta(#gamma)",-1,-1,true,100);

 overLay("leta_CDF",pythiaF,"h_letaa",wgnloF,"h16",
	 wgF,"h_letaa","#eta(electron)",-1,-1,true,100);
 
 overLay("dR_CDF",pythiaF,"h_dR2",wgnloF,"h27",
	 wgF,"h_dR2","#DeltaR(e-#gamma)",-1,-1,true,10   );
 



}
