{

  std::string pythiaF = "egamma_pthat0_zhijun_masshisto.root";
  std::string wgnloF  = "wgamma_nlo_wpg_ATLAS_CTEQ5L.root";
  std::string wgF     = "wpgamma_ATLAS_baur.root";

 gROOT->ProcessLine(".L overLay.C");


 overLay("mw_ATLAS",pythiaF,"h_mW2",wgnloF,"h300",
	 wgF,"h_mW2","M(W) [GeV]",-1,-1,false,1.1);

 overLay("mwg_ATLAS",pythiaF,"h_mWg2",wgnloF,"h23",
	 wgF,"h_mWg2","M(W#gamma) [GeV]",50,250,false,0.35);

 overLay("ctm_ATLAS",pythiaF,"hmT3",wgnloF,"h24",
	 wgF,"hmT3","Cluster mass [GeV]",-1,-1,false,0.1);

 overLay("gpt_ATLAS",pythiaF,"h_gpta",wgnloF,"h12",
	 wgF,"h_gpta","p_{T}(#gamma) [GeV]",0,150,true);

 overLay("ept_ATLAS",pythiaF,"h_lpta",wgnloF,"h13",
	 wgF,"h_lpta","p_{T}(electron) [GeV]",-1,-1,true);

 overLay("npt_ATLAS",pythiaF,"h_npta",wgnloF,"h14",
	 wgF,"h_npta","p_{T}(#nu) [GeV]",0,200,true);

 
 overLay("geta_ATLAS",pythiaF,"h_getaa",wgnloF,"h17",
	 wgF,"h_getaa","#eta(#gamma)",-1,-1,true,0.2);

 overLay("leta_ATLAS",pythiaF,"h_letaa",wgnloF,"h16",
	 wgF,"h_letaa","#eta(electron)",-1,-1,true,0.2);
 
 overLay("dR_ATLAS",pythiaF,"h_dR2",wgnloF,"h27",
	 wgF,"h_dR2","#DeltaR(e-#gamma)",-1,-1,true,10   );
 



}
