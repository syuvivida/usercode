#include "compare3.C"

void runAll()
{

  compare3("h_mh_parton");
  compare3("h_mll_parton");
  compare3("h_mjj_parton");

  compare3("h_mh_stable");
  compare3("h_mll_stable");
  compare3("h_mjj_stable");

  compare3("h_mh_rec");
  compare3("h_mll_rec");
  compare3("h_mjj_rec");

  compare3("h_dR_ll");
  compare3("h_dR_jj");

  const double dRArray[]={0.0,1.0,2.0,3.5};
//   const double dRArray[]={0.5,1.0,1.5,2.0,2.5,3.0,3.5};
  const int NBINS = sizeof(dRArray)/sizeof(dRArray[0])-1;

  std::string jetName[2]={"leading","subleading"};

  for(int j=0; j<2; j++){

    compare3(Form("h_jec_%s",jetName[j].data()));

    for(int i=0; i< NBINS; i++){

      compare3(Form("h_mll_recdR%d",i));
      
      compare3(Form("h_mjj_recdR%d",i));

      compare3(Form("h_jec_dR%d_%s",i,jetName[j].data()));


    } // end of loop over  bins

  } // end of loop over jets


}
