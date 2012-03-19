#include "compareDataMC.C"

void call_compareDataMC(std::string dataFile, std::string mcFile){



  compareDataMC(dataFile.data(),mcFile.data(),"h_zmass_ID",76,106);


  for(int ip=0; ip < 5; ip++){
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_cost_COM3D_%d",ip));
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_cost_COMZ_%d",ip));
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_cost_sumJet_%d",ip));

    compareDataMC(dataFile.data(),mcFile.data(),Form("h_ystar_%d",ip),-2.0,2.0);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_ystar_COM3D_%d",ip),-2.0,2.0);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_ystar_COMZ_%d",ip),-2.0,2.0);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_ystar_sumJet_%d",ip),-2.0,2.0);

    compareDataMC(dataFile.data(),mcFile.data(),Form("h_zpt_%d",ip),0,200);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_zy_%d",ip),-2.0,2.0);

    compareDataMC(dataFile.data(),mcFile.data(),Form("h_leadingjetpt_%d",ip),30,200);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_leadingjety_%d",ip),-2.4,2.4);

    compareDataMC(dataFile.data(),mcFile.data(),Form("h_sumjetpt_%d",ip),30,200);
    compareDataMC(dataFile.data(),mcFile.data(),Form("h_sumjety_%d",ip),-2.4,2.4);

  }

}
