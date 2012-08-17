#include "compareDressed_darko.C"

void call_compareDressed_darko(std::string File){

  compareDressed_darko(File.data(),"h_ystar");
  compareDressed_darko(File.data(),"h_yB");
  compareDressed_darko(File.data(),"h_zy");
  compareDressed_darko(File.data(),"h_jety");
  compareDressed_darko(File.data(),"h_zpt");
  compareDressed_darko(File.data(),"h_jetpt");
  compareDressed_darko(File.data(),"h_mZ",70,110,0.5,2.0);


}
