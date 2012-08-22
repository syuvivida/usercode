#include "compare900GeV_random.C"
void call_compare900GeV_random()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
 
  compare900GeV_random("h_mjj_stable_event");
  compare900GeV_random("h_mjj_parton_event");

  compare900GeV_random("h_mh_stable_event");
  compare900GeV_random("h_mh_parton_event");

  compare900GeV_random("h_mjj_stable_random");
  compare900GeV_random("h_mjj_parton_random");

  compare900GeV_random("h_mh_stable_random");
  compare900GeV_random("h_mh_parton_random");


  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
