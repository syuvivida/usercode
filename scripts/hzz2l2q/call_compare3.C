#include "compare3.C"
void call_compare3()
{
  TStopwatch* myWatch = new TStopwatch();
  myWatch->Start();
  
  compare3("h_mh_parton_mother");
  compare3("h_mh_parton_mother",true);
  compare3("h_mll_parton_daughter",false,70,110);
  compare3("h_mjj_parton_daughter",false,70,110);


  compare3("h_mjj_parton_event0",false,60,130);
  compare3("h_mjj_parton_event1",false,60,130);
  compare3("h_mjj_parton_event2",false,60,130);
  compare3("h_mjj_parton_event3",false,60,130);
  compare3("h_mjj_parton_event4",false,60,130);
  compare3("h_mjj_parton_event5",false,60,130);

  compare3("h_mjj_stable_event0",false,60,130);
  compare3("h_mjj_stable_event1");
  compare3("h_mjj_stable_event2");
  compare3("h_mjj_stable_event3");
  compare3("h_mjj_stable_event4");
  compare3("h_mjj_stable_event5");

  compare3("h_mjj_parton_random0",false,60,130);
  compare3("h_mjj_parton_random1",false,60,130);
  compare3("h_mjj_parton_random2",false,60,130);
  compare3("h_mjj_parton_random3",false,60,130);
  compare3("h_mjj_parton_random4",false,60,130);
  compare3("h_mjj_parton_random5",false,60,130);

  compare3("h_mll_stable_random0",false,60,130);

  compare3("h_mjj_stable_random0",false,60,130);
  compare3("h_mjj_stable_random1");
  compare3("h_mjj_stable_random2");
  compare3("h_mjj_stable_random3");
  compare3("h_mjj_stable_random4");
  compare3("h_mjj_stable_random5");

  compare3("h_mjj_rec_truth0",false,70,110);
  compare3("h_mjj_rec_truth1",false,70,110);
  compare3("h_mjj_rec_truth2",false,70,110);
  compare3("h_mjj_rec_truth3",false,70,110);
  compare3("h_mjj_rec_truth4",false,70,110);
  compare3("h_mjj_rec_truth5",false,70,110);

  compare3("h_mll_rec_random0",false,60,130);

  compare3("h_mjj_rec_random0",false,60,130);
  compare3("h_mjj_rec_random1",false,60,130);
  compare3("h_mjj_rec_random2",false);
  compare3("h_mjj_rec_random3",false,60,130);
  compare3("h_mjj_rec_random4",false,60,130);
  compare3("h_mjj_rec_random5",false,60,130);



  myWatch->Stop();
  cout << myWatch->RealTime() << " seconds has passed! " << endl; 

}
