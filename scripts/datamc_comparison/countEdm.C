{

  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  const int NFILES=96;
  std::string a[NFILES];
  a[0]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_88.root";
  a[1]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_89.root";
  a[2]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_90.root";
  a[3]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_91.root";
  a[4]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_92.root";
  a[5]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_93.root";
  a[6]="BSCFilter_Jan29_v8_MinBiasFilter_Run123596_94.root";
  a[7]="BSCFilter_Jan29_v8_MinBiasFilter_Run123615_86.root";
  a[8]="BSCFilter_Jan29_v8_MinBiasFilter_Run123615_87.root";
  a[9]="BSCFilter_Jan29_v8_MinBiasFilter_Run123732_81.root";
  a[10]="BSCFilter_Jan29_v8_MinBiasFilter_Run123732_82.root";
  a[11]="BSCFilter_Jan29_v8_MinBiasFilter_Run123732_83.root";
  a[12]="BSCFilter_Jan29_v8_MinBiasFilter_Run123732_84.root";
  a[13]="BSCFilter_Jan29_v8_MinBiasFilter_Run123732_85.root";
  a[14]="BSCFilter_Jan29_v8_MinBiasFilter_Run123815_80.root";
  a[15]="BSCFilter_Jan29_v8_MinBiasFilter_Run123818_77.root";
  a[16]="BSCFilter_Jan29_v8_MinBiasFilter_Run123818_78.root";
  a[17]="BSCFilter_Jan29_v8_MinBiasFilter_Run123818_79.root";
  a[18]="BSCFilter_Jan29_v8_MinBiasFilter_Run123906_76.root";
  a[19]="BSCFilter_Jan29_v8_MinBiasFilter_Run123908_75.root";
  a[20]="BSCFilter_Jan29_v8_MinBiasFilter_Run123909_73.root";
  a[21]="BSCFilter_Jan29_v8_MinBiasFilter_Run123909_74.root";
  a[22]="BSCFilter_Jan29_v8_MinBiasFilter_Run123970_71.root";
  a[23]="BSCFilter_Jan29_v8_MinBiasFilter_Run123970_72.root";
  a[24]="BSCFilter_Jan29_v8_MinBiasFilter_Run123977_66.root";
  a[25]="BSCFilter_Jan29_v8_MinBiasFilter_Run123977_67.root";
  a[26]="BSCFilter_Jan29_v8_MinBiasFilter_Run123977_68.root";
  a[27]="BSCFilter_Jan29_v8_MinBiasFilter_Run123977_69.root";
  a[28]="BSCFilter_Jan29_v8_MinBiasFilter_Run123977_70.root";
  a[29]="BSCFilter_Jan29_v8_MinBiasFilter_Run123978_63.root";
  a[30]="BSCFilter_Jan29_v8_MinBiasFilter_Run123978_64.root";
  a[31]="BSCFilter_Jan29_v8_MinBiasFilter_Run123978_65.root";
  a[32]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_47.root";
  a[33]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_48.root";
  a[34]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_49.root";
  a[35]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_50.root";
  a[36]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_51.root";
  a[37]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_52.root";
  a[38]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_53.root";
  a[39]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_54.root";
  a[40]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_55.root";
  a[41]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_56.root";
  a[42]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_57.root";
  a[43]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_58.root";
  a[44]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_59.root";
  a[45]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_60.root";
  a[46]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_61.root";
  a[47]="BSCFilter_Jan29_v8_MinBiasFilter_Run123987_62.root";
  a[48]="BSCFilter_Jan29_v8_MinBiasFilter_Run124006_46.root";
  a[49]="BSCFilter_Jan29_v8_MinBiasFilter_Run124008_45.root";
  a[50]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_35.root";
  a[51]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_36.root";
  a[52]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_37.root";
  a[53]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_38.root";
  a[54]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_39.root";
  a[55]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_40.root";
  a[56]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_41.root";
  a[57]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_42.root";
  a[58]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_43.root";
  a[59]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_44.root";
  a[60]="BSCFilter_Jan29_v8_MinBiasFilter_Run124009_99.root";
  a[61]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_29.root";
  a[62]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_30.root";
  a[63]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_31.root";
  a[64]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_32.root";
  a[65]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_33.root";
  a[66]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_34.root";
  a[67]="BSCFilter_Jan29_v8_MinBiasFilter_Run124020_95.root";
  a[68]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_21.root";
  a[69]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_22.root";
  a[70]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_23.root";
  a[71]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_24.root";
  a[72]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_25.root";
  a[73]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_26.root";
  a[74]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_27.root";
  a[75]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_28.root";
  a[76]="BSCFilter_Jan29_v8_MinBiasFilter_Run124022_96.root";
  a[77]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_16.root";
  a[78]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_17.root";
  a[79]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_18.root";
  a[80]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_19.root";
  a[81]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_20.root";
  a[82]="BSCFilter_Jan29_v8_MinBiasFilter_Run124023_97.root";
  a[83]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_10.root";
  a[84]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_11.root";
  a[85]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_12.root";
  a[86]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_13.root";
  a[87]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_14.root";
  a[88]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_15.root";
  a[89]="BSCFilter_Jan29_v8_MinBiasFilter_Run124024_98.root";
  a[90]="BSCFilter_Jan29_v8_MinBiasFilter_Run124025_9.root";
  a[91]="BSCFilter_Jan29_v8_MinBiasFilter_Run124027_7.root";
  a[92]="BSCFilter_Jan29_v8_MinBiasFilter_Run124027_8.root";
  a[93]="BSCFilter_Jan29_v8_MinBiasFilter_Run124030_4.root";
  a[94]="BSCFilter_Jan29_v8_MinBiasFilter_Run124030_5.root";
  a[95]="BSCFilter_Jan29_v8_MinBiasFilter_Run124030_6.root";
//   a[96]="BSCFilter_Jan29_v8_MinBiasFilter_Run124120_2.root";
//   a[97]="BSCFilter_Jan29_v8_MinBiasFilter_Run124120_3.root";
//   a[98]="BSCFilter_Jan29_v8_MinBiasFilter_Run124275_1.root";


 int runnumber[]={
   123596,
   123615,
   123818,
   123906,
   124006,
   124008,
   124009,
   124020,
   124022,
   124023,
   124024,
   124025,
   124027,
   124030};

 const int nruns = sizeof(runnumber)/sizeof(runnumber[0]);
 cout << "There are " << nruns << " runs in total" <<  endl;
 Long64_t total_event = 0;
 Long64_t run_event[nruns];
 for(int i=0; i< nruns; i++)run_event[i] =0;
 TFile *temp;  
 for(int ifile=0; ifile< NFILES; ifile++)
   {
     temp= TFile::Open(a[ifile].data());
     cout << "File" << a[ifile].data() << " opend " << endl;
     int nevent = 
       Events->GetEntries("EventAuxiliary.run()");
     total_event += nevent;
     for(int irun = 0; irun < nruns; irun++)
       {
 	 char cutname[4000];
 	 sprintf(cutname,"EventAuxiliary.run()==%d",runnumber[irun]);
 	 cout << cutname << endl;
 	 int runEve = 
 	   Events->Draw("EventAuxiliary.run()",cutname);
	 run_event[irun] += runEve;
       } // finishing looping over runs
    
   } // finishing looping over files

 cout << "total_event = " << total_event << endl;
 for(int irun=0; irun < nruns; irun++)
   cout << "Run " << runnumber[irun] << " : " << run_event[irun] << endl;
}
