{
const int NFILES=3;
float npass[NFILES];
float ntotal[NFILES];
float npass_real[NFILES];
float ntotal_real[NFILES];
float xsec[NFILES];
float filt[NFILES];
float ngen[NFILES];
ifstream fin;

fin.open("event.dat");
for(int i=0; i < NFILES;i++)
  fin >> ntotal[i] >> npass[i] >> ntotal_real[i] >> npass_real[i];
fin.close();


fin.clear();
fin.open("xsec.dat");
for(int i=0; i< NFILES; i++)
  fin >> xsec[i] >> filt[i] >> ngen[i];

float signal_pass = 0;
float signal_total= 0;
float signal_pass_real = 0;
float signal_total_real= 0;


for(int i=0; i<3;i++)
{
  signal_total += xsec[i]*filt[i]*ntotal[i]/ngen[i];
  signal_pass += xsec[i]*filt[i]*npass[i]/ngen[i];

  signal_total_real += xsec[i]*filt[i]*ntotal_real[i]/ngen[i];
  signal_pass_real  += xsec[i]*filt[i]*npass_real[i]/ngen[i];
}


float eff = (signal_pass)/(signal_total);
float eff_real = (signal_pass_real)/(signal_total_real);

 cout << "eff before matching = " << eff << endl;
 cout << "eff after matching = " << eff_real << endl;


}
