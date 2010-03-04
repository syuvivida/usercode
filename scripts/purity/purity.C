{
const int NFILES=5;
float npass[NFILES];
float ntotal[NFILES];
float xsec[NFILES];
float filt[NFILES];
float ngen[NFILES];
ifstream fin;

fin.open("event.dat");
for(int i=0; i < NFILES;i++)
  fin >> ntotal[i] >> npass[i];
fin.close();


fin.clear();
fin.open("xsec.dat");
for(int i=0; i< NFILES; i++)
  fin >> xsec[i] >> filt[i] >> ngen[i];

float signal_pass = 0;
float signal_total= 0;

for(int i=0; i<3;i++)
{
  signal_total += xsec[i]*filt[i]*ntotal[i]/ngen[i];
  signal_pass += xsec[i]*filt[i]*npass[i]/ngen[i];
}

float background_pass=0;
float background_total=0;

for(int i=3; i<5;i++)
{
  background_total += xsec[i]*filt[i]*ntotal[i]/ngen[i];
  background_pass += xsec[i]*filt[i]*npass[i]/ngen[i];
}

float purity = (signal_pass+background_pass)/(signal_total+background_total);

cout << "purity after cut = " << purity << endl;


}
