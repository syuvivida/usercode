// c++ -o addBCMassInLHE `root-config --glibs --cflags` -lm addBCMassInLHE.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"

using namespace std ;


int main (int argc, char ** argv)
{
  if(argc < 3)
    {
      cout << "Usage: " << argv[0]
           << " input.lhe output.lhe" << endl ;
      return -1;
    }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  ofstream outputStream (argv[2]) ;
  LHEF::Writer writer (outputStream) ;

  writer.headerBlock () << reader.headerBlock ;
  writer.initComments () << reader.initComments ;
  writer.heprup = reader.heprup ;
  writer.init () ;

//PG mu massless in phantom
  const double mb = 5.0;
  const double mc = 1.5;
  const double mb2 = mb*mb;
  const double mc2 = mc*mc;

  int count = 0 ;
  //PG loop over input events
  while (reader.readEvent ())
    {
      ++count ;
      if ( reader.outsideBlock.length ()) std::cout << reader.outsideBlock;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
           // outgoing particles
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
	       int PID = abs (reader.hepeup.IDUP.at (iPart));
               if (PID == 4 || PID == 5)
                 {
                   TLorentzVector dummy
                     (
                       reader.hepeup.PUP.at (iPart).at (0), // px
                       reader.hepeup.PUP.at (iPart).at (1), // py
                       reader.hepeup.PUP.at (iPart).at (2), // pz
                       reader.hepeup.PUP.at (iPart).at (3) // E
                     ) ;
		   double m2    = pow(reader.hepeup.PUP.at (iPart).at (4),2);
                   double pVec2 = dummy.Vect ().Mag2 () ;

		   double k2 = PID==4? (m2-mc2): (m2-mb2);
                   double scale = sqrt (1 + k2 / pVec2) ;
                   if (pVec2 < (-1 * k2))
                     {
                       cout << "warning: pVec2 is smaller than the mass difference " << pVec2 << endl ;
                       scale = 1 ;
                     }
                   reader.hepeup.PUP.at (iPart).at (0) *= scale ; // px
                   reader.hepeup.PUP.at (iPart).at (1) *= scale ; // py
                   reader.hepeup.PUP.at (iPart).at (2) *= scale ; // pz
		   reader.hepeup.PUP.at (iPart).at (4)  = PID==4? mc: mb;
                 }

             } // outgoing particles
        } // loop over particles in the event
      writer.eventComments () << reader.eventComments ;
      writer.hepeup = reader.hepeup ;
      bool written = writer.writeEvent () ;
      if (!written)
        {
          cout << "warning: event " << count << " not written" << endl ;
        }

    } //PG loop over input events

  cout << "end loop over " << count << " events" << endl ;
  return 0 ;
}

