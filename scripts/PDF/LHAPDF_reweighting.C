#include "LHAPDF/LHAPDF.h"

class MyPDF {
    public:
      MyPDF(const char* pdfname) {
            //LHAPDF::setVerbosity(LHAPDF::SILENT);
            //LHAPDF::initPDFSet(1,"cteq6ll.LHpdf");
            LHAPDF::initPDFSet(1,"cteq6ll.LHgrid");
            // Typical pdfnames to reweight to are:
            //    CT10.LHgrid (this gives 90% CL uncertainties)
            //    MSTW2008nlo68cl.LHgrid, MSTW2008nlo90cl.LHgrid
            //    NNPDF21_100.LHgrid (this gives 68% CL uncertainties)
            LHAPDF::initPDFSet(2,pdfname);
      }
      virtual ~MyPDF(){};

      // This reweights the event to reference to set "pdfname", member "member"
      //
      // For NNPDF2x_100 one should just use all 100 members and calculate 
      // the r.m.s. of the 100 values of the physical observable obtained
      // via reweighting
      //
      // For CTEQ or MSTW, in order to assign systematics one should 
      // use the master formulae. See for instance:
      // http://cmsdoc.cern.ch/cms/PRS/gentools/www/pdfuncert/uncert.html
      double weight(GenInfo& genInfo, unsigned int member) {
            int id1 = genInfo.id1;
            int id2 = genInfo.id2;
            double x1 = genInfo.x1;
            double x2 = genInfo.x2;
            double Q = genInfo.scalePDF;
            LHAPDF::usePDFMember(1,0);
            double pdf1 = LHAPDF::xfx(1,x1,Q,id1)*LHAPDF::xfx(1,x2,Q,id2);
            LHAPDF::usePDFMember(2,member);
            double pdf2 = LHAPDF::xfx(2,x1,Q,id1)*LHAPDF::xfx(2,x2,Q,id2);
            if (pdf1>0) {
                  return pdf2/pdf1;
            } else {
                  //printf ("pdf1 = %e, pdf2 = %e\n", pdf1, pdf2);
                  return 1.;
            }
      }
};

int main(int argc, char** argv){

  MyPDF* pdf = new MyPDF("MSTW2008nlo68cl.LHgrid"); 
  //MyPDF* pdf = new MyPDF("CT10.LHgrid"); 
  //MyPDF* pdf = new MyPDF("NNPDF21_100.LHgrid"); 

  // Loop on events for the current sample
  for (int iEvent=0; iEvent<nevents; iEvent++) {

      // Calculate event pdf weight (this depends on your Root structure)
      GenInfo geninfo = ....YourGenInfoInThisEvent();
      int pdfmember = 0; // 0 is the central PDF member of the set
      double weight = pdf->weight(geninfo,pdfmember); // weight for pdfmember
      ...
  }

  ...
}
