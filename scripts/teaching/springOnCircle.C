#include <TMatrixD.h>
#include <TMatrixDEigen.h>

void springOnCircle(int ndim)
{
 TMatrixD a(ndim,ndim);
 for(int i=0; i<ndim; i++)
 { 	 
   a(i,i)=-2;
   int left  = (i-1)<0     ? ndim-1: i-1;
   int right = (i+1)>ndim-1?      0: i+1; 	 
   for(int j=0;j<ndim; j++)
     {
       if(j==i)continue;
       else if(j==left) 
	 a(i,left)=1;
       else if(j==right) 
	 a(i,right)=1;
       else
	 a(i,j)=0;
     }
 }
 a.Print(); 
 TMatrixDEigen b(a);
 TMatrixD eigen=b.GetEigenValues();
 eigen.Print();
 TMatrixD eigenV=b.GetEigenVectors();
 eigenV.Print();

}
