//
//  main.cpp
//  CE 33 Project 2 Curve-Fitting
//  Bontogon & Violeta 2019
//

#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

//FUNCTION DECLARATIONS:
long double *GetSk(int numpts, long double x[], long double y[]);    
long double *GetDif(int numpts, long double array[]);
long double GetSumXDifSq(int numpts, long double xdif[]);
long double GetSumDifSk(int numpts, long double dif[],
  long double sk[]);
long double GetSumSkSq(int numpts, long double sk[]);
long double GetSumYDifXDif(int numpts, long double xdif[],
  long double ydif[]);
long double acquirematvalue(long double secondMatrixA[2][2],
  long double secondMatrixB[2][1], int matnum);
long double acquireThetaSum(long double xArray[], long double c,
  int numpoints);
long double acquireThetaSquaredSum(long double xArray[], long double c,
  int numpoints);
long double acquireYSum(long double yArray[], int numpoints);
long double acquireYThetaSum(long double xArray[],
  long double yArray[], long double c, int numpoints);
long double residual(int numpts, long double y[], long double x[],
  long double a, long double b, long double c);

int main(int argc, const char * argv[]) {
     
  //CHECK FOR AVAILABILITY OF FILE AND OPEN IT
  string CFFilename;
  cout << "Filename: ";
  cin >> CFFilename;
  fstream CFFileStr(CFFilename.c_str(), ios::in);
    
  //CHECK IF FILE IS VALID===========================================
  if(CFFileStr.fail())
  {
    cerr << "Unable to open " << CFFilename << endl;
    system("pause");
    exit(0);
  }

  //CHECK IF FILE IS EMPTY================================================
  if (CFFileStr.peek() == EOF)
  {
    cout << "The file is empty. Exiting..." << endl;
    system("pause");
    return (0);
  }
    
  //CHECK IF FILE HAS INTEGERS IN IT========================================
  int numpoints = 0;
  long double* xArray = NULL;
  long double* yArray = NULL;
  if (isdigit(CFFileStr.peek()))
  {
    CFFileStr >> numpoints;
    xArray = new long double[numpoints-1];
    yArray = new long double[numpoints-1];
    cout << "Number of points: " << numpoints << endl;
        
    for (int xory=0; xory<=(2*numpoints); xory++)
    {
      if (xory%2==0)
      {
        CFFileStr >> xArray[xory/2];
      }
      if (xory%2==1)
      {
        CFFileStr >> yArray[xory/2];
      }          
    }
  }

  //CHECK IF FILE HAS NON-NUMERIC CHARACTERS IN IT=======================
  else
  {
    cerr << "You inputted an invalid file. Enter a .txt file containing numbers only." << endl;
    system("pause");
    exit(0);
  }
    
  //DISPLAY X and Y Values
  cout<<"Pt. No."<<setw(5)<<"x"<<setw(19)<<"y"<<endl;
  cout<<"-----------------------------------------------------------------\n";
  for (int i=0;i<numpoints;i++)
  cout<<i+1<<"."<<setw(10)<<xArray[i]<<setw(19)<<yArray[i]<<endl;
    
  //SORT X VALUES FROM LEAST TO GREATEST
  for(int i=0;i<numpoints;i++)
  {
    for(int j=0;j<numpoints-i;j++)
    {
      long double a = xArray[j];
      long double b = xArray[j+1];
      if(xArray[j]>xArray[j+1])
      {
        long double tempx, tempy;
        tempx=xArray[j];
        tempy=yArray[j];
        xArray[j]=xArray[j+1];
        yArray[j]=yArray[j+1];
        xArray[j+1]=tempx;
        yArray[j+1]=tempy;
      }
    }
  }
    
  //Get Inputs for Matrix Values
  long double *Sk = GetSk(numpoints,xArray,yArray);
  long double *XDif = GetDif(numpoints,xArray);
  long double *YDif = GetDif(numpoints,yArray);
    
  /*//DISPLAY SORTED X and Y Values
  cout<<"Pt. No."<<setw(5)<<"x"<<setw(8)<<"y"<<setw(8)<<endl;
  cout<<"----------------------------------------------------------------------------------------------------------\n";
  for (int i=0;i<numpoints;i++)
  {
    cout<<i+1<<"."<<setw(10)<<xArray[i]<<setw(8)<<yArray[i]<<endl;
  }

  cout<<"Sum of Sk^2: " << SumSkSq << endl;
  cout<<"Sum of xdif^2: " << SumXDifSq << endl;
  cout<<"Sum of xdifsk: " << SumXDifSk << endl;
  cout<<"Sum of ydifsk: " << SumYDifSk << endl;
  cout<<"Sum of ydifxdif: " << SumYDifXDif << endl;
  */

  long double firstMatrixA[2][2], firstMatrixB[2][1];
  long double secondMatrixA[2][2], secondMatrixB[2][1];    
  firstMatrixA[0][0] = GetSumXDifSq(numpoints, XDif);
  firstMatrixA[0][1] = GetSumDifSk(numpoints, XDif, Sk);
  firstMatrixA[1][0] = GetSumDifSk(numpoints, XDif, Sk);
  firstMatrixA[1][1] = GetSumSkSq(numpoints, Sk);
  firstMatrixB[0][0] = GetSumYDifXDif(numpoints, XDif, YDif);
  firstMatrixB[1][0] = GetSumDifSk(numpoints, YDif, Sk);

  long double c = acquirematvalue(firstMatrixA, firstMatrixB, 1);
  secondMatrixA[0][0] = numpoints;
  secondMatrixA[0][1] = acquireThetaSum(xArray, c, numpoints);
  secondMatrixA[1][0] = acquireThetaSum(xArray, c, numpoints);
  secondMatrixA[1][1] = acquireThetaSquaredSum(xArray, c, numpoints);
  secondMatrixB[0][0] = acquireYSum(yArray, numpoints);
  secondMatrixB[1][0] = acquireYThetaSum(xArray, yArray, c, numpoints);

  /*
  cout << "c is equal to: " << c << endl;
  cout << "theta is equal to: " << theta << endl;
  cout << "theta squared is equal to: " << thetaSquared << endl;
  cout << "y sum is equal to: " << ySum << endl;
  cout << "y theta sum  is equal to: " << yThetaSum << endl;
  */
  long double a = acquirematvalue(secondMatrixA, secondMatrixB, 0);
  long double b = acquirematvalue(secondMatrixA, secondMatrixB, 1);
    
  cout << "The best fit equation for the set of data is given as:\n";
  cout << "y = " << a << "+" << b << "e^ (" << c << "x)" << endl;  
  long double res = residual(numpoints, yArray, xArray, a, b, c);
  cout << "S_r is: " << res << endl;
    
  system("pause");
  return 0;
}


long double acquirematvalue(long double firstMatrixA[2][2],
  long double firstMatrixB[2][1], int matnum)
{
  long double det, denom;
  denom = (firstMatrixA[0][0]*firstMatrixA[1][1])
    - (firstMatrixA[0][1]*firstMatrixA[1][0]);
  det = pow(denom, -1);
    
  long double inverseMatrix[2][2];
  inverseMatrix[0][0] = firstMatrixA[1][1];
  inverseMatrix[0][1] = firstMatrixA[0][1]*(-1);
  inverseMatrix[1][0] = firstMatrixA[1][0]*(-1);
  inverseMatrix[1][1] = firstMatrixA[0][0];
  
  //multiply by the determinant
  for(int i=0; i<=1; i++){
    for(int j=0; j<=1; j++){
      inverseMatrix[i][j] = inverseMatrix[i][j]*det;
    }
  }
    
  long double cMatrix[2][1];
  cMatrix[0][0] = (inverseMatrix[0][0]*firstMatrixB[0][0])
    +(inverseMatrix[0][1]*firstMatrixB[1][0]);
  cMatrix[1][0] = (inverseMatrix[1][0]*firstMatrixB[0][0])
    +(inverseMatrix[1][1]*firstMatrixB[1][0]);
  return cMatrix[matnum][0];
}

long double acquireThetaSum(long double xArray[], long double c,
  int numpoints)
{
  long double theta = 0;
  for(int i=0; i<numpoints; i++)
  {
    theta += exp(c*xArray[i]);
  }
  return theta;
}

long double acquireThetaSquaredSum(long double xArray[], long double c,
  int numpoints)
{
  long double thetaSquared = 0;
  for(int i=0; i<numpoints; i++)
  {
    thetaSquared += pow(exp(c*xArray[i]),2);
  }
  return thetaSquared;
}

long double acquireYSum(long double yArray[], int numpoints)
{
  long double ySum = 0;
  for(int i=0; i<numpoints; i++)
  {
    ySum += yArray[i];
  }
  return ySum;
}

long double acquireYThetaSum(long double xArray[], long double yArray[],
  long double c, int numpoints)
{
  long double yThetaSum = 0;
  long double theta = 0;
  for(int i=0; i<numpoints; i++)
  {
    theta = exp(c*xArray[i]);
    yThetaSum += yArray[i]*theta;
  }
  return yThetaSum;  
}

long double *GetSk(int numpts, long double x[], long double y[])
{
  long double* sk = NULL;
  sk = new long double[numpts];
  sk[0]=0;
  for(int i=1;i<=numpts;i++)
  {
    sk[i]=sk[i-1]+((0.5)*(y[i]+y[i-1])*(x[i]-x[i-1]));
  }
  return sk;
}

long double *GetDif(int numpts, long double array[])
{
  long double* difarr = NULL;
  difarr = new long double[numpts];
  for(int i=0;i<=numpts;i++)
  {
    difarr[i]=array[i]-array[0];
  }
  return difarr;
}

long double GetSumXDifSq(int numpts, long double xdif[])
{
  long double sumxdifsq=0;
  for(int i=0;i<numpts;i++)
  {
    sumxdifsq+=pow(xdif[i],2);
  }
  return sumxdifsq;
}

long double GetSumDifSk(int numpts,long double dif[],long double sk[])
{
  long double sumdifsk = 0;
  for(int i=0; i<numpts;i++)
  {
    sumdifsk+=dif[i]*sk[i];
  }
  return sumdifsk;
}

long double GetSumSkSq(int numpts, long double sk[])
{
  long double sumsksq=0;
  for(int i=0;i<numpts;i++)
  {
    sumsksq+=pow(sk[i],2);
  }
  return sumsksq;
}

long double GetSumYDifXDif(int numpts, long double xdif[],
  long double ydif[])
{
  long double sumydifxdif=0;
  for(int i=0; i<numpts;i++)
  {
    sumydifxdif+=xdif[i]*ydif[i];
  }
  return sumydifxdif;
}

long double residual(int numpts, long double y[], long double x[],
  long double a, long double b, long double c)
{
  long double Residual = 0;
  for(int i=0; i<numpts;i++)
  {
    Residual+=pow((y[i]-(a+b*exp(c*x[i]))),2);
  }
  return Residual;  
}    
