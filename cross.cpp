#ifndef __cross__section__cpp
#define __cross__section__cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <TROOT.h>

using namespace std;

class CS {

public:

  Int_t val = 0;
  std::vector <Double_t> Energy;
  std::vector <Double_t> Cross;
  std::vector <Double_t> cross_dx;

  Double_t maxCrossSection;
  Double_t maxCSEnergy;
 
 
  CS(){;}

  //Int_t ReadFile(const char*);
  //template<class T> Double_t Thickness(T*); 
  //template<class T> Double_t CrossSection(Double_t, T*); 




  Int_t ReadFile(const char* CrossSectionFile) {

    ifstream cross_file;
    string line;
    Int_t line_counter;
   
    cross_file.open(CrossSectionFile);
    if(!cross_file.is_open())
      cout << " cross section file: " << CrossSectionFile << " is NOT open " << endl;
    else {

      cout << " cross section file: " << CrossSectionFile << " is being read " << endl;
      line_counter = 0;

      do {
	getline(cross_file, line);
	if(!line.empty())
	  line_counter++;
      } while(!cross_file.eof());

      cross_file.close();

      val = line_counter-1;
      Energy.resize(val,0);
      Cross.resize(val,0);

      cout << val << endl;

      cross_file.open(CrossSectionFile);
      getline(cross_file,line);
      cout << line << endl;

      for(Int_t i=0; i<val; i++)
	cross_file >> Energy[i] >> Cross[i];

      maxCrossSection = *max_element(Cross.begin(),Cross.end());

      Int_t index = distance(Cross.begin(),max_element(Cross.begin(),Cross.end()));
      maxCSEnergy = Energy[index];

    }

    return 1;

  }


  //function Thickness is made so that it calculates the product of 
  //cross-section and target thickness dx to be used for the Nout probability calculation
  template<class T> 
  Double_t Thickness(T* S)
  {
    ofstream out("cross_thick_product.txt");
    
    Int_t length = Energy.size();
   
    cross_dx.resize(length,0.0);

    vector <Double_t> thickness(length,0.0);
    vector <Double_t> Enlab(length,0.0);
    
    Double_t maxProduct = 0.0;
    
    for(Int_t i=0; i<length; i++)
      {
	Enlab[i] = Energy[i]*22./4;
      }

    
    for(Int_t i=0; i<length-1; i++)
      {
	thickness[i] = S->IonInGas->GetDistance_new(Enlab[i+1],Enlab[i],0.001);
	cross_dx[i] = Cross[i]*thickness[i];
	
	out << " Energy: " << Energy[i] 
	    << " dx: " << thickness[i] << " cross: " << Cross[i] << " product: " << cross_dx[i] << endl;
	
      }
    
    maxProduct = *std::max_element(cross_dx.begin(),cross_dx.end());
    Int_t index = std::distance(cross_dx.begin(),std::max_element(cross_dx.begin(),cross_dx.end()));
    
    //cout << " maxProduct: " << maxProduct << " index: " << index << " Energy: " << Energy[index] << endl;
    
    return maxProduct;
   
  }
 

  /* must put template class definitions here for linkage purposes */
  template<class T> 
  Double_t CrossSection(Double_t En, T* S) {

    Double_t density = (1.1617e-4)*(6.022e23)*0.96/4.0;
    Double_t Nin = 9.18e8;

    Double_t Nout = 0.0;
    Double_t cross = 0.0;
    Double_t crossDx_product = 0.0;

    Double_t maxNorm = Thickness(S);

    En = En*(4./22);
    
    if(En > Energy[Energy.size()-1])
      return 0;
    
    Int_t index = distance(Energy.begin(),lower_bound(Energy.begin(),Energy.end(),En));    
 
    cross = Cross[index];
    crossDx_product = cross_dx[index];

    // here calculate probability of outgoing particles 
    Nout = (crossDx_product)/maxNorm;

    if(Nout >= 1.0)
      cout << " Nout: " << Nout << " maxNorm: " << maxNorm << " En: " << En << " Energy: " << Energy[index]
	   << " Cross[index]:  " << cross << " Cross[index-1]: " << Cross[index-1] << " product: " << crossDx_product 
	   << " maxCSEnergy: " << maxCSEnergy << " maxCS: " << maxCrossSection << endl;


    return Nout;
  }

};

#endif
