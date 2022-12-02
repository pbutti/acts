#ifdef __CINT__

#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 

#pragma link C++ nestedclasses;

// This is to create the dictionary for stl containers

#pragma link C++ class vector<vector<int>>   +;
//#pragma link C++ class vector<TObject>     +;
//#pragma link C++ class vector<TObject*>    +;
#endif
