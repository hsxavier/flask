#include <iostream>
#include "ParameterList.hpp"
#include "Utilities.hpp"

std::ofstream debugfile;

/*** Main Program ***/

int main (int argc, char *argv[]) {
  using std::cout; using std::endl; 
  using std::string; using std::ofstream;
  using namespace ParDef;
  ParameterList config; 

  /*** Opening debug file for dumping information about the program ***/
  debugfile.open("debug.log");
  if (!debugfile.is_open()) error("main: cannot open debug file.");

  /*** Loading config file ***/
  
  if (argc<=1) { cout << "You must supply a config file." << endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl; 

  /*** End of the program ***/

  debugfile.close(); // Close debug file.

  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}

