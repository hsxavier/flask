/***************************************************************
2012-09-18: function definitions for the class 'ParameterList' defined 
in ParameterList.hpp.
***************************************************************/

#include <iostream>  // standard I/O
#include <fstream>   // file I/O
#include <cstring>   // to use strcpy();

#include "ParameterList.hpp"
#include "Utilities.hpp"

// Default constructor
ParameterList::ParameterList () {
  using namespace ParDef;
  int i;
  parloaded = 0;
  for (i=0; i<MAXPARS; i++) {
    list[i].name[0]='\0';
    list[i].value.inum=0;
    list[i].type=i1;
  }
}
// Constructor with input
ParameterList::ParameterList (const char *filename) {
  ParameterList();
  load(filename);
}

// Finds on the ParDef namespace the parameter mentioned in file.
int ParameterList::findpar (std::string word) const {
  using namespace ParDef;
  void warning(const std::string message);
  int index=0;
  word=word.substr(0,word.size()-1);
  while (index<npars && par_name[index]!=word) index++;
  if (index>=npars) { warning("ParameterList::findpar: "+word+" is not listed in code!"); return -1;}
  return index;
}

// Returns int value of parameter (at pos if array).
int ParameterList::readi(int index, int pos) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50];
  bool flag_noElement=0;
  switch (list[index].type) {
  case i1:
    if (pos==-1) return list[index].value.inum;
    else flag_noElement=1;
    break;
  case i2:
    if (pos>=0 && pos<2) return list[index].value.ivec[pos];
    else flag_noElement=1;
    break;
  case i3:
    if (pos>=0 && pos<3) return list[index].value.ivec[pos];
    else flag_noElement=1;
    break;
  case d1:
  case d2:
  case d3:
  case c:
  case s:
    error("ParameterList::readi<id>: not suited for type "+typelabel[list[index].type]);
    break;
  default:
    error("ParameterList::readi<id>: not prepared for type "+typelabel[list[index].type]);
  }   
  if (flag_noElement==1) 
    {sprintf(message,"ParameterList::readi<id>: element %s[%d] inexists.",list[index].name,pos); error (message);}
}
int ParameterList::readi(std::string name, int pos) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::readi<name>: cannot find parameter "+name+".");
  return readi(index, pos);
}
// Returns double value of parameter (at pos if array).
double ParameterList::readd(int index, int pos) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50];
  bool flag_noElement=0;
  switch (list[index].type) {
  case d1:
    if (pos==-1) return list[index].value.dnum;
    else flag_noElement=1;
    break;
  case d2:
    if (pos>=0 && pos<2) return list[index].value.dvec[pos];
    else flag_noElement=1;
    break;
  case d3:
    if (pos>=0 && pos<3) return list[index].value.dvec[pos];
    else flag_noElement=1;
    break;
  case i1:
  case i2:
  case i3:
  case c:
  case s:
    error("ParameterList::readd<id>: not suited for type "+typelabel[list[index].type]);
    break;
  default:
    error("ParameterList::readd<id>: not prepared for type "+typelabel[list[index].type]);
  }   
  if (flag_noElement==1) 
    {sprintf(message,"ParameterList::readd<id>: element %s[%d] inexists.",list[index].name,pos); error (message);}
}
double ParameterList::readd(std::string name, int pos) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::readd<name>: cannot find parameter "+name+".");
  return readd(index, pos);
}
// Returns char value of parameter.
char ParameterList::readc(int index, int pos) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50];
  bool flag_noElement=0;
  switch (list[index].type) {
  case c:
    if (pos==-1) return list[index].value.cnum;
    else flag_noElement=1;
    break;
  case s:
    if (pos>=0) return list[index].value.cvec[pos];
    else error("ParameterList::readc<id>: especify string element or try reads().");
    break;
  case i1:
  case i2:
  case i3:
  case d1:
  case d2:
  case d3:
    error("ParameterList::readc<id>: not suited for type "+typelabel[list[index].type]);
    break;
  default:
    error("ParameterList::readc<id>: not prepared for type "+typelabel[list[index].type]);
  }   
  if (flag_noElement==1) 
    {sprintf(message,"ParameterList::readc<id>: element %s[%d] inexists.",list[index].name,pos); error (message);}
}
char ParameterList::readc(std::string name, int pos) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::readc<name>: cannot find parameter "+name+".");
  return readc(index, pos);
}
// Returns string value of parameter.
std::string ParameterList::reads(int index) const {
  using namespace ParDef;
  void error (const std::string message);
  if (list[index].type==s) return list[index].value.cvec;
  if (list[index].type==ph) return list[index].value.cvec;
  else error("ParameterList::reads<id>: not suited for type "+typelabel[list[index].type]);
}
std::string ParameterList::reads(std::string name) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::reads<name>: cannot find parameter "+name+".");
  return reads(index);
}


// Copy parameter int value to int variable.
void ParameterList::copy (int index, long *value) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50];
  if (index >=npars) // Error in case index does not exist.
    {sprintf(message,"ParameterList::copy<id,int>: unkown index %d.",index); error (message);}
  switch (list[index].type) {
  case i1:
    (*value)=list[index].value.inum;
    break;
  case i2:
    value[0]=list[index].value.ivec[0];
    value[1]=list[index].value.ivec[1];
    break;
  case i3:
    value[0]=list[index].value.ivec[0];
    value[1]=list[index].value.ivec[1];
    value[2]=list[index].value.ivec[2];
    break;
  default:
    error("ParameterList::copy<id,int>: not prepared for type "+typelabel[list[index].type]);
  }
}
void ParameterList::copy (std::string name, long *value) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::copy<name,int>: cannot find parameter "+name+".");
  copy(index, value);
}
// Copy parameter double value to double variable.
void ParameterList::copy (int index, double *value) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50]; 
  if (index >=npars) // Error in case index does not exist.
    {sprintf(message,"ParameterList::copy<id,double>: unkown index %d.",index); error (message);}
  switch (list[index].type) {
  case d1:
    (*value)=list[index].value.dnum;
    break;
  case d2:
    value[0]=list[index].value.dvec[0];
    value[1]=list[index].value.dvec[1];
    break;
  case d3:
    value[0]=list[index].value.dvec[0];
    value[1]=list[index].value.dvec[1];
    value[2]=list[index].value.dvec[2];
    break;
  default:
    error("ParameterList::copy<id,double>: not prepared for type "+typelabel[list[index].type]);
  }
}
void ParameterList::copy (std::string name, double *value) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::copy<name,double>: cannot find parameter "+name+".");
  copy(index, value);
}


// Copy parameter char value to char variable.
void ParameterList::copy (int index, char *value) const {
  using namespace ParDef;
  void error (const std::string message);
  char message[50]; 
  if (index >=npars) // Error in case index does not exist.
    {sprintf(message,"ParameterList::copy<id,char>: unkown index %d.",index); error (message);} 
  switch (list[index].type) {
  case c:
    (*value)=list[index].value.cnum;
    break;
  case s:
    int i;
    for(i=0; list[index].value.cvec[i]!='\0'; i++) value[i]=list[index].value.cvec[i];  
    value[i]='\0';
    break;
  default:
    error("ParameterList::copy<id,char>: not prepared for type "+typelabel[list[index].type]);
  }
}
void ParameterList::copy (std::string name, char *value) const {
  void error (const std::string message);
  int index;
  index = findpar(name+":");
  if (index < 0) error ("ParameterList::copy<name,char>: cannot find parameter "+name+".");
  copy(index, value);
}


// Read parameters from file.
void ParameterList::load (const char *filename) {
  using std::ifstream;
  using std::string;
  using namespace ParDef;
  void error (const std::string message);
  void warning (const std::string message);
  ifstream parfile;
  string word;
  int index;
  parfile.open(filename);
  if (!parfile.is_open()) error("ParameterList::load: cannot open file.");
  
  parloaded=0;
  while (parfile >> word) {
    if (word[word.size()-1] == ':') {
      index=findpar(word);          // Lookup parameter in 'ParDef' namespace and get its index.
      if (index>=0) {
	parloaded++;
	if (parloaded>MAXPARS) error("ParameterList::load: too many parameters. Increase MAXPARS.");
	if (par_name[index].size() > MAXPARNAME-1) error("ParameterList::load: name "+par_name[index]+" is too big. Increase MAXPARNAME.");
	list[index].type=par_type[index];                   // Copy parameter type to parameter list.
	strcpy(list[index].name, par_name[index].c_str());  // Copy parameter name to parameter list.
	switch (par_type[index]) {                          // Copy parameter value to parameter list.
	case i1:
	  parfile >> list[index].value.inum;
	  break;
	case i2:
	  parfile >> list[index].value.ivec[0];
	  parfile >> list[index].value.ivec[1];
	  break;
	case i3:
	  parfile >> list[index].value.ivec[0];
	  parfile >> list[index].value.ivec[1];  
	  parfile >> list[index].value.ivec[2];  
	  break;
	case d1:
	  parfile >> list[index].value.dnum;
	  break;
	case d2:
	  parfile >> list[index].value.dvec[0];
	  parfile >> list[index].value.dvec[1];
	  break;
	case d3:
	  parfile >> list[index].value.dvec[0];
	  parfile >> list[index].value.dvec[1];  
	  parfile >> list[index].value.dvec[2];  
	  break; 
	case c: 
	  parfile >> list[index].value.cnum;
	  break;
	case s: 
	  parfile >> list[index].value.cvec;
	  break;
	case ph:
	  parfile.getline(list[index].value.cvec, STRSIZE);
	  if (strlen(list[index].value.cvec) == STRSIZE-1) 
	    warning("ParameterList::load: "+par_name[index]+" too big, increase STRSIZE in code.");
	  break;
	default:
	  error("ParameterList::load: not prepared for type "+typelabel[par_type[index]]);
	}
      }
    }
  }
  if (parloaded < npars) warning("ParameterList::load: some parameters are missing. Check ParDef namespace.");
  parfile.close();
}


// Read parameters from command line.
// It assumes Parameterlist::load was ran before it.
void ParameterList::lineload (int argc, char *argv[]) {
  using std::string;
  using namespace ParDef;
  void error (const std::string message);
  void warning (const std::string message);
  char word[MAXPARNAME+2];
  int index, n=2, parupdated=0;
  
  while (n<argc) {
    n++;
    if (argv[n-1][strlen(argv[n-1])-1]==':') sprintf(word,"%s",argv[n-1]);
    else sprintf(word,"%s:",argv[n-1]);
    index=findpar(word);          // Lookup parameter in 'ParDef' namespace and get its index.
    if (index>=0) {
      if (par_name[index].size() > MAXPARNAME-1) error("ParameterList::lineload: name "+par_name[index]+" is too big. Increase MAXPARNAME.");
      if (strcmp(list[index].name, "")!=0) parupdated++;  // Parameter was loaded before and is being updated.
      else {
	parloaded++;
	if (parloaded>MAXPARS) error("ParameterList::lineload: too many parameters. Increase MAXPARS.");
	list[index].type=par_type[index];                   // Copy parameter type to parameter list.
	strcpy(list[index].name, par_name[index].c_str());  // Copy parameter name to parameter list.
      }
      switch (par_type[index]) {                          // Copy parameter value to parameter list.
      case i1:
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.inum);
	break;
      case i2:
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.ivec[0]);
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.ivec[1]);
	break;
      case i3:
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.ivec[0]);
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.ivec[1]);
	n++;
	sscanf(argv[n-1],"%ld", &list[index].value.ivec[2]);
	break;
      case d1:
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dnum);
	break;
      case d2:
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dvec[0]);
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dvec[1]);
	break;
      case d3:
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dvec[0]);
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dvec[1]);
	n++;
	sscanf(argv[n-1],"%lf", &list[index].value.dvec[2]);  
	break; 
      case c: 
	n++;
	sscanf(argv[n-1],"%c", &list[index].value.cnum);
	break;
      case s:
	n++;
	sscanf(argv[n-1],"%s", list[index].value.cvec);
	break;
      default:
	error("ParameterList::lineload: not prepared for type "+typelabel[par_type[index]]);
      }
    }
  }
  if (parloaded < npars) warning("ParameterList::lineload: some parameters are missing. Check ParDef namespace.");
  std::cout<<"   Updated "<<parupdated<<" parameters from command line.\n";
}




// Displays on screen the parameters loaded. 
void ParameterList::show (std::ostream * output) const {
  using namespace ParDef;
  using std::endl;
  void warning (const std::string message);
  int i;
  const int numWidth=8, typeWidth=10;
  
  for (i=0; i<parloaded; i++) {
    (*output).width(MAXPARNAME);
    *output << list[i].name << ": ";    // Print parameter name.
    (*output).width(numWidth);
    switch (list[i].type) {             // Print parameter value.
    case i1:
      *output << list[i].value.inum;           (*output).width(2*(numWidth+1)+typeWidth);
      break;
    case i2:
      *output << list[i].value.ivec[0] << " "; (*output).width(numWidth);
      *output << list[i].value.ivec[1];        (*output).width(numWidth+1+typeWidth);
      break;
    case i3:
      *output << list[i].value.ivec[0] << " "; (*output).width(numWidth);
      *output << list[i].value.ivec[1] << " "; (*output).width(numWidth);
      *output << list[i].value.ivec[2];        (*output).width(typeWidth);
      break;
    case d1:
      *output << list[i].value.dnum;           (*output).width(2*(numWidth+1)+typeWidth);
      break;
    case d2:
      *output << list[i].value.dvec[0] << " "; (*output).width(numWidth);
      *output << list[i].value.dvec[1];        (*output).width(numWidth+1+typeWidth);
      break;
    case d3:
      *output << list[i].value.dvec[0] << " "; (*output).width(numWidth);
      *output << list[i].value.dvec[1] << " "; (*output).width(numWidth);
      *output << list[i].value.dvec[2];        (*output).width(typeWidth);
      break;
    case c:
      *output << list[i].value.cnum;           (*output).width(2*(numWidth+1)+typeWidth);
      break;
    case s:
      (*output).width(3*numWidth+2);
      *output << list[i].value.cvec;           (*output).width(typeWidth);        
      break;
    case ph:
      (*output).width(3*numWidth+2);
      *output << list[i].value.cvec;           (*output).width(typeWidth); 
      break;
    default:
      warning("ParameterList::show: not prepared for type "+typelabel[list[i].type]);
    }
    *output <<"  (type: ";                             // Print parameter type.
    (*output).width(2);
    *output << typelabel[list[i].type] << ")" << endl;
  }
}
