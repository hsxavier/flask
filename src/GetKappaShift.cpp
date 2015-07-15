

void ClassColName2ftypez(std::string word, int *ft1, int *z1, int *ft2, int *z2) {
  int start1, start2, end1, end2;

  // Get the type of the first element:
  start1 = word.find(":dens[");
  if (start1!=word.npos) *ft1=1;
  else {
    start1 = word.find(":lens[");
    if (start1!=word.npos) *ft1=2;
  }
  // Set to -1 if unknown:
  else *ft1=-1;

  // Get the type of the second element:
  start2 = word.find("-dens[");
  if (start2!=word.npos) *ft2=1;
  else {
    start2 = word.find(":lens[");
    if (start2!=word.npos) *ft2=2;
  }
  // Set to -1 if unknown:
  else *ft2=-1;


  if(start1!=word.npos && start2!=word.npos) {
    end1 = word.find("]", start1);
    end2 = word.find("]", start2);
    *i   = atoi(word.substr(start1+6,end1-1).c_str());
    *j   = atoi(word.substr(start2+6,end2-1).c_str());
  }
  else {
    *i   = -1;
    *j   = -1;
  }
}


double **LoadClassCls(std::string filename, long *nr, long *nc, int offset=0, int verbose=0) {
  using std::ifstream;
  using std::string;
  using std::istringstream;
  using std::ostringstream;
  long nrows=0, ncols=0, i, j, nheaders=0;
  ifstream file;
  istringstream inputline; ostringstream outputline;
  string word, phrase, ColumnNames;
  type **table;
  int datapos=0, start1, start2, end1, end2, max;
  
  // Open file
  file.open(filename.c_str());
  if (!file.is_open()) error("LoadTable: cannot open file "+filename);
  
  // Detect headers (must start with # or empty lines):
  getline(file,phrase);
  while(!file.eof() && (phrase[0]=='#' || phrase.length()==0)) {
    if (phrase.substr(0,4)=="# 1:l") ColumnNames = phrase;       // Save column names for later.
    datapos     = file.tellg(); 
    nheaders++; 
    getline(file,phrase);
  }

  // Get position of covariance matrix elements: 
  inputline.str(ColumnNames);
  max = 0;
  while (inputline >> word) {

    start1 = word.find(":dens[");
    start2 = word.find("-dens[");
    
    if(start1!=word.npos && start2!=word.npos) {
      end1   = word.find("]", start1);
      end2   = word.find("]", start2);
      i      = atoi(word.substr(start1+6,end1-1).c_str());
      j      = atoi(word.substr(start2+6,end2-1).c_str());
      if (i>max) max=i;
      if (j>max) max=j; 
    }
  }

  // Count number of columns (using first data row):
  outputline << phrase;
  inputline.str(outputline.str());
  while (inputline >> word) ncols++;
  // Count number of rows (ignoring comments and empty spaces):
  while(!file.eof() && phrase[0]!='#' && phrase.length()!=0) {getline(file,phrase); nrows++;}
  if (phrase.length()!=0 && phrase[0]!='#') nrows++;
  // Inform number of rows and columns if requested:
  if (verbose!=0) 
    std::cout<<"LoadTable will allocate "<<nrows<<" lines and "<<ncols<<" columns for file "<<filename<<std::endl;
  
  // Loading values to table:
  file.clear();
  file.seekg(datapos);
  table = matrix<type>(offset,nrows+offset-1,offset,ncols+offset-1);
  for (i=offset; i<nrows+offset; i++)
    for (j=offset; j<ncols+offset; j++) 
      if (!(file >> table[i][j])) error("LoadTable: more data expected in file "+filename); // DO NOT put comments in the end of the file!
  if(file >> word && word[0]!='#') error("LoadTable: data was ignored in "+filename);
  *nr=nrows; *nc=ncols;

  file.close();
  return table;
}
