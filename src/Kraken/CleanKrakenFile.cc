#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



string Append(const string & s, const string &a) {
  char tmp[256];
  string out;
  strcpy(tmp, s.c_str());
  tmp[strlen(tmp)-2] = 0;
  out = tmp;
  out += a;
  out += "\";";
  return out;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Checks a GTF file for consistency.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  string last;
  string chr;
  string ori;
  int startLast = -1;
  string suffix;
  int i;

  while (parser.ParseLine()) {
    if(parser.GetItemCount() < 7)
      continue;
    string s;
    int transIndex = -1;
    for (i=0; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "transcript_id") {
	s = parser.AsString(i+1);
	transIndex = i+1;
      }
    }
    if (s != last) {
      suffix = "";
    }
    const string & c = parser.AsString(0);
    const string & o = parser.AsString(6);

    int start = parser.AsInt(3);

    if (parser.AsInt(4)-start > 50000)
      continue;

    if (parser.AsInt(4) < parser.AsInt(3)) {
      chr = c;
      ori = o;
      continue;
    }

    //if (s == "\"ENST00000395088\";") {
    //  cout << parser.Line() << endl;
    //  cout << last << " " << s << " " << o << " " << ori << endl;
    //}
    
    if (s == last) {
      if (start - startLast < 0 || start - startLast > 2000000)
	suffix += "X";
	startLast = start;
      if (c != chr)
	continue;
      if (o != ori)
	continue;
    } else {
      last = s;
      suffix = ""; 
    }
    startLast = start;

    chr = c;
    ori = o;
    for (i=0; i<8; i++) {
      if (i == 3) {
	int ss = parser.AsInt(i);
	if (ss < 0)
	  ss = 0;
	cout << ss << "\t";	
      } else {
	cout << parser.AsString(i) << "\t";
      }
    }

    for (i=8; i<parser.GetItemCount(); i++) {
      
      if (i == transIndex)
	cout << Append(parser.AsString(i), suffix);
      else
	cout << parser.AsString(i);
      if (i+1 < parser.GetItemCount())
	cout << " ";
    }
    cout << endl;
    //cout << parser.Line() << endl;
  }

  return 0;
}
