#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "src/DNAVector.h"

void Fill(char * p) {
  int l = strlen(p);
  for (int i=0; i<l; i++) {
    if (p[i] == ' ')
      p[i] = '0';
  }
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fixCmmd("-s","species suffix");
  commandArg<string> faCmmd("-f","fasta file", "");
  commandLineParser P(argc,argv);
  P.SetDescription("Assigns RUM gene and transcript IDs.");
  P.registerArg(fileCmmd);
  P.registerArg(fixCmmd);
  P.registerArg(faCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string ss = P.GetStringValueFor(fixCmmd);
  string fa = P.GetStringValueFor(faCmmd);
  
  vecDNAVector dna;
  if (fa != "")
    dna.Read(fa);


  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  int k=0;
  int l=0;

  string lastGene, lastTrans;

  char tmp[1024];

  string rumGene, rumTrans;
  int i;
  int lastStart = 0;
  string lastOri;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    for (i=0; i<8; i++) {
      if (i == 3) {
	if (parser.AsInt(3) < 1) 
	  cout << "1\t";
	else
	  if (parser.AsInt(3) < parser.AsInt(4)) {
	    cout << parser.AsString(i) << "\t";
	  } else {
	    cout << parser.AsString(4) << "\t";
	  }
	continue;
      }
      if (dna.isize() > 0 && i == 4) {
	if (parser.AsInt(4)+1 >= dna(parser.AsString(0)).isize()) {
	  cout << dna(parser.AsString(0)).isize()-1 << "\t";
	} else  {
	  if (parser.AsInt(3) < parser.AsInt(4)) {
	    cout << parser.AsString(i) << "\t";
	  } else{
	    cout << parser.AsString(3) << "\t";
	  }
	}
	continue;
      }
       
      cout << parser.AsString(i) << "\t";
    }
    string gene, trans;
    int start = parser.AsInt(3);
    const string & ori = parser.AsString(6);
    for (i=8; i<parser.GetItemCount(); i++) {
      if (parser.AsString(i) == "gene_id")
	gene = parser.AsString(i+1);
      if (parser.AsString(i) == "transcript_id")
	trans = parser.AsString(i+1);
    }
    bool bForce = false;
    if (start - lastStart > 500000 || ori != lastOri || lastStart - start > 500000) {
      bForce = true;
    }
      
    if (trans != lastTrans || bForce) {
      k++;
    }
    sprintf(tmp, "\"RUM%sT%10d.%s\";", ss.c_str(), k, parser.AsString(0).c_str());
    Fill(tmp);
    rumTrans = tmp;
    if (gene != lastGene || bForce) {
      l++;
    }
    sprintf(tmp, "\"RUM%sG%10d.%s\";", ss.c_str(), l, parser.AsString(0).c_str());
    Fill(tmp);
    rumGene = tmp;

    lastStart = start;
    lastOri = ori;

    lastGene = gene;
    lastTrans = trans;
    cout << "gene_id " << rumGene << " transcript_id " << rumTrans;

    for (i=8; i<parser.GetItemCount(); i++) {
      cout << " ";
      if (parser.AsString(i) == "gene_id") {
	cout << "rum_oref_gene_id";
	continue;
      }
      if (parser.AsString(i) == "transcript_id") {
	cout << "rum_oref_transcript_id";
	continue;
      }
      cout << parser.AsString(i);
    }
    cout << endl;
  }



  return 0;
}
