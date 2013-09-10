#include "src/Kraken/KrakenConfig.h"
#include "base/FileParser.h"
#include "src/Kraken/KrakenXMFA.h"


KrakenConfig::KrakenConfig(Kraken * p)
{
  m_pKraken = p;
}

enum K_SECTION
{
  K_SECTION_NONE,
  K_SECTION_GENOME,
  K_SECTION_MAP,
  K_SECTION_XMFA
};


bool KrakenConfig::Configure(const string & fileName)
{
  FlatFileParser parser;
  
  parser.Open(fileName);

  K_SECTION s = K_SECTION_NONE;


  svec<string> genome, file;
  svec<string> kmap, g1, g2;

  int i;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsString(0) == "//" || parser.AsString(0) == "#") {
      continue;
    }

    if (parser.AsString(0) == "[genomes]") {
      s = K_SECTION_GENOME;
      continue;
    }
    if (parser.AsString(0) == "[pairwise-maps]") {      
      s = K_SECTION_MAP;
      continue;
    }
    if (parser.AsString(0) == "[multiple-alignments]") {
      s = K_SECTION_XMFA;
      continue;
    }    
    switch(s) {
    case K_SECTION_NONE:
      break;
    case K_SECTION_GENOME:
      genome.push_back(parser.AsString(0));
      file.push_back(parser.AsString(1));     
      break;
    case K_SECTION_MAP:
      kmap.push_back(parser.AsString(2));
      g1.push_back(parser.AsString(0));     
      g2.push_back(parser.AsString(1));     
      break;
    case K_SECTION_XMFA:
//      cout << "Creating directory " << output << endl;
      system("mkdir Kraken_temp");
 
      string xmfaFile = parser.AsString(0);
      KrakenXMFA handler;
      svec<string> outFiles;
      handler.convertXMFA(parser.AsString(0), "Kraken_temp", outFiles);
      for(svec<string>::iterator iter=outFiles.begin(); iter!=outFiles.end(); ++iter) {
        StringParser sParser;
        sParser.SetLine(*iter, "\t");
        kmap.push_back(sParser.AsString(2));
        g1.push_back(sParser.AsString(0));     
        g2.push_back(sParser.AsString(1));     
      }
      break;
    }

  }

  for (i=0; i<kmap.isize(); i++) {
    m_pKraken->Allocate(g1[i], g2[i]);
  }

  m_pKraken->DoneAlloc();
  FILE_LOG(logDEBUG) << "Done allocating genomes.";  
  for (i=0; i<kmap.isize(); i++) {
    FILE_LOG(logDEBUG) << "Reading map: " << kmap[i];  
    m_pKraken->ReadMap(kmap[i], g1[i], g2[i]);
  }

  for (i=0; i<genome.isize(); i++) {
    FILE_LOG(logDEBUG) << "Reading genome: " << genome[i] << "\t" << file[i]; 
    m_pKraken->ReadGenome(file[i], genome[i]);
  }
  FILE_LOG(logDEBUG) << "Done reading!";
  return true;
}
