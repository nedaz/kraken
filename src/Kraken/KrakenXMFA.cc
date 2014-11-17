#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Kraken/KrakenXMFA.h"


//======================================================
void KrakenXMFA::convertXMFA(const string& inFile, const string& outDir, svec<string>& outFiles) {
  FlatFileParser parser;
  parser.Open(inFile);
  map<string, Coordinate> coords;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string startStr = parser.AsString(0);
    if (startStr.at(0) == '>') { //Alignment info
      int posColon  = startStr.find(':'); 
      string specie = startStr.substr(1, posColon-1); 
      int posDash   = startStr.find('-'); 
      int start     = atol(startStr.substr(posColon+1, posDash-1).c_str()); 
      int stop      = atol(startStr.substr(posDash+1).c_str()); 
      bool orient   = (parser.AsString(1)=="+")?true:false;
      string chr    = parser.AsString(2);
      coords[specie] = Coordinate(chr, orient, start, stop);
    }

    if (startStr.at(0) == '=') { //End of set
      map<string, string> outStrings;
      outSatsumaBlocks(coords, outStrings);
      for(map<string, string>::iterator iter=outStrings.begin(); iter!=outStrings.end(); ++iter) {
        ofstream fout;
        string pair = iter->first;
        replace(pair.begin(), pair.end(), '\t', '_');
        string file = outDir+"/"+"synteny_"+pair+".chained"; //For other OS, need to change 
        fout.open(file.c_str()); 
        fout<<iter->second;
        fout.close();
        outFiles.push_back(iter->first+'\t'+file);
      }
      coords.clear();
    }
  }
}

void KrakenXMFA::outSatsumaBlocks(const map<string, Coordinate>& coords,
                                  map<string, string>& outStrings) {
  for(map<string, Coordinate>::const_iterator iter1=coords.begin(); iter1!=coords.end(); ++iter1) {
    for(map<string, Coordinate>::const_iterator iter2=iter1; iter2!=coords.end(); ++iter2) {
      // Map is sorted on the key in the same way so no need to check for the reverse order
      string spPair = iter1->first + '\t' + iter2->first;
      outStrings[spPair] += (iter1->second.toString_noOrient('\t') + '\t' + 
                             iter2->second.toString_noOrient('\t') + '\t' + 
                             (iter1->second.isSameOrient(iter2->second)?"+":"-") + '\n'); 
    }
  }
}

