#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "MultiAlignParser.h"


//======================================================
void MultiAlignParser::convertXMFA(const string& inFile, const string& outDir, svec<string>& outFiles, const string& onlyOneSpecie) {
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
      outSatsumaBlocks(coords, outStrings, onlyOneSpecie);
      for(map<string, string>::iterator iter=outStrings.begin(); iter!=outStrings.end(); ++iter) {
        ofstream fout;
        string pair = iter->first;
        replace(pair.begin(), pair.end(), '\t', '_');
        string file = outDir+"/"+"synteny_"+pair+".chained"; //For other OS, need to change 
        fout.open(file.c_str(), std::ofstream::app); 
        fout<<iter->second;
        fout.close();
        outFiles.push_back(iter->first+'\t'+file);
      }
      coords.clear();
    }
  }
}

void MultiAlignParser::convertMAF(const string& inFile, const string& outDir, svec<string>& outFiles, const string& onlyOneSpecie) {
  FlatFileParser parser;
  parser.Open(inFile);
  map<string, Coordinate> coords;
  while (parser.ParseLine()) {
//    if (parser.GetItemCount() == 0)
//      continue;
    string startStr = parser.AsString(0);
    if (startStr.at(0) == 's') { //Alignment info
      string ident  = parser.AsString(1);
      int posColon  = ident.find('.'); 
      string specie = ident.substr(0, posColon); 
      string chr    = ident.substr(posColon+1, ident.length()); 
      int start     = parser.AsInt(2);
      int stop      = start + parser.AsInt(3);
      bool orient   = (parser.AsString(4)=="+")?true:false;
      coords[specie] = Coordinate(chr, orient, start, stop);
    }

//    if (startStr.at(0) == '\n') { //End of set
    if (parser.GetItemCount()==0) { //Emprty line : End of set
      map<string, string> outStrings;
      outSatsumaBlocks(coords, outStrings, onlyOneSpecie);
      for(map<string, string>::iterator iter=outStrings.begin(); iter!=outStrings.end(); ++iter) {
        ofstream fout;
        string pair = iter->first;
        replace(pair.begin(), pair.end(), '\t', '_');
        string file = outDir+"/"+"synteny_"+pair+".chained"; //For other OS, need to change 
        fout.open(file.c_str(), std::ofstream::app); 
        fout<<iter->second;
        fout.close();
        outFiles.push_back(iter->first+'\t'+file);
      }
      coords.clear();
    }
  }
}

void MultiAlignParser::outSatsumaBlocks(const map<string, Coordinate>& coords,
                                  map<string, string>& outStrings, const string& onlyOneSpecie) {
  for(map<string, Coordinate>::const_iterator iter1=coords.begin(); iter1!=coords.end(); ++iter1) {
    for(map<string, Coordinate>::const_iterator iter2=coords.begin(); iter2!=coords.end(); ++iter2) {
      // Map is sorted on the key in the same way so no need to check for the reverse order
      if(iter1->first==iter2->first)                      { continue; }
      if(onlyOneSpecie!="" && iter1->first!=onlyOneSpecie) { continue; }
      string spPair = iter1->first + '\t' + iter2->first;
      outStrings[spPair] += (iter1->second.toString_noOrient('\t') + '\t' + 
                             iter2->second.toString_noOrient('\t') + '\t' + 
                             (iter1->second.isSameOrient(iter2->second)?"+":"-") + '\n'); 
    }
  }
}

