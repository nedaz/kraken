#ifndef _KRAKEN_XMFA_H_
#define _KRAKEN_XMFA_H_

#include <map>
#include "src/AnnotationQuery/AnnotationQuery.h"

class KrakenXMFA
{ 
public:
  /** Default Ctor */
  KrakenXMFA() {}

  /** */
  void convertXMFA(const string& inFile, const string& outDir, svec<string>& outFiles); 

private:
  void outSatsumaBlocks(const map<string, Coordinate>& coords, map<string, string>& outStreams); 
};

#endif //_KRAKEN_XMFA_H_
