#ifndef _KRAKENMULTIALIGN_H_
#define _KRAKENMULTIALIGN_H_

#include <map>
#include "../annotationQuery/AnnotationQuery.h"

class KrakenMultiAlign
{ 
public:
  /** Default Ctor */
  KrakenMultiAlign() {}

  /** */
  void convertXMFA(const string& inFile, const string& outDir, svec<string>& outFiles); 

private:
  void outSatsumaBlocks(const map<string, Coordinate>& coords, map<string, string>& outStreams); 
};

#endif //_KRAKENMULTIALIGN_H_
