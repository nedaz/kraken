#ifndef _MULTIALIGNPARSE_H_
#define _MULTIALIGNPARSE_H_

#include <map>
#include <string>
#include <sstream>
#include "ryggrad/src/base/SVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/AlignmentBlock.h"
#include "ryggrad/src/general/Coordinate.h"



class MultiAlignParser
{ 
public:
  /** Default Ctor */
  MultiAlignParser() {}

  /** */
  void convertXMFA(const string& inFile, const string& outDir, svec<string>& outFiles, const string& onlyOneSpecie=""); 
  void convertMAF(const string& inFile, const string& outDir, svec<string>& outFiles, const string& onlyOneSpecie=""); 

private:
  void outSatsumaBlocks(const map<string, Coordinate>& coords, map<string, string>& outStreams, const string& onlyOneSpecie=""); 
};

#endif //_MULTIALIGNPARSE_H_
