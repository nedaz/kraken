#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include<fstream>
#include "src/AlignmentBlock.h"

//=====================================================================

string AlignmentBlock::toString() const{
  stringstream outStream;
  outStream << "Target:  "  << getTargetChrom() << ":" << getTargetStart() << "-" << getTargetStop() 
            <<  " Query:  " << getQueryChrom() << ":" << getQueryStart() << "-" << getQueryStop();
  return outStream.str();
}

void AlignmentBlock::print() const { cout << toString() << endl; } 

bool AlignmentBlock::parse(FlatFileParser& parser, bool flip) {
  if(!parser.ParseLine()) { return false; }
  //cout << "Next block: " << parser.Line() << endl;

  if (flip) {
    queryChrom  = parser.AsString(0);
    targetChrom = parser.AsString(3);
    queryStart  = parser.AsInt(1);
    queryStop   = parser.AsInt(2);
    targetStart = parser.AsInt(4);
    targetStop  = parser.AsInt(5);
    orientation = parser.AsChar(7);
  } else {
    targetChrom = parser.AsString(0);
    queryChrom  = parser.AsString(3);
    targetStart = parser.AsInt(1);
    targetStop  = parser.AsInt(2);
    queryStart  = parser.AsInt(4);
    queryStop   = parser.AsInt(5);
    orientation = parser.AsChar(7);
  }

  if (targetStop - targetStart > 10000)
    targetStop = targetStart + 100;
  if (queryStop - queryStart > 10000)
    queryStop = queryStart + 100;

  return true;
}

bool AlignmentBlock::merge(const AlignmentBlock& toMerge) {
  // Make sure the alignment blocks are from the same chromosomes and orientation
  if(!isCompatible(toMerge)) { return false; }

  if (orientation == '-' && isWithinBound(targetStop, toMerge.targetStart) 
      && isWithinBound(queryStop, toMerge.queryStart)) {
    targetStop = toMerge.targetStop;
    queryStop  = toMerge.queryStop;
    return true;
  } 
  //  If we are dealing with reverse orientation, query should be used in opposite direction
  if (orientation == '-' && isWithinBound(targetStop, toMerge.targetStart) 
      && isWithinBound(toMerge.queryStop, queryStart)) {
    targetStop = toMerge.targetStop;
    queryStart = toMerge.queryStart;
    return true;
  }
  // If the past two conditions haven't passed, merge cannot be performed
  return false;
}

bool AlignmentBlock::isWithinBound(int num1, int num2) { 
  return (num1 <= num2 + mergeBoundary && num1 >= num2 - mergeBoundary); 
}

bool AlignmentBlock::isCompatible(const AlignmentBlock& otherBlock)const {
  if (targetChrom == otherBlock.getTargetChrom()  
      && queryChrom == otherBlock.queryChrom  
      && orientation == otherBlock.orientation) {
    return true;
  } else { 
    return false;
  }
}
