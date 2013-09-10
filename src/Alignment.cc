#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <fstream>
#include <iomanip>
#include <cmath>
#include "src/Alignment.h"

//=====================================================================

void Alignment::print(int outputType, double pValLimit, ostream& sout, int screenWidth) const{
  switch(outputType) {
    case 0: printFull(pValLimit, sout, screenWidth); break;
    case 1: printInfoCSV(sout);
  }
}

string Alignment::toString(int screenWidth) const{
  stringstream sout;
  sout << "**********************************************"              << endl
       << "Target sequence size:             " << getTargetLength()        << endl
       << "Query sequence size:              " << getQueryLength()         << endl
       << "Target offset:                    " << getTargetOffset()     << endl
       << "Query offset:                     " << getQueryOffset()      << endl 
       << "Target aligned basepairs:         " << getTargetBaseAligned()<< endl
       << "Query aligned basepairs:          " << getQueryBaseAligned() << endl
       << "Raw Score:                        " << getRawScore()         << endl 
       << "Identity score:                   " << calcIdentityScore()   << endl
       << "Total Edit Count                  " << calcEditCount()       << endl
       << "Mean Contiguity length            " << calcMeanContig()      << endl
       << "Mod-Smith-waterman score:         " << calcModSWScore()      << endl
       << "Significance P-value:             " << calcPVal()            << endl
       << "***********************************************"             << endl;
  for( int i=0; i<=(int)(matchesStr.size()/screenWidth); i++) {
    string query = queryStr.substr(i*screenWidth,screenWidth); 
    sout << "Query: " << setw(6) << getQueryOffset()+(i*screenWidth) << " "
         << query << " " 
         << getQueryOffset()+(i*screenWidth)+min(screenWidth,(int)query.size()) << endl
         << setw(14) << " " << matchesStr.substr(i*screenWidth,screenWidth) << endl
         << "Sbjct: " << setw(6) << getTargetOffset()+(i*screenWidth) << " " 
         << targetStr.substr(i*screenWidth,screenWidth) << " " 
         << getTargetOffset()+(i*screenWidth)+min(screenWidth,(int)query.size()) << endl
         << endl << endl << endl;
  }
  return sout.str();
}

void Alignment::printFull(double pValLimit, ostream& sout, int screenWidth) const{
  if(calcPVal()>pValLimit) { return; } // Significance is less than min required
  sout << toString(screenWidth);
}

void Alignment::printInfoCSV( ostream& sout) const{
  sout<< getTargetOffset()     <<","
      << getQueryOffset()      <<","
      << getTargetBaseAligned()<<","
      << getQueryBaseAligned() <<","
      << calcEditCount()       <<","
      << calcMeanContig()      <<","
      << getRawScore()         <<","
      << calcIdentityScore()   <<","
      << getSWScore()          <<","
      << calcPVal()            <<","
      << getRuntime()          <<","
      << getRuntimeCoef()      <<endl;
}

char Alignment::getQueryAlignCharForTarget(int targetIndex) const { 
  int qIdx = targetIdxsInQuery[targetIndex];
  if( qIdx >= 0 ) { 
    return querySeq[qIdx]; 
  } else {
    return GAP_CHAR; 
  }
}

int Alignment::getQueryAlignIndexForTarget(int targetIndex) const { 
  return targetIdxsInQuery[targetIndex];
}

char Alignment::getTargetAlignCharForQuery(int queryIndex) const { 
  int tIdx = queryIdxsInTarget[queryIndex];
  if( tIdx >= 0 ) { 
    return targetSeq[tIdx]; 
  } else {
    return GAP_CHAR; 
  }
}

int Alignment::getTargetAlignIndexForQuery(int targetIndex) const { 
  return targetIdxsInQuery[targetIndex];
}

double Alignment::calcMeanContig() const {
  int lenCount = 0;
  int totCount = 0;
  int totSum   = 0;
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(matchesStr[i]==MATCH_CHAR) { lenCount++; }
    else { 
      totSum += lenCount; 
      if(lenCount != 0 ) { totCount++; }   
      lenCount  = 0; //reset contig count  
    }
  }  
  totSum += lenCount; 
  if(lenCount !=0 || totCount==0 ) { totCount++; }   
  return double(totSum)/totCount;
}

int Alignment::calcEditCount() const {
  int numOfEdits   = 0;  
  int gapCnt       = 0;
  // Note that one indel (regardless of its length) is considered as one edit
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(targetStr[i]==GAP_CHAR || queryStr[i]==GAP_CHAR) { gapCnt++;}
    else {
      if(gapCnt != 0) { 
        numOfEdits++; // gap ended
        gapCnt = 0;   //reset 
      }
      if(targetStr[i] != queryStr[i]) { numOfEdits++; }  //Mismatch
    }
  }
  return numOfEdits;
}

int Alignment::calcModSWScore() const {
  int totScore = 0;  
  int gapCnt   = 0;
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(targetStr[i]==GAP_CHAR || queryStr[i]==GAP_CHAR) { gapCnt++;}
    else {
      if(gapCnt != 0) { 
        int maxPenalty = min(gapCnt, 10); // Limit gap penalty to 10
        totScore -= maxPenalty; 
        gapCnt = 0;   //reset - gap ended
      }
      if(targetStr[i] != queryStr[i]) { totScore--; }  //Mismatch
      else { totScore++; } //Match
    }
  }
  return totScore;
}

double Alignment::calcPVal() const {
  double lambda = 0.51;
  double mu     = 15.0;
  double x      = calcModSWScore();
  double p      = 1 - exp(-exp(-lambda*(x-mu)));
  return p;
}

void Alignment::resetContainers() {
  targetStr  = "";
  queryStr   = "";
  matchesStr = "";
  targetIdxsInQuery.clear();
  queryIdxsInTarget.clear();
}

//=====================================================================
void AlignmentInfo::resetCalculations(int tos, int qos) {
  tBaseAligned       = 0;
  qBaseAligned       = 0;
  baseMatched        = 0;
  alignmentLen       = 0;
  smithWatermanScore = 0;

  // Set offset values
  tOffset            = tos;
  qOffset            = qos;  
}

double AlignmentInfo::calcIdentity() const {
  return double(baseMatched)/getMaxBaseAligned();
}

