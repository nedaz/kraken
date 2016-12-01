#ifndef KRAKENMAP_H
#define KRAKENMAP_H

#include "ryggrad/src/base/Logger.h"
#include "ryggrad/src/general/DNAVector.h"
#include "ryggrad/src/base/FileParser.h"
#include "ryggrad/src/general/MultXCorr.h"
#include "ryggrad/src/general/AlignmentBlock.h"
#include "../annotationQuery/AnnotationQuery.h"
#include "KrakenParams.h"

class GenomeWideMap
{
public:
  GenomeWideMap() {
    m_distance = 0.5;
  }

  void Set(const string & source, const string & target, double distance = 0.5) {
    m_source = source;
    m_target = target;
    m_distance = distance;
  }

  void Read(const string & fileName, const string & source, const string & target, bool flip, double distance = 0.5);
  bool Map(const Coordinate& lookup, svec<Coordinate>& results);

  bool operator < (const GenomeWideMap & m) const {
    if (m_source != m.m_source) {
      return (m_source < m.m_source);
    }
    return (m_target < m.m_target);
  }

  void Print() const {
    FILE_LOG(logDEBUG) << "T=" << m_source << " Q=" << m_target << endl;
  }

  int GetBlockCount() const {return m_blocks.isize();}
  const AlignmentBlock& GetBlock(int i) const {return  m_blocks[i];}

  const string & Destination() const {return m_target;}
  const string & Origin() const {return m_source;}
  double Distance() const {return m_distance;}
private:
  void SetAnchors(const Coordinate & lookup, const AlignmentBlock& beginBlock, 
                  const AlignmentBlock& endBlock, int startExtend, int stopExtend,
                  Coordinate & result); 
  void MergeBlocks();


  string m_source;
  string m_target;
  double m_distance;
  svec<AlignmentBlock> m_blocks;
};



class GenomeSeq
{
public:
  GenomeSeq() {}
  GenomeSeq(const string &n) {m_name = n;}

  void Read(const string & fileName, const string & genome) {
    m_dna.Read(fileName);
    m_name = genome;
  }

  void SetName(const string & name) {
    m_name = name;
  }

  bool operator < (const GenomeSeq & s) const {
    return (m_name < s.m_name);
  }

  const vecDNAVector & DNA() const {return m_dna;}
  const string & Name() const {return m_name;}
  
private:
  vecDNAVector m_dna;
  string m_name;
};


//=========================================================

class Route
{
public:
  Route() {
    m_dist = 0.;
    m_invalid = false;
  }

  int GetCount() const {return m_source.isize();}
  const string & Origin(int i) const {return m_source[i];}
  const string & Destination(int i)  const {return m_target[i];}

  bool IsInvalid() const {return m_invalid;}
  void SetInvalid() {m_invalid = true;}

  void Add(const string & source, const string & target, double dist = 0.) {
    m_source.push_back(source);
    m_target.push_back(target);
    m_dist += dist;
  }

private:
 
  svec<string> m_source;
  svec<string> m_target;
  double m_dist;
  bool m_invalid;


};


class Kraken;

class RouteFinder
{
public:
  RouteFinder() {}
  
  void SetNumGenomes(int genomes) {
    m_routes.resize(genomes * genomes);
  }

  bool FindRoute(Route & out, const string & source, const string & target, Kraken & rum);

private:
  bool FindRecursive(svec<string> & final, svec<string> & path, const string & target, Kraken & rum) const;
  
  svec<Route> m_routes;
};


//=========================================================


class Kraken
{
friend class RouteFinder;
public:
  //Default ctor
  Kraken():m_seq(), m_maps(), m_xc(), m_router(), m_params() {}
  //Ctor to set custom params
  Kraken(const KrakenParams& params):m_seq(), m_maps(), m_xc(), m_router(), m_params(params) {}

  void    setLocalAlignAdjust(bool laa)    { m_params.setLocalAlignAdjust(laa);   }
  void    setOverflowAdjust(bool ofa)      { m_params.setOverflowAdjust(ofa);     } 
  void    setTransSizeLimit(int tsl)       { m_params.setTransSizeLimit(tsl);     }
  void    setPValThresh(double pvt)        { m_params.setPValThresh(pvt);         }
  void    setMinIdent(double mi)           { m_params.setMinIdent(mi);            }
  void    setMinAlignCover( double mac)    { m_params.setMinAlignCover(mac);      } 

  void Allocate(const string & source, const string & target, double distance = 0.5);
  void DoneAlloc();

  void ReadMap(const string & fileName, const string & source, const string & target, double distance = 0.5);
  void ReadGenome(const string & fileName, const string & name);

  int GenomeCount() const                   {return m_seq.isize();  }
  const string & GenomeName(int i) const    {return m_seq[i].Name();}
  const svec<GenomeSeq>& GetGenomes() const {return m_seq;          }
  const GenomeWideMap & GetMap(const string & source) const;
  
  bool Find(const Coordinate & lookup, 
	    const string & source, 
	    const string & target,
            Coordinate& result);

  bool FindWithEdges(const Coordinate& lookup, const string & source,
                     const string & target,
                     int edgeLength, Coordinate& result);

private:
  bool RoughMap(const Coordinate& lookup, const string& source, const string& target,
                DNAVector& sourceSeq, DNAVector& targetSeq, int& maxPos,
                float& maxVal, int& len, Coordinate& result); 
  bool SetSequence(const vecDNAVector& genome, Coordinate& coords, DNAVector& resultSeq);
  bool RoughAlign(DNAVector& target, DNAVector& source, int& maxPos, float& maxVal, int& len, Coordinate& result); 
  void Ccorrelate(const DNAVector& q, const DNAVector& t, double size, float& maxValOut, int& maxPosOut); 
  bool ExhaustAlign(DNAVector& trueDestination, DNAVector& source, int slack, Coordinate& result);
  int  Index(const string & source, const string & target);
  int  Genome(const string & name);
  
  bool MapThroughRoute(const Route & route, svec<Coordinate>& results, const Coordinate & lookup);


  svec<GenomeSeq> m_seq;
  svec<GenomeWideMap> m_maps;

  MultiSizeXCorr m_xc;
  RouteFinder m_router;

  KrakenParams m_params;  /// Set of parameters containing, Threshold of pValue, min ident value and trans limit size among other items
};

#endif // KRAKENMAP_H

