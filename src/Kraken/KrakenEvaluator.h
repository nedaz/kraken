
#ifndef _EVALUATE_RUM_H_
#define _EVALUTATE_RUM_H_

#include <map>
#include "src/Kraken/KrakenMap.h"
#include "src/Kraken/KrakenConfig.h"

/**
 * 
 */
class KrakenEvaluator 
{ 
public:
  /** Default Ctor */
  KrakenEvaluator() {}

  /** 
   * Given two different synteny graphs (i.e. two separatly configured Kraken objects)
   * This function will compare a given percentage of the genomic blocks mapped via the
   * different synteny graphs & report the results.
   */ 
  void compareAllMapGenomes(Kraken& mapper1, Kraken& mapper2, int percentage, int blockSize); 

  /** 
   * Given two different synteny graphs (i.e. two separatly configured Kraken objects)
   * This function will compare the annotation items (e.g. exon, transcript..) across 
   * all other genomes in the config file genomes from a 
   * given annotation file mapped via the synteny graphs & report the results.
   */ 
  void compareAnnotationAllMapGenomes(Kraken& mapper1, Kraken& mapper2, const Annotation& origAnnot, 
                                      const string& origin, AnnotField qFieldType); 
  /** 
   * Given two different synteny graphs (i.e. two separatly configured Kraken objects)
   * This function will compare the annotation items (e.g. exon, transcript..) from a 
   * given annotation file mapped via the synteny graphs & report the results.
   */ 
  void compareAnnotation(Kraken& mapper1, Kraken& mapper2, const Annotation& origAnnot, 
                         const string& origin, const string& destin, AnnotField qFieldType); 
  
private:
  /** 
   * Given two different synteny graphs (i.e. two separatly configured Kraken objects)
   * This function will compare a given origCoords mapping from origin to destin genomes.
   * The operation can be made faster by using the limited edge mapping mode, which only
   * maps the edges of the given region (edge size can be provided or default to 200). 
   */ 
  void compareMapOutputs(Kraken& mapper1, Kraken& mapper2, const Coordinate& origCoords,
                         const string& origin, const string& destin, 
                         bool limitToEdgeMap=false, int edgeSize=200);
};

#endif //_EVALUATE_RUM_H_
