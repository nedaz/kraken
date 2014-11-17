
#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"
#include "src/Kraken/KrakenEvaluator.h"

int main(int argc,char** argv) {
  commandArg<string> aStringCmmd("-config1","first configuration file");
  commandArg<string> bStringCmmd("-config2","second configuration file");
  commandArg<string> cStringCmmd("-oAnnot","origin annotation GTF file","");
  commandArg<string> dStringCmmd("-oName","origin genome id","");
  commandArg<string> eStringCmmd("-dName","destination genome id","");
  commandArg<string> fStringCmmd("-log","Application logging file","application.log");
 
  commandLineParser P(argc,argv);
  P.SetDescription("Given two rum config files, this program can be run in two differen modes\
                    1. maps blocks from one genome to all others via both configurations and compares.\
                    2. if an annotation file and origin/destination genome ids are given, then the annotation\
                       is mapped to the destination genome via the different mappers and compared");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd); 
  P.registerArg(cStringCmmd);
  P.registerArg(dStringCmmd);
  P.registerArg(eStringCmmd);
  P.registerArg(fStringCmmd);
  P.parse();
  string configFile1     = P.GetStringValueFor(aStringCmmd);
  string configFile2     = P.GetStringValueFor(bStringCmmd);
  string oAnnotFile      = P.GetStringValueFor(cStringCmmd);
  string origId          = P.GetStringValueFor(dStringCmmd);
  string destId          = P.GetStringValueFor(eStringCmmd);
  string applicationFile = P.GetStringValueFor(fStringCmmd);
 
  FILE* pFile = fopen(applicationFile.c_str(), "w");
  Output2FILE::Stream() = pFile;
  FILELog::ReportingLevel() = logINFO; //logDEBUG3; //logINFO;
 
  KrakenEvaluator rEval;
  Kraken mapper1, mapper2;
  KrakenConfig config1(&mapper1);
  config1.Configure(configFile1);
  KrakenConfig config2(&mapper2);
  config2.Configure(configFile2);
  if(oAnnotFile =="") {
    rEval.compareAllMapGenomes(mapper1, mapper2, 10, 200);
  } else if(destId=="") { 
    Annotation origAnnot = Annotation(oAnnotFile, origId);
    rEval.compareAnnotationAllMapGenomes(mapper1, mapper2, origAnnot, origId, AITEM);
  } else {
    Annotation origAnnot = Annotation(oAnnotFile, origId);
    rEval.compareAnnotation(mapper1, mapper2, origAnnot, origId, destId, AITEM);
  }
}

//======================================================

void KrakenEvaluator::compareAllMapGenomes(Kraken& mapper1, Kraken& mapper2, int percentage, int blockSize) {
  // Go through all the genomes in mapper1, split into blocks
  // and map over to all other genomes via both mappers
  svec<GenomeSeq> genomes = mapper1.GetGenomes();
  for(int i=0; i<genomes.isize(); i++) {
    string origin           = genomes[i].Name();
    vecDNAVector origGenome = genomes[i].DNA();
    for(int j=i+1; j<genomes.isize(); j++) {
      string destin           = genomes[j].Name();
      FILE_LOG(logINFO)  << "Mapping blocks from genome: " << origin << " to: " << destin; 
      for(int x=0; x<origGenome.isize(); x++) {
        if(origGenome[x].size() < 5000000) { continue; }  // Only choose blocks from larger scaffolds
        for(int y=0; y<origGenome[x].isize()-blockSize*2; y+=blockSize) {
          if(rand()%(100)>=percentage) { continue; } //Only analyse a given percentage of the dataset
          string chrName = origGenome[x].getName();
          chrName.erase(std::remove(chrName.begin(), chrName.end(), '>'), chrName.end());
          Coordinate origCoords(chrName, true, y, y+blockSize-1);
          if(origCoords.findLength()>1000) {
            compareMapOutputs(mapper1, mapper2, origCoords, origin, destin, true, 100);
          } else {
            compareMapOutputs(mapper1, mapper2, origCoords, origin, destin, false);
          }
        }
      }
    }
  }
}

void KrakenEvaluator::compareAnnotationAllMapGenomes(Kraken& mapper1, Kraken& mapper2, const Annotation& origAnnot, 
                                     const string& origin, AnnotField qFieldType) {
  svec<GenomeSeq> genomes = mapper1.GetGenomes();
  for(int i=0; i<genomes.isize(); i++) {
    string destin = genomes[i].Name();
    compareAnnotation(mapper1, mapper2, origAnnot, origin, destin, qFieldType); 
  }
}


void KrakenEvaluator::compareAnnotation(Kraken& mapper1, Kraken& mapper2, const Annotation& origAnnot, 
                                     const string& origin, const string& destin, AnnotField qFieldType) {
  // Translate the mappings of the origin annotaion into destination annotation space
  svec<AnnotItemBase*> annotItems = origAnnot.getDataByCoord(qFieldType);
  FILE_LOG(logINFO) <<">>Size of entire set to be mapped: "<<annotItems.size()<<endl;
  FILE_LOG(logINFO) << "Mapping blocks from genome: " << origin << " to: " << destin; 
  for (int i=0; i<annotItems.isize(); i++) {
    if(rand()%100!=1 || !annotItems[i]->isCodingExon()) { continue; } //Only analyse 1% of the entire set
    //if(annotItems[i]->getCategory()!="exon") { continue; } //Only analyse 1% of the entire set
    FILE_LOG(logINFO) << annotItems[i]->getCoords().toString(',');  
    if(annotItems[i]->findLength()>1000) {
      compareMapOutputs(mapper1, mapper2, annotItems[i]->getCoords(), origin, destin, true, 100);
    } else {
      compareMapOutputs(mapper1, mapper2, annotItems[i]->getCoords(), origin, destin, false);
    }
  }
}

void KrakenEvaluator::compareMapOutputs(Kraken& mapper1, Kraken& mapper2, const Coordinate& origCoords,
                                     const string& origin, const string& destin, 
                                     bool limitToEdgeMap, int edgeSize) {
  Coordinate result1, result2;
  bool found1, found2;
  if(limitToEdgeMap) {
    found1 = mapper1.FindWithEdges(origCoords, origin, destin, false, edgeSize, result1);
    found2 = mapper2.FindWithEdges(origCoords, origin, destin, false, edgeSize, result2);
  } else {
    found1 = mapper1.Find(origCoords, origin, destin, true, result1);
    found2 = mapper2.Find(origCoords, origin, destin, true, result2);
  }

  //if((found1 && destGenome(result1.getChr()).size()<1000000) || 
  //   (found2 && destGenome(result2.getChr()).size()<1000000)) { 
  //  FILE_LOG(logINFO)  << "RESULT: MAPPED-SMALL-SCAFFOLD"; 
  //  return;
  //}

  if(!found1 || !found2) {
    FILE_LOG(logDEBUG1) << "Origin Coords: " << origCoords.toString('\t'); 
    if(found1) { 
      FILE_LOG(logDEBUG1) << "Mapping 1: " << result1.toString('\t'); 
      FILE_LOG(logINFO)  << "RESULT: FIRST-MAPPED-ONLY"; 
    } else if(found2) { 
      FILE_LOG(logDEBUG1) << "Mapping 2: " << result2.toString('\t'); 
      FILE_LOG(logINFO)  << "RESULT: SECOND-MAPPED-ONLY"; 
    } else {
      FILE_LOG(logINFO)  << "RESULT: NONE-MAPPED"; 
    }
  } else { //find1 && find2
    // output results
    FILE_LOG(logDEBUG1) << "Origin Coords: " << origCoords.toString('\t') 
                        << " Mapped on to: " << result1.toString('\t') 
                        << " and " << result2.toString('\t');
    if(result1.isSameCoords(result2)) {
      FILE_LOG(logINFO)  << "RESULT: SAME";
    } else if(result1.hasOverlap(result2)) {
      FILE_LOG(logINFO)  << "RESULT: OVERLAP-TRUE";
    } else {
      FILE_LOG(logINFO)  << "RESULT: OVERLAP-FALSE";
    } 
  }
}
