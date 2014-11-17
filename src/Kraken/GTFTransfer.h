#ifndef _GTF_TRANSFER_H_
#define _GTF_TRANSFER_H_

#include <map>
#include "src/Kraken/KrakenMap.h"
#include "src/Kraken/KrakenConfig.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


class TransAnnotation: public Annotation 
{ 
public:
  /** Construct by reading from fileName and setting speciesId */
  TransAnnotation(const string& fileName, const string& specie):Annotation(fileName, specie), translateSpace(specie) {}

  /** Construct from an Annotation object */
  TransAnnotation(const Annotation& annot):Annotation(annot) {
    translateSpace = annot.getSpecieId();
  }

  /** Note: Copying this object is time consuming - to be used only when absolutely necessary */   
  TransAnnotation(const TransAnnotation& other)
    :Annotation(static_cast<Annotation>(other)) { // Use the inherited copy function to copy the rest of the fields
    translateSpace = other.translateSpace;
  }

  void operator=(const Annotation& other) {
    clear();
    translateSpace = other.getSpecieId();
    // Use the inherited copy function to copy the rest of the fields
    copy(static_cast<Annotation>(other)); 
  }

  /** The coordinates of this annotation are translated into the space of the given destination specie destSpecieId*/
  const string& getTranslateSpace()  { return translateSpace;                }
  bool hasBeenTranslated()           { return (translateSpace != speciesId); } 

  void translateCoordinates(const string& destSpecieId, Kraken& mapper, bool localAlign); 
  virtual void writeGTF(ostream& sout); 
  
private:
  /**   */
  void updateAnnotItem(AnnotItemBase* item, const Coordinate& tCoord); 

  string translateSpace; /// The genomeId of the annotaion to which this  has been translated to 
};

/**
 * This class is based on GTFCompare and is used for comparing two annotations
 * The difference is that one should use this class when the annotations belong
 * to two different genomes and the origin annotation needs to be first translated into
 * the the space of the destination annotation before any comparison can be done. Hence,
 * the class uses an instance of Kraken, for which the config is used in object construction.
 */
class GTFTransfer:public GTFCompare
{
public:
  GTFTransfer(const string& configFile):mapper() {
    KrakenConfig config(&mapper);
    config.Configure(configFile);
  } 

  void SetPValThresh(double p)     { mapper.SetPValThresh(p);     }
  void SetMinIdent(double d)       { mapper.SetMinIdent(d);       }
  void SetTransSizeLimit(double l) { mapper.SetTransSizeLimit(l); }

  virtual void reportAllOverlaps(const Annotation& qA, const Annotation& tA, 
                         AnnotField qFieldType, ostream& sout);

  void translate(TransAnnotation& origAnnot, const string& destId, bool localAlign);

private:
  Kraken mapper;
};  


#endif //_GTF_TRANSFER_H_
