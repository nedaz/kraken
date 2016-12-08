#ifndef _GTF_TRANSFER_H_
#define _GTF_TRANSFER_H_

#include <map>
#include "KrakenMap.h"
#include "KrakenConfig.h"
#include "../annotationQuery/AnnotationQuery.h"


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

  void translateCoordinates(const string& destSpecieId, Kraken& m_mapper); 
  virtual void writeGTF(ostream& sout, bool outputAll); 
  
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
  GTFTransfer(const string& configFile):m_mapper() {
    KrakenConfig config(&m_mapper);
    config.Configure(configFile);
  } 

  void    setLocalAlignAdjust(bool laa)    { m_mapper.setLocalAlignAdjust(laa);   }
  void    setOverflowAdjust(bool ofa)      { m_mapper.setOverflowAdjust(ofa);     } 
  void    setTransSizeLimit(int tsl)       { m_mapper.setTransSizeLimit(tsl);     }
  void    setPValThresh(double pvt)        { m_mapper.setPValThresh(pvt);         }
  void    setMinIdent(double mi)           { m_mapper.setMinIdent(mi);            }
  void    setMinAlignCover( double mac)    { m_mapper.setMinAlignCover(mac);      } 

  virtual void reportAllOverlaps(const Annotation& qA, const Annotation& tA, 
                         AnnotField qFieldType, ostream& sout);

  void translate(TransAnnotation& origAnnot, const string& destId);

private:
  Kraken m_mapper;
};  


#endif //_GTF_TRANSFER_H_
