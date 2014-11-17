#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "extern/logger/log.h"
#include "src/Kraken/GTFTransfer.h"

//======================================================

void TransAnnotation::translateCoordinates( const string& targetSpecieId, 
                                           Kraken& mapper, bool localAlign) { 
  if(hasBeenTranslated()) { 
    FILE_LOG(logERROR) <<"Cannot translate coordinates of this annotation \
                          as they have already been translated once"; 
    return;
  }
  const svec<AnnotItemBase*>& annotItems = getDataByCoord(AITEM);
  FILE_LOG(logDEBUG) << "Total annotation items to translate: " << annotItems.size();  
  for (int i=0; i<annotItems.isize(); i++) {
    FILE_LOG(logDEBUG)  << "Translating annotation item: " << i;  
    FILE_LOG(logDEBUG1) << annotItems[i]->toString('\t');  
    Coordinate res;
    bool bOK = mapper.Find(annotItems[i]->getCoords(), this->getTranslateSpace(),
                           targetSpecieId, localAlign, res);
    if (bOK) {
      FILE_LOG(logDEBUG) << "Item was Translated: " << res.toString('\t');  
      updateAnnotItem(annotItems[i], res);  
    }    
  }
  //TODO temp - until sublist has been decoupled
  for(int mode=0; mode<4; mode++) {
    for(svec<AnnotItemBase*>::const_iterator it = getDataByCoord(AnnotField(mode)).begin();
        it != getDataByCoord(AnnotField(mode)).end() ; ++it) {
      (*it)->setSublist(-1);
    }
  }
  FILE_LOG(logDEBUG2) << "Start to Sort lists";
  // Sort and Rearrange the data structures now that the coordinates have been transformed
  sortSetNCLists(); 
  // Set the Id of the space that this annotaion has been translated into
  translateSpace = targetSpecieId; 
}

void TransAnnotation::updateAnnotItem(AnnotItemBase* item, const Coordinate& tCoord) { 
  // Update parent transcript and parent gene
  AnnotItemBase* parentTrans = item->getParentNC();
  AnnotItemBase* parentGene  = parentTrans->getParentNC();
  item->transCoords(tCoord);
  if(parentTrans->transCoords(item->getCoords())) {
    parentGene->transCoords(parentTrans->getCoords()); 
  }
}

void TransAnnotation::writeGTF(ostream& sout) {
  svec<AnnotItemBase*> aItems(getDataByCoord(AITEM)); //Need a copy
  // Sort items based on primarily geneId and secondly transcriptId
  sort(aItems.begin(), aItems.end(), CompareGeneTransIdLess());
  for(svec<AnnotItemBase*>::iterator it = aItems.begin(); 
    it != aItems.end(); ++it) {
    if((*it)->getTransferred()) {
      sout << (*it)->toString('\t') << endl;
    }
  }
}


//======================================================
void GTFTransfer::translate(TransAnnotation& sourceAnnot, const string& targetId, bool localAlign) {
  // Translate the mappings of the source annotaion into targetination annotation space
  sourceAnnot.translateCoordinates(targetId, mapper, localAlign); 
}

void GTFTransfer::reportAllOverlaps(const Annotation& qA, const Annotation& tA, 
                                   AnnotField qFieldType, ostream& sout) {
  const svec<AnnotItemBase*>& annotItems = qA.getDataByCoord(qFieldType);
  FILE_LOG(logDEBUG)<<"---Size of entire set to be mapped: "<<annotItems.size()<<endl;
  for (int i=0; i<annotItems.isize(); i++) {
    FILE_LOG(logDEBUG)<<">Index: "<<i<<" - "<<annotItems[i]->toString('\t');
    if(annotItems[i]->getTransferred()) {
      FILE_LOG(logDEBUG)<<"WAS TRANSLATED";
      annotItems[i]->reportOverlaps(tA, sout);
    }
  }
}

