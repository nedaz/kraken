#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "extern/logger/log.h"
#include "src/AnnotationQuery/AnnotationQuery.h"



//======================================================

string AnnotItemBase::toString(char sep) const {
  stringstream outStream;
  outStream << getId() << sep << getChr() << sep << getStart() << sep 
            << getStop() << sep << getOrient() << sep << getCategory() 
            << sep << getBioType();
  return outStream.str();
} 

bool AnnotItemBase::transCoords(const Coordinate& transCoords) {
  if(!getTransferred()) { //First time being transferred
    setCoords(transCoords);
    setTransferred(true);
    return true;
  } else {
    //TODO The size limit is arbitrary, this should be thought through
    if(transCoords.getChr()==getChr() && (abs(transCoords.getStop()-getStop())<500000)) {
      // Extend the coordinates if necessary
      if(transCoords.getStart() < getStart())  { setStart(transCoords.getStart()); }
      if(transCoords.getStop()  > getStop())   { setStop(transCoords.getStop());   }
      return true;
    } else {
      FILE_LOG(logDEBUG2) << toString('\t') << ": was not updated";
      FILE_LOG(logDEBUG2) << transCoords.toString('\t') << ": Chromosomes did not match or bases were too far apart";
    }
  }
  return false;
}
 
bool AnnotItemBase::isIntronic(const AnnotItemBase& other) const {
  // This function cannot be used for AITEM type as it has no children & intronic region is undefined
  if(other.getType()==AITEM || getChildren().isize()==0) { return false; }
  Coordinate rollingCoords(other.getChildren()[0]->getCoords());
  for(int i=0; i<other.getChildren().isize()-1; i++) {
    rollingCoords.setStart(other.getChildren()[i]->getStop()+1);
    rollingCoords.setStop(other.getChildren()[i+1]->getStart()-1);
    if(rollingCoords.contains(getCoords())) { return true; }
  }
  return false;
}

//======================================================
/** Used in GTF writer, important to keep correct formatting */
string AnnotItem::toString(char sep) const {
  stringstream outStream;
  string bioType  = getParent()->getBioType();
  string transId  = getParentTransId();
  string geneId   = getParentGeneId();
  // Coordinates are published in 1-based (GTF format)
  char fill = '.';
  outStream << getChr() << sep << bioType << sep << getCategory() << sep 
            << getStart()+1 << sep << getStop()+1
            << sep << fill << sep << getOrient()  << sep << fill << sep 
            << "gene_id \"" << geneId << "\"; "
            << "transcript_id \"" << transId << "\"; "  << aux.toString(); 
  return outStream.str();
} 

void AnnotItem::reportOverlaps(const Annotation& tA, ostream& sout)const {
  // Identify six different types of overlap
  svec<AnnotItemBase*> overlaps; 
  int overlapCount = tA.getAnyOverlapping(this->getCoords(), getType(), overlaps);
  char delim = '\t';
  sout<< ">>"<< delim << toString('\t') << delim;
  string tag       = "NONE";
  if(overlapCount>0) { tag = "OVERLAP-PARTIAL"; }
  else { sout << endl; return; }
  int  maxOverlap = 0;
  int  minDiff    = 100000000;
  int  olIndex    = 0; //Overlapping sequence index saved to find the ratio later
  FILE_LOG(logDEBUG2) << "Analysing item: " << toString('\t');
  for (int ol=0; ol<overlapCount; ol++) {
    FILE_LOG(logDEBUG3) << "Testing overlap: " << overlaps[ol]->toString('\t');
    if(overlaps[ol]->getCoords().isSameCoords(getCoords())) { 
      tag = "OVERLAP-EXACT"; 
      maxOverlap = findLength();
      olIndex = ol; 
      break;
    }
    int totDiff = abs(getStart()-overlaps[ol]->getStart())
                + abs(overlaps[ol]->getStop()-getStop());  
    int overlap = getCoords().findOverlapCnt(overlaps[ol]->getCoords()); 
    if(totDiff < minDiff) { 
      minDiff    = totDiff;
      maxOverlap = overlap; 
      olIndex    = ol; 
    }
    FILE_LOG(logDEBUG3) << "Testing overlap: " << totDiff << " " << overlap << " " << olIndex << endl; 
  }
  FILE_LOG(logDEBUG1) << "Chosen: " << overlaps[olIndex]->toString('\t');
  sout << delim << tag << delim << overlaps[olIndex]->getParentTransId()
       << delim << overlaps[olIndex]->getCategory() << delim
       << overlaps[olIndex]->getStart() - getStart() << delim 
       << getStop() - overlaps[olIndex]->getStop() << endl;
}

//======================================================
void Transcript::reportOverlaps(const Annotation& tA, ostream& sout)const {
  // Identify six different types of overlap
  svec<AnnotItemBase*> fullSense; 
  svec<AnnotItemBase*> fullAnti; 
  svec<AnnotItemBase*> partialSense;
  svec<AnnotItemBase*> partialAnti;
  svec<AnnotItemBase*> intronicSense; 
  svec<AnnotItemBase*> intronicAnti; 
  svec<AnnotItemBase*> otherSense; 
  svec<AnnotItemBase*> otherAnti; 

  svec<AnnotItemBase*> anyOverlaps;
  tA.getAnyOverlapping(getCoords(), getType(), anyOverlaps); 
  // 1. get all children and keep counts based on which parent they come from
  // This needs to be done as nodes have records of their parents but not their children
  map<const AnnotItemBase*, int> counts;
  //For all of the children in subject find the children/parent of destination genome that they overlap with.
  //For items that have exons as children only use exons
  svec<AnnotItemBase*> subjectChildren = getExons();
  if(subjectChildren.size()==0) { subjectChildren = getChildren(); }
  for (int ch=0; ch<subjectChildren.isize(); ch++) {
    svec<AnnotItemBase*> oChildren;
    // Code relies on AnnotField orders and AITEM (child of) TRANS (child of) GENE...
    tA.getAnyOverlapping(subjectChildren[ch]->getCoords(), static_cast<AnnotField>(getType()-1), oChildren); 
    for (int j=0; j<oChildren.isize(); j++) {
      // If AnnotationItems are AITEM don't count if not an exon
      if(oChildren[j]->getType()==AITEM && !oChildren[j]->isCodingExon()) { continue; } 
      counts[oChildren[j]->getParent()]++;
    }
  }
  // 2. compare the counts for each parent with all the number of children of the 
  // potential items. If this number is the same all children are overlapping
  for (int k=0; k<anyOverlaps.isize(); k++) {
    map<const AnnotItemBase*, int>::iterator it = counts.find(anyOverlaps[k]);
    if(it == counts.end()) {
      if(!isIntronic(*anyOverlaps[k])) {
        if(isSameOrient(anyOverlaps[k])) { otherSense.push_back(anyOverlaps[k]); }
        else { otherAnti.push_back(anyOverlaps[k]); }
      } else {
        if(isSameOrient(anyOverlaps[k])) { intronicSense.push_back(anyOverlaps[k]); }
        else { intronicAnti.push_back(anyOverlaps[k]); }
      }
    } else { //it!=counts.end()
      if(it->second == anyOverlaps[k]->getChildren().isize() && it->second == subjectChildren.isize()) {
        if(isSameOrient(anyOverlaps[k])) { fullSense.push_back(anyOverlaps[k]); }
        else { fullAnti.push_back(anyOverlaps[k]); }
      } else {
        if(isSameOrient(anyOverlaps[k])) { partialSense.push_back(anyOverlaps[k]); }
        else { partialAnti.push_back(anyOverlaps[k]); }
      }
    }
  }

  char delim = '\t';
  sout << "## " << toString(delim) << endl;
  string parentId = ".";
  if(getParent()) { parentId = getParent()->getId(); }
  sout<< ">> " << getId() << delim << parentId << delim
      << fullSense.size() << delim << fullAnti.size() << delim << partialSense.size() << delim 
      << partialAnti.size() << delim << intronicSense.size() << delim << intronicAnti.size() 
      << delim << otherSense.size() << delim << otherAnti.size() <<endl; 

  sout<<"FULL_SENSE"<<delim;
  for (int fs=0; fs<fullSense.isize(); fs++) {
    sout<<fullSense[fs]->getId()<<delim;
  }
  sout<<endl<<"FULL_ANTI"<<delim;
  for (int fa=0; fa<fullAnti.isize(); fa++) {
    sout<<fullAnti[fa]->getId()<<delim;
  }
  sout<<endl<<"PARTIAL_SENSE"<<delim;
  for (int ps=0; ps<partialSense.isize(); ps++) {
    sout<<partialSense[ps]->getId()<<delim;
  }
  sout<<endl<<"PARTIAL_ANTI"<<delim;
  for (int pa=0; pa<partialAnti.isize(); pa++) {
    sout<<partialAnti[pa]->getId()<<delim;
  }
  sout<<endl<<"INTRONIC_SENSE"<<delim;
  for (int is=0; is<intronicSense.isize(); is++) {
    sout<<intronicSense[is]->getId()<<delim;
  }
  sout<<endl<<"INTRONIC_ANTI"<<delim;
  for (int ia=0; ia<intronicAnti.isize(); ia++) {
    sout<<intronicAnti[ia]->getId()<<delim;
  }
  sout<<endl<<"OTHER_SENSE"<<delim;
  for (int nes=0; nes<otherSense.isize(); nes++) {
    sout<<otherSense[nes]->getId()<<delim;
  }
  sout<<endl<<"OTHER_ANTI"<<delim;
  for (int nea=0; nea<otherAnti.isize(); nea++) {
    sout<<otherAnti[nea]->getId()<<delim;
  }

  sout<<endl<<endl;
}

//======================================================
void Gene::reportOverlaps(const Annotation& tA, ostream& sout)const {
  svec<AnnotItemBase*> overlaps; 
  int overlapCount = tA.getAnyOverlapping(this->getCoords(), getType(), overlaps);
  char delim = '\t';
  sout<< ">>"<< delim << toString('\t') << delim;
  string tag       = "NONE";
  if(overlapCount>0) { tag = "OVERLAP"; }
  sout << tag << endl;  
}

//======================================================
bool AIAux::getNext(string& key, string& value) { 
  iter++;
  if(hasMore()) {
    key   = iter->first;
    value = iter->second; 
    return true;
  } else {
    return false;
  }
}

string AIAux::toString() const{
  stringstream outStream;
  string key, value;
  map<string, string>::const_iterator iter = data.begin();
  for(; iter!=data.end(); iter++) {
    outStream << iter->first <<  " \"" << iter->second << "\"; "; 
  }
  return outStream.str();
}

//======================================================

void Annotation::readGTF(const string& fileName, const string& specie) {

  speciesId = specie; //Set the specie name/Id
  
  FlatFileParser parser;
  parser.Open(fileName);
  parser.ParseLine();

  // Use to track new gene/transcript
  Gene curr_gene          = Gene();
  Transcript curr_trans   = Transcript(); 
  bool orient             = false;
  int start               = 0;
   int stop               = 0;
  AnnotItemBase* aItem    = NULL;
  string chr, category, geneId="", transId="";
  while (parser.ParseLine()) {
    int itemCount = parser.GetItemCount();
    if(itemCount<9 ) { 
     if(itemCount!=0) {
       FILE_LOG(logWARNING)<<"GTF file inconsistency"; 
     }
     continue;
    } 
    orient    = (parser.AsString(6)=="+")?true:false;
    start     = parser.AsInt(3) - 1; //1-based GTF, internal 0-based
    stop      = parser.AsInt(4) - 1;
    chr       = parser.AsString(0);
    category  = parser.AsString(2);
    AIAux aux; // Get all the key value pairs starting from index 9 
    string key, value;
    for (int i=8; i+1 < itemCount; i+=2) {
      key   = parser.AsString(i);
      value = parser.AsString(i+1);
      cleanKeyValue(value); // Remove extra characters
      if(key=="transcript_id") { transId = value; }
      else if(key=="gene_id")  { geneId  = value; }
      else { aux.add(key, value); }
    }

    // 1. New Annotation Item
    Coordinate crds = Coordinate(chr, orient, start, stop);
    aItem         = addAnnotItem(AnnotItem(crds, category, transId, geneId, aux));
    if(!aItem) { 
      FILE_LOG(logWARNING) << "Could not allocate annotation item: " << crds.toString('\t');
      continue; 
    }

    // 2. New Transcript
    if(curr_trans.getId() != transId) { 
      if(curr_trans.getId()!="") { // Don't add if curr_trans has not been set yet
        AnnotItemBase* i = addTranscript(curr_trans);
        curr_gene.addNode(i);
      }

      const string& bioType = parser.AsString(1);
      curr_trans = Transcript(Coordinate(chr, orient, start, stop),
                              bioType, transId);
      curr_trans.addNode(aItem); 
    } else {
      curr_trans.addNode(aItem); 
    }

    // 3. New Gene
    if(curr_gene.getId() != geneId) { 
      if(curr_gene.getId() != "") { // Don't add if curr_gene has not been set yet
        addGene(curr_gene);
      }
      curr_gene = Gene(Coordinate(chr, orient, start, stop), geneId);
    } 
  }

  // Last gene and transcript after loop ended need to be added
  curr_gene.addNode(addTranscript(curr_trans));
  addGene(curr_gene);

  // Sort the vectors and construct nested containment lists
  sortSetNCLists();
}

void Annotation::writeGTF(ostream& sout) {
  svec<AnnotItemBase*> aItems(getDataByCoord(AITEM)); //Need a copy
  // Sort items based on primarily geneId and secondly transcriptId
  sort(aItems.begin(), aItems.end(), CompareGeneTransIdLess());
  for(svec<AnnotItemBase*>::iterator it = aItems.begin(); 
      it != aItems.end(); ++it) {
    sout << (*it)->toString('\t') << endl;
  }
}

void Annotation::copy(const Annotation& annot) {
  speciesId = annot.speciesId;
  // Copy dynamically allocated memory
  for(svec<AnnotItemBase*>::const_iterator it = annot.annotsByCoord.begin(); 
  it != annot.annotsByCoord.end(); ++it) {
    addAnnotItem(*dynamic_cast<AnnotItem*>(*it));
  }
  for(svec<AnnotItemBase*>::const_iterator it = annot.transByCoord.begin(); 
  it != annot.transByCoord.end(); ++it) {
    addTranscript(*dynamic_cast<Transcript*>(*it));
  }
  for(svec<AnnotItemBase*>::const_iterator it = annot.genesByCoord.begin(); 
  it != annot.genesByCoord.end(); ++it) {
    addGene(*dynamic_cast<Gene*>(*it));
  }
  for(svec<AnnotItemBase*>::const_iterator it = annot.lociByCoord.begin(); 
  it != annot.lociByCoord.end(); ++it) {
    addLocus(*dynamic_cast<Locus*>(*it));

  }
  setNCLists();
}

// Sort the vectors and construct nested containment lists
void Annotation::sortSetNCLists() {
  sortAll();
  setNCLists();
}
 
void Annotation::sortAll() {
  sort(annotsByCoord.begin(), annotsByCoord.end(), AnnotItemBase());
  sort(transByCoord.begin(),  transByCoord.end(), AnnotItemBase());
  sort(genesByCoord.begin(),  genesByCoord.end(), AnnotItemBase());
}

void Annotation::setNCLists() {
  annotsNCList.constructSublists(annotsByCoord);
  transNCList.constructSublists(transByCoord);
  genesNCList.constructSublists(genesByCoord);
}
 
const svec<AnnotItemBase*>& Annotation::setLoci() {
  map<AnnotItemBase*, int> usedItems; // Track genes used in loci and number of times
  const svec<AnnotItemBase*>& genes = getDataByCoord(GENE);
  for(svec<AnnotItemBase*>::const_iterator it = genes.begin(); it != genes.end() ; ++it) {
    //Has not been used in any prior loci before
    if(usedItems.find(*it) == usedItems.end()) { 
      Locus loc((*it)->getCoords());
      // Find all the overlapping genes and iterate through them to add to the locus and flag them as used
      svec<AnnotItemBase*> overlapGenes;
      getAnyOverlapping(*it, GENE, overlapGenes);
      for(svec<AnnotItemBase*>::iterator it2 = overlapGenes.begin(); it2 != overlapGenes.end() ; ++it2) {
        loc.addNode(*it2);
        if(usedItems.find(*it2) == usedItems.end()) { usedItems[*it2] = 1; }  
        else { usedItems[*it2]++ ; } 
      }
      addLocus(loc);
    }
  }
  lociNCList.constructSublists(lociByCoord);
  return lociByCoord;
}

int Annotation::getAnyOverlapping(AnnotItemBase* subject,
AnnotField mode, svec<AnnotItemBase*>& results) const {
  getNCListByCoord(mode).getAnyOverlapping(subject, results); // Return results
  return results.isize();
}

int Annotation::getAnyOverlapping(const Coordinate& subject, 
AnnotField mode, svec<AnnotItemBase*>& results) const {
  AnnotItemBase tempItem(subject); // Make an object with the coordinates to pass onto NClist
  getNCListByCoord(mode).getAnyOverlapping(&tempItem, results); // Return results
  return results.isize();
}

int Annotation::getFullyContained(AnnotItemBase* subject, 
AnnotField mode, svec<AnnotItemBase*>& results) const {
  svec<AnnotItemBase*>::const_iterator lBound = lower_bound( getDataByCoord(mode).begin(), 
                                                             getDataByCoord(mode).end(), 
                                                             subject, AnnotItemBase());
  // Scan forward from the lBound to find the overlapping items 
  for(svec<AnnotItemBase*>::const_iterator it = lBound; it != getDataByCoord(mode).end() ; ++it) {
    if(subject->getCoords().contains((*it)->getCoords())) { results.push_back(*it); }
    else { break; } // Gone beyond any possible overlap
  }
  // Scan backward from the lBound to find the overlapping items
  svec<AnnotItemBase*>::const_reverse_iterator rlBound(lBound);
  for(svec<AnnotItemBase*>::const_reverse_iterator rit = rlBound; 
  rit != getDataByCoord(mode).rend() ; ++rit) {
    if(subject->getCoords().contains((*rit)->getCoords())) { results.push_back(*rit); }
    else { break; } // Gone beyond any possible overlap
  }
  return results.isize();
}

int Annotation::getFullyContained(const Coordinate& subject, 
AnnotField mode, svec<AnnotItemBase*>& results) const {
  AnnotItemBase tempItem(subject); // Make an object with the coordinates to pass onto NClist
  return getFullyContained(&tempItem, mode, results); // Return size of results 
}

void Annotation::clear() {
  //Dealocate dynamically allocated memory
  for(int mode=0; mode<4; mode++) {
    for(svec<AnnotItemBase*>::const_iterator it = getDataByCoord(AnnotField(mode)).begin();
        it != getDataByCoord(AnnotField(mode)).end() ; ++it) {
      delete (*it);
    }
  }
  annotsByCoord.clear(); 
  annotsNCList.clear();
  transByCoord.clear();
  transNCList.clear();  
  genesByCoord.clear(); 
  genesNCList.clear();
  lociByCoord.clear();
  lociNCList.clear(); 
}

void Annotation::cleanKeyValue(string& value) {
    value.erase(std::remove(value.begin(), value.end(), ';'), value.end());
    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
}

//======================================================
void GTFCompare::reportAllOverlaps(const Annotation& qA, const Annotation& tA, 
                                   AnnotField qFieldType, ostream& sout) {
  const svec<AnnotItemBase*>& annotItems = qA.getDataByCoord(qFieldType);
  FILE_LOG(logDEBUG)<<"---Size of entire set to be mapped: "<<annotItems.size()<<endl;
  for (int i=0; i<annotItems.isize(); i++) {
    FILE_LOG(logDEBUG)<<">Index: "<<i<<" - "<<annotItems[i]->toString('\t');
    annotItems[i]->reportOverlaps(tA, sout);
  }
}

