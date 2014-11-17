#ifndef _ANNOTATION_QUERY_H_
#define _ANNOTATION_QUERY_H_

#include <map>
#include <string>
#include <sstream>
#include "base/SVector.h"
#include "base/FileParser.h"
#include "src/AlignmentBlock.h"
#include "src/AnnotationQuery/NCList.h"
#include "src/Coordinate.h"

// Forward declaration 
class Annotation; 

//======================================================
/** Annotation Item's Auxilliary data (key-value pairs) */
class AIAux {
public:
  AIAux(): data(), iter() {}
  AIAux(map<string, string> d): data(d), iter(d.begin()) {}

  int getSize() const { return data.size(); }
  /** given a key and value, add to the dataset */
  void add(const string& key, const string& value) { data[key] = value; }
  /** Get the relevant value for a given key. If doesnt exist return empty string */
  const string& getValue(const string& key) { return data[key]; }
  /** Get the next item */
  bool getNext(string& key, string& value); 
  /** Reset the item iterator */
  void resetIter() { iter = data.begin(); }
  /** Check if iterator has more */
  bool hasMore() const { return iter!=data.end(); }
  /** Return all key/values as a GTF string */
  string toString() const;

private:
  map<string, string> data; /// key-value pairs
  map<string, string>::iterator iter;
};


//======================================================
/**
 * Enumeration specifying the different types of annotation fields
 * AnnotationItem (specified by AITEM:0) 
 * Transcript     (specified by TRANS:1)
 * Gene           (specified by GENE: 2)
 * Locus          (specified by LOCUS:3)
 */
enum AnnotField { AITEM, TRANS, GENE, LOCUS, NONE };


//======================================================
/** Base class for AnnotItem (e.g. each row in a GTF) */
class AnnotItemBase
{
public:
  AnnotItemBase(): coords(), children(), exons(), parent(NULL), transferredCoords(false), sublistIndex(-1) {}
  AnnotItemBase(const Coordinate& crds): coords(crds), children(), exons(), parent(NULL), transferredCoords(false), sublistIndex(-1) {}

  virtual ~AnnotItemBase() {}

  const string & getChr() const                         { return coords.getChr();     }
  void setChr(const string& chr)                        { coords.setChr(chr);         }
  int getStart() const                                  { return coords.getStart();   }
  void setStart(int start)                              { coords.setStart(start);     }
  int getStop() const                                   { return coords.getStop();    }
  void setStop(int stop)                                { coords.setStop(stop);       }
  char getOrient() const                                { return coords.getOrient();  }
  void setOrient(bool orient)                           { coords.setOrient(orient);   }
  const Coordinate& getCoords() const                     { return coords;              }
  void setCoords(const Coordinate& crds)                  { coords = crds;              }
  const svec<AnnotItemBase*>& getChildren() const       { return children;            }
  bool isReversed() const                               { return coords.isReversed(); }
  /** For Transcripts & Genes, this function returns exons belonging to them */
  const svec<AnnotItemBase*>& getExons() const          { return exons;               }
  const AnnotItemBase* getParent() const                { return parent;              }
  /** Note: Needed for modifying parent use the constant version otherwise */
  AnnotItemBase* getParentNC() const                    { return parent;              }   
  void setParentNode(AnnotItemBase* pNode)              { parent = pNode;             }
  bool getTransferred() const                           { return transferredCoords;   }
  void setTransferred(bool flag)                        { transferredCoords = flag;   }
  int  getSublist() const                               { return sublistIndex;        }
  void setSublist(int sli)                              { sublistIndex = sli;         }
  bool hasSublist() const                               { return (sublistIndex!=-1);  }
  virtual bool isSameOrient(AnnotItemBase* other)const  { return coords.isSameOrient(other->getCoords());   } 
  bool contains(const AnnotItemBase& other) const       { return( getCoords().contains(other.getCoords())); }
  /** Checks if the whole of the given item falls within one intronic/intergenic region */
  bool isIntronic(const AnnotItemBase& other) const;
  /** Return the length of the coordinate by subtracting start from stop (absolute value) */
  int findLength() const { return coords.findLength(); }

  // Following functions useful for children
  virtual bool isCodingExon()const                            { return false; }
  virtual AnnotField getType()const                     { return NONE;  } 
  virtual const string getCategory()const               { return "";    }
  virtual const string getParentTransId()const          { return "";    }
  virtual const string getParentGeneId()const           { return "";    }
  virtual const string getId()const                     { return "";    }
  virtual const string getBioType()const                { return "";    }
  /** Given an annotation, the item is compared to it and a type specific report is produced in sout */
  virtual void reportOverlaps(const Annotation& tA, ostream& sout)const { return; }

  bool operator < (const AnnotItemBase& i) const {
    return (getCoords() < i.getCoords()); 
  }

  /** To use object as functor for Annotation data sorting and other searches */
  bool operator()(const AnnotItemBase* a, const AnnotItemBase* b) {
    return (*a < *b);
  }

  /** 
   * Add children that this object is constituted from.
   * Note this is only applicable for the Annotation and
   * Gene objects and not for the AnnotItem object
   */
  void addNode(AnnotItemBase* item) { 
    // Extend the coordinates if necessary
    if(item->getStart() < getStart())  { coords.setStart(item->getStart()); }
    if(item->getStop()  > getStop())   { coords.setStop(item->getStop());   }
    // Add the node 
    children.push_back(item); 
    // Add to exon list if necessary 
    if(item->isCodingExon()) { 
      exons.push_back(item); 
    }
  }

/** 
   * Update the items coordinates based on the given coords
   * Note this is only applicable for items that require translation
   * as part of the TransAnnotation class
   * Returns flag signifying whether the new coordinates were accepted
   */
  bool transCoords(const Coordinate& transCoords);
 
  /** 
   * Creates a string containing the details of the object
   * A separator character is provided to use for separating the fields
   * Note that the id, category, or bioType might be empty, one can 
   * overload function in child classess if need be.
   */
  virtual string toString(char sep) const;
                
 protected:
  Coordinate coords;                 /// Coordinates (start/stop location, orientation, and chromosome)
  svec<AnnotItemBase*> children;   /// list of pointers to Transcripts or AnnotationItems (used only for Gene and Transcript)
  svec<AnnotItemBase*> exons;      /// list of pointers to Exons that this item constitutes of (only applies to Transcripts) 
  AnnotItemBase* parent;           /// Pointer to parent (either a gene for transcripts or a transcript for annoItems)
  bool transferredCoords;          /// Flag identifies if the coordinates of this item  have been transferred into another space (annotation)

private:
  int sublistIndex;                /// The index of a sublist associated with this interval (if any, if not -1)
};

//======================================================
/** 
 * Each annotation item (e.g. each row in a GTF)
 * is represented with one AnnotationItem object
 */
class AnnotItem: public AnnotItemBase
{
public:
  AnnotItem(): AnnotItemBase(), category(""), parentTransId(""), parentGeneId(""), aux() {}

  // Ctor used for creating dummy object for comparison of coordinates
  AnnotItem(Coordinate crds):AnnotItemBase(crds), category(), aux() {}

  AnnotItem(Coordinate crds, const string& ctgry, const string& tId,
            const string& gId, const AIAux& ax)
           :AnnotItemBase(crds), category(ctgry), parentTransId(tId),
            parentGeneId(gId), aux(ax) {}

  virtual const string getCategory()const      { return category;      }
  virtual const string getParentTransId()const { return parentTransId; }
  virtual const string getParentGeneId()const  { return parentGeneId;  }

  void set(const string& ch, const string& cat, const string& tId, 
           const string& gId, int str, int stp, bool  ori) {
    category      = cat;
    parentTransId = tId; 
    parentGeneId  = gId; 
    coords        = Coordinate(ch, ori, str, stp);
  }
  virtual AnnotField getType()const      { return AITEM;} 
  virtual bool isCodingExon()const             { return (getCategory()=="CDS"); }

  /** 
   * Creates a string containing the details of the object
   * in the GTF format (the separator should be a tab for 
   * standard GTF)
   */
  virtual string toString(char sep) const;
  virtual void reportOverlaps(const Annotation& tA, ostream& sout)const; 

 private:
  string   category;       /// Type category (i.e. exon, CDS, etc.)
  string   parentTransId;  /// The Id of the Transcript to which this item belongs to
  string   parentGeneId;   /// The Id of the Gene to which this item belongs to
  AIAux    aux;          /// Key-value pairs for free form data
};

//======================================================
/** 
 * Each transcript which is a collection of GTF 
 * Annotation Items is represented with a Transcript object
 */
class Transcript: public AnnotItemBase
{
public:
  Transcript(): AnnotItemBase(), bioType(), transId() {} 

  Transcript(Coordinate crds, const string& bType,
             const string& tId)
             : AnnotItemBase(crds), bioType(bType),
               transId(tId) {} 

  virtual const string getBioType() const { return bioType; }
  virtual const string getId() const      { return transId; }
  virtual AnnotField getType()const       { return TRANS;   } 

  void set(const string& ch, const string& bType, int str, 
	   int stp, bool  ori, const string& tId) {
    bioType   = bType;
    transId   = tId;
    coords    = Coordinate(ch, ori, str, stp);
  }

  virtual void reportOverlaps(const Annotation& tA, ostream& sout)const;

private:
  string   bioType;      /// Biological class of annotation
  string   transId;      /// Transcript ID
};

//======================================================
/** 
 * Each gene which is a collection of 
 * GTF Items is represented with a Gene object
 */
class Gene: public AnnotItemBase
{
public:
  Gene(): geneId("") {}

  Gene(const Coordinate& crds, const string& gId)
      : AnnotItemBase(crds), geneId(gId) {} 

  const string getId()const          { return geneId; }
  virtual AnnotField getType()const  { return GENE;   } 
  void set(const string& ch, int str, int stp,
           bool  ori, const string& gId) {
    coords   = Coordinate(ch, ori, str, stp);
    geneId   = gId;
  }
  virtual void reportOverlaps(const Annotation& tA, ostream& sout)const;

private:
  string   geneId;       /// Gene ID
};

//======================================================
/** 
  * Each Locus which is a collection of 
  * genes is represented with a Locus object
  */
class Locus: public AnnotItemBase
{
public:
  Locus() {}

  Locus(const Coordinate& crds)
      : AnnotItemBase(crds) {} 

  virtual AnnotField getType()const  { return LOCUS; } 

  void set(const string& ch, int str, int stp,
           bool  ori) {
    coords   = Coordinate(ch, ori, str, stp);
  }
};


//======================================================
/** To use as functor for sorting AITEMS based on geneId, and secondly transcriptId */
struct CompareGeneTransIdLess { 
  bool operator() (const AnnotItemBase* a, const AnnotItemBase* b) const {
    if((*a).getParentGeneId() != (*b).getParentGeneId()) {
      return ((*a).getParentGeneId() < (*b).getParentGeneId());
    } else {
      if((*a).getParentTransId() != (*b).getParentTransId()) {
        return ((*a).getParentTransId() < (*b).getParentTransId());
      } else {
        return (*a < *b); //If both gene/trans id are the same sort on coordinates
      }
    }
  }
};

//======================================================
/**
 * The Annotation class can be used to read a GTF file and create structures 
 * that enable one to query for and establish correspondece among exons, 
 * transcripts, genes, and loci. 
 */
class Annotation {
public:
  Annotation(const string& fileName, const string& specie) {
    readGTF(fileName, specie);
  }

  /** Note: Copying this object is time consuming - to be used only when absolutely necessary */   
  Annotation(const Annotation& other) {
    copy(other);
  }
  void operator=(const Annotation& other) {
    clear();
    copy(other);
  }

  ~Annotation() { clear(); }
  
  const string& getSpecieId() const { return speciesId; }

  /** Mode 0: Annotation Items, Mode 1: Transcripts, Mode 2: Genes, Mode 3: Loci */
  const svec<AnnotItemBase*>& getDataByCoord(AnnotField type) const { 
    switch(type) {
      case AITEM:
        return annotsByCoord;
      case TRANS:
        return transByCoord;
      case GENE:
        return genesByCoord;
      case LOCUS:
      default:
        return lociByCoord;
    }
  }

  /** Only set the loci if needed through calling this function */
  const svec<AnnotItemBase*>& setLoci();

  /** 
   * Return objects overlapping with the given coordinates in the results object passed in by user.
   * By choosing the right mode you can obtain of the following:
   * AITEM:  annotationItems, TRANS: transcripts, GENE: genes, or LOCUS: Loci
   * Returns integer represnting the number of overlapping items found
   */
  int getAnyOverlapping(AnnotItemBase* subject, AnnotField mode, svec<AnnotItemBase*>& results) const;     
  /** Same as other overload but with AICoord as input instead of AnnotItem */
  int getAnyOverlapping(const Coordinate& subject, AnnotField mode, svec<AnnotItemBase*>& results) const;     

  /** 
   * Return objects fully containing the given coordinates in the results object passed in by user.
   * By choosing the right mode you can obtain of the following:
   * AITEM:  annotationItems, TRANS: transcripts, GENE: genes, or LOCUS: Loci
   * Returns integer represnting the number of items found
   */
  int getFullyContained(AnnotItemBase* subject, AnnotField mode, svec<AnnotItemBase*>& results) const;     
  /** Same as other overload but with AICoord as input instead of AnnotItem */
  int getFullyContained(const Coordinate& subject, AnnotField mode, svec<AnnotItemBase*>& results) const;     
 
  /** Returns number of items of the same type as the query where all or some
   * of  their child nodes have overlaps. The overlapMode parameter determines
   * whether those with all children overlapping(overlapMode:3) or only those 
   * with part of children overlapping (overlapMode:2), or alternatively those 
   * with none overlapping as a whole contained within the target -intornic 
   * (overlapMode:1) or not contained within the target - none-exonic
   * (overlapMode:0) should be returned. Also the sameOrient flag chooses
   * whether the given results should be the sameOrient or antisameOrient ones.
   */
  int getExonicOverlaps(Transcript* subject, int overlapMode, bool sameOrient, svec<AnnotItemBase*>& results) const;

  /** Clear all the data content and reset the object */
  void clear();

  /** Write from an annotation object into gtf format */
  virtual void writeGTF(ostream& sout); 

protected:
  /** Function used for copy constructor & assignment operator only
   * Note: Time consuming - to be used only when absolutely necessary 
   */
  void copy(const Annotation& annot);

  /** Read the give GTF file into the rlevant annotation/transcript/gene structures */
  void readGTF(const string& fileName, const string& specie);

  /** Sort the vectors and construct nested containment lists - used by the read or update function */
  void sortSetNCLists();
  void sortAll(); 
  void setNCLists();

  /** 
   * Returns the pointer to object that has been added
   * so that it can be used for links to transcripts
   */
  AnnotItemBase* addAnnotItem(const AnnotItem& i) { 
    AnnotItemBase* newItem = new AnnotItem(i);
    annotsByCoord.push_back(newItem);
    return newItem;
  }

  /** 
   * Returns the pointer to object that has been added
   * so that it can be used for links to genes 
   */
  AnnotItemBase* addTranscript(const Transcript& t) { 
    AnnotItemBase* newItem = new Transcript(t);
    transByCoord.push_back(newItem);
    for(int i=0; i<t.getChildren().isize(); i++) {
      t.getChildren()[i]->setParentNode(newItem);
    }
    return newItem;
  }

  /** 
   * Returns the pointer to object that has been added
   */
  AnnotItemBase* addGene(const Gene& g) { 
    AnnotItemBase* newItem = new Gene(g);
    genesByCoord.push_back(newItem);
    for(int i=0; i<g.getChildren().isize(); i++) {
      g.getChildren()[i]->setParentNode(newItem);
    }
    return newItem;
  }

  /** 
   * Returns the pointer to object that has been added
   */
  AnnotItemBase* addLocus(const Locus& l) { 
    AnnotItemBase* newItem = new Locus(l);
    lociByCoord.push_back(newItem);
    return newItem;
  }

  /** Mode 0: Annotation Items, Mode 1: Transcripts, Mode 2: Genes, Mode 3: Loci */
  const NCList<AnnotItemBase>& getNCListByCoord(AnnotField type) const { 
    switch(type) {
      case AITEM:
        return annotsNCList;
        break;
      case TRANS:
        return transNCList;
        break;
      case GENE:
        return genesNCList;
        break;
      case LOCUS:
      default:
        return lociNCList;
    }
  }

  /** Accessory function used for removing unwanted characters from GTF key-value terms */
  void cleanKeyValue(string& value);


  string speciesId;                    /// The specie to which the annotation belongs 
  svec<AnnotItemBase*>  annotsByCoord; /// Annotation items sorted by the coordinates 
  NCList<AnnotItemBase> annotsNCList;  /// Annotation items None-containment Lists 
  svec<AnnotItemBase*>  transByCoord;  /// Transcripts sorted by the coordinates
  NCList<AnnotItemBase> transNCList;   /// Transcripts None-containment Lists
  svec<AnnotItemBase*>  genesByCoord;  /// Genes sorted by the coordinates
  NCList<AnnotItemBase> genesNCList;   /// Genes None-containment lists
  svec<AnnotItemBase*>  lociByCoord;   /// Loci sorted by the coordinates
  NCList<AnnotItemBase> lociNCList;    /// Loci None-containment lists
};

/**
 */
class GTFCompare
{
public:
  GTFCompare() {} 
  /**Used to report all overlapping items of the qFiledType and 
   * determine whether they have full/partial/none child overlapping 
   * Sense and antisense categories based on whether transcripts share
   * the same orientation and the fullOverlap, partialOverlap of exons.
   * Full overlap is where all exons have some kind of overlap with eachother
   * partial is when some exons have overlaps. Intronic is when a transcript
   * is completely contained within the intronic region of another.
   */
  virtual void reportAllOverlaps(const Annotation& qA, const Annotation& tA, 
                         AnnotField qFieldType, ostream& sout);
};  

#endif //_ANNOTATION_QUERY_H_
