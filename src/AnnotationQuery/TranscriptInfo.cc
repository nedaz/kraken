#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "src/AnnotationQuery/AnnotationQuery.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input GTF file");
  commandLineParser P(argc,argv);
  P.SetDescription("Prints info for transcripts, exons, loci etc.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  
  
  Annotation a(fileName, "");


  const svec<AnnotItemBase*> & d = a.getDataByCoord(TRANS);

  int i, j;
  for (i=0; i<d.isize(); i++) {
    AnnotItemBase* p = d[i];
    const string & bt = p->getBioType();
    const string & id = p->getId();

    const svec<AnnotItemBase*> & exons = p->getChildren();
    
    cout << id << "\t" << bt << "\texons: " << exons.isize();
    const Coordinate& coord = p->getCoords();
    svec<AnnotItemBase*> locus;
    a.getAnyOverlapping(p, GENE, locus);     
    //if (locus.isize() == 1) {
    //  cout << "\tstandalone" << endl;
    //  continue;
    //}
    int contains = 0;
    int contained = 0;
    int antisense = 0;
    int overlap = 0;

    cout << "\ttranscripts in locus: " << locus.isize();
    for (j=0; j<locus.isize(); j++) {
      // Don't do self-matches!
      const svec<AnnotItemBase*> & trans = locus[j]->getChildren();
      bool bSelf = false;
      for (int k=0; k<trans.isize(); k++) {
	if (trans[k] == p)
	  bSelf = true;
      }
      if (bSelf) {
	//cout << "Found myself!!" << endl;
	continue;
      }


      bool bContained = false;
      bool bContainer = false;
      bool bOverlap = false;
      bool bAntisense = false;
      const AnnotItemBase * t = locus[j];
      const Coordinate& c = t->getCoords();
      if (coord.getStart() < c.getStart() && coord.getStop() > c.getStop())
	bContainer = true;
      if (coord.getStart() > c.getStart() && coord.getStop() < c.getStop())
	bContained = true;
      if (coord.getOrient() != c.getOrient())
	bAntisense = true;
      if (!bContainer && !bContained)
	bOverlap = true;

      if (bContained)
	contained++;
      if (bContainer)
	contains++;
      if (bAntisense)
	antisense++;
      if (bOverlap)
	overlap++;
    }
    cout << "\tcontains: " << contains;
    cout << "\tcontained by: " << contained;
    cout << "\tantisense to: " << antisense;
    cout << "\toverlaps with: " << overlap;
    cout << endl;
  }


  return 0;
}
