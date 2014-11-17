#include <string>
#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-q","query GTF File");
  commandArg<string> bStringCmd("-t","target GTF File");

  commandLineParser P(argc,argv);
  P.SetDescription("Parse and read in GTF entries");
  P.registerArg(aStringCmd);
  P.registerArg(bStringCmd);

  P.parse();
  string aString = P.GetStringValueFor(aStringCmd);
  string bString = P.GetStringValueFor(bStringCmd);
 
  cout << "Reading qyery" << endl;
  Annotation annotA(aString, ""); 
  cout << "Reading target" << endl;
  Annotation annotB(bString, "");

  int cnt = 0;
  for (int m = 0; m <3; m++) {
    AnnotField mode = AnnotField(m); // choose AITEM for annotationItems, TRANS for transcripts, and GENE for genes
    svec<AnnotItemBase*> out = annotA.getDataByCoord(mode);
    int n = 0;
    int lap = 0;
    int exact = 0;
    for(svec<AnnotItemBase*>::iterator it = out.begin(); it != out.end() ; ++it) {
      
      svec<AnnotItemBase*> match;
      annotB.getAnyOverlapping(*it, mode, match);
      int l = 0;
      int e = 0;
      for (int i=0; i<match.isize(); i++) {
	AnnotItemBase* b = match[i];
	if (b->getOrient() == (*it)->getCoords().getOrient()) {
	  l++;
	  if ((*it)->getCoords().getStart() == b->getStart() && (*it)->getCoords().getStop() == b->getStop()) {
	    e++;
	  }
	}
      }
      n++;
      if (l > 0)
	lap++;
      if (e > 0)
	exact++;
    }
    cout << "MODE: " << (int)mode << endl;
    cout << "Total: " << n << " Overlap: " << lap << " Exact: " << exact << endl;
  }
  return 0;
}

