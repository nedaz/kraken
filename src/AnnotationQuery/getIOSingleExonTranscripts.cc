#include <string>
#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-f","GTF File");

  commandLineParser P(argc,argv);
  P.SetDescription("Parse and read in GTF entries, output transcript ids of entries that \
                    are single exon and entirely overlap with other transcript introns");
  P.registerArg(aStringCmd);

  P.parse();
  string aString = P.GetStringValueFor(aStringCmd);

  Annotation annot(aString, "");

  int count = 0; // Count number of output entries
  svec<AnnotItemBase*> allTrans = annot.getDataByCoord(TRANS);
  for(svec<AnnotItemBase*>::iterator it = allTrans.begin(); it != allTrans.end() ; ++it) {
    if((*it)->getChildren().size()==1) { // Single exon transcript
      // Find all the ones that do not overlap with any exons & are not intergeneic (i.e. Overlap with a gene)
      svec<AnnotItemBase*> tempOut1;
      svec<AnnotItemBase*> tempOut2;
      if(annot.getAnyOverlapping(*it, AITEM, tempOut1)-1 == 0 && annot.getAnyOverlapping(*it, GENE, tempOut2)-1 > 0)  {
        cout<<(*it)->getId()<<endl; 
        count++; 
      }
    }
  }
  cout<<count<<" Out of a total of "<<allTrans.size()<<" transcripts have been detected."<<endl;
  return 0;
}

