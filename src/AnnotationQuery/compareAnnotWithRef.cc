#include <string>
#include <vector>
#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-q","query GTF File");
  commandArg<string> bStringCmd("-t","target GTF File");
  commandArg<int>    fieldTypeCmd("-f","The field type for what will be compared, i.e. 1: transcripts, 2: genes, 3:loci", 1);

  commandLineParser P(argc,argv);
  P.SetDescription("Compare the fields of a given annotation file (query) with a reference annotation");
  P.registerArg(aStringCmd);
  P.registerArg(bStringCmd);
  P.registerArg(fieldTypeCmd);

  P.parse();
  string aString   = P.GetStringValueFor(aStringCmd);
  string bString   = P.GetStringValueFor(bStringCmd);
  int    fieldType = P.GetIntValueFor(fieldTypeCmd);
 
  Annotation annotA(aString, "");
  Annotation annotB(bString, "");

  AnnotField mode = AnnotField(fieldType); // choose AITEM for annotationItems, TRANS for transcripts, and GENE for genes
  svec<AnnotItemBase*> out = annotA.getDataByCoord(mode);
  for(svec<AnnotItemBase*>::iterator it = out.begin(); it != out.end() ; ++it) {
    string detInfo; // Detailed exon, transcript, and gene info
    svec<AnnotItemBase*> genes;
    annotB.getAnyOverlapping(*it, GENE, genes);
    int aItemCount = 0;
    int transCount = 0;
    int geneCount  = 0;
    for (int i=0; i<genes.isize(); i++) {
      if (genes[i]->getOrient() == (*it)->getCoords().getOrient()) {
        geneCount++;
        detInfo.append(genes[i]->toString('\t')); // Use tab for delimiter
        detInfo.append("\n");
        // Get Transcripts & for every transcript get the annotation items
        svec<AnnotItemBase*> trans = genes[i]->getChildren();
        for (int j=0; j<trans.isize(); j++) {
          transCount++;
          detInfo.append("\t");
          detInfo.append(trans[j]->toString('\t')); // Use tab for delimiter
          detInfo.append("\n");
          svec<AnnotItemBase*> aItems = trans[j]->getChildren();
          for (int k=0; k<aItems.isize(); k++) {
            aItemCount++;
            detInfo.append("\t\t");
            detInfo.append(aItems[k]->toString('\t')); // Use tab for delimiter
            detInfo.append("\n");
          }
        }
      }
    }
    cout << endl <<">\t" << (*it)->toString('\t')
         << "\t" << aItemCount << "\t" << transCount 
         << "\t" << geneCount << "\n" << detInfo;
  }
  return 0;
}

