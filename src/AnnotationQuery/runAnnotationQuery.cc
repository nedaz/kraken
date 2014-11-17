#include <string>
#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-f","GTF File");

  commandLineParser P(argc,argv);
  P.SetDescription("Parse and read in GTF entries");
  P.registerArg(aStringCmd);

  P.parse();
  string aString = P.GetStringValueFor(aStringCmd);

  Annotation annot(aString, "");

   cout<<"Number of annotation items: "<<annot.getDataByCoord(AITEM).size()<<endl;
   cout<<"Number of transcripts: "     <<annot.getDataByCoord(TRANS).size()<<endl;
   cout<<"Number of genes: "           <<annot.getDataByCoord(GENE).size()<<endl;

/*
  int cnt = 0;
  AnnotField mode = GENE; // choose AITEM for annotationItems, TRANS for transcripts, and GENE for genes
  svec<AnnotItemBase*> out = annot.getDataByCoord(mode);
  for(svec<AnnotItemBase*>::iterator it = out.begin(); it != out.end() ; ++it) {
    cout<<cnt++<<": overlapping count: ";
    svec<AnnotItemBase*> tempOut;
    cout << annot.getAnyOverlapping((*it), mode, tempOut) << endl;
  }
*/
  //Only set loci if you need them 
  cout<<"Number of Loci: "<<(annot.setLoci()).size()<<endl;

  return 0;
}

