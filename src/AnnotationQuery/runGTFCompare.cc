#include <string>

#include "base/CommandLineParser.h"
#include "src/AnnotationQuery/AnnotationQuery.h"


int main(int argc,char** argv)
{
  commandArg<string> bStringCmmd("-sourceAnnot","Source annotation GTF file");
  commandArg<string> cStringCmmd("-targetAnnot","Target GTF file");
  commandArg<int>    fieldTypeCmd("-f","The field type for what will be compared, i.e. 0:exons  1: transcripts, 2: genes", 0);
  commandLineParser P(argc,argv);
  P.SetDescription("Compare all transcripts/exons/genes from origin GTF to dest GTF and report overlaps.");
  P.registerArg(bStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(fieldTypeCmd);
  P.parse();
  string origAnnotFile  = P.GetStringValueFor(bStringCmmd);
  string destAnnotFile  = P.GetStringValueFor(cStringCmmd);
  int    fieldType      = P.GetIntValueFor(fieldTypeCmd);
  
  FILELog::ReportingLevel() = logINFO; 

  GTFCompare comparer;
  Annotation origAnnot = Annotation(origAnnotFile, "Source");
  Annotation destAnnot = Annotation(destAnnotFile, "Target");
  
  AnnotField qFieldType = AnnotField(fieldType); // choose AITEM for annotationItems, TRANS for transcripts, and GENE for genes
  comparer.reportAllOverlaps(origAnnot, destAnnot, qFieldType, cout);
  return 0;
}
  
