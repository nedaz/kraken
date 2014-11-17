#include <string>

#include "base/CommandLineParser.h"
#include "extern/logger/log.h"
#include "src/AnnotationQuery/AnnotationQuery.h"
#include "src/Kraken/GTFTransfer.h"


int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-c", "Configuration file");
  commandArg<string> bStringCmmd("-s", "Source GTF file");
  commandArg<string> dStringCmmd("-S", "Source genome id");
  commandArg<string> eStringCmmd("-T", "Target genome id");
  commandArg<string> cStringCmmd("-t", "Target annotation GTF file", "");
  commandArg<string> fStringCmmd("-o", "Mapping output in GTF format", "mapped.gtf");
  commandArg<string> gStringCmmd("-O", "Output from comparing mapped transcripts to targetAnnot", "compared");
  commandArg<string> hStringCmmd("-l", "Application logging file","application.log");
  commandArg<bool>   iStringCmmd("-L", "Choose if final boundaries should be set by local alignment", false);
  commandArg<double> jStringCmmd("-m", "Limit on the maximum size of a region that will be translated", 200000);
  commandArg<double> kStringCmmd("-p", "P-value threshold for acceptable alignment of translated region", 0.0001);
  commandArg<double> lStringCmmd("-i", "Minimum sequence identity acceptable for a translated region", 0.0);
  commandLineParser P(argc,argv);
  P.SetDescription("Batch mode GTF transfer/comparison from an source to target genome.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(dStringCmmd);
  P.registerArg(eStringCmmd);
  P.registerArg(cStringCmmd);
  P.registerArg(fStringCmmd);
  P.registerArg(gStringCmmd);
  P.registerArg(hStringCmmd);
  P.registerArg(iStringCmmd);
  P.registerArg(jStringCmmd);
  P.registerArg(kStringCmmd);
  P.registerArg(lStringCmmd);
  P.parse();
  string rumConfigFile    = P.GetStringValueFor(aStringCmmd);
  string sourceAnnotFile  = P.GetStringValueFor(bStringCmmd);
  string targetAnnotFile  = P.GetStringValueFor(cStringCmmd);
  string sourceGenomeId   = P.GetStringValueFor(dStringCmmd);
  string targetGenomeId   = P.GetStringValueFor(eStringCmmd);
  string outputGTFFileStr = P.GetStringValueFor(fStringCmmd);
  string outputCmpFileStr = P.GetStringValueFor(gStringCmmd);
  string applicationFile  = P.GetStringValueFor(hStringCmmd);
  bool   lAlign           = P.GetBoolValueFor(iStringCmmd);
  double transSizeLimit   = P.GetDoubleValueFor(jStringCmmd);
  double pValThreshold    = P.GetDoubleValueFor(kStringCmmd);
  double minIdent         = P.GetDoubleValueFor(lStringCmmd);
 
  FILE* pFile = fopen(applicationFile.c_str(), "w");
  Output2FILE::Stream()     = pFile;
  FILELog::ReportingLevel() = logINFO; 
  FILE_LOG(logINFO) <<"Running GTF transfer";
 
  GTFTransfer transer(rumConfigFile);
  transer.SetTransSizeLimit(transSizeLimit);
  transer.SetMinIdent(minIdent);
  transer.SetPValThresh(pValThreshold);
  TransAnnotation sourceAnnot = TransAnnotation(sourceAnnotFile, sourceGenomeId);
  
  // Map Transcripts onto corresponding exons and infer corresponding 
  // transcripts based on number of overlapping exons - nonoverlapping exons
  transer.translate(sourceAnnot, targetGenomeId, lAlign); 
  ofstream outGTFStream;
  outGTFStream.open(outputGTFFileStr.c_str(), ios_base::out);
  sourceAnnot.writeGTF(outGTFStream);
  FILE_LOG(logINFO) <<"Done - writing GTF output";
 
  if(targetAnnotFile != "") {
    Annotation targetAnnot = Annotation(targetAnnotFile, targetGenomeId);
    ofstream sout_exon, sout_trans, sout_gene;
    sout_exon.open((outputCmpFileStr + ".exons").c_str(), ios_base::out);
    sout_trans.open((outputCmpFileStr + ".transcripts").c_str(), ios_base::out);
    sout_gene.open((outputCmpFileStr + ".genes").c_str(), ios_base::out);
    transer.reportAllOverlaps(sourceAnnot, targetAnnot, AITEM, sout_exon); 
    transer.reportAllOverlaps(sourceAnnot, targetAnnot, TRANS, sout_trans); 
    transer.reportAllOverlaps(sourceAnnot, targetAnnot, GENE, sout_gene); 
  }
  return 0;
}
  


