#include <string>

#include "ryggrad/src/base/CommandLineParser.h"
#include "MultiAlignParser.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmmd("-i" ,"MAF format input file");
  commandArg<string> bStringCmmd("-o","Output directory name - make sure directory exists", "Current Directory");
  commandArg<string> cStringCmmd("-s","Only produce pairwise syntenies for given specie", "Produce all pairwise syntenies");

  commandLineParser P(argc,argv);

  P.SetDescription("Parser that will produce pairwise alignment files from a given Maf file.");
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(cStringCmmd);
  P.parse();
  string inputFile     = P.GetStringValueFor(aStringCmmd);
  string outputName    = P.GetStringValueFor(bStringCmmd);
  string onlyOneSpecie = P.GetStringValueFor(cStringCmmd);
  
  FILELog::ReportingLevel() = logINFO; 

  if(outputName=="") { outputName = "."; }

  MultiAlignParser multiParser;
  svec<string> outFiles;
  multiParser.convertMAF(inputFile, outputName, outFiles, onlyOneSpecie);
  return 0;
}
  
