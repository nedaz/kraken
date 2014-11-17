#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Kraken/KrakenMap.h"
#include "base/FileParser.h"
#include "src/Cola/Cola.h"

void GenomeWideMap::Read(const string & fileName, 
			 const string & source, 
			 const string & target, 
			 bool flip,
			 double distance)
{
  m_source = source;
  m_target = target;

  FlatFileParser parser;
  parser.Open(fileName);
  
  int i = -1;
  do {
    i++;
    if (i >= m_blocks.isize())
      m_blocks.resize(i + 1000000);
  } while (m_blocks[i].parse(parser, flip));
  m_blocks.resize(i+1);

  Sort(m_blocks);
  MergeBlocks();
  FILE_LOG(logDEBUG) << "Syntenic blocks: " << m_blocks.isize();
}

void GenomeWideMap::MergeBlocks()
{
  
  int i, j;
  svec<AlignmentBlock> tmp;
  tmp.resize(m_blocks.isize());
  j = 0;
  int k = 0;

  for (i=0; i<m_blocks.isize(); i++) {
    if (m_blocks[i].getOrient() == tmp[k].getOrient() &&
	m_blocks[i].getTargetChrom() == tmp[k].getTargetChrom() &&
	m_blocks[i].getQueryChrom() == tmp[k].getQueryChrom() &&	
	m_blocks[i].getTargetStart() < tmp[k].getTargetStop()) {      
      /*      if (m_blocks[i].getTargetStop() > tmp[k].getTargetStop())
	tmp[k].setStopT(m_blocks[i].getTargetStop());
      if (m_blocks[i].getQueryStop() > tmp[k].getQueryStop())
	tmp[k].setStopQ(m_blocks[i].getQueryStop());
      if (m_blocks[i].getQueryStart() < tmp[k].getQueryStart())
      tmp[k].setStartQ(m_blocks[i].getQueryStart());*/
      continue;
    }
    tmp[k] = m_blocks[i];
    k++;
  }

  tmp.resize(k);
  m_blocks = tmp;
  
}

//TODO Kraken functions are not const where they should be, underlying FuzzySearch.. functions need to be fixed first
bool GenomeWideMap::Map(const Coordinate & lookup, svec<Coordinate>& results)
{
  AlignmentBlock tmp;
  tmp.set(lookup.getChr(), lookup.getStart(), lookup.getStart());

  long long index = BinSearchFuzzy(m_blocks, tmp);
  if ((index == 0) || (index >= m_blocks.isize())) {
    FILE_LOG(logDEBUG3) << "Initial target region not found for lookup start - code2";
    return false;
  }

  AlignmentBlock begin = m_blocks[index-1];

  FILE_LOG(logDEBUG3) << "Index=" << index;
  FILE_LOG(logDEBUG3) << tmp.toString();
  FILE_LOG(logDEBUG3) << begin.toString();
  if (begin.getTargetChrom() != lookup.getChr()) {
    FILE_LOG(logDEBUG1) << "Not found synteny for start of source - Code1"; 
    FILE_LOG(logDEBUG3) << "Found block source chrom: " << begin.getTargetChrom() 
                        << "  look up chromosome: "<< lookup.getChr();
    return false;
  }

  tmp.set(lookup.getChr(), lookup.getStop(), lookup.getStop());
  index = BinSearchFuzzy(m_blocks, tmp);
  if ((index == 0) || (index >= m_blocks.isize())) {
    FILE_LOG(logDEBUG1) << "Initial target region not found for lookup stop - code3";
    return false;
  }

  AlignmentBlock end = m_blocks[index-1];
  // If lookup region is not covered then extend
  if ((end.getTargetStop() < lookup.getStop()) 
      && m_blocks[index].getTargetChrom()==lookup.getChr()) {  
    end = m_blocks[index]; 
  }

  if (end.getTargetChrom() != lookup.getChr()) {
    FILE_LOG(logDEBUG1) << "Not found synteny for end of source - code4"; 
    return false;
  }

  FILE_LOG(logDEBUG3) << "Index=" << index;
  FILE_LOG(logDEBUG3) << tmp.toString();
  FILE_LOG(logDEBUG3) << end.toString();

  bool split = false;
  if (begin.getQueryChrom() != end.getQueryChrom()) {
    FILE_LOG(logDEBUG2) << "Start & stop of initial target region not on the same Chromosome - Code5";
    FILE_LOG(logDEBUG3) << "Start: " << begin.toString() << " Stop: "<< end.toString(); 
    split = true; 
  }
  int threshold = max(100000, 10*lookup.findLength());
  if (abs(end.getQueryStop()-begin.getQueryStart()) > threshold 
     || abs(begin.getQueryStop() - end.getQueryStart()) > threshold) { //TODO parametrize
    FILE_LOG(logDEBUG1) << "Initial target region too big - code6: "
                        << end.getQueryStop() - begin.getQueryStart() << " "
                        << begin.getQueryStop() - end.getQueryStart();
    split = true;
  }
  if(split) {
    int size = lookup.getStop() - lookup.getStart();
    Coordinate result;
    SetAnchors(lookup, begin, end, 0, size, result);
    results.push_back(result);
    SetAnchors(lookup, begin, end, size, 0, result);
    results.push_back(result);
    return true;
  } 

  Coordinate result;
  SetAnchors(lookup, begin, end, 0, 0, result);
  results.push_back(result);
  
  return true;
}

void GenomeWideMap::SetAnchors(const Coordinate & lookup, const AlignmentBlock& begin, 
                               const AlignmentBlock& end, int startExtend, int stopExtend,
                               Coordinate & result) {

  AlignmentBlock beginAdj(begin), endAdj(end); //Blocks adjusted by considering reversed start/stop
  if(begin.isReversed() && end.isReversed()) {
    beginAdj = end;
    endAdj   = begin;
  }

  result.setChr(begin.getQueryChrom());
  result.setOrient(begin.getOrient());

  if(startExtend==0 && stopExtend==0) { 
    result.setStart(beginAdj.getQueryStart());
    result.setStop(endAdj.getQueryStop());
  } else if(startExtend!=0) {
    result.setStart(begin.getQueryStop() - startExtend);
    result.setStop(begin.getQueryStop() + startExtend);
  } else { //stopExtend!=0
    result.setChr(end.getQueryChrom()); //chromosome name should be end.chr
    result.setOrient(end.getOrient()); //chromosome orientation should be end.orient
    result.setStart(end.getQueryStop() - stopExtend);
    result.setStop(end.getQueryStop() + stopExtend);
  }
  
  if (result.getStop() < result.getStart()) {
    int tt = result.getStart();
    result.setStart(result.getStop());
    result.setStop(tt);
  }

  if (lookup.isReversed()) {
    if (result.isReversed())
      result.setOrient(true);  // Positive orientation
    else
      result.setOrient(false); // Negative orientation
  }

  FILE_LOG(logDEBUG2) << "Mapped "
                      << lookup.toString('\t')
                      << " LEN: " << lookup.getStop() - lookup.getStart() <<" "
                      << " to "
                      << result.toString('\t')
                      << " LEN: " << result.getStop() - result.getStart();
}

//==================================================

void Kraken::Allocate(const string & source, const string & target, double distance)
{
  GenomeWideMap tmp;
  tmp.Set(source, target, distance);
  m_maps.push_back(tmp);

  tmp.Set(target, source, distance);
  m_maps.push_back(tmp);

  m_seq.push_back(GenomeSeq(source));
  m_seq.push_back(GenomeSeq(target));
}

void Kraken::DoneAlloc()
{
  Sort(m_maps);
  UniqueSort(m_seq);
  for (int i=0; i<m_maps.isize(); i++)
    m_maps[i].Print();

  m_router.SetNumGenomes(m_seq.isize());
  int i, j;
  for (i=0; i<m_seq.isize(); i++) {
    for (j=0; j<m_seq.isize(); j++) {
      if (i == j)
	continue;
      Route out;
      // Find the route between genomes
      m_router.FindRoute(out, m_seq[i].Name(), m_seq[j].Name(), *this);
    }
  }
}

void Kraken::ReadMap(const string & fileName, const string & source, const string & target, double distance)
{
  int index = Index(source, target);
  bool bSort = false;
  if (index == -1) {
    FILE_LOG(logWARNING) << "Warning: allocating new map!";
    bSort = true;
    index = m_maps.isize();
    m_maps.resize(index+1);
  }
  
  m_maps[index].Read(fileName, source, target, false, distance);

  index = Index(target, source);
  if (index == -1) {
    FILE_LOG(logWARNING) << "Warning: allocating new map!";
    bSort = true;
    index = m_maps.isize();
    m_maps.resize(index+1);
  }

  m_maps[index].Read(fileName, target, source, true, distance);
  
  if (bSort)
    Sort(m_maps);
}
 
void Kraken::ReadGenome(const string & fileName, const string & name)
{
  UniqueSort(m_seq);
  int i = Genome(name);
  bool bSort = false;
  if (i == -1) {
    bSort = true;
    i = m_seq.isize();
    m_seq.resize(i+1);
    FILE_LOG(logWARNING) << "Warning: genome " << name 
                         << " has not been pre-allocated!!";
  }

  m_seq[i].Read(fileName, name);

  if (bSort)
    UniqueSort(m_seq);
}
 
bool Kraken::MapThroughRoute(const Route & route, svec<Coordinate>& results, const Coordinate & lookup)
{
  Coordinate tempLookup = lookup;
  FILE_LOG(logDEBUG4) << "Route: count=" << route.GetCount();
  int i;
  for (i=0; i<route.GetCount(); i++) {
    const string & source = route.Origin(i);
    const string & target = route.Destination(i);
    FILE_LOG(logDEBUG3) << source << " -> " << target;
  }

  for (i=0; i<route.GetCount(); i++) {
    results.clear();
    const string & source = route.Origin(i);
    const string & target = route.Destination(i);
    FILE_LOG(logDEBUG3) << "Mapping " << source << " and " << target;
    int index = Index(source, target);
    Coordinate tmp;
    if (!m_maps[index].Map(tempLookup, results)) {
      FILE_LOG(logDEBUG2) << "No map found between " << source << " and " << target;
      return false;
    }
    tempLookup = results[0]; //TODO all results might need to be looked at eveyr stage
  }
  
  FILE_LOG(logDEBUG2) << "Found " << results.size() << " set/sets of possible mappings";
  return (results.size()>0);
}

const GenomeWideMap & Kraken::GetMap(const string & name) const 
{
  double dist = 9999.;
  int index = -1;
  for (int i=0; i<m_maps.isize(); i++) {
    if (m_maps[i].Origin() == name) {
      if (m_maps[i].Distance() < dist) {
	dist = m_maps[i].Distance();
	index = i;
      }
    }
  }
  return m_maps[index];
}

bool Kraken::FindWithEdges(const Coordinate& lookup, const string & source,
                   const string & target, bool lAlign, int edgeLength,
                   Coordinate& result) 
{
    int from = lookup.getStart();
    int to   = lookup.getStop();
    // Don't continue if lookup is smaller or equal to 3 bases
    if (to-from <= 3)
      return false;
    // Reverse to/from if coordinates are given in reverse order
    if (to < from) {
      int tmp = to;
      to = from;
      from = tmp;
    }

    FILE_LOG(logDEBUG2) << "Mapping (left): ";
    Coordinate s, resLeft;
    s.set(lookup.getChr(), true, from, min(from + edgeLength, to));
    bool bLeft = Find(s, source, target, lAlign, resLeft);

    FILE_LOG(logDEBUG3) << "Mapping (right): ";
    Coordinate resRight;
    s.set(lookup.getChr(), true,  max(to - edgeLength, from), to);
    bool bRight = Find(s, source, target, lAlign, resRight);

    if (!bLeft && !bRight) {
      FILE_LOG(logDEBUG2) << "failed left and right mapping ";
      return false;
    }

    FILE_LOG(logDEBUG3) << "Right result: " << resRight.toString('\t');
    FILE_LOG(logDEBUG3) << "Left result: "  << resLeft.toString('\t');

    if (bRight && bLeft && resRight.getChr() != resLeft.getChr()) {
      FILE_LOG(logDEBUG3) << "left and right chromosome don't match!";
      if(resRight.findLength() >= resLeft.findLength()) { bLeft=false; }
      else { bRight=false; }
    }

    if (bLeft) {
      result.setChr(resLeft.getChr());
      if (!(resLeft.isReversed())) {
	result.setStart(resLeft.getStart());
      } else {
	result.setStop(resLeft.getStop());
	result.setOrient(false);
      }
    }
    if (bRight) {
      result.setChr(resRight.getChr());
      if (!resRight.isReversed()) {
	result.setStop(resRight.getStop());
      } else {
	result.setStart(resRight.getStart());
	result.setOrient(false);
      }
    }
    if (result.getStart() == -1) {
      result.setStart(result.getStop() - (to - from));
    }
    if (result.getStop() == -1) {
      result.setStop(result.getStart() + (to - from));
    }

    return true;
}

bool Kraken::Find(const Coordinate & lookup, 
               const string & source, const string & target,
               bool lAlign, Coordinate & result)
{
  Route route;
  if(! m_router.FindRoute(route, source, target, *this)) {
    FILE_LOG(logDEBUG2) << "NO route!";
    return false;
  }

  svec<Coordinate> results;
  if(!MapThroughRoute(route, results, lookup)) {
    FILE_LOG(logDEBUG2) << "Mapping failed!!!";
    return false;
  }
  
  int bestMaxPos=0, bestLen=0;
  float bestMaxVal=0;
  DNAVector sourceSeq, bestDestSeq; 
  for(int i=0; i<results.isize(); i++) {
    DNAVector targetSeq; 
    int maxPos=-1, len=-1;
    float maxVal=.0;
    if(RoughMap(lookup, source, target, 
                 sourceSeq, targetSeq, maxPos, 
                 maxVal, len, results[i])) {
      if(maxVal>bestMaxVal) {
        bestMaxVal  = maxVal;
        bestMaxPos  = maxPos;
        bestLen     = len;
        bestDestSeq = targetSeq;
        result      = results[i];
      }
    }
  }
  if(bestLen==0) { return false; } 

  int slack=12;
  if (bestMaxPos-slack < 0) { slack = bestMaxPos; }
  FILE_LOG(logDEBUG3) << "Slack for finer alignment: " << slack;
  DNAVector trueDestination;
  trueDestination.SetToSubOf(bestDestSeq, bestMaxPos-slack, bestLen+2*slack);
  return ExhaustAlign(trueDestination, sourceSeq, slack, 0.3, lAlign, result);
}
 
bool Kraken::RoughMap(const Coordinate& lookup, const string& source,
                   const string& target, DNAVector& sourceSeq, DNAVector& targetSeq, 
                   int& maxPos, float& maxVal, int& len, Coordinate& result) {
  int sourceIndex = Genome(source);
  int targetIndex = Genome(target);
  const vecDNAVector & sourceGenome = m_seq[sourceIndex].DNA();
  const vecDNAVector & targetGenome = m_seq[targetIndex].DNA();
  
  FILE_LOG(logDEBUG3) << sourceGenome.Name(0);
  FILE_LOG(logDEBUG3) << targetGenome.Name(0);
  FILE_LOG(logDEBUG3) << "Check."; 
  FILE_LOG(logDEBUG3) << "Raw Origin: "
                      << lookup.toString('\t')
                      << "Raw Destination:  " 
                      << result.toString('\t');  
  if(!sourceGenome.SetSequence(lookup, sourceSeq)) { return false; }
  result.setStart(result.getStart() - 5000);
  result.setStop(result.getStop() + 5000);
  if(!SetSequence(targetGenome, result, targetSeq)) { return false; }
  bool successAlign = RoughAlign(targetSeq, sourceSeq, maxPos, maxVal, len, result);  
  if(!successAlign) { return false; }
  FILE_LOG(logDEBUG2) << "Final Origin: " 
                      << lookup.toString('\t')
                      << "Final Destination:  "
                      << result.toString('\t');
  if ((maxPos + len >= targetSeq.isize())) {
    FILE_LOG(logDEBUG2) << "Out of bound sequence - Code9";
  }
  return true;
}

bool Kraken::SetSequence(const vecDNAVector& genome, Coordinate& coords, DNAVector& resultSeq) {
  if(!genome.HasChromosome(coords.getChr())) { 
    FILE_LOG(logWARNING) << "Check Genome data! - Chromosome: "  << coords.getChr() 
                         << " was not found in the given fasta file";
    return false; 
  }   
  int size = abs(coords.getStop() - coords.getStart());
  if(coords.getStart()<0) { coords.setStart(0); }
  if(coords.getStop()<0)  { coords.setStop(0); }
  if(coords.getStart()+1 > genome(coords.getChr()).isize()) {
    FILE_LOG(logDEBUG3) << "Limiting initial start to fit in with chromosome."; 
    coords.setStart(genome(coords.getChr()).isize()-1);
  }
  if(coords.getStop()+1 > genome(coords.getChr()).isize()) {
    FILE_LOG(logDEBUG3) << "Limiting initial stop to fit in with chromosome."; 
    coords.setStop(genome(coords.getChr()).isize()-1);
  }
  bool set = resultSeq.SetToSubOf(genome(coords.getChr()), coords.getStart(), coords.getStop()-coords.getStart()+1);
  if (coords.isReversed())
    resultSeq.ReverseComplement();
  return set;
}

bool Kraken::RoughAlign(DNAVector& q, DNAVector& t, 
                      int& maxPos, float& maxVal, int& len, Coordinate& result) {

  const int BLOCK_LIMIT   = 163840;
  const int BLOCK_OVERLAP = BLOCK_LIMIT/10;
  if(t.isize() > m_transSizeLimit_p) {  
      FILE_LOG(logWARNING) << "Requested region to be mapped: " 
                           << t.isize()<<" is too large";
      return false;
  }

  DNAVector qBlock  = q;
  int currStart     = 0;
  int currLen       = 0;
  do {
    currLen = min(BLOCK_LIMIT, q.isize()-currStart);
    qBlock.SetToSubOf(q, currStart, currLen);
    int size = m_xc.Size(t.isize(), qBlock.isize());
    float maxVal_temp;
    int maxPos_temp;
    Ccorrelate(qBlock, t, size, maxVal_temp, maxPos_temp);
    if(maxVal_temp > maxVal) { 
      maxVal = maxVal_temp;
      maxPos = currStart + maxPos_temp;
    }
    currStart += BLOCK_LIMIT - min(t.isize(), BLOCK_OVERLAP);
    FILE_LOG(logDEBUG2) << "q=" << qBlock.isize() << " t=" << t.isize() << " size=" << size;
    FILE_LOG(logDEBUG2) << "Remaining number of bases= " << q.isize()-currStart;
  } while(currStart < q.isize());

  if (maxPos >= q.isize()) {
    FILE_LOG(logDEBUG2) << "maxPos>q.size => No alignment.";
    return false;
  }

  if ( maxPos < 0) {
    FILE_LOG(logDEBUG3) << "Zero-ing maxPos " << maxPos;
    maxPos = 0;
  }

  len = t.isize();
  if (maxPos + t.isize() >= q.isize())
    len = q.isize() - maxPos;

  if (t.isize() >= q.isize()) {
    len = q.isize();
    FILE_LOG(logDEBUG3) << "Zero-ing maxPos (2): " << maxPos;
    maxPos = 0;
  }

  if (result.isReversed()) {
    result.setStart(result.getStart() + q.isize() - (maxPos + len));
    result.setStop(result.getStart() +len-1);    
  } else {
    result.setStart(result.getStart() + maxPos);
    result.setStop(result.getStart() +len-1);
  }

  FILE_LOG(logDEBUG2) << "best pos=" << maxPos << " max=" << maxVal 
                      << " tsize=" << t.isize() << " qsize=" << q.isize() << " len=" << len;

  if((float) maxVal/t.size() < 0.05) {
    FILE_LOG(logDEBUG2) << "Rejecting as cross-correlation maximum not significant - Code8";
    return false;
  } 

  return true;
}

void Kraken::Ccorrelate(const DNAVector& q, const DNAVector& t, double size, 
                        float& maxValOut, int& maxPosOut) {

  CCSignal sigsource, sigtarget;      
  sigsource.SetSequence(t, size);
  sigtarget.SetSequence(q, size);
  svec<float> signal;
  m_xc.CrossCorrelate(signal, sigtarget, sigsource);
  
  svec<float>::iterator it = max_element(signal.begin(), signal.end());
  maxValOut = *it;
  maxPosOut = distance(signal.begin(), it);
  maxPosOut = signal.isize()/2 - maxPosOut;
  
}

bool Kraken::ExhaustAlign(DNAVector& trueDestination, DNAVector& source,
                          int slack, float alignedRatio, 
                          bool localAlign, Coordinate& result) {
  Cola aligner;
  int bound;
  // Optimal align with a band of slack+5% of the query sequence size using Smithwaterman-gap-affine
  bound=min(source.size()/20, 20)+slack; 
  AlignmentCola pAlign =  aligner.createAlignment(source, trueDestination, 
                             AlignerParams(bound, SWGA, -5, -2, -1));

  FILE_LOG(logDEBUG3) << endl << pAlign.toString(100);
     
  int adjustBegin, adjustEnd; 
  if(localAlign) {
    adjustBegin = pAlign.getQueryOffset() - slack;
    adjustEnd   = pAlign.getQueryLength() - 2*slack - adjustBegin - pAlign.getQueryBaseAligned();
  } else {
    adjustBegin = pAlign.getQueryOffset() - slack - pAlign.getTargetOffset();
//  Full formula  adjustEnd   = (pAlign.getQueryLength() -2*slack - pAlign.getQueryOffset() + slack - pAlign.getAlignmentLen()) 
//                  - (pAlign.getTargetLength() - p.Align.getTargetOffset() - pAlign.getAlignmentLen());  
    adjustEnd   = pAlign.getQueryLength() - 2*slack - adjustBegin - max(pAlign.getTargetLength(), pAlign.getAlignmentLen());
  } 

  FILE_LOG(logDEBUG2) << "Adjust begin: " << adjustBegin << " Adjust end: " << adjustEnd;
  if (result.isReversed()){
    result.setStart(result.getStart() + adjustEnd);
    result.setStop(result.getStop() - adjustBegin);    
  } else {
    result.setStart(result.getStart() + adjustBegin);
    result.setStop(result.getStop() - adjustEnd);
  }

  double ratio = (double)pAlign.getTargetBaseAligned()/(double)source.isize();
  if (ratio<alignedRatio || pAlign.calcPVal()>m_pValThresh_p || pAlign.calcIdentityScore()<m_minIdent_p) {
    FILE_LOG(logDEBUG1) << "Rejecting...based on exhaustive alignment - Code7";
    return false;
  }

  return true;
}

int Kraken::Index(const string & source, const string & target)
{
  GenomeWideMap tmp;
  tmp.Set(source, target, 0.);

  int i = (int)BinSearch(m_maps, tmp);
  return i;
}

int Kraken::Genome(const string & name)
{
  GenomeSeq tmp(name);

  int i = (int)BinSearch(m_seq, tmp);
  return i;
}

bool RouteFinder::FindRoute(Route & out, const string & source, const string & target,  Kraken & rum)
{
  int iT = rum.Genome(source);
  int iQ = rum.Genome(target);
  int index = iQ + rum.GenomeCount() * iT;
  int i;

  Route & r = m_routes[index];
  if (r.IsInvalid()) {
    FILE_LOG(logDEBUG2) << "NO possible route from " << source << " to " << " target";
    return false;
  }
  if (r.GetCount() > 0) {
    out = r;
    return true;
  }
  svec<string> path;
  svec<string> final;
  path.reserve(100);

  FILE_LOG(logDEBUG1) << "Finding path from " << source << " to " << target;
  path.push_back(source);
  FindRecursive(final, path, target, rum);
  
  if (final.isize() == 1) {
    FILE_LOG(logDEBUG1) << "NO possible route from " << source << " to " << " target";
    r.SetInvalid();
    return false;
  }
  
  for (i=1; i<final.isize(); i++) {
    FILE_LOG(logDEBUG1) << "path: " << final[i-1] << " -> " << final[i];
    r.Add(final[i-1], final[i]);
  }
  
  return true;
}

bool RouteFinder::FindRecursive(svec<string> & final, svec<string> & path, const string & target, Kraken & rum) const
{
  string last = path[path.isize()-1];
  
  if (last == target)
    return true;

  int index = rum.Index(last, target);
  if (index != -1) {
    path.push_back(target);
    FILE_LOG(logDEBUG2) << "**** Found it!!";
    final = path;
    return true;
  }
  bool b = false;
  for (int i=0; i<rum.GenomeCount(); i++) {
    const string & sourceGenome = rum.GenomeName(i);
    int index = rum.Index(last, sourceGenome);
    if (index == -1)
      continue;
	    
    bool bBad = false;
    for (int j=0; j<path.isize(); j++) {
      if (path[j] == sourceGenome) {
	bBad = true;
	break;
      }
    }
    if (bBad)
      continue;
    
    path.push_back(sourceGenome);
    b = FindRecursive(final, path, target, rum);
    if (!b) {
      path.resize(path.isize()-1);
    }
    
  }
  return b;
}
