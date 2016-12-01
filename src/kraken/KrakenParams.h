#ifndef _KRAKENPARAMS_H_
#define _KRAKENPARAMS_H_
//======================================================
class KrakenParams
{
public:
  KrakenParams(bool laAdjust=false, bool ofAdjust=true, int transSizeLimit=200000,
               double pValThreshold=0.001, double minIdent=0.2, double minAlignCover=0.3,
               int alignmentBound=3, int minBasePerScaf=200)
              :m_laAdjust(laAdjust), m_ofAdjust(ofAdjust), m_transSizeLimit(transSizeLimit),
               m_pValThreshold(pValThreshold), m_minIdent(minIdent), m_minAlignCover(minAlignCover) {}
 
    bool    isLocalAlignAdjust() const  { return m_laAdjust;       }
    bool    isOverflowAdjust() const    { return m_ofAdjust;       } 
    int     getTransSizeLimit() const   { return m_transSizeLimit; }
    double  getPValThresh() const       { return m_pValThreshold;  }
    double  getMinIdent() const         { return m_minIdent;       }
    double  getMinAlignCover() const    { return m_minAlignCover;  }

    void    setLocalAlignAdjust(bool laa)    { m_laAdjust = laa;       }
    void    setOverflowAdjust(bool ofa)      { m_ofAdjust = ofa;       } 
    void    setTransSizeLimit(int tsl)       { m_transSizeLimit = tsl; }
    void    setPValThresh(double pvt)        { m_pValThreshold = pvt;  }
    void    setMinIdent(double mi)           { m_minIdent =  mi;       }
    void    setMinAlignCover( double mac)    { m_minAlignCover = mac;  }

private: 
  bool   m_laAdjust;          /// Choose if mapped region boundaries should be adjusted/limited with local alignment values
  bool   m_ofAdjust;          /// Choose if mapped regions should be adjusted to remove overflow/negative indexes
  int    m_transSizeLimit;    /// Limit on the maximum size of a region that will be translated
  double m_pValThreshold;     /// P-value threshold for acceptable alignment of translated region
  double m_minIdent;          /// Minimum alignment sequence identity acceptable for a translated region
  double m_minAlignCover;     /// Minimum acceptable portion of sequence covered by exhasustive alignment
 
};
//======================================================

#endif // _ASSEMBLYPARAMS_H_
